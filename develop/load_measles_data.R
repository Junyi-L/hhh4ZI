#!/usr/bin/env Rscript
## this script generates data/measles.rda and data/Coverage.rda

library("surveillance")
library("openxlsx")
library("here")

measles <- read.xlsx(here("develop", "data.xlsx"))
measles$Unknown <- NULL
measles[is.na(measles)] <- 0
names(measles)[1] <- "date"
measles$year <- as.numeric(substr(measles$date, 1, 4))
#measles$week <- as.numeric(substr(measles$date, 7, 8))

start_year <- which(measles$year == 2005)[1]
# year 2005 - 2018
measles <- measles[start_year : (start_year + 14*52 - 1),]
measles[ ,c("date", "year")] <- NULL

load(system.file("shapes", "districtsD.RData", package = "surveillance"))
#plot(districtsD)

statesD <- rgeos::gUnaryUnion(districtsD,
                              id = factor(substr(districtsD$KEY, 1, 2)))
# summary(statesD)
# plot(statesD)
# text(coordinates(statesD), labels = row.names(statesD))

row.names(statesD) <- namesDE <- c(
  "Schleswig-Holstein",
  "Hamburg",
  "Lower Saxony",
  "Bremen",
  "North Rhine-Westphalia",
  "Hesse",
  "Rhineland-Palatinate",
  "Baden-Wuerttemberg",
  "Bavaria",
  "Saarland",
  "Berlin",
  "Brandenburg",
  "Mecklenburg-Western Pomerania",
  "Saxony",
  "Saxony-Anhalt",
  "Thuringia"
)

## add population as map data
namesDE <- sort(namesDE)
Population <- read.xlsx(here("develop", "population.xlsx"),
                        sheet = 3, startRow = 2, rowNames = TRUE)
colnames(Population) <- namesDE[order(colnames(Population))]
Population <- Population[, row.names(statesD)]
rownames(Population) <- paste0("POP", rownames(Population))
statesD <- SpatialPolygonsDataFrame(statesD, as.data.frame(t(Population)))

## compute matrix of adjacency orders
adjmat <- poly2adjmat(statesD)
adjord <- nbOrder(adjmat, maxlag = Inf)

## sync columns
colnames(measles) <- namesDE[order(colnames(measles))]
measles <- measles[, colnames(adjord)] # reorder

## create "sts" object
measles <- sts(observed = measles, start = c(2005, 1), frequency = 52,
               map = statesD, neighbourhood = adjord)
# aggregate to bi-weekly counts
measles <- aggregate(measles, by = "time", nfreq = 26)

## add population fractions to sts object
Population <- Population[rep(1 : nrow(Population), each = 26),]
Population <- Population/rowSums(Population)
measles@populationFrac <- as.matrix(Population, rownames = FALSE)

## done with measles data
save(measles, file = here("data", "measles.rda"))


### load vaccination coverage

TotalKids <- read.xlsx(here("develop", "coverage.xlsx"),
                       sheet = 1, startRow = 2, rowNames = TRUE, sep.names = " ")
TotalKids <- TotalKids[, colnames(measles)]

VacPass <- read.xlsx(here("develop", "coverage.xlsx"),
                     sheet = 2, startRow = 2, rowNames = TRUE, sep.names = " ")
VacPass <- VacPass[, colnames(measles)]

Dosis1 <- read.xlsx(here("develop", "coverage.xlsx"),
                    sheet = 3, startRow = 2, rowNames = TRUE, sep.names = " ")
Dosis1 <- Dosis1[, colnames(measles)]

Dosis2 <- read.xlsx(here("develop", "coverage.xlsx"),
                    sheet = 4, startRow = 2, rowNames = TRUE, sep.names = " ")
Dosis2 <- Dosis2[, colnames(measles)]


save(list = c("TotalKids", "VacPass",
              "Dosis1", "Dosis2"), file = here("./data/Coverage.rda"))
