library("surveillance")
library("xlsx")
library("here")

#measles <- read.csv(file = "measles.csv", header = TRUE)

measles <- read.xlsx(file = here("./develop/data.xlsx"),sheetIndex = 1, header = TRUE)
measles$Unknown <- NULL
measles[is.na(measles)] <- 0
names(measles)[1] <- "date"
measles$year <- as.numeric(substr(measles$date, 1, 4))
measles$week <- as.numeric(substr(measles$date, 7, 8))

start_year <- which(measles$year == 2005)[1]
# year 2005 - 2018
measles <- measles[start_year : (start_year + 14*52 - 1),]
measles[ ,c("date", "year")] <- NULL

load(system.file("shapes/districtsD.RData", package = "surveillance"))
#plot(districtsD)

statesD <- rgeos::gUnaryUnion(districtsD,
                              id = factor(substr(districtsD$KEY, 1, 2)))
# summary(statesD)
# plot(statesD)
# text(coordinates(statesD), labels = row.names(statesD))

row.names(statesD) <- namesD <- c(
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

adjmat <- poly2adjmat(statesD)
adjord <- nbOrder(adjmat, maxlag = Inf)

namesDE <- colnames(adjord)
namesDE <- namesDE[order(namesDE)]
colnames(measles) <- namesDE[order(colnames(measles))]
measles <- measles[, colnames(adjord)] # reorder

measles <- sts(observed = measles, start = c(2005, 1),
                 neighbourhood = adjord)
# aggregate data when using
measles <- aggregate(measles, by = "time", nfreq = 26)

measles@map <- statesD
#measles@neighbourhood <- adjord

Population <- read.xlsx(here("./develop/population.xlsx"),
                        sheetIndex = 3, colIndex = 2:17, startRow = 2)
colnames(Population) <- namesDE[order(colnames(Population))]
Population <- Population[, colnames(adjord)]
rownames(Population) <- 2005:2018
Population <- Population[rep(1 : nrow(Population), each = 26),]
Population <- Population/rowSums(Population)
measles@populationFrac <- as.matrix(Population)

save(measles, file = here("./data/measles.rda"))


TotalKids <- read.xlsx(here("./develop/coverage.xlsx"),
                sheetIndex = 1, colIndex = 2:17, startRow = 2)
colnames(TotalKids) <- namesDE[order(colnames(TotalKids))]
TotalKids <- TotalKids[, colnames(adjord)]
rownames(TotalKids) <- 2005:2018

VacPass <- read.xlsx(here("./develop/coverage.xlsx"),
                       sheetIndex = 2, colIndex = 2:17, startRow = 2)
colnames(VacPass) <- namesDE[order(colnames(VacPass))]
VacPass <- VacPass[, colnames(adjord)]
rownames(VacPass) <- 2005:2018

Dosis1 <- read.xlsx(here("./develop/coverage.xlsx"),
                     sheetIndex = 3, colIndex = 2:17, startRow = 2)
colnames(Dosis1) <- namesDE[order(colnames(Dosis1))]
Dosis1 <- Dosis1[, colnames(adjord)]
rownames(Dosis1) <- 2005:2018

Dosis2 <- read.xlsx(here("./develop/coverage.xlsx"),
                     sheetIndex = 4, colIndex = 2:17, startRow = 2)
colnames(Dosis2) <- namesDE[order(colnames(Dosis2))]
Dosis2 <- Dosis2[, colnames(adjord)]
rownames(Dosis2) <- 2005:2018


save(list = c("TotalKids", "VacPass",
              "Dosis1", "Dosis2"), file = here("./data/Coverage.rda"))


