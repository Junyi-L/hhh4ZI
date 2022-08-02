#' @title Measles in the 16 states of Germany
#' @description Bi-weekly number of measles cases in the 16 states (Bundeslaender) of Germany
#' for the years 2005 to 2018.
#' @details
#' The \code{population} slot contains the population fractions
#' of each state from 2005 to 2018, obtained from the Federal Statistical Office of Germany.
#' @source SurvStat@RKI 2.0 (https://survstat.rki.de),
#' Robert Koch Institute; Queried on 26 March 2021.
#'
#' Statistisches Bundesamt (https://www.destatis.de/DE/Home/_inhalt.html);
#' Queried on 10 April 2021.
"measles"

#' @title MMR coverage levels in the 16 states of Germany
#' @description Coverage levels at school entry for the first and second dose
#' of the combined measles-mumps-rubella (MMR) vaccine for the years 2005 to 2018,
#' estimated from children presenting vaccination documents at school entry examinations.
#' @name Coverage
#' @format Four data frames, each with 14 rows (years) and 16 columns (states):
#' \itemize{
#' \item{TotalKids} Number of children examined.
#' \item{VacPass} Number of children who presented vaccination documents.
#' \item{Dosis1} Percentage of children with vaccination documents,
#' who received at least 1 dose of MMR vaccine.
#' \item{Dosis2} Percentage of children with vaccination documents,
#'  who received at least 2 doses of MMR vaccine.
#' }
#'
#' Coverage levels were derived from vaccination documents presented at medical examinations,
#' which are conducted by local health authorities at school entry each year.
#' Records include information about the receipt of 1st and 2nd doses of MMR,
#' from year 2005 to 2018. Note that information from children who did not
#' present a vaccination document on the day of the medical examination,
#' is not included in the estimated coverage.
#' @source https://www.gbe-bund.de,
#' Gesundheitsbericherstattung des Bundes; Queried on 23 February 2021.
NULL

#' @rdname Coverage
#' @format NULL
"TotalKids"

#' @rdname Coverage
#' @format NULL
"VacPass"

#' @rdname Coverage
#' @format NULL
"Dosis1"

#' @rdname Coverage
#' @format NULL
"Dosis2"
