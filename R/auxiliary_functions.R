#' @title Measles in the 16 states of Germany
#' @description Weekly number of measles cases in the 16 states (Bundeslaender) of Germany
#' for years 2005 to 2019.
#'
#' @name measles
#' @docType data
#' @usage data(measles)
#' @format An \code{sts} object containing \eqn{749\times 16}{749 x 16} observations starting
#' from week 1 in 2005.The \code{population} slot contains the population fractions
#' of each state at 31.12.2006, obtained from the Federal Statistical Office of Germany.
#' @source SurvStat@RKI 2.0 (https://survstat.rki.de),
#' Robert Koch Institute; Queried on 10 April 2021.
#' @keywords datasets
NULL

#' @title MMR coverage levels in the 16 states of Germany
#' @description Coverage levels at school entry for the first and second dose
#' of the combined measles-mumps-rubella (MMR) vaccine in 2005 to 2018,
#' estimated from children presenting vaccination documents at school entry examinations.
#'
#' @name Coverage
#' @docType data
#' @usage data(Coverage)
#' @format Four \code{data.frame},
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
#' @keywords datasets
NULL
