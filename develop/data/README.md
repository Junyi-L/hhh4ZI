# Data preparation
## Description
This folder contains code and raw data to prepare data in folder `hhh4ZI/data`.

* `survstat/` contains downloaded files from SurvStat@RKI 2.0, Robert Koch Institute (https://survstat.rki.de, accessed 26 March 2021).
    * `survstat/Info.pdf` contains the search query to download measles data.
    * `survstat/Data.csv` is the downloaded raw data of measles cases.
* `data.xlsx` is the converted version of `survstat/Data.csv`.
* `Impfquoten.xlsx` is the raw data of vaccination coverage downloaded from the Information System of the Federal Health Monitoring (https://www.gbe-bund.de/, accessed 23 February 2021).
* `coverage.xlsx` is the cleaned version of `Impfquoten.xlsx`.
* `population1.xlsx` is the raw data of state population downloaded from the Federal Statistical Office of Germany (Statistisches Bundesamt, https://www.destatis.de/, accessed 10 April 2021).
* `population.xlsx` is the cleaned version of `population1.xlsx`.
* `load_measles_data.R` converts above data to `Coverage.rda` and `measles.rda` in `hhh4ZI/data`.

