# hhh4ZI

The [R](https://www.R-project.org/) package **hhh4ZI** implements the
statistical model introduced by:

> Lu J, Meyer S (2023).
> "A zero-inflated endemic-epidemic model with an application to measles time series in Germany."
> *Biometrical Journal*, **65**(8), 2100408.
> <https://doi.org/10.1002/bimj.202100408>

## Installation

To install the **hhh4ZI** package from the GitHub repository, use:

```R
## install.packages("remotes")
remotes::install_github("Junyi-L/hhh4ZI")
```

## Example

A [demo script](demo/measles.R) of `hhh4ZI()` modelling of measles in
Germany (`data("measles")`) is available as

```R
system.file("demo", "measles.R", package = "hhh4ZI")
```

in the installed package.
To run the *whole* script (which takes around 1 hour), you could use:

```R
library("hhh4ZI")
demo(measles, ask = FALSE)
```
