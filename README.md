# hhh4ZI

The [R](https://www.R-project.org/) package **hhh4ZI** implements the
statistical model introduced in the manuscript
*A zero-inflated endemic-epidemic model
with an application to measles time series in Germany*
by Junyi Lu and Sebastian Meyer.

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
