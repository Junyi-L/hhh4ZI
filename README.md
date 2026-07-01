# hhh4ZI

The [R](https://www.R-project.org/) package **hhh4ZI** implements
an extension of the endemic-epidemic model from
[`surveillance::hhh4()`](https://CRAN.R-project.org/package=surveillance/vignettes/hhh4_spacetime.pdf),
adding a zero-inflation component as described in `citation("hhh4ZI")`:

> Lu J, Meyer S (2023).
> "A zero-inflated endemic-epidemic model with an application to measles time series in Germany."
> *Biometrical Journal*, **65**(8), 2100408.
> <https://doi.org/10.1002/bimj.202100408>

## Installation

- Directly from the Git repository at Codeberg, using the **remotes** package:

  ```R
  ## install.packages("remotes")
  remotes::install_git("https://codeberg.org/EE-hub/hhh4ZI")
  ```

- Or using the automated builds from <https://ee-lib.r-universe.dev/hhh4ZI>:

  ```R
  setRepositories(name = "CRAN", addURLs = "https://ee-lib.r-universe.dev")
  install.packages("hhh4ZI")
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
