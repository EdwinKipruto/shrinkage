# Post-estimation shrinkage methods
## Edwin Kipruto, Willi Sauerbrei
### Maintained by Edwin Kipruto
This project implements post-estimation shrinkage methods: global shrinkage, parameter-wise shrinkage (PWS), nonnegative PWS (NPWS), and quadratic PWS (QPWS), as described in Kipruto and Sauerbrei (2024). It also includes tools for running simulations to compare post-estimation shrinkage approaches with penalized regression methods (lasso and ridge).

For our discussion paper on "Post-Estimation Shrinkage in Full and Selected Linear Regression Models in Low-Dimensional Data Revisited," see: [provide link].

### Install the R package

To install the shrinkage R package and view the vignette run the
following codes in R: 

```{r}
devtools::install_github(repo="EdwinKipruto/shrinkage", build_vignettes = TRUE)
vignette("shrinkage")
```
