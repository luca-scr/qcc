# qcc 3.0 (NOT ON CRAN)

- Modification of rules for out-of-control points. A subset of Western Eletric Rules (WER) have been implemented. See `qccRules()`.
- All functions in `qcc` now return an object with associated `print`, `summary`, and `plot` methods.
- Several modifications to plot and print appearances. In particular graphs are produced using `ggplot2` package, with layout obtained using `patchworks` packakge, and print/summary uses `cli` package.
- Added function `describe()` to statistically describe variables in a data.frame according to the type.
- `qccGroups()` now has first argument `data` to extract variables from 
a data frame.
- `sd.p()` calculation now takes into sizes.
- Removed usage of `attach/detach` from all examples and documentation.
- `plot.qcc()` allows to use `as.Date()` objects along x-axis
- Add option to remove the title from all plots using `title = NULL`.
- Plotting single data values does not issue an error.

# qcc 2.7 

- Created an html vignette entitled "A quick tour of qcc".
- Moved R News paper to documentation.
- Removed demos.
- Improved appearance of graphs.
- Control limits for p and np charts computed based on binomial quantiles (and not on normal approximation).
- Control limits for c chart computed based on Poisson quantiles (and not on normal approximation).
- Added head.start argument to cusum function.
- Implemented OC curves for R and S charts.
- Make xbar OC curves when `confidence.level` is used.
- Fix a bug in `sd.u`
- Fix typos in the documentation.
- Fix a bug in cex labels.
- Add a link in the `qcc-package.Rd` help page to yhat blog post describing how to implement the Western Eletric Rules (WER).
- `pareto.chart()` function now returns an object of class `paretoChart` with associated print and plot method.
- Included a link in the `qcc-package.Rd` to R News article.

# qcc 2.6

- Allow one sided specification limit in `process.capability()`. 
- Add an example in qcc man page showing how to add warning limits.

# qcc 2.5 

- Bug fix for plotting np-chart with variable sample size.
- Pdf files (previously included as vignettes) moved to inst/doc with corresponding index.html.
    
# qcc 2.4

- Added a demo showing how to extend the package by defining a new control chart, i.e. the standardized p chart.
- Fix a bug in `qcc()` to allow a user defined control chart.

# qcc 2.3

- Fix a bug in pareto.chart when compute percentages based on quantiles for the right axis. Added the argument plot which if sets to FALSE won't produce the chart but only return a frequency summary table.
- Added R-news article as vignette.
- Modified the `qcc.option()` function so changing a parameter in `.qcc.options` inside a function will make the modification persistent even outside the function environment.
- Bug fix in cusum function if a vector is provided as input with `sizes` > 1 but no `std.dev`.      

# qcc 2.2

- Removed code producing plots on process capability.

# qcc 2.1

- Some modifications and bug fixes to `pareto.chart()` function.

# qcc 2.0

- Multivariate control charts (T^2, T^2 for individual observations, Ellipse chart for bivariate data) have been included.
- `sd.xbar.one` now allows to estimate the standard deviation by using the moving range (default as in the previous versions) or the method suggested by Ryan (2000) using the scaled std deviation of the observations.
- `cusum()` and `ewma()` have been rewritten to be called directly with data and arguments. No need to create a `qcc` object before.
- Added a `NAMESPACE`, exporting all functions with names not starting with a dot `.`. 

# qcc 1.3

- Redefined `std.dev` argument in `qcc()` function. It now allows to give a numerical value or a sting identifying a method for estimating the standard deviation of a continuous process variable. Thus, these methods are only available to "xbar", "R", and "S" charts. For details see `help(qcc)`.
Functions involved in such change are `std.xbar()`, `std.R()`, `std.S()`, whereas the other `sd.*` functions only have "..." argument added.
- If xbar chart is needed and the number of observations is larger than 25 use "RMSDF" method for computing process standard deviation.
- Fixed bug in highlighting violating runs and out of control points when `chart.all = FALSE`.
- Changed background color for Control Charts.
- Changed position of control limit labels in Shewhat and CusSum charts.

# qcc 1.2

- Added functions to plot Shewhart g chart (geometric distribution): `stats.g()`, `sd.g()`, `limits.g()`.  
Contributed by Greg Snow (greg.snow@ihc.com).     
- Bug fix in `violating.runs()`
- Changes in `qcc.options()`:
  - `run.length` set by default at 7 (it was 5 previously);
  - `font.stats` and `cex.stats` control font and character expansion used to draw text at the bottom of control charts.
- Added qcc_Rnews.pdf paper in doc directory
   
# qcc 1.1

- Fixed some minor bugs
- Reworked on par settings to allow multiple figures
- Corrected typos in `*.Rd` files
   
# qcc 1.0 

- Package release on CRAN.
