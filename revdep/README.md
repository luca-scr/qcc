# Setup

## Platform

|setting  |value                        |
|:--------|:----------------------------|
|version  |R version 3.4.1 (2017-06-30) |
|system   |x86_64, darwin15.6.0         |
|ui       |X11                          |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|tz       |Europe/Rome                  |
|date     |2017-07-10                   |

## Packages

|package |*  |version |date       |source                  |
|:-------|:--|:-------|:----------|:-----------------------|
|qcc     |   |2.7     |2017-07-10 |local (luca-scr/qcc@NA) |

# Check results

9 packages

|package           |version | errors| warnings| notes|
|:-----------------|:-------|------:|--------:|-----:|
|EnvStats          |2.2.1   |      0|        0|     1|
|IPSUR             |1.5     |      0|        0|     3|
|IQCC              |0.6     |      0|        0|     2|
|mistat            |1.0-4   |      0|        0|     1|
|qcr               |1.0     |      0|        0|     0|
|RcmdrPlugin.IPSUR |0.2-1   |      0|        0|     2|
|RcmdrPlugin.qual  |2.2.6   |      0|        0|     3|
|rSARP             |1.0.0   |      0|        0|     0|
|SixSigma          |0.9-4   |      0|        0|     0|

## EnvStats (2.2.1)
Maintainer: Steven P. Millard <EnvStats@ProbStatInfo.com>

0 errors | 0 warnings | 1 note 

```
checking installed package size ... NOTE
  installed size is  7.8Mb
  sub-directories of 1Mb or more:
    doc    3.2Mb
    help   3.3Mb
```

## IPSUR (1.5)
Maintainer: G. Jay Kerns <gkerns@ysu.edu>

0 errors | 0 warnings | 3 notes

```
checking dependencies in R code ... NOTE
'library' or 'require' call to ‘HH’ in package code.
  Please use :: or requireNamespace() instead.
  See section 'Suggested packages' in the 'Writing R Extensions' manual.

checking S3 generic/method consistency ... NOTE
Found the following apparent S3 methods exported but not registered:
  plot.htest
See section ‘Registering S3 methods’ in the ‘Writing R Extensions’
manual.

checking R code for possible problems ... NOTE
clt1: no visible global function definition for ‘graphics.off’
clt1: no visible global function definition for ‘curve’
clt1 : <anonymous>: no visible global function definition for ‘dt’
clt1: no visible global function definition for ‘abline’
clt1: no visible global function definition for ‘text’
clt1: no visible global function definition for ‘dt’
clt1: no visible global function definition for ‘var’
clt1: no visible global function definition for ‘dev.new’
clt1: no visible global function definition for ‘dev.set’
... 39 lines ...
read: no visible global function definition for ‘vignette’
Undefined global functions or variables:
  abline curve dev.new dev.set dgamma dnorm dt dunif graphics.off hist
  locator normal.and.t.dist optimize sd text var vignette
Consider adding
  importFrom("graphics", "abline", "curve", "hist", "locator", "text")
  importFrom("grDevices", "dev.new", "dev.set", "graphics.off")
  importFrom("stats", "dgamma", "dnorm", "dt", "dunif", "optimize", "sd",
             "var")
  importFrom("utils", "vignette")
to your NAMESPACE file.
```

## IQCC (0.6)
Maintainer: Emanuel P. Barbosa <emanuel@ime.unicamp.br>

0 errors | 0 warnings | 2 notes

```
checking package namespace information ... NOTE
  Namespace with empty importFrom: ‘miscTools’

checking R code for possible problems ... NOTE
add.data: no visible global function definition for ‘points’
add.data: no visible global function definition for ‘lines’
alpha.risk : risco: no visible global function definition for ‘ptukey’
cchart.R: no visible global function definition for ‘mtext’
cchart.R: no visible global function definition for ‘qtukey’
cchart.S: no visible global function definition for ‘qchisq’
cchart.T2.1: no visible global function definition for ‘qbeta’
cchart.T2.1: no visible global function definition for ‘plot’
cchart.T2.1: no visible global function definition for ‘qf’
... 18 lines ...
table.qtukey : g: no visible global function definition for ‘qtukey’
table.qtukey : d: no visible global function definition for ‘qtukey’
Undefined global functions or variables:
  abline axis integrate lines mtext par plot points ptukey qbeta qchisq
  qf qtukey title
Consider adding
  importFrom("graphics", "abline", "axis", "lines", "mtext", "par",
             "plot", "points", "title")
  importFrom("stats", "integrate", "ptukey", "qbeta", "qchisq", "qf",
             "qtukey")
to your NAMESPACE file.
```

## mistat (1.0-4)
Maintainer: Daniele Amberti <daniele.amberti@gmail.com>

0 errors | 0 warnings | 1 note 

```
checking top-level files ... NOTE
Non-standard file/directory found at top level:
  ‘README.html’
```

## qcr (1.0)
Maintainer: Miguel Flores <ma.flores@outlook.com>

0 errors | 0 warnings | 0 notes

## RcmdrPlugin.IPSUR (0.2-1)
Maintainer: G. Jay Kerns <gkerns@ysu.edu>

0 errors | 0 warnings | 2 notes

```
checking dependencies in R code ... NOTE
'library' or 'require' calls in package code:
  ‘abind’ ‘e1071’ ‘qcc’
  Please use :: or requireNamespace() instead.
  See section 'Suggested packages' in the 'Writing R Extensions' manual.

checking S3 generic/method consistency ... NOTE
Found the following apparent S3 methods exported but not registered:
  print.numSummaryIPSUR
See section ‘Registering S3 methods’ in the ‘Writing R Extensions’
manual.
```

## RcmdrPlugin.qual (2.2.6)
Maintainer: Erin Hodgess <hodgesse@uhd.edu>

0 errors | 0 warnings | 3 notes

```
checking package namespace information ... NOTE
  Namespaces with empty importFrom:
  ‘stats’ ‘tcltk’

checking DESCRIPTION meta-information ... NOTE
Malformed Description field: should contain one or more complete sentences.

checking R code for possible problems ... NOTE
.onAttach: no visible global function definition for ‘packageVersion’
.onAttach: no visible global function definition for ‘.gettext’
cchart: no visible binding for global variable ‘top’
cchart: no visible binding for global variable ‘buttonsFrame’
cschart: no visible binding for global variable ‘top’
cschart : onOK: no visible binding for global variable ‘top’
cschart: no visible binding for global variable ‘buttonsFrame’
cusumMod: no visible binding for global variable ‘top’
cusumMod : onOK: no visible binding for global variable ‘top’
... 44 lines ...
xbarnew: no visible binding for global variable ‘top’
xbarnew : onOK: no visible binding for global variable ‘top’
xbarnew: no visible binding for global variable ‘buttonsFrame’
Undefined global functions or variables:
  .gettext abline barplot buttonsFrame dweibull fitdistr ks.test
  newDataSet packageVersion plot qcc.options sd top
Consider adding
  importFrom("graphics", "abline", "barplot", "plot")
  importFrom("stats", "dweibull", "ks.test", "sd")
  importFrom("utils", "packageVersion")
to your NAMESPACE file.
```

## rSARP (1.0.0)
Maintainer: John Hutcheson <jacknx8a@gmail.com>

0 errors | 0 warnings | 0 notes

## SixSigma (0.9-4)
Maintainer: Emilio L. Cano <emilio.lcano@uclm.es>

0 errors | 0 warnings | 0 notes

