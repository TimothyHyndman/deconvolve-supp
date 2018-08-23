deconvolve-supp
==========

Supplementary files for the R package [*deconvolve*](https://github.com/timothyhyndman/deconvolve).

MATLAB Symmetric Error
------------

MATLAB files to perform deconvolution of data $W = X + U$, where $U$ is unknown and assumed symmetric. The R code in deconvolve can be unreliable so this has been offered as an alternative. Methods used are as indicated in the documentation for deconvolve.

### Syntax
```matlab
[y, x, Q] = decon_err_sym(W)
[y, x, Q] = decon_err_sym(W, xx)
[y, x, Q] = decon_err_sym(W, xx, m)
[y, x, Q] = decon_err_sym(W, xx, m, pmf)
[y, x, Q] = decon_err_sym(W, xx, m, pmf, bw)
[y, x, Q] = decon_err_sym(W, xx, m, pmf, bw, show_diagnostics)
```

### Arguments

All arguments except for `W` will use default values when they are set to `[]`.

* `W` - A vector of the univariate contaminated data. This argument must be supplied, however all other arguments are optional.
* `xx` - A vector of x values on which to compute the density. Defaults to `linspace(min(W), max(W), 100)`.
* `m` - The number of point masses used to estimate the distribution of $X$. Defaults to 10.
* `pmf` - If `true`, does not calculate `y`, instead only calculates `Q`. Defaults to `false`.
* `bw` - The bandwidth to use. Defaults to a bandwidth calculated using an appropriate plug-in estimator.
* `show_diagnostic` - If `true`, prints messages to screen showing the exitflags and minimum values of each of the optimizations performed. Intended to be used for development only.


### Examples

Typical usage
```matlab
load('framingham.mat')
[y, x, Q] = decon_err_sym(W)
```

Use $m=20$ and show diagnostic messages

```matlab
load('framingham.mat')
[y, x, Q] = decon_err_sym(W, [], 20, [], [], true)
```