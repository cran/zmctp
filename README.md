# zmctp: Zero-Modified Complex Triparametric Parametric Distribution

[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)](https://github.com/roladoja/zmctp/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

## Overview

The **zmctp** package extends the Complex Triparametric Pearson (CTP) distribution with zero-modified versions for handling overdispersed count data. It addresses limitations in existing implementations when the parameter *b* approaches zero.

### Key Features

- 🎯 **Robust CTP fitting** - Handles cases where *b* → 0
- ⭐ **Zero-Modified CTP** - Unique feature for zero-inflated/deflated data
- 📊 **Complete S3 methods** - print, summary, and plot functions
- 📈 **Diagnostic tools** - Goodness-of-fit statistics and visualizations
- 📚 **Comprehensive documentation** - Vignette with examples

## Installation
```r
# Install from GitHub (when available)
# devtools::install_github("yourusername/zmctp")

# Or install from source
devtools::install_local("path/to/zmctp")
```

## Quick Start
```r
library(zmctp)

# Generate data
x <- rctp(200, a = 1, b = 0.5, gama = 5)

# Fit CTP model
fit <- ctp.fit(x)
print(fit)
plot(fit)

# Fit Zero-Modified CTP
x_zi <- rzictp(200, a = 1, b = 0.5, gama = 5, omega = 0.3)
fit_zi <- zictp.fit(x_zi)
plot(fit_zi)
```

## Why zmctp?

Existing implementations (e.g., the `cpd` package) struggle when *b* ≈ 0, often estimating *b* = 0 which reduces model flexibility. The `zmctp` package solves this through:

1. **Reparameterization** - Ensures variance constraint is always satisfied
2. **Zero-Modified variant** - Explicitly models zero-inflation/deflation
3. **Better optimization** - Robust default starting values

## Documentation

- **Vignette**: Run `vignette("introduction", package = "zmctp")`
- **Help**: `?ctp.fit`, `?zictp.fit`, `?dctp`

## Example: Comparison with cpd
```r
library(cpd)
library(zmctp)

# Data where cpd estimates b ≈ 0
x <- rzictp(200, a = 1, b = 0.001, gama = 8, omega = 0.2)

# cpd may fail
fit_cpd <- cpd::fitCTP(x)
# b estimate ≈ 0

# zmctp handles it better
fit_zmctp <- zictp.fit(x)
# Recovers both b and omega
```

## Citation

If you use this package, please cite:
```
@Manual{zmctp,
  title = {zmctp: Zero-Modified Complex Triparametric Pearson Distribution},
  author = {Rasheedat Oladoja},
  year = {2025},
  note = {R package version 0.1.0},
}
```

And the original CTP paper:
```
@article{rodriguez2003,
  title={A new class of discrete distributions with complex parameters},
  author={Rodríguez-Avi, J and Conde-Sánchez, A and Sáez-Castillo, AJ},
  journal={Statistical Papers},
  volume={44},
  pages={67--88},
  year={2003},
  doi={10.1007/s00362-002-0134-7}
}
```

## License

GPL-3

## Author

Rasheedat Oladoja - roladoja@ttu.edu

