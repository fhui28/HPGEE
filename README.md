# HPGEEs (Simultaneous homogeneity pursuit and variable selection in regression models for multivariate abundance data)

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

<!-- badges: end -->

This repository code contain (template) code associated with the manuscript "Simultaneous homogeneity pursuit and variable selection in regression models for multivariate abundance data" by [Hui](https://francishui.netlify.app/), [Maestrini](https://sites.google.com/view/lucamaestrini), and [Welsh](https://cbe.anu.edu.au/about/staff-directory/professor-alan-welsh).

# Getting started

There are currently three directories in this repository

-   `code`, which contains `R` scripts implementing the proposed HPGEE method along with other scripts for simulating multivariate abundance data and fitting adaptive lasso GEEs. Many of the functions in those scripts contain pseudo-help files.

-   `simulations/setting1`, which contains template scripts to implement the simulation study in the manuscript.

-   `application_GBR`, which contains scripts applying the proposed method to the Great Barrier Reef data in the manuscript. Note the original data sources is from [Pichler et al., 2007](http://www.frdc.com.au/Archived-Reports/FRDC%20Projects/2003-021-DLD.pdf). We provide an example dataset which possess the same structure as the actual data used for analysis in the manuscript, but the values are altered to mask their original values.


# If you find any bugs and issues...

If you find something that looks like a bug/issue, you can make use of the Github issues and post it up there. As much as possible, please include in the issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem e.g., using the [reprex](https://cran.r-project.org/web/packages/reprex/index.html) package. If you also have an idea of how to fix the problem (Francis tends to make a lot of mistakes in my code, so some may be easy amendments!), then that is also much appreciated.
3.  Required data files etc...

Alternatively, please contact the corresponding author at [francis.hui\@anu.edu.au](mailto:francis.hui@anu.edu.au){.email}
