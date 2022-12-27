Deconvolution in 1D
================

-   <a href="#quick-start" id="toc-quick-start">Quick start</a>
-   <a href="#model" id="toc-model">Model</a>
-   <a href="#preprocessing-dna-dna-matrix"
    id="toc-preprocessing-dna-dna-matrix">Preprocessing: DNA-DNA matrix</a>
    -   <a href="#normalization-and-smoothing"
        id="toc-normalization-and-smoothing">Normalization and smoothing</a>
    -   <a href="#setting-gamma-the-amount-of-deconvolution"
        id="toc-setting-gamma-the-amount-of-deconvolution">Setting gamma, the
        amount of deconvolution</a>
-   <a href="#simulated-data" id="toc-simulated-data">Simulated data</a>
    -   <a href="#gamma-and-degrees-of-freedom"
        id="toc-gamma-and-degrees-of-freedom">Gamma and degrees of freedom</a>
        -   <a href="#using-true-df" id="toc-using-true-df">Using true df</a>
        -   <a href="#multiple-peaks" id="toc-multiple-peaks">Multiple peaks</a>
    -   <a href="#plots-for-show" id="toc-plots-for-show">Plots for show</a>

This vignette walks through deconvolution in 1D.

``` r
library(data.table)
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:data.table':

        between, first, last

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
library(tidyr)
library(ggplot2)
library(dcon)
set.seed(13579)
```

# Quick start

We first provide a quick example to show how the method works. For more
details, see everything below the quick start section.

Here, we use a small sample of the RD-SPRITE data
([GSE151515](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151515)).

We focus on the lncRNA *Airn* which is a part of the *Igf2r* imprinted
cluster on chromosome 17 (mm10) and is thought to orchestrate silencing
of the gene *Slc22a3*.

``` r
airn <- fread('/rafalab/lzou/rdsprite/count_windows/Airn_RNA_count_windows_10000_clusters_2-1000_chr17_129S1_SvImJ.csv') |>
  filter(chrom=='chr17')
```

``` r
airn |>
  filter(start <= 30e6) |>
  ggplot(aes(x = start/1e6, y = count)) +
  geom_rect(xmin = 12739311/1e6, xmax = 12861884/1e6,
            ymin = 0, ymax = max(airn$count),
            fill = 'lightpink') +
  geom_rect(xmin = 12417972/1e6, xmax = 12509704/1e6,
            ymin = 0, ymax = max(airn$count),
            fill = 'lightgreen') +
  geom_point(size=0.1) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  theme_classic() +
  xlab('position along chr17 (mb)') +
  ylab('Airn lncRNA raw count')
```

![](demo_1d_files/figure-gfm/unnamed-chunk-2-1.png)

# Model

The model for deconvolution is

$$Y \sim \text{Poisson}(D\lambda)$$

$$\text{log}(\lambda) = B\alpha$$

where $D$ represents the DNA-DNA contact matrix and $\lambda$ the true,
underlying signal. $B$ are Gaussian basis functions and $\alpha$ are the
coefficients to be estimated.

# Preprocessing: DNA-DNA matrix

Before we deconvolve, it is important to inspect the DNA-DNA contact
matrix and preprocess it. The default parameters should work in many
cases, but it is also good to manually inspect what is happening.

Here, we use an example DNA matrix, taken from the *Kcnq1ot1* imprinted
region (a 1mb region around the *Kcnq1ot1* locus) at 10kb bin
resolution.

``` r
dd <- readRDS(system.file("extdata", "rdsprite_kcnq1ot1_domain_DD.rds", package="dcon"))
dd <- as.matrix(dd)
plot_mat(dd, main = expression(paste('Raw DNA-DNA contacts at the ', italic('Kcnq1ot1'), ' locus (1mb)')))
```

![](demo_1d_files/figure-gfm/unnamed-chunk-3-1.png)

## Normalization and smoothing

There are a few important tuning parameters to normalize and smooth the
DNA-DNA contact matrix with `normalize_hic`.

1.  `threshold` - for setting a threshold for the total number of
    contacts necessary per column to be considered for deconvolution;
    otherwise it can be quite noisy if there are gaps. For example, note
    that the center cross where there is no signal is because of the
    repetitive LINE element located inside of the Kcnq1ot1 locus, which
    is difficult to map reads to. The default threshold is set to be the
    bottom 0.5% quantile of the `colSums` of the matrix, but this can
    also be manually adjusted with `threshold=`.
2.  `smooth` - `TRUE` or `FALSE`, this uses a Gaussian blur to smooth
    the matrix. This takes care of instances like above where such
    mapping discrepancies can result in ‘holes’ in the matrix.
3.  `sigma` - the SD of the Gaussian blur which controls the amount of
    smoothing to apply. Default is `sigma=1`.

Finally, after thresholding and smoothing, the matrix is
column-normalized.

We illustrate the importance of these steps below:

With no thresholding, we see that there are indeed some spurious reads
in those ‘empty’ gaps, and if we just column-normalize, we end up
getting a messy signal:

``` r
ddnorm <- normalize_hic(dd, threshold=0)
plot_mat(ddnorm, main = 'Column-normalized, with no thresholding or smoothing')
```

![](demo_1d_files/figure-gfm/unnamed-chunk-4-1.png)

The default amount of thresholding clears out the signal in the lowest
0.5% of bins and sets all the signal to be on the diagonal. Smoothing
helps essentially impute nearby values that are lost. Note that we
choose not to make any assumptions about how the ‘empty’ region
interacts with distant regions; the blur effectively smooths over
contacts with nearby regions.

``` r
ddnorm <- normalize_hic(dd, smooth = T)
plot_mat(ddnorm, main = 'Column-normalized, with thresholding and smoothing')
```

![](demo_1d_files/figure-gfm/unnamed-chunk-5-1.png)

## Setting gamma, the amount of deconvolution

In practice, even after thresholding and smoothing, the DNA-DNA contact
matrix can still be quite noisy. The `gamma` parameter is important to
set because it controls a penalty on the off-diagonal of the DNA-DNA
contact matrix, essentially controlling how much of the original signal
we want to preserve vs. deconvolve. We have $$D=(1-\gamma)I+\gamma O,$$
where $O$ is the smoothed, column-normalized matrix we created above. We
will explore the effects of different values of `gamma` below.

# Simulated data

## Gamma and degrees of freedom

We now generate a random true and observed signal assuming different
`gamma` values, and assuming we use different `gamma` values to
deconvolve.

Interpretation of the true `gamma`:

1.  When the true `gamma` is low, then effectively there is no
    convolution happening. The signal has less noise, and the observed
    data more accurately reflects the true signal.
2.  When the true `gamma` is high, then the signal is heavily convolved.
    The observed data appears much more spread out than the true signal.

Interpretation of the `gamma` used to fit:

1.  In general, higher `gamma` results in effectively more of the signal
    being placed at the center peak.
2.  In general, lower `gamma` deconvolves less and the fitted signal is
    close to the observed data.

And in general, fitting with a lower degrees of freedom, set using `df`,
results in fewer spurious peaks while still capturing the main true
signal.

### Using true df

``` r
par(mfrow=c(2,2))
for (g in c(0.25, 0.5, 0.75, 0.9)) {
  ddnorm <- normalize_hic(dd, gamma = g, smooth = T)
  sim <- dcon:::simulate_y(len = 100, df = 10, D = ddnorm, npeaks=1)
  fit0.25 <- fit_decon(sim$y, normalize_hic(dd, gamma = 0.25, smooth = T), df = 10)
  fit0.5 <- fit_decon(sim$y, normalize_hic(dd, gamma = 0.5, smooth = T), df = 10)
  fit0.75 <- fit_decon(sim$y, normalize_hic(dd, gamma = 0.75, smooth = T), df = 10)
  fit0.9 <- fit_decon(sim$y, normalize_hic(dd, gamma = 0.9, smooth = T), df = 10,
                      hessian = T)
  plot(1:100, sim$y, ylim=c(0,200), main = paste0('gamma=',g))
  lines(1:100, exp(sim$B%*%sim$a), col='red')
  lines(1:100, exp(fit0.25$est), col='#8cb0fb')
  lines(1:100, exp(fit0.5$est), col='#748cd0')
  lines(1:100, exp(fit0.75$est), col='#5b6aa6')
  lines(1:100, exp(fit0.9$est), col='#424b7e')
}
```

![](demo_1d_files/figure-gfm/unnamed-chunk-6-1.png)

### Multiple peaks

``` r
par(mfrow=c(2,2))
for (g in c(0.25, 0.5, 0.75, 0.95)) {
  ddnorm <- normalize_hic(dd, gamma = g, smooth = T)
  sim <- dcon:::simulate_y(len = 100, df = 10, D = ddnorm, npeaks=2)
  fit0.25 <- fit_decon(sim$y, normalize_hic(dd, gamma = 0.25, smooth = T), df = 10)
  fit0.5 <- fit_decon(sim$y, normalize_hic(dd, gamma = 0.5, smooth = T), df = 10)
  fit0.75 <- fit_decon(sim$y, normalize_hic(dd, gamma = 0.75, smooth = T), df = 10)
  fit0.95 <- fit_decon(sim$y, normalize_hic(dd, gamma = 0.95, smooth = T), df = 10,
                      hessian = T)
  plot(1:100, sim$y, ylim=c(0,200), main = paste0('gamma=',g))
  lines(1:100, exp(sim$B%*%sim$a), col='red')
  lines(1:100, exp(fit0.25$est), col='#8cb0fb')
  lines(1:100, exp(fit0.5$est), col='#748cd0')
  lines(1:100, exp(fit0.75$est), col='#5b6aa6')
  lines(1:100, exp(fit0.95$est), col='#424b7e')
}
```

![](demo_1d_files/figure-gfm/unnamed-chunk-7-1.png)

## Plots for show

``` r
data.frame(
  pos = 1:100,
  observed = sim$y,
  true = exp(sim$B%*%sim$a),
  fit = exp(fit0.95$est)
)
```

        pos observed       true          fit
    1     1        9   1.057117 7.052660e-01
    2     2       13   1.078882 1.088398e+00
    3     3        7   1.108073 1.586907e+00
    4     4       15   1.147007 2.188727e+00
    5     5       21   1.198726 2.865803e+00
    6     6       18   1.267270 3.581626e+00
    7     7       31   1.358074 4.302264e+00
    8     8       28   1.478552 5.006553e+00
    9     9       23   1.638979 5.692223e+00
    10   10       18   1.853797 6.377354e+00
    11   11       27   2.143573 7.098750e+00
    12   12       23   2.537939 7.909735e+00
    13   13       21   3.079937 8.879391e+00
    14   14       28   3.832410 1.009436e+01
    15   15       35   4.887124 1.166339e+01
    16   16       37   6.377314 1.372407e+01
    17   17       40   8.493769 1.645023e+01
    18   18       37  11.503101 2.005732e+01
    19   19       38  15.763540 2.480072e+01
    20   20       43  21.727973 3.095922e+01
    21   21       44  29.916184 3.879322e+01
    22   22       42  40.831618 4.846788e+01
    23   23       28  54.800263 5.993895e+01
    24   24       44  71.732045 7.281851e+01
    25   25       52  90.856302 8.626614e+01
    26   26       44 110.550105 9.897671e+01
    27   27       48 128.416544 1.093327e+02
    28   28       48 141.717450 1.157355e+02
    29   29       44 148.098790 1.170355e+02
    30   30       54 146.343038 1.128954e+02
    31   31       53 136.794693 1.039185e+02
    32   32       43 121.237379 9.147233e+01
    33   33       57 102.287848 7.728838e+01
    34   34       44  82.612568 6.301630e+01
    35   35       53  64.306760 4.990156e+01
    36   36       50  48.618478 3.866225e+01
    37   37       43  35.996385 2.953842e+01
    38   38       42  26.319715 2.243359e+01
    39   39       37  19.162175 1.706963e+01
    40   40       38  13.999804 1.310841e+01
    41   41       31  10.336614 1.022668e+01
    42   42       43   7.760715 8.151270e+00
    43   43       35   5.956122 6.668132e+00
    44   44       43   4.692571 5.617788e+00
    45   45       39   3.807893 4.885620e+00
    46   46       28   3.190521 4.391622e+00
    47   47       18   2.765073 4.081639e+00
    48   48       25   2.481617 3.920671e+00
    49   49       38   2.308180 3.888147e+00
    50   50       35   2.225851 3.974856e+00
    51   51       44   2.225851 4.181220e+00
    52   52       40   2.308180 4.516679e+00
    53   53       45   2.481617 5.000025e+00
    54   54       42   2.765073 5.660640e+00
    55   55       37   3.190521 6.540640e+00
    56   56       30   3.807893 7.697934e+00
    57   57       45   4.692571 9.210212e+00
    58   58       36   5.956122 1.117970e+01
    59   59       50   7.760715 1.373822e+01
    60   60       40  10.336614 1.705167e+01
    61   61       33  13.999804 2.132199e+01
    62   62       38  19.162175 2.678380e+01
    63   63       45  26.319715 3.369147e+01
    64   64       48  35.996385 4.229132e+01
    65   65       42  48.618478 5.277403e+01
    66   66       41  64.306760 6.520472e+01
    67   67       37  82.612568 7.943431e+01
    68   68       41 102.287848 9.500596e+01
    69   69       56 121.237379 1.110824e+02
    70   70       55 136.794693 1.264293e+02
    71   71       49 146.343038 1.394915e+02
    72   72       55 148.098790 1.485814e+02
    73   73       52 141.717450 1.521715e+02
    74   74       54 128.416544 1.492378e+02
    75   75       42 110.550105 1.395665e+02
    76   76       47  90.856302 1.239214e+02
    77   77       49  71.732045 1.039884e+02
    78   78       33  54.800263 8.207580e+01
    79   79       42  40.831618 6.063103e+01
    80   80       37  29.916184 4.171516e+01
    81   81       40  21.727973 2.660771e+01
    82   82       34  15.763540 1.567161e+01
    83   83       23  11.503101 8.498996e+00
    84   84       32   8.493769 4.238600e+00
    85   85       34   6.377314 1.945632e+00
    86   86       30   4.887124 8.249328e-01
    87   87       21   3.832410 3.252208e-01
    88   88       17   3.079937 1.204230e-01
    89   89       11   2.537939 4.245557e-02
    90   90       18   2.143573 1.449799e-02
    91   91       14   1.853797 4.894253e-03
    92   92       16   1.638979 1.671487e-03
    93   93       14   1.478552 5.921520e-04
    94   94       26   1.358074 2.233455e-04
    95   95       23   1.267270 9.204027e-05
    96   96       31   1.198726 4.247171e-05
    97   97       32   1.147007 2.243399e-05
    98   98       22   1.108073 1.381559e-05
    99   99       26   1.078882 1.005642e-05
    100 100       22   1.057117 8.724422e-06

``` r
sessionInfo()
```

    R version 4.2.2 (2022-10-31)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Rocky Linux 8.7 (Green Obsidian)

    Matrix products: default
    BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.15.so

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
     [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] dcon_0.0.0.9000   ggplot2_3.4.0     tidyr_1.2.1       dplyr_1.0.10     
    [5] data.table_1.14.6

    loaded via a namespace (and not attached):
     [1] Rcpp_1.0.9             lattice_0.20-45        deldir_1.0-6          
     [4] assertthat_0.2.1       digest_0.6.31          utf8_1.2.2            
     [7] spatstat.core_2.4-4    R6_2.5.1               GenomeInfoDb_1.34.4   
    [10] plyr_1.8.8             stats4_4.2.2           evaluate_0.19         
    [13] tensor_1.5             pillar_1.8.1           zlibbioc_1.44.0       
    [16] rlang_1.0.6            rstudioapi_0.14        S4Vectors_0.36.1      
    [19] rpart_4.1.19           Matrix_1.5-3           goftest_1.2-3         
    [22] rmarkdown_2.19         labeling_0.4.2         splines_4.2.2         
    [25] stringr_1.5.0          RCurl_1.98-1.9         polyclip_1.10-4       
    [28] munsell_0.5.0          spatstat.data_3.0-0    compiler_4.2.2        
    [31] xfun_0.35              pkgconfig_2.0.3        BiocGenerics_0.44.0   
    [34] mgcv_1.8-41            htmltools_0.5.4        tidyselect_1.2.0      
    [37] spatstat.random_3.0-1  tibble_3.1.8           GenomeInfoDbData_1.2.9
    [40] IRanges_2.32.0         fansi_1.0.3            withr_2.5.0           
    [43] bitops_1.0-7           grid_4.2.2             nlme_3.1-161          
    [46] jsonlite_1.8.4         gtable_0.3.1           lifecycle_1.0.3       
    [49] DBI_1.1.3              magrittr_2.0.3         scales_1.2.1          
    [52] cli_3.4.1              stringi_1.7.8          farver_2.1.1          
    [55] XVector_0.38.0         reshape2_1.4.4         generics_0.1.3        
    [58] vctrs_0.5.1            spatstat.utils_3.0-1   RColorBrewer_1.1-3    
    [61] tools_4.2.2            glue_1.6.2             purrr_0.3.5           
    [64] abind_1.4-5            fastmap_1.1.0          yaml_2.3.6            
    [67] spatstat.sparse_3.0-0  colorspace_2.0-3       GenomicRanges_1.50.1  
    [70] spatstat.geom_3.0-3    knitr_1.41            
