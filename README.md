
-   <a
    href="#dcon-probabilistic-deconvolution-of-chromatin-interaction-signal-in-3d-genome-data"
    id="toc-dcon-probabilistic-deconvolution-of-chromatin-interaction-signal-in-3d-genome-data">dcon:
    Probabilistic deconvolution of chromatin interaction signal in 3D genome
    data</a>
    -   <a href="#quick-start" id="toc-quick-start">Quick Start</a>
        -   <a href="#installation" id="toc-installation">Installation</a>
        -   <a href="#normalize-dna-dna-contact-matrix"
            id="toc-normalize-dna-dna-contact-matrix">Normalize DNA-DNA contact
            matrix</a>
        -   <a href="#fit-deconvolution-in-1d" id="toc-fit-deconvolution-in-1d">Fit
            deconvolution in 1D</a>
    -   <a href="#vignettes" id="toc-vignettes">Vignettes</a>
    -   <a href="#citation" id="toc-citation">Citation</a>
    -   <a href="#figures-from-the-manuscript"
        id="toc-figures-from-the-manuscript">Figures from the manuscript</a>

# dcon: Probabilistic deconvolution of chromatin interaction signal in 3D genome data

## Quick Start

### Installation

``` r
# install.packages("devtools")
devtools::install_github("lulizou/dcon", build_vignettes = FALSE)
```

### Normalize DNA-DNA contact matrix

``` r
dd <- normalize_hic(raw_matrix)
```

### Fit deconvolution in 1D

``` r
y_hat <- fit_decon(y, dd)
```

For details, see the vignettes below.

## Vignettes

-   [Deconvolution in
    1D](https://github.com/lulizou/dcon/blob/main/vignettes/demo_1d.md)
-   [Loop
    deconvolution](https://github.com/lulizou/dcon/blob/main/vignettes/demo_2d.md)
-   [Simulating
    data](https://github.com/lulizou/dcon/blob/main/vignettes/simulations.md)

## Citation

## Figures from the manuscript

-   [RD-SPRITE
    tracks](https://github.com/lulizou/dcon/blob/main/vignettes/rdsprite_track_plots.md)
