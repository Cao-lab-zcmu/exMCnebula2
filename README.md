# exMCnebula2

Add-on package to the 'MCnebula2' package for illustrative analysis and mapping.

## Installation

Please install package of MCnebula2 first: <https://github.com/Cao-lab-zcmu/MCnebula2>

In general, the following code will install the dependent or imported packages
as well as the exMCnebula2 package:

```r
remotes::install_github("Cao-lab-zcmu/exMCnebula2")
```

If you do not want to download the additional files ('inst/extdata/*') and just
want to use the functions in the package, install them using the following
code:

```r
remotes::install_github("Cao-lab-zcmu/exMCnebula2", "light")
```

However, packages that do not exist in CRAN but exist in Bioconductor may not be installed
automatically, in which case they may need to be installed first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("BiocParallel", "FELLA"))
```

## Analysis in the article of MCnebula

Load the packages:

```r
library(MCnebula2)
library(exMCnebula2)
```

The scripts or pdf documents or figures were as following:

```r
dir <- system.file("inst", "extdata", "scripts_evaluation", package = "exMCnebula2")
list.files(dir, recursive = T)
```

```r
#  [1] "eucommia_workflow/eucommia_workflow.R"     "eucommia_workflow/report.pdf"             
#  [3] "evaluation_workflow/evaluation_workflow.R" "evaluation_workflow/report.pdf"           
#  [5] "mcn_principle/a_elements.R"                "mcn_principle/b_gather.R"                 
#  [7] "mcn_principle/figure_mech.pdf"             "mcn_structure/a_project.R"                
#  [9] "mcn_structure/b_mcn_dataset.R"             "mcn_structure/c_nebulae.R"                
# [11] "mcn_structure/d_across.R"                  "mcn_structure/data_stream.pdf"            
# [13] "serum_workflow/report.pdf"                 "serum_workflow/serum_workflow.R"          
```

For evaluation of MCnebula:

```r
envir <- new.env()
## It may be time-consumed. The output files would be saved in `tempdir()`
source(paste0(dir, "/evaluation_workflow/evaluation_workflow.R"), local = envir)
```

For analysis of herbal dataset (_E. ulmoides_):

```r
envir <- new.env()
source(paste0(dir, "/eucommia_workflow/eucommia_workflow.R"), local = envir)
```

For analysis of serum dataset:

```r
envir <- new.env()
source(paste0(dir, "/serum_workflow/serum_workflow.R"), local = envir)
```

For drawing 'figure_mech.pdf':

```r
envir <- new.env()
lapply(c("a_elements.R", "b_gather.R"),
  function(script) {
    source(paste0(dir, "/mcn_principle/", script), local = envir)
  })
```

For drawing 'data_stream.pdf':

```r
envir <- new.env()
lapply(c("a_project.R", "b_mcn_dataset.R", "c_nebulae.R", "d_across.R"),
  function(script) {
    source(paste0(dir, "/mcn_structure/", script), local = envir)
  })
```

The codes of these scripts were run successfully in:

```r
Sys.info()
```

```r
#                                                            sysname 
#                                                            "Linux" 
#                                                            release 
#                                         "5.17.15-76051715-generic" 
#                                                            version 
# "#202206141358~1655919116~22.04~1db9e34 SMP PREEMPT Wed Jun 22 19" 
#                                                           nodename 
#                                                           "pop-os" 
#                                                            machine 
#                                                           "x86_64" 
#                                                              login 
#                                                             "echo" 
#                                                               user 
#                                                             "echo" 
#                                                     effective_user 
#                                                             "echo" 
```


