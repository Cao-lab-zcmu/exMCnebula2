% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_EIC_stack.R
\name{plot_EIC_stack}
\alias{plot_EIC_stack}
\alias{bioc.par}
\title{Draw extracted ions chromatography for 'features'}
\usage{
plot_EIC_stack(
  idset,
  metadata,
  quant.path,
  mzml.path,
  palette = (ggsci::pal_npg())(10),
  mz.tol = 0.01,
  rt.tol = 0.1,
  cl = NULL,
  data = NULL
)

bioc.par(cl = 4)
}
\description{
Use quantification table (with peak start time and end time)
exported by MCmine to draw EIC plot.

\code{plot_EIC_stack}: ...

\code{bioc.par}: ...
}
