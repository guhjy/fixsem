% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/engineV1.R
\name{fixsem}
\alias{fixsem}
\title{fixsem}
\usage{
fixsem(model, data, group = NULL, growth = F,
  remote = getOption("kaefaServers"), fastrun = F, method = if (fastrun) {
      ((c("ML", "MLM", "MLMV", "MLMVS", "MLF", "MLR"))) } else {     ((c("WLS",
  "DWLS", "WLSM", "WLSMV", "WLSMVS", "ML", "MLM", "MLMV", "MLMVS", "MLF", "MLR",
  "ULS", "ULSM", "ULSMV", "ULSMVS"))) })
}
\arguments{
\item{model}{a lavaan-style syntax}

\item{data}{data for calibration}

\item{group}{specify the group if you want}

\item{growth}{logical; your model is growth model? default is FALSE}

\item{remote}{see future::plan()}

\item{fastrun}{logical; default is FALSE. If you turn on to TRUE, it will be run with ML family.}

\item{method}{calibration method what can lavaan run}
}
\value{
calibrated lavaan models
}
\description{
fixsem
}
\examples{
\dontrun{
geocode("3817 Spruce St, Philadelphia, PA 19104")
geocode("Philadelphia, PA")
dat <- data.frame(value=runif(3),address=c("3817 Spruce St, Philadelphia, PA 19104","Philadelphia, PA","Neverneverland"))
geocode(dat)
}
}
