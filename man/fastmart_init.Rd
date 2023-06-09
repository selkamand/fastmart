% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cache.R
\name{fastmart_init}
\alias{fastmart_init}
\title{Initialise fastmart}
\usage{
fastmart_init(cache_dir = fastmart_default_cache())
}
\arguments{
\item{cache_dir}{directory to store cache in}
}
\value{
path to the newly created fastmart cache
}
\description{
Creates a folder (by default at ~/.fastmart) where any database caches will be stored
}
