% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataspace.R
\name{pa_initialize_dataspace}
\alias{pa_initialize_dataspace}
\title{Register the Data Space credentials to the R
  environment}
\usage{
pa_initialize_dataspace(username, password, verbose = TRUE)
}
\arguments{
\item{username}{username used to authenticate the HTTP
request}

\item{password}{password used to authenticate the HTTP
request}

\item{verbose}{whether to print information about this operation}
}
\value{
No return value, called for side effects
}
\description{
Register the Data Space credentials to the R
  environment
}
\details{
`pa_initialize_dataspace()` registers the username
  and password to the machine's R environment. All the
  other functions that rely on authentication will search
  for the username and password in the R environment. Do
  not share your R environment with others, as they will
  be able to read your username and password. You can
  register at \url{https://dataspace.copernicus.eu/}.
}
\examples{
\dontrun{
pa_initialize_dataspace('my-username', 'my-password')
}

}
\author{
Caio dos Santos and Fernando Miguez
}
