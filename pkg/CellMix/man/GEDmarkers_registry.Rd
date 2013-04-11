\docType{data}
\name{GEDmarkers_registry}
\alias{GEDmarkers_registry}
\title{Internal Registry for Marker Lists}
\format{List of 8
 $ get_field             : NULL
 $ get_fields            : NULL
 $ get_field_names       : NULL
 $ set_field             : NULL
 $ has_entry             : NULL
 $ get_entry             : NULL
 $ get_entries           : NULL
 $ get_entry_names       : NULL
 - attr(*, "class")= chr [1:4] "GEDmarkers_registry" "markers_package_registry" "package_registry" "registry"
 - attr(*, "package")= chr "CellMix"}
\usage{
  GEDmarkers_registry
}
\description{
  Registry object of class
  \code{\link[registry:regobj]{registry}} for storing data
  on marker lists. The data are stored as plain lists but
  are wrapped into \code{\linkS4class{MarkerList}} objects
  when retrieved with the factory function
  \code{\link{MarkerList}}.
}
\keyword{datasets}
\keyword{internal}

