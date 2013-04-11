\docType{data}
\name{GEDalgorithm_registry}
\alias{GEDalgorithm_registry}
\title{Registry Object for GED Algorithms}
\format{List of 15
 $ get_field             : NULL
 $ get_fields            : NULL
 $ get_field_names       : NULL
 $ set_field             : NULL
 $ has_entry             : NULL
 $ get_entry             : NULL
 $ get_entries           : NULL
 $ get_entry_names       : NULL
 $ grep_entries          : NULL
 $ set_entry             : NULL
 $ modify_entry          : NULL
 $ delete_entry          : NULL
 $ n_of_entries          : NULL
 $ get_field_entries     : NULL
 $ get_permissions       : NULL
 - attr(*, "class")= chr [1:4] "GEDalgorithm_registry" "algorithms_package_registry" "package_registry" "registry"
 - attr(*, "package")= chr "CellMix"}
\usage{
  GEDalgorithm_registry
}
\description{
  Registry object of class
  \code{\link[registry:regobj]{registry}} for storing data
  on gene expression deconvolution algorithms. The data are
  stored as plain lists but are wrapped into
  \code{\linkS4class{GEDStrategy}} objects when retrieved
  with the factory function \code{\link{GEDStrategy}}.
}
\keyword{datasets}
\keyword{internal}

