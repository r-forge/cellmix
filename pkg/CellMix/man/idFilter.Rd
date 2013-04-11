\name{idFilter}
\alias{idFilter}
\alias{idFilterAffy}
\alias{idFilterAll}
\alias{idFilterAuto}
\alias{idFilterFirstN}
\alias{idFilterInjective}
\alias{idFilterMAuto}
\alias{idFilterOneToMany}
\alias{idFilterOneToOne}
\title{Gene Identifier Filtering Strategies}
\usage{
  idFilter(name, ..., .wrap = FALSE)

  idFilterAll(map, ...)

  idFilterFirstN(map, n = 1L, exact.first = TRUE, ...)

  idFilterOneToOne(map, ...)

  idFilterInjective(map, strict = TRUE, ..., .last)

  idFilterOneToMany(map, ..., .last)

  idFilterAffy(map, secondary = NULL, ...)

  idFilterAuto(map, ..., .last)

  idFilterMAuto(map, ..., .last)
}
\arguments{
  \item{name}{name of the strategy}

  \item{...}{in \code{idFilter}, these are extra arguments
  that are used to pre-build a call to the strategy, which
  is wrapped into a function. These arguments are only used
  when \code{.wrap=TRUE}.

  For the different filtering stratgies, e.g.,
  idFilterAuto, these arguments are either passed to
  internal calls to other strategies, or not used, but
  required so that filters can be composed.}

  \item{.wrap}{logical that indicates if the function call
  should be build and wrapped into a caller function, which
  takes only two arguments, \code{map} and \code{.last},
  that will receive the current mapping and a logical that
  is \code{TRUE} only when \code{map} is the last mapping
  in the whole sequence of mappings.

  The value of \code{.last} will only be passed to the
  original filter function if this one has an argument of
  the same name.}

  \item{map}{mapping \code{list} object that maps source
  identifiers to another type of identifier.}

  \item{n}{maximum number of matches to keep. The default
  is 1, which only keeps the first match}

  \item{exact.first}{logical that indicates if exact
  matches should be prioritized, and appear as first
  match.}

  \item{strict}{logical that indicates if all intermediate
  mappings are required to be \emph{injective} one to one,
  or only the last mapping step.}

  \item{.last}{logical that indicates if \code{map}
  constitutes the last mapping or an intermediate mapping.}

  \item{secondary}{specifies what should be done with
  secondary probes. If \code{secondary=TRUE} then secondary
  probes are kept in the mapping but primary probes are
  prioritized, in the sense that they will appear in front
  of secondary probes in each element of the mapping. If
  \code{secondary=FALSE} then secondary probes are removed
  from the mapping, meaning that identifiers that only
  mapped to secondary probes are then unmapped. If
  \code{secondary=NULL} (default), then the filter is
  applied as when \code{secondary=TRUE}, \strong{except}
  when source identifiers are also Affymetrix probes and do
  not contain any secondary probes, in which case the the
  filters is applied as when \code{secondary=FALSE}.

  Any other value for \code{secondary} makes the filter to
  behave as a pass through, i.e. the map is returned
  unchanged.}
}
\value{
  a function
}
\description{
  \code{idFilter} provides access to a small registry that
  contains a set of filtering strategy functions, which can
  be used in calls to the \code{\link{mapIDs}} function.
}
\details{
  Filter strategies are functions that take a mapping as
  its first argument and returns it after performing some
  filtering. They should always also have argument
  \code{...}, so that they can be combined in sequences.
  Filter strategies may optionally have an argument
  \code{.last}, whose value will be filled by
  \code{mapIDs}: \code{TRUE} if and only if the mapping
  constitutes the last step in the identifier conversion.

  The following filter strategies are currently defined:

  \code{idFilterAll} [key \sQuote{all}]: does not filter at
  all, and simply return the input map unchanged.

  \code{idFilterFirstN} [key \sQuote{firstN}]: keeps the
  first \code{n} matches.

  \code{idFilterOneToOne} [key \sQuote{1:1}]: only keeps
  single matches.

  \code{idFilterBiOnetoOne} [key \sQuote{1:1}]: only keeps
  bidirectional one to one mapping, i.e. that the reverse
  mapping is also one to one. In other words it return the
  largest injective mapping.

  \code{idFilterOneToMany} [key \sQuote{1:*}]: only keeps
  single matches for intermediate mappings, but keeps all
  matches (i.e. single and multiple) when mapping to the
  last identifier type.

  \code{idFilterAffy} [key \sQuote{affy}]: enables to
  filter out secondary affy probes, i.e., that do not match
  the primary pattern "^[0-9]+_at$". This filter only
  applies when the destination identifiers are suitable
  Affymetrix probes.

  \code{idFilterAuto} [key \sQuote{auto}]: this filter is a
  composed filter, that successively applies filters
  \sQuote{affy} and \sQuote{firstN}. It is suitable when
  one wants to map 1 to 1 as many identifiers as possible.

  \code{idFilterMAuto} [key \sQuote{mauto}]: this filter is
  a composed filter, that successively applies filters
  \sQuote{affy} and \sQuote{1:*}.
}
\examples{
# list of all available filters
idFilter()

# show a given filter function
idFilter('1:1')

# composed filter
f <- idFilter(c('1:1', 'affy'))
}

