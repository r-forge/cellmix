\name{lverbose}
\alias{is.verbose}
\alias{lmessage}
\alias{log_append}
\alias{lsilent}
\alias{lverbose}
\alias{vmessage}
\title{Logging Feature}
\usage{
  lverbose(val)

  lsilent()

  is.verbose()

  lmessage(level, ..., appendLF = TRUE, sep = "",
    force = FALSE)

  vmessage(...)

  log_append(...)
}
\arguments{
  \item{val}{logical/numeric value that sets the verbosity
  level.}

  \item{...}{parts of a character message that are
  concatenated and passed to the current logger's write
  function.}

  \item{level}{verbosity level threshold (numeric value)
  above which the message should be printed out. This
  threshold is compared with the current verbosity level as
  returned by \code{lverbose}.}

  \item{appendLF}{logical indicating if an endline
  character should be appended at the end of the message.}

  \item{sep}{separation character, used when concatenating
  all arguments in \code{...}.}

  \item{force}{logical that indicates if one should output
  messages or return a non null logger, even if the verbose
  mode is not high enough.}
}
\value{
  the old verbose level
}
\description{
  \code{lverbose} returns/sets the current verbosity level.

  \code{lsilent} tells if all verbose messages are
  silenced.

  \code{is.verbose} tells if verbosity is on, i.e. at level
  greater than 0.

  \code{lmessage} prints out a message (on sdtout) if the
  verbosity level is greater than a given value.

  \code{vmessage} prints out a log message (at level 1)
  using the current logger, typically on stdout. It is a
  shortcut for \code{lmessage(1L, ...)}.

  \code{log_append} directly appends some message to the
  current log line.
}
\keyword{internal}

