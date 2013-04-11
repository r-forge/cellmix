\name{MySQLtoSQLite}
\alias{MySQLtoSQLite}
\title{Load a MySQL Dump for Using with SQLite}
\usage{
  MySQLtoSQLite(file, conn = NULL)
}
\arguments{
  \item{file}{MySQL dump file}

  \item{conn}{SQLite connection object}
}
\value{
  the statements as a character vector if \code{conn=NULL}
  or nothing
}
\description{
  Converts a MySQL dump file into SQLite compatible SQL
  statements.
}
\note{
  This function was developed and used to load the VeryGene
  database, in order to create the \code{\link{VeryGene}}
  marker list.
}

