\name{read_gasex_file}

\alias{read_gasex_file}

\title{Reading a gas exchange log file}

\description{
  Tool for reading log files created by gas exchange measurement instruments and
  storing their contents in \code{\link{exdf}} objects.
}

\usage{
  read_gasex_file(
    file_name,
    timestamp_colname = NA,
    posix_options = list(),
    file_type = 'AUTO',
    instrument_type = 'AUTO',
    standardize_columns = TRUE,
    remove_NA_rows = TRUE,
    \dots
  )
}

\arguments{
  \item{file_name}{
    A relative or absolute path to a log file containing gas exchange data.
  }

  \item{timestamp_colname}{
    The name of the column that contains the timestamp of each measurement;
    typically, this is something like \code{'time'} or \code{'TIMESTAMP'}.
  }

  \item{posix_options}{
    Optional arguments to pass to \code{\link{as.POSIXlt}}; must be formatted as
    a list of named elements. See details below for more information.
  }

  \item{file_type}{
    The type of file to be loaded. If \code{file_type} is \code{'AUTO'}, then
    the file type will be automatically determined from the extension of
    \code{file_name}. The other supported options are \code{'plaintext'},
    \code{'Excel'}, and \code{'data'}.
  }

  \item{instrument_type}{
    The type of measurement instrument that produced the log file. If
    \code{instrument_type} is \code{'AUTO'}, then the instrument type will be
    determined from the \code{file_type}. The other supported options are
    \code{'Licor LI-6800'} and \code{'CR3000'}.
  }

  \item{standardize_columns}{
    A logical value indicating whether to standardize columns; see details
    below.
  }

  \item{remove_NA_rows}{
    A logical value indicating whether to remove any rows whose values are all
    \code{NA}; this argument will be passed to the specialized reading
    functions; see below for more details.
  }

  \item{\dots}{
    Additional arguments to be passed to specialized reading functions; see
    below for more details.
  }
}

\details{
  Some log files contain Unicode characters in some column names and units, but
  these characters cannot be represented properly in R. To address this, Unicode
  characters are replaced with reasonable alternatives; for example, the
  character for the capital Greek letter delta is replaced with the word
  \code{Delta}. The replacement rules are stored in a data frame that can be
  accessed via \code{PhotoGEA:::UNICODE_REPLACEMENTS}, and more information can
  be found in the source code (\code{R/unicode_replacements.R}).

  Sometimes it is useful to "standardize" the names, units, or categories of
  columns in instrument log files. This can be helpful in several situations:
  \itemize{
    \item An instrument may not be consistent with the name of a column; for
          example, Licor LI-6800s may may have a \code{PhiPs2} or \code{PhiPS2}
          column depending on the version of the operating system running on the
          machine.
    \item An instrument may not specify the units of a column; for example,
          Licor LI-6800s do not specify that \code{PhiPS2} has units of
          \code{dimensionless}.
    \item An instrument may use different names or different units than another
          instrument for the same measured quantity.
  }
  To deal with these situations, it is possible to "standardize" the column
  names, units, and categories when reading an instrument file. A list of
  definitions for all standardizations can be accessed from an R session by
  typing \code{View(PhotoGEA:::gasex_column_conversions)}.

  When reading a log file, it can be useful to identify the timestamp column so
  its values can be properly interpreted as \code{\link{POSIXlt}} objects. If
  \code{timestamp_colname} is \code{NA}, this conversion will be skipped. By
  default, \code{read_gasex_file} calls \code{\link{as.POSIXlt}} with
  \code{origin = '1970-01-01'} and \code{tz = ''}. With these options, any
  numeric timestamps (such as \code{1692386305.5}) will be interpreted as the
  number of seconds since January 1, 1970 (the UNIX standard) and the time will
  be expressed using the local system time. This works well in many situations.
  However, if a log file was created in a different time zone than the local
  one, it may be necessary to specify the time zone. This can be done via the
  \code{posix_options} argument. For example, to interpret the timestamp as a
  time in US Central time, set \code{posix_options = list(tz = 'US/Central')}.
  This may be necessary when using \code{\link{pair_gasex_and_tdl}} to match
  timestamps between different log files.

  When automatically determining the file type from its extension, the following
  rules are used:
  \itemize{
    \item A \code{.xlsx} extension corresponds to \code{file_type = 'Excel'}.
    \item A \code{.dat} extension corresponds to \code{file_type = 'data'}.
    \item A \code{.txt} extension or a file with no extension corresponds to
          \code{file_type = 'plaintext'}.
  }

  When automatically determining the instrument type from the file type, the
  following rules are used:
  \itemize{
    \item File types of \code{'Excel'} and \code{'plaintext'} correspond to
          \code{instrument_type = 'Licor LI-6800'}.
    \item A file type of \code{'data'} corresponds to \code{instrument_type =
          'CR3000'}.
  }

  Internally, this function calls one of several other (non-exported) functions
  depending on the values of \code{instrument_type} and \code{file_type}:
  \itemize{
    \item \code{\link{read_licor_6800_plaintext}} (for \code{instrument_type =
          'LI-6800'} and \code{file_type = 'plaintext'})
    \item \code{\link{read_licor_6800_Excel}} (for \code{instrument_type =
          'LI-6800'} and \code{file_type = 'Excel'})
    \item \code{\link{read_cr3000}} (for \code{instrument_type =
          'CR3000'} and \code{file_type = 'data'})
  }
  Any additional arguments specified via \code{\dots} will be passed to these
  functions, along with the value of \code{remove_NA_rows}.

  \strong{IMPORTANT NOTE ABOUT LICOR EXCEL FILES}: by default, Licor Excel files
  do not "calculate" formula values. This causes a problem when reading them in
  R, since any data entry determined from a formula will be read as 0. To fix
  this issue for a Licor Excel file, open it in in Excel, go to the
  \code{Formulas} menu, and choose \code{Calculate Now}. (Alternatively, press
  F9.) Then save the file and close it. See \code{\link{read_licor_6800_Excel}}
  for more details.
}

\value{
  An \code{exdf} object that fully includes all the data from the log file. In
  addition to the required elements of an \code{exdf} object, the following
  "extra" elements are also included:
  \itemize{
    \item \code{file_name}: A copy of the input argument with the same name.
    \item \code{instrument_type}: A copy of the input argument with the same
          name.
    \item \code{file_type}: A copy of the input argument with the same name,
          unless it was set to \code{'AUTO'}; in that case, the file type that
          was determined from the file's extension.
    \item \code{timestamp_colname}: A copy of the input argument with the same
          name, unless it was set to \code{'AUTO'}; in that case, the instrument
          type that was determined from the file type.
  }
}

\examples{
# Example: Eeading a Licor Excel file that is included with the PhotoGEA
# package. Here we specify 'time' as the name of the timestamp column.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx'),
  'time'
)

licor_file$file_name     # A record of where the data came from
str(licor_file)          # View the contents of the exdf object's main_data
str(licor_file$preamble) # View the Licor file's preamble data
}

\concept{exdf}
