\name{error_function_c3_aci}

\alias{error_function_c3_aci}

\title{Generate an error function for C3 A-Ci curve fitting}

\description{
  Creates a function that returns an error value (the negative of the natural
  logarithm of the likelihood) representing the amount of agreement between
  modeled and measured \code{An} values. When this function is minimized, the
  likelihood is maximized.

  Internally, this function uses \code{\link{apply_gm}} to calculate \code{Cc},
  and then uses \code{link{calculate_c3_assimilation}} to calculate assimilation
  rate values that are compared to the measured ones.
}

\usage{
  error_function_c3_aci(
    replicate_exdf,
    fit_options = list(),
    sd_A = 1,
    Wj_coef_C = 4.0,
    Wj_coef_Gamma_star = 8.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    gamma_star_norm_column_name = 'Gamma_star_norm',
    gmc_norm_column_name = 'gmc_norm',
    j_norm_column_name = 'J_norm',
    kc_norm_column_name = 'Kc_norm',
    ko_norm_column_name = 'Ko_norm',
    oxygen_column_name = 'oxygen',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    tp_norm_column_name = 'Tp_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    hard_constraints = 0,
    debug_mode = FALSE,
    \dots
  )
}

\arguments{
  \item{replicate_exdf}{
    An \code{exdf} object representing one CO2 response curve.
  }

  \item{fit_options}{
    A list of named elements representing fit options to use for each parameter.
    Values supplied here override the default values (see details below). Each
    element must be \code{'fit'}, \code{'column'}, or a numeric value. A value
    of \code{'fit'} means that the parameter will be fit; a value of
    \code{'column'} means that the value of the parameter will be taken from a
    column in \code{replicate_exdf} of the same name; and a numeric value means
    that the parameter will be set to that value. For example,
    \code{fit_options = list(alpha_g = 0, Vcmax_at_25 = 'fit', Tp_at_25 = 'column')}
    means that \code{alpha_g} will be set to 0, \code{Vcmax_at_25} will be fit,
    and \code{Tp_at_25} will be set to the values in the \code{Tp_at_25} column
    of \code{replicate_exdf}.
  }

  \item{sd_A}{
    The standard deviation of the measured values of the net CO2 assimilation
    rate, expressed in units of \code{micromol m^(-2) s^(-1)}. If \code{sd_A} is
    not a number, then there must be a column in \code{replicate_exdf} called
    \code{sd_A} with appropriate units. A numeric value supplied here will
    overwrite the values in the \code{sd_A} column of \code{replicate_exdf} if
    it exists.
  }

  \item{Wj_coef_C}{
    A coefficient in the equation for RuBP-regeneration-limited carboxylation,
    whose value depends on assumptions about the NADPH and ATP requirements of
    RuBP regeneration; see \code{\link{calculate_c3_assimilation}} for more
    information.
  }

  \item{Wj_coef_Gamma_star}{
    A coefficient in the equation for RuBP-regeneration-limited carboxylation,
    whose value depends on assumptions about the NADPH and ATP requirements of
    RuBP regeneration; see \code{\link{calculate_c3_assimilation}} for more
    information.
  }

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    intercellular CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{gamma_star_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Gamma_star} values (with units of
    \code{normalized to Gamma_star at 25 degrees C}).
  }

  \item{gmc_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    mesophyll conductance values (with units of
    \code{normalized to gmc at 25 degrees C}).
  }

  \item{j_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{J} values (with units of \code{normalized to J at 25 degrees C}).
  }

  \item{kc_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Kc} values (with units of \code{normalized to Kc at 25 degrees C}).
  }

  \item{ko_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Ko} values (with units of \code{normalized to Ko at 25 degrees C}).
  }

  \item{oxygen_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    concentration of O2 in the ambient air, expressed as a percentage (commonly
    21\% or 2\%); the units must be \code{percent}.
  }

  \item{rl_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{RL} values (with units of \code{normalized to RL at 25 degrees C}).
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{replicate_exdf} that contains the total
    pressure in \code{bar}.
  }

  \item{tp_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Tp} values (with units of \code{normalized to Tp at 25 degrees C}).
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{cj_crossover_min}{
    The minimum value of \code{Cc} (in ppm) where \code{Aj} is allowed to become
    the overall rate-limiting factor. If \code{cj_crossover_min} is set to
    \code{NA}, this restriction will not be applied.
  }

  \item{cj_crossover_max}{
    The maximim value of \code{Cc} (in ppm) where \code{Wj} is allowed to be
    smaller than \code{Wc}. If \code{cj_crossover_max} is set to \code{NA}, this
    restriction will not be applied.
  }

  \item{hard_constraints}{
    To be passed to \code{\link{calculate_c3_assimilation}}; see that function
    for more details.
  }

  \item{debug_mode}{
    A logical (\code{TRUE} or \code{FALSE}) variable indicating whether to
    operate in debug mode. In debug mode, information about the \code{guess} is
    printed each time the error function is called; this can be helpful when
    troubleshooting issues with an optimizer.
  }

  \item{\dots}{
    Additional arguments to be passed to \code{\link{calculate_c3_assimilation}}.
  }
}

\details{
  When fitting A-Ci curves using a maximum likelihood approach, it is necessary
  to define a function that calculates the likelihood of a given set of
  \code{alpha_g}, \code{alpha_old}, \code{alpha_s}, \code{alpha_t},
  \code{Gamma_star_at_25}, \code{gmc_at_25}, \code{J_at_25}, \code{Kc_at_25},
  \code{Ko_at_25},  \code{RL_at_25}, \code{Tp_at_25}, and \code{Vcmax_at_25}
  values by comparing a model prediction to a measured curve. This function will
  be passed to an optimization algorithm which will determine the values that
  produce the largest likelihood.

  The \code{error_function_c3_aci} returns such a function, which is based on a
  particular A-Ci curve and a set of fitting options. It is possible to just fit
  a subset of the available fitting parameters; by default, the fitting
  parameters are \code{alpha_old}, \code{J_at_25}, \code{RL_at_25},
  \code{Tp_at_25}, and \code{Vcmax_at_25}. This behavior can be changed via the
  \code{fit_options} argument.

  For practical reasons, the function actually returns values of \code{-ln(L)},
  where \code{L} is the likelihood. The logarithm of \code{L} is simpler to
  calculate than \code{L} itself, and the minus sign converts the problem from
  a maximization to a minimization, which is important because most optimizers
  are designed to minimize a value.

  Sometimes an optimizer will choose biologically unreasonable parameter values
  that nevertheless produce good fits to the supplied assimilation values. A
  common problem is that the fit result may not indicate Ac-limited assimilation
  at low CO2 values, which should be the case for any A-Ci curves measured at
  saturating light. In this case, the optional \code{cj_crossover_min} and
  \code{cj_crossover_max} can be used to constrain the range of \code{Cc} values
  (in ppm) where \code{Aj} is allowed to be the overall rate limiting factor.
  If the crossover from Rubisco-limited to RuBP-regeneration limited
  assimilation occurs outside these bounds (when they are supplied), a heavy
  penalty will be added to the error function, preventing the optimizer from
  choosing those parameter values.

  A penalty is also added for any parameter combination where \code{An} is not a
  number, or where \code{\link{calculate_c3_assimilation}} produces an error.
}

\value{
  A function with one input argument \code{guess}, which should be a numeric
  vector representing values of the parameters to be fitted (which are specified
  by the \code{fit_options} input argument.) Each element of \code{guess} is the
  value of one parameter (arranged in alphabetical order.) For example, with the
  default settings, \code{guess} should contain values of \code{alpha_old},
  \code{J_at_25}, \code{RL_at_25}, \code{Tp_at_25}, and \code{Vcmax_at_25} (in
  that order).
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c3_aci_1.xlsx')
)

# Define a new column that uniquely identifies each curve
licor_file[, 'species_plot'] <-
  paste(licor_file[, 'species'], '-', licor_file[, 'plot'] )

# Organize the data
licor_file <- organize_response_curve_data(
    licor_file,
    'species_plot',
    c(9, 10, 16),
    'CO2_r_sp'
)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_file <- calculate_temperature_response(licor_file, c3_temperature_param_bernacchi)

# Define an error function for one curve from the set
error_fcn <- error_function_c3_aci(
  licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE]
)

# Evaluate the error for:
#  alpha_old = 0
#  J_at_25 = 236
#  RL_at_25 = 4e-8
#  Tp_at_25 = 22.7
#  Vcmax_at_25 = 147
error_fcn(c(0, 236, 4e-8, 22.7, 147))

# Make a plot of likelihood vs. Vcmax when other parameters are fixed to the
# values above.
vcmax_error_fcn <- function(Vcmax) {error_fcn(c(0, 236, 4e-8, 22.7, Vcmax))}
vcmax_seq <- seq(135, 152, length.out = 41)

lattice::xyplot(
  exp(-sapply(vcmax_seq, vcmax_error_fcn)) ~ vcmax_seq,
  type = 'b',
  xlab = 'Vcmax_at_25 (micromol / m^2 / s)',
  ylab = 'Negative log likelihood (dimensionless)'
)
}

\concept{exdf}
