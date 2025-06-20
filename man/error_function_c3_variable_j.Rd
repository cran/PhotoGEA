\name{error_function_c3_variable_j}

\alias{error_function_c3_variable_j}

\title{Generate an error function for C3 Variable J curve fitting}

\description{
  Creates a function that returns an error value (the negative of the natural
  logarithm of the likelihood) representing the amount of agreement between
  modeled and measured \code{An} values. When this function is minimized, the
  likelihood is maximized.

  Internally, this function uses \code{link{calculate_c3_variable_j}} and
  \code{link{calculate_c3_assimilation}} to calculate assimilation rate values
  that are compared to the measured ones.
}

\usage{
  error_function_c3_variable_j(
    replicate_exdf,
    fit_options = list(),
    sd_A = 1,
    Wj_coef_C = 4.0,
    Wj_coef_Gamma_star = 8.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    gamma_star_norm_column_name = 'Gamma_star_norm',
    j_norm_column_name = 'J_norm',
    kc_norm_column_name = 'Kc_norm',
    ko_norm_column_name = 'Ko_norm',
    oxygen_column_name = 'oxygen',
    phips2_column_name = 'PhiPS2',
    qin_column_name = 'Qin',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    tp_norm_column_name = 'Tp_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    hard_constraints = 0,
    require_positive_gmc = 'positive_a',
    gmc_max = Inf,
    check_j = TRUE,
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

  \item{phips2_column_name}{
    The name of the column in \code{replicate_exdf} that contains values of the
    operating efficiency of photosystem II (dimensionless).
  }

  \item{qin_column_name}{
    The name of the column in \code{replicate_exdf} that contains values of the
    incident photosynthetically active flux density in
    \code{micromol m^(-2) s^(-1)}.
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
    To be passed to \code{\link{calculate_c3_assimilation}} and
    \code{\link{calculate_c3_variable_j}}; see those functions for more details.
  }

  \item{require_positive_gmc}{
    A character string specifying the method to be used for requiring positive
    values of mesophyll conductance. Can be \code{'none'}, \code{'all'}, or
    \code{'positive_a'}. See below for more details.
  }

  \item{gmc_max}{
    The maximum value of mesophyll conductance that should be considered to be
    acceptable. See below for more details.
  }

  \item{check_j}{
    A logical (TRUE/FALSE) value indicating whether to check whether
    \code{J_F > J_tl}. See below for more details.
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
  When fitting A-Ci + chlorophyll fluorescence curves using the Variable J
  method, it is necessary to define a function that calculates the likelihood
  of a given set of \code{alpha_g}, \code{alpha_old}, \code{alpha_s},
  \code{alpha_t}, \code{Gamma_star}, \code{J_at_25}, \code{RL_at_25},
  \code{tau}, \code{Tp_at_25}, and \code{Vcmax_at_25} values by comparing a
  model prediction to a measured curve. This function will be passed to an
  optimization algorithm which will determine the values that produce the
  smallest error.

  The \code{error_function_c3_variable_j} returns such a function, which is
  based on a particular replicate and a set of fitting options. It is possible
  to just fit a subset of the available fitting parameters; by default, the
  fitting parameters are \code{alpha_old}, \code{J_at_25}, \code{RL_at_25},
  \code{Tp_at_25}, \code{tau}, and \code{Vcmax_at_25}. This behavior can be
  changed via the \code{fit_options} argument.

  For practical reasons, the function actually returns values of \code{-ln(L)},
  where \code{L} is the likelihood. The logarithm of \code{L} is simpler to
  calculate than \code{L} itself, and the minus sign converts the problem from
  a maximization to a minimization, which is important because most optimizers
  are designed to minimize a value.

  Sometimes an optimizer will choose biologically unreasonable parameter values
  that nevertheless produce good fits to the supplied assimilation values. There
  are several options for preventing an optimizer from choosing such parameter
  values:
  \itemize{
    \item \code{Enforcing Rubisco limitations}: A common problem is that the fit
          result may not indicate Rubisc-limited assimilation at low CO2 values,
          which should be the case for any A-Ci curves measured at saturating
          light. In this case, the optional \code{cj_crossover_min} and
          \code{cj_crossover_max} can be used to constrain the range of
          \code{Cc} values (in ppm) where \code{Wj} is allowed to be the overall
          rate limiting factor. If the crossover from Rubisco-limited to
          RuBP-regeneration limited carboxylation occurs outside these bounds
          (when they are supplied), a heavy penalty will be added to the error
          function, preventing the optimizer from choosing those parameter
          values.
    \item \code{Requiring positive gmc}: The Variable J method sometimes
          predicts negative values of the mesophyll conductance (\code{gmc}).
          Such values are probably not realistic, especially when \code{Cc} is
          above the CO2 compensation point. The \code{require_positive_gmc}
          input argument can be used to penalize negative values of \code{gmc}.
          When \code{require_positive_gmc} is \code{'all'}, a heavy penalty will
          be added to the error function if there are any negative \code{gmc}
          values. When \code{require_positive_gmc} is \code{'positive_a'}, a
          heavy penalty will be added to the error function if there are any
          negative \code{gmc} values when \code{A} is positive; negative
          \code{gmc} for negative \code{A} will be allowed. When
          \code{require_positive_gmc} is \code{'none'}, these restrictions are
          disabled and no penalties are added for negative \code{gmc}.
    \item \code{Preventing large values of gmc}: The Variable J method sometimes
          produces unreasonably high values of \code{gmc}. The \code{gmc_max}
          argument can be used to address this. If any predicted \code{gmc}
          values exceed \code{gmc_max} when \code{A} is positive, a heavy
          penalty will be added to the error function.
    \item \code{Enforcing consistent RuBP regeneration rates}: In principle, the
          actual RuBP regeneration rate (\code{J_F}) should always be less than
          or equal to its maximum value for a given \code{Qin} and leaf
          temperature (\code{J_tl}), with equality only occuring when
          assimilation is RuBP-regeneration-limited. When \code{check_j} is
          \code{TRUE}, a heavy penalty will be added to the error function for
          any parameter values where \code{J_F} is greater than \code{J_tl} at
          any point in the curve.
  }

  A penalty is also added for any parameter combination where \code{An} is not a
  number, or where \code{\link{calculate_c3_variable_j}} or
  \code{\link{calculate_c3_assimilation}} produce an error.
}

\value{
  A function with one input argument \code{guess}, which should be a numeric
  vector representing values of the parameters to be fitted (which are specified
  by the \code{fit_options} input argument.) Each element of \code{guess} is the
  value of one parameter (arranged in alphabetical order.) For example, with the
  default settings, \code{guess} should contain values of \code{alpha_old},
  \code{J_at_25}, \code{RL_at_25}, \code{tau}, \code{Tp_at_25}, and
  \code{Vcmax_at_25} (in that order).
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
error_fcn <- error_function_c3_variable_j(
  licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE]
)

# Evaluate the error for:
#  alpha_old = 1.9
#  J_at_25 = 270
#  RL_at_25 = 1.9
#  tau = 0.42
#  Tp_at_25 = 8.7
#  Vcmax_at_25 = 258
error_fcn(c(1.9, 270, 1.9, 0.42, 8.7, 258))

# Make a plot of error vs. Tp_at_25 when the other parameters are fixed to the
# values above. As Tp_at_25 increases, it eventually stops limiting the
# assimilation rate and its value stops influencing the error.
tpu_error_fcn <- function(Tp_at_25) {error_fcn(c(1.9, 270, 1.9, 0.42, Tp_at_25, 258))}
tpu_seq <- seq(5, 12, by = 0.25)

lattice::xyplot(
  sapply(tpu_seq, tpu_error_fcn) ~ tpu_seq,
  type = 'b',
  xlab = 'Tp at 25 degrees C (micromol / m^2 / s)',
  ylab = 'Negative log likelihood (dimensionless)'
)
}

\concept{exdf}
