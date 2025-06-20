\name{fit_c4_aci}

\alias{fit_c4_aci}

\title{Fits a C4 assimilation model to an A-Ci curve}

\description{
  Fits the von Caemmerer model to an experimentally measured C4 A-Ci curve.

  It is possible to fit the following parameters: \code{alpha_psii}, \code{gbs},
  \code{gmc_at_25}, \code{J_at_25}, \code{RL_at_25}, \code{Rm_frac},
  \code{Vcmax_at_25}, \code{Vpmax_at_25}, and \code{Vpr}.

  By default, only a subset of these parameters are actually fit:
  \code{RL_at_25}, \code{Vcmax_at_25}, and \code{Vpmax_at_25}. This can be
  altered using the \code{fit_options} argument, as described below.

  Best-fit parameters are found using maximum likelihood fitting, where the
  optimizer (\code{optim_fun}) is used to minimize the error function (defined
  by \code{\link{error_function_c4_aci}}).

  Once best-fit parameters are found, confidence intervals are calculated
  using \code{\link{confidence_intervals_c4_aci}}, and unreliable parameter
  estimates are removed.

  For temperature-dependent parameters, best-fit values and confidence intervals
  are returned at 25 degrees C and at leaf temperature.

  See below for more details.
}

\usage{
  fit_c4_aci(
    replicate_exdf,
    Ca_atmospheric = NA,
    ao_column_name = 'ao',
    a_column_name = 'A',
    ca_column_name = 'Ca',
    ci_column_name = 'Ci',
    gamma_star_column_name = 'gamma_star',
    gmc_norm_column_name = 'gmc_norm',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    kp_column_name = 'Kp',
    oxygen_column_name = 'oxygen',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm',
    sd_A = 'RMSE',
    x_etr = 0.4,
    optim_fun = optimizer_deoptim(200),
    lower = list(),
    upper = list(),
    fit_options = list(),
    relative_likelihood_threshold = 0.147,
    hard_constraints = 0,
    calculate_confidence_intervals = TRUE,
    remove_unreliable_param = 2,
    debug_mode = FALSE
  )
}

\arguments{
  \item{replicate_exdf}{
    An \code{exdf} object representing one CO2 response curve.
  }

  \item{Ca_atmospheric}{
    The atmospheric CO2 concentration (with units of \code{micromol mol^(-1)});
    this will be used by \code{\link{estimate_operating_point}} to estimate the
    operating point. A value of \code{NA} disables this feature.
  }

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ao_column_name}{
    The name of the column in \code{exdf_obj} that contains the dimensionless
    ratio of solubility and diffusivity of O2 to CO2.
  }

  \item{ca_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    ambient CO2 concentration in \code{micromol mol^(-1)}. If values of
    \code{Ca} are not available, they can be set to \code{NA}. In this case, it
    will not be possible to estimate the operating point, and
    \code{\link{apply_gm}} will not be able to calculate the CO2 drawdown across
    the stomata.
  }

  \item{ci_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    intercellular CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{gamma_star_column_name}{
    The name of the column in \code{exdf_obj} that contains the dimensionless
    \code{gamma_star} values.
  }

  \item{gmc_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    mesophyll conductance values (with units of
    \code{normalized to gmc at 25 degrees C}).
  }

  \item{j_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{J} values (with units of
    \code{normalized to J at 25 degrees C}).
  }

  \item{kc_column_name}{
    The name of the column in \code{exdf_obj} that contains the Michaelis-Menten
    constant for rubisco carboxylation in \code{microbar}.
  }

  \item{ko_column_name}{
    The name of the column in \code{exdf_obj} that contains the Michaelis-Menten
    constant for rubisco oxygenation in \code{mbar}.
  }

  \item{kp_column_name}{
    The name of the column in \code{exdf_obj} that contains the Michaelis-Menten
    constant for PEP carboxylase carboxylation in \code{microbar}.
  }

  \item{oxygen_column_name}{
    The name of the column in \code{exdf_obj} that contains the concentration
    of O2 in the ambient air, expressed as a percentage (commonly 21\% or 2\%);
    the units must be \code{percent}.
  }

  \item{rl_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{RL} values (with units of \code{normalized to RL at 25 degrees C}).
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{exdf_obj} that contains the total pressure
    in \code{bar}.
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{vpmax_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Vpmax} values (with units of
    \code{normalized to Vpmax at 25 degrees C}).
  }

  \item{sd_A}{
    A value of the standard deviation of measured \code{A} values, or the name
    of a method for determining the deviation; currently, the only supported
    option is \code{'RMSE'}.
  }

  \item{x_etr}{
    The fraction of whole-chain electron transport occurring in the mesophyll
    (dimensionless). See Equation 29 from S. von Caemmerer (2021).
  }

  \item{optim_fun}{
    An optimization function that accepts the following input arguments: an
    initial guess, an error function, lower bounds, and upper bounds. It should
    return a list with the following elements: \code{par}, \code{convergence},
    \code{feval}, and \code{convergence_msg}. See \code{\link{optimizers}}
    for a list of available options.
  }

  \item{lower}{
    A list of named numeric elements representing lower bounds to use when
    fitting. Values supplied here override the default values (see details
    below). For example, \code{lower = list(Vcmax_at_25 = 10)} sets the lower
    limit for \code{Vcmax_at_25} to 10 micromol / m^2 / s.
  }

  \item{upper}{
    A list of named numeric elements representing upper bounds to use when
    fitting. Values supplied here override the default values (see details
    below). For example, \code{upper = list(Vcmax_at_25 = 200)} sets the upper
    limit for \code{Vcmax_at_25} to 200 micromol / m^2 / s.
  }

  \item{fit_options}{
    A list of named elements representing fit options to use for each parameter.
    Values supplied here override the default values (see details below). Each
    element must be \code{'fit'}, \code{'column'}, or a numeric value. A value
    of \code{'fit'} means that the parameter will be fit; a value of
    \code{'column'} means that the value of the parameter will be taken from a
    column in \code{exdf_obj} of the same name; and a numeric value means that
    the parameter will be set to that value. For example, \code{fit_options =
    list(RL_at_25 = 0, Vcmax_at_25 = 'fit', Vpr = 'column')} means that
    \code{RL_at_25} will be set to 0, \code{Vcmax_at_25} will be fit, and
    \code{Vpr} will be set to the values in the \code{Vpr} column of
    \code{exdf_obj}.
  }

  \item{relative_likelihood_threshold}{
    To be passed to \code{\link{confidence_intervals_c4_aci}} when
    \code{calculate_confidence_intervals} is \code{TRUE}.
  }

  \item{hard_constraints}{
    To be passed to \code{\link{calculate_c4_assimilation}}; see that function
    for more details.
  }

  \item{calculate_confidence_intervals}{
    A logical value indicating whether or not to estimate confidence intervals
    for the fitting parameters using \code{\link{confidence_intervals_c4_aci}}.
  }

  \item{remove_unreliable_param}{
    An integer value indicating the rules to use when identifying and removing
    unreliable parameter estimates. A value of 2 is the most conservative
    option. A value of 0 disables this feature, which is not typically
    recommended. It is also possible to directly specify the trust values to
    remove; for example, \code{'unreliable (process never limiting)'} is
    equivalent to 1. See below for more details.
  }

  \item{debug_mode}{
    A logical (\code{TRUE} or \code{FALSE}) variable indicating whether to
    operate in debug mode. In debug mode, information about
    \code{replicate_exdf}, the initial guess, each guess supplied from the
    optimizer, and the final outcome is printed; this can be helpful when
    troubleshooting issues with a particular curve.
  }
}

\details{
  This function calls \code{\link{calculate_c4_assimilation}} to calculate
  values of net assimilation. The user-supplied optimization function is used to
  vary the values of \code{alpha_psii}, \code{gbs}, \code{gmc_at_25},
  \code{J_at_25}, \code{RL_at_25}, \code{Rm_frac}, \code{Vcmax_at_25},
  \code{Vpmax_at_25}, and \code{Vpr} to find ones that best reproduce the
  experimentally measured values of net assimilation. By default, the following
  options are used for the fits:
  \itemize{
    \item \code{alpha_psii}: lower = -1, upper = 10, fit_option = 0
    \item \code{gbs}: lower = -1, upper = 10, fit_option = 0.003
    \item \code{gmc_at_25}: lower = -1, upper = 10, fit_option = 1
    \item \code{J_at_25}: lower = -50, upper = 1000, fit_option = 1000
    \item \code{RL_at_25}: lower = -10, upper = 100, fit_option = \code{'fit'}
    \item \code{Rm_frac}: lower = -10, upper = 10, fit_option = 0.5
    \item \code{Vcmax_at_25}: lower = -50, upper = 1000, fit_option = \code{'fit'}
    \item \code{Vpmax_at_25}: lower = -50, upper = 1000, fit_option = \code{'fit'}
    \item \code{Vpr}: lower = -50, upper = 1000, fit_option = 1000
  }
  With these settings, \code{J_at_25} and \code{Vpr} are set to 1000 (so net
  assimilation is essentially never limited by light or PEP carboxylase
  regeneration), \code{alpha_psii}, \code{gbs}, \code{gmc_at_25}, and
  \code{Rm_frac} are set to default values used in von Caemmerer (2000), and the
  other parameters are fit during the process (see \code{fit_options} above).
  The bounds are chosen liberally to avoid any bias.

  An initial guess for the parameters is generated by calling
  \code{\link{initial_guess_c4_aci}} as follows:
  \itemize{
    \item \code{pcm_threshold_rlm} is set to 40 microbar.
    \item If \code{alpha_psii} is being fit, the \code{alpha_psii} argument of
          \code{initial_guess_c4_aci} is set to 0.1; otherwise, the argument is
          set to the value specified by the fit options.
    \item If \code{gbs} is being fit, the \code{gbs} argument of
          \code{initial_guess_c4_aci} is set to 0.003; otherwise, the argument
          is set to the value specified by the fit options.
    \item If \code{gmc_at_25} is being fit, the \code{gmc_at_25} argument of
          \code{initial_guess_c4_aci} is set to 1; otherwise, the argument is
          set to the value specified by the fit options.
    \item If \code{Rm_frac} is being fit, the \code{Rm_frac} argument of
          \code{initial_guess_c4_aci} is set to 0.5; otherwise, the argument is
          set to the value specified by the fit options.
  }
  Note that any fixed values specified in the fit options will override the
  values returned by the guessing function.

  The fit is made by creating an error function using
  \code{\link{error_function_c4_aci}} and minimizing its value using
  \code{optim_fun}, starting from the initial guess described above.  The
  \code{\link{optimizer_deoptim}} optimizer is used by default since it has been
  found to reliably return great fits. However, it is a slow optimizer. If speed
  is important, consider reducing the number of generations or using
  \code{\link{optimizer_nmkb}}, but be aware that this optimizer is more likely
  to get stuck in a local minimum.

  The photosynthesis model represented by \code{calculate_c4_assimilation} is
  not smooth in the sense that small changes in the input parameters do not
  necessarily cause changes in its outputs. This is related to the calculation
  of the PEP carboxylase activity \code{Vp}, which is taken to be the minimum of
  \code{Vpr} and \code{Vpc}. For example, if \code{Vpr} is high and \code{Vp =
  Vpc} at all points along the curve, modifying \code{Vpr} by a small amount
  will not change the model's outputs. Similar issues can occur when calculating
  \code{An} as the minimum of \code{Ac} and \code{Aj}. Because of this,
  derivative-based optimizers tend to struggle when fitting C4 A-Ci curves. Best
  results are obtained using derivative-free methods.

  Sometimes the optimizer may choose a set of parameter values where one of the
  potential limiting rates \code{Vpc} or \code{Vpr} is never the smallest rate.
  In this case, the corresponding parameter estimates (\code{Vpmax} or
  \code{Vpr}) will be severely unreliable. Likewise, it may happen that one of
  \code{Ac} or \code{Aj} is never the smallest rate. In this case the
  corresponding parameter estimates (\code{Vpmax}, \code{Vpr}, and \code{Vcmax},
  or \code{J}) will be severely unreliable. This will be indicated by a value
  of \code{'unreliable (process never limiting)'} in the corresponding trust
  column (for example, \code{Vcmax_trust}). If \code{remove_unreliable_param} is
  \code{1} or larger, then such parameter estimates (and the corresponding
  rates) will be replaced by \code{NA} in the fitting results.

  It is also possible that the upper limit of the confidence interval for a
  parameter is infinity; this indicates a potentially unreliable parameter
  estimate. This will be indicated by a value of
  \code{'unreliable (infinite upper limit)'} in the corresponding trust column
  (for example, \code{Vcmax_trust}). If \code{remove_unreliable_param} is
  \code{2} or larger, then such parameter estimates (but not the corresponding
  rates) will be replaced by \code{NA} in the fitting results.

  The trust value for fully reliable parameter estimates is set to
  \code{'reliable'} and they will never be replaced by \code{NA}.

  Once the best-fit parameters have been determined, this function also
  estimates the operating value of `\code{PCm} from the atmospheric CO2
  concentration \code{atmospheric_ca} using
  \code{\link{estimate_operating_point}}, and then uses that value to estimate
  the modeled \code{An} at the operating point via
  \code{\link{calculate_c4_assimilation}}. It also estimates the
  \href{https://en.wikipedia.org/wiki/Akaike_information_criterion}{Akaike information criterion (AIC)}.

  This function assumes that \code{replicate_exdf} represents a single
  C4 A-Ci curve. To fit multiple curves at once, this function is often used
  along with \code{\link{by.exdf}} and \code{\link{consolidate}}.
}

\value{
  A list with two elements:
  \itemize{
    \item \code{fits}: An \code{exdf} object including the original contents of
          \code{replicate_exdf} along with several new columns:
          \itemize{
            \item The fitted values of net assimilation will be stored in a
                  column whose name is determined by appending \code{'_fit'} to
                  the end of \code{a_column_name}; typically, this will be
                  \code{'A_fit'}.
            \item Residuals (measured - fitted) will be stored in a column whose
                  name is determined by appending \code{'_residuals'} to the end
                  of \code{a_column_name}; typically, this will be
                  \code{'A_residuals'}.
            \item Values of fitting parameters at 25 degrees C will be stored in
                  the \code{gmc_at_25}, \code{J_at_25}, \code{RL_at_25},
                  \code{Vcmax_at_25}, \code{Vpmax_at_25}, and \code{Vpr}
                  columns.
            \item The other outputs from \code{\link{calculate_c4_assimilation}}
                  will be stored in columns with the usual names:
                  \code{alpha_psii}, \code{gbs}, \code{gmc_tl}, \code{Rm_Frac},
                  \code{Vcmax_tl}, \code{Vpmax_tl}, \code{RL_tl}, \code{RLm_tl},
                  \code{Vp}, \code{Apc}, \code{Apr}, \code{Ap}, \code{Ar},
                  \code{Ajm}, \code{Ajbs}, \code{Ac}, and \code{Aj}.
          }
    \item \code{fits_interpolated}: An \code{exdf} object including the
          calculated assimilation rates at a fine spacing of \code{Ci} values
          (step size of 1 \code{micromol mol^(-1)}).
    \item \code{parameters}: An \code{exdf} object including the identifiers,
          fitting parameters, and convergence information for the A-Ci curve:
          \itemize{
            \item The number of points where \code{Vpc} and \code{Vpr} are each
                  the smallest potential carboxylation rate are stored in the
                  \code{n_Vpc_smallest} and \code{n_Vpr_smallest} columns.
            \item The best-fit values are stored in the \code{alpha_psii},
                  \code{gbs}, \code{gmc_at_25}, \code{J_at_25},
                  \code{RL_at_25}, \code{Rm_frac}, \code{Vcmax_at_25},
                  \code{Vpmax_at_25}, and \code{Vpr} columns. If
                  \code{calculate_confidence_intervals} is \code{TRUE}, upper
                  and lower limits for each of these parameters will also be
                  included.
            \item For parameters that depend on leaf temperature, the average
                  leaf-temperature-dependent values are stored in
                  \code{X_tl_avg} columns: \code{gmc_tl_avg}, \code{J_tl_avg},
                  \code{Jmax_tl_avg}, \code{RL_tl_avg}, \code{Vcmax_tl_avg}, and
                  \code{Vpmax_tl_avg}.
            \item The average leaf temperature is also stored in the
                  \code{Tleaf_avg} column.
            \item Information about the operating point is stored in
                  \code{operating_PCm}, \code{operating_Ci},
                  \code{operating_An}, and \code{operating_An_model}.
            \item The \code{convergence} column indicates whether the fit was
                  successful (\code{==0}) or if the optimizer encountered a
                  problem (\code{!=0}).
            \item The \code{feval} column indicates how many cost function
                  evaluations were required while finding the optimal parameter
                  values.
            \item The residual stats as returned by \code{\link{residual_stats}}
                  are included as columns with the default names: \code{dof},
                  \code{RSS}, \code{RMSE}, etc.
            \item The Akaike information criterion is included in the \code{AIC}
                  column.
          }
  }
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c4_aci_1.xlsx')
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

# Calculate temperature-dependent values of C4 photosynthetic parameters
licor_file <- calculate_temperature_response(licor_file, c4_temperature_param_vc)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# For these examples, we will use a faster (but sometimes less reliable)
# optimizer so they run faster
optimizer <- optimizer_nmkb(1e-7)

# Fit just one curve from the data set (it is rare to do this).
one_result <- fit_c4_aci(
  licor_file[licor_file[, 'species_plot'] == 'maize - 5', , TRUE],
  Ca_atmospheric = 420,
  optim_fun = optimizer
)

# Fit all curves in the data set (it is more common to do this)
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c4_aci,
  Ca_atmospheric = 420,
  optim_fun = optimizer
))

# View the fitting parameters for each species / plot
col_to_keep <- c(
  'species', 'plot',                                       # identifiers
  'RL_at_25', 'Vcmax_at_25', 'Vpmax_at_25', 'Vpr',         # parameters scaled to 25 degrees C
  'RL_tl_avg', 'Vcmax_tl_avg', 'Vpmax_tl_avg',             # average temperature-dependent values
  'operating_Ci', 'operating_An', 'operating_An_model',    # operating point info
  'dof', 'RSS', 'MSE', 'RMSE', 'RSE',                      # residual stats
  'convergence', 'convergence_msg', 'feval', 'optimum_val' # convergence info
)

aci_results$parameters[ , col_to_keep, TRUE]

# View the fits for each species / plot
plot_c4_aci_fit(aci_results, 'species_plot', 'Ci', ylim = c(0, 100))

# View the residuals for each species / plot
lattice::xyplot(
  A_residuals ~ Ci | species_plot,
  data = aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste('Intercellular CO2 concentration [', aci_results$fits$units$Ci, ']'),
  ylab = paste('Assimilation rate residuals [', aci_results$fits$units$A_residuals, ']')
)
}

\concept{exdf}
