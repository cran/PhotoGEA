\name{initial_guess_c3_aci}

\alias{initial_guess_c3_aci}

\title{Make an initial guess of FvCB model parameter values for one curve}

\description{
  Creates a function that makes an initial guess of FvCB model parameter values
  for one curve. This function is used internally by \code{\link{fit_c3_aci}}.

  Values estimated by this guessing function should be considered inaccurate,
  and should always be improved upon by an optimizer.
}

\usage{
  initial_guess_c3_aci(
    alpha_g,
    alpha_old,
    alpha_s,
    alpha_t,
    Gamma_star_at_25,
    gmc_at_25,
    Kc_at_25,
    Ko_at_25,
    cc_threshold_rl = 100,
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
    debug_mode = FALSE
  )
}

\arguments{
  \item{alpha_g}{
    A dimensionless parameter where \code{0 <= alpha_g <= 1}, representing the
    proportion of glycolate carbon taken out of the photorespiratory pathway as
    glycine. \code{alpha_g} is often assumed to be 0. If \code{alpha_g} is not a
    number, then there must be a column in \code{rc_exdf} called \code{alpha_g}
    with appropriate units. A numeric value supplied here will overwrite the
    values in the \code{alpha_g} column of \code{rc_exdf} if it exists.
  }

  \item{alpha_old}{
    A dimensionless parameter where \code{0 <= alpha_old <= 1}, representing the
    fraction of remaining glycolate carbon not returned to the chloroplast
    after accounting for carbon released as CO2. \code{alpha_old} is often
    assumed to be 0. If \code{alpha_old} is not a number, then there must be a
    column in \code{rc_exdf} called \code{alpha_old} with appropriate units.
    A numeric value supplied here will overwrite the values in the
    \code{alpha_old} column of \code{rc_exdf} if it exists.
  }

  \item{alpha_s}{
    A dimensionless parameter where \code{0 <= alpha_s <= 0.75 * (1 - alpha_g)}
    representing the proportion of glycolate carbon taken out of the
    photorespiratory pathway as serine. \code{alpha_s} is often assumed to be 0.
    If \code{alpha_s} is not a number, then there must be a column in
    \code{rc_exdf} called \code{alpha_s} with appropriate units. A numeric
    value supplied here will overwrite the values in the \code{alpha_s} column
    of \code{rc_exdf} if it exists.
  }

  \item{alpha_t}{
    A dimensionless parameter where \code{0 <= alpha_t <= 1} representing the
    proportion of glycolate carbon taken out of the photorespiratory pathway as
    CH2-THF. \code{alpha_t} is often assumed to be 0. If \code{alpha_t} is not a
    number, then there must be a column in \code{rc_exdf} called \code{alpha_t}
    with appropriate units. A numeric value supplied here will overwrite the
    values in the \code{alpha_t} column of \code{rc_exdf} if it exists.
  }

  \item{Gamma_star_at_25}{
    The chloroplastic CO2 concentration at which CO2 gains from Rubisco
    carboxylation are exactly balanced by CO2 losses from Rubisco oxygenation,
    at 25 degrees C, expressed in \code{micromol mol^(-1)}. If
    \code{Gamma_star_at_25} is not a number, then there must be a column in
    \code{rc_exdf} called \code{Gamma_star_at_25} with appropriate units. A
    numeric value supplied here will overwrite the values in the
    \code{Gamma_star_at_25} column of \code{rc_exdf} if it exists.
  }

  \item{gmc_at_25}{
    The mesophyll conductance to CO2 diffusion at 25 degrees C, expressed in
    \code{mol m^(-2) s^(-1) bar^(-1)}. In the absence of other reliable
    information, \code{gmc_at_25} is often assumed to be infinitely large. If
    \code{gmc_at_25} is not a number, then there must be a column in
    \code{rc_exdf} called \code{gmc_at_25} with appropriate units. A numeric
    value supplied here will overwrite the values in the \code{gmc_at_25}
    column of \code{rc_exdf} if it exists.
  }

  \item{Kc_at_25}{
    The Michaelis-Menten constant for Rubisco carboxylation at 25 degrees C,
    expressed in \code{micromol mol^(-1)}. If \code{Kc_at_25} is not a number,
    then there must be a column in \code{rc_exdf} called \code{Kc_at_25} with
    appropriate units. A numeric value supplied here will overwrite the values
    in the \code{Kc_at_25} column of \code{rc_exdf} if it exists.
  }

  \item{Ko_at_25}{
    The Michaelis-Menten constant for Rubisco oxygenation at 25 degrees C,
    expressed in \code{mmol mol^(-1)}. If \code{Ko_at_25} is not a number, then
    there must be a column in \code{rc_exdf} called \code{Ko_at_25} with
    appropriate units. A numeric value supplied here will overwrite the values
    in the \code{Ko_at_25} column of \code{rc_exdf} if it exists.
  }

  \item{cc_threshold_rl}{
    An upper cutoff value for the chloroplast CO2 concentration in
    \code{micromol mol^(-1)} to be used when estimating \code{RL}.
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
    The name of the column in \code{rc_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{rc_exdf} that contains the
    intercellular CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{gamma_star_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{Gamma_star} values (with units of
    \code{normalized to Gamma_star at 25 degrees C}).
  }

  \item{gmc_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    mesophyll conductance values (with units of
    \code{normalized to gmc at 25 degrees C}).
  }

  \item{j_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{J} values (with units of \code{normalized to J at 25 degrees C}).
  }

  \item{kc_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{Kc} values (with units of \code{normalized to Kc at 25 degrees C}).
  }

  \item{ko_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{Ko} values (with units of \code{normalized to Ko at 25 degrees C}).
  }

  \item{oxygen_column_name}{
    The name of the column in \code{rc_exdf} that contains the concentration
    of O2 in the ambient air, expressed as a percentage (commonly 21\% or 2\%);
    the units must be \code{percent}.
  }

  \item{rl_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{RL} values (with units of \code{normalized to RL at 25 degrees C}).
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{rc_exdf} that contains the total
    pressure in \code{bar}.
  }

  \item{tp_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{Tp} values (with units of \code{normalized to Tp at 25 degrees C}).
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{debug_mode}{
    A logical (\code{TRUE} or \code{FALSE}) variable indicating whether to
    operate in debug mode. In debug mode, information about the linear fit used
    to estimate \code{RL} is printed; this can be helpful when troubleshooting
    issues with a particular curve.
  }
}

\details{
  Here we estimate values of \code{J_at_25}, \code{RL_at_25}, \code{Tp_at_25},
  and \code{Vcmax_at_25} from a measured C3 CO2 response curve. It is difficult
  to estimate values of \code{alpha_g}, \code{alpha_old}, \code{alpha_s},
  \code{alpha_t}, \code{Gamma_star_at_25}, \code{gmc_at_25}, \code{Kc_at_25},
  \code{Ko_at_25} from a curve, so they must be supplied beforehand. For more
  information about these parameters, see the documentation for
  \code{\link{calculate_c3_assimilation}}.

  \itemize{
    \item \strong{Estimating RL}: Regardless of which process is limiting at low
          \code{Cc}, it is always true that \code{An = -RL} when
          \code{Cc = Gamma_star_agt}. Here we make a linear fit of the measured
          \code{An} vs. \code{Cc} values where \code{Cc} is below
          \code{cc_threshold_rl}, and evaluate it at at
          \code{Cc = Gamma_star_agt} to estimate \code{RL}. If there are fewer
          than two points with \code{Cc <= cc_threshold_rl}, the fit cannot be
          made, and we use a typical value instead
          (1.0 \code{micromol m^(-2) s^(-1)}). Likewise, if the linear fit
          predicts a negative or \code{NA} value for \code{RL}, we use the same
          typical value instead.
    \item \strong{Estimating Vc}: Once an estimate for \code{RL} has been found,
          the RuBP carboxylation rate \code{Vc} can be estimated using
          \code{Vc = (An + RL) / (1 - Gamma_star_agt / Cc)}. This is useful for
          the remaining parameter estimates.
    \item \strong{Estimating Vcmax}: An estimate for \code{Vcmax} can be
          obtained by solving the equation for \code{Wc} for \code{Vcmax}, and
          evaluating it with \code{Wc = Vc} as estimated above. In the
          rubisco-limited part of the curve, \code{Vc = Wc} and the estimated
          values of \code{Vcmax} should be reasonable. In other parts of the
          curve, \code{Wc} is not the limiting rate, so \code{Vc < Wc}.
          Consequently, the estimated values of \code{Vcmax} in these parts of
          the curve will be smaller. So, to make an overall estimate, we choose
          the the largest estimated \code{Vcmax} value.
    \item \strong{Estimating J and Tp}: Estimates for these parameters can be
          made using the equations for \code{Wj} and \code{Wp}, similar to the
          approach followed for \code{Vcmax}.
  }

  For the parameter values estimated above, the values of \code{RL_norm},
  \code{Vcmax_norm}, and \code{J_norm} are used to convert the values at leaf
  temperature to the values at 25 degrees C.
}

\value{
  A function with one input argument \code{rc_exdf}, which should be an
  \code{exdf} object representing one C3 CO2 response curve. The return value of
  this function will be a numeric vector with twelve elements, representing the
  values of \code{alpha_g}, \code{alpha_old}, \code{alpha_s}, \code{alpha_t},
  \code{Gamma_star_at_25}, \code{gmc_at_25}, \code{J_at_25}, \code{Kc_at_25},
  \code{Ko_at_25}, \code{RL_at_25}, \code{Tp_at_25}, and \code{Vcmax_at_25} (in
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

# Create the guessing function; here we set:
# - All alpha values to 0
# - Gamma_star_at_25 to 40 micromol / mol
# - gmc to infinity
# - Kc_at_25 to 400 micromol / mol
# - Ko_at_25 to 275 mmol / mol
guessing_func <- initial_guess_c3_aci(
  alpha_g = 0,
  alpha_old = 0,
  alpha_s = 0,
  alpha_t = 0,
  Gamma_star = 40,
  gmc_at_25 = Inf,
  Kc_at_25 = 400,
  Ko_at_25 = 275
)

# Apply it and see the initial guesses for each curve
print(by(licor_file, licor_file[, 'species_plot'], guessing_func))

# A simple way to visualize the guesses is to "fit" the curves using the null
# optimizer, which simply returns the initial guess
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c3_aci,
  fit_options = list(alpha_old = 0),
  optim_fun = optimizer_null(),
  remove_unreliable_param = 0
))

plot_c3_aci_fit(aci_results, 'species_plot', 'Ci')
}

\concept{exdf}
