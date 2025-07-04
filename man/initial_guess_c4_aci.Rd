\name{initial_guess_c4_aci}

\alias{initial_guess_c4_aci}

\title{Make an initial guess of C4 photosynthesis parameter values for one curve}

\description{
  Creates a function that makes an initial guess of C4 photosynthesis model
  parameter values for one curve. This function is used internally by
  \code{\link{fit_c4_aci}}.

  Values estimated by this guessing function should be considered inaccurate,
  and should always be improved upon by an optimizer.
}

\usage{
  initial_guess_c4_aci(
    alpha_psii,
    gbs,
    gmc_at_25,
    Rm_frac,
    pcm_threshold_rlm = 40,
    x_etr = 0.4,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    gmc_norm_column_name = 'gmc_norm',
    j_norm_column_name = 'J_norm',
    kp_column_name = 'Kp',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm',
    debug_mode = FALSE
  )
}

\arguments{
  \item{alpha_psii}{
    The fraction of photosystem II activity in the bundle sheath
    (\code{dimensionless}). If \code{alpha_psii} is not a number, then there
    must be a column in \code{rc_exdf} called \code{alpha_psii} with
    appropriate units. A numeric value supplied here will overwrite the values
    in the \code{alpha_psii} column of \code{rc_exdf} if it exists.
  }

  \item{gbs}{
    The bundle sheath conductance to CO2 in \code{mol m^(-2) s^(-1) bar^(-1)}.
    If \code{gbs} is not a number, then there must be a column in
    \code{rc_exdf} called \code{gbs} with appropriate units. A numeric value
    supplied here will overwrite the values in the \code{gbs} column of
    \code{rc_exdf} if it exists.
  }

  \item{gmc_at_25}{
    The mesophyll conductance to CO2 diffusion at 25 degrees C, expressed in
    \code{mol m^(-2) s^(-1) bar^(-1)}. If \code{gmc_at_25} is not a number, then
    there must be a column in \code{rc_exdf} called \code{gmc_at_25} with
    appropriate units. A numeric value supplied here will overwrite the values
    in the \code{gmc_at_25} column of \code{rc_exdf} if it exists.
  }

  \item{Rm_frac}{
    The fraction of the total mitochondrial respiration that occurs in the
    mesophyll. If \code{Rm_frac} is not a number, then there must be a column in
    \code{rc_exdf} called \code{Rm_frac} with appropriate units. A numeric
    value supplied here will overwrite the values in the \code{Rm_frac} column
    of \code{rc_exdf} if it exists.
  }

  \item{pcm_threshold_rlm}{
    An upper cutoff value for the partial pressure of CO2 in the mesophyll (in
    \code{microbar}) to be used when estimating \code{RLm}.
  }

  \item{x_etr}{
    The fraction of whole-chain electron transport occurring in the mesophyll
    (dimensionless). See Equation 29 from S. von Caemmerer (2021).
  }

  \item{a_column_name}{
    The name of the column in \code{rc_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{rc_exdf} that contains the
    intercellular CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{gmc_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    mesophyll conductance values (with units of
    \code{normalized to gmc at 25 degrees C}).
  }

  \item{j_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{J} values (with units of
    \code{normalized to J at 25 degrees C}).
  }

  \item{kp_column_name}{
    The name of the column in \code{rc_exdf} that contains the Michaelis-Menten
    constant for PEP carboxylase carboxylation in \code{microbar}.
  }

  \item{rl_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{RL} values (with units of \code{normalized to RL at 25 degrees C}).
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{rc_exdf} that contains the total pressure
    in \code{bar}.
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{vpmax_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{Vpmax} values (with units of
    \code{normalized to Vpmax at 25 degrees C}).
  }

  \item{debug_mode}{
    A logical (\code{TRUE} or \code{FALSE}) variable indicating whether to
    operate in debug mode. In debug mode, information about the linear fit used
    to estimate \code{RL} is printed; this can be helpful when troubleshooting
    issues with a particular curve.
  }
}

\details{
  Here we estimate values of \code{J_at_25}, \code{RL_at_25},
  \code{Vcmax_at_25}, \code{Vpmax_at_25}, and \code{Vpr} from a measured C4 CO2
  response curve. It is difficult to estimate values of \code{alpha_psii},
  \code{gbs}, \code{gmc_at_25}, and \code{Rm_frac} from a curve, so they must be
  supplied beforehand. For more information about these parameters, see the
  documentation for \code{\link{calculate_c4_assimilation}}. To estimate these
  parameter values, we use several equations from S. von Caemmerer, "Biochemical
  Models of Leaf Photosynthesis" (CSIRO Publishing, 2000)
  [\doi{10.1071/9780643103405}]. Any equation numbers referenced below are from
  this book.

  \itemize{
    \item \strong{Estimating RL}: An estimate for \code{RLm} can be obtained
          using Equation 4.26, which applies for low values of \code{PCm}. In
          this situation, \code{PCm + Kp} can be approximated by \code{Kp}, and
          Equation 4.26 simplifies to a linear relationship between the net
          assimilation \code{An} and \code{PCm}: \code{An =
          (gbs + Vpmax / kP) * PCm - RLm}. So, to estimate \code{RLm}, we make a
          linear fit of \code{An} vs. \code{PCm} in the low \code{PCm} range
          (\code{PCm <= pcm_threshold_rlm}) where this equation is expected to
          be valid. Then \code{RLm} is given by the negative of the intercept
          from the fit. In the C4 assimilation model, we assume that
          \code{RLm = Rm_frac * RL}, so we can also estimate
          \code{RL = RLm / Rm_frac} from this value.

          If there are fewer than two points with
          \code{PCm <= pcm_threshold_rlm}, the fit cannot be made, and we use a
          typical value instead (0.5 \code{micromol m^(-2) s^(-1)}). Likewise,
          if the linear fit predicts a negative or \code{NA} value for
          \code{RLm}, we use the same typical value instead.
    \item \strong{Estimating Vpmax}: An estimate for \code{Vpmax} can also be
          obtained from Equation 4.26. In this case, we simply solve the
          equation for \code{Vpmax} and use it to calculate a value of
          \code{Vpmax} at each point in the curve from the measured values of
          \code{An} and \code{PCm}, the input value of \code{gbs}, and the value
          of \code{RLm} estimated above. In the PEP-carboxylation-limited range,
          the estimated values of \code{Vpmax} should be reasonable. In other
          parts of the curve, the assimilation rate is limited by other factors,
          so \code{An} will be smaller than the PEP-carboxylation-limited
          values, causing the estimated values of \code{Vpmax} to be smaller.
          So, to make an overall estimate, we choose the largest estimated
          \code{Vpmax} value.
    \item \strong{Estimating Vcmax}: An estimate for \code{Vcmax} can be
          obtained by solving \code{An = Vcmax - RL} for \code{Vcmax}, similar
          to the method used to estimate \code{Vpmax}.
    \item \strong{Estimating Vpr}: An estimate for \code{Vpr} can be
          obtained by solving \code{An = Vpr + gbs * PCm - RLm} for \code{Vpr},
          similar to the method used to estimate \code{Vpmax}.
    \item \strong{Estimating J}: First, an estimate for \code{J} can be
          obtained by solving \code{An = (1 - x_etr) * J / 3 - RL} for \code{J}.
          Then, estimates of \code{J} can be made from \code{J} and
          \code{Qin}. The largest value of \code{J / J_norm} is chosen as
          the best estimate for \code{J_at_25}.
  }

  Note that a key assumption underlying this approach is that the net
  assimilation can be reasonably approximated by
  \code{An = min(Apc, Apr, Ar, Ajm)} (Equations 4.19, 4.25, 4.45, and 4.47
  combined). While this approximation seems to work well for low values of
  \code{PCm}, it tends to deviate significantly from the more accurate version
  at higher values of \code{PCm}, predicting values that are noticably smaller.
  Thus, the values of \code{Vcmax} and \code{Vpr} estimated using this procedure
  are unlikely to be accurate. This is not a problem; instead it simply
  highlights the importance of improving this initial guess using an optimizer,
  which can be accomplished via \code{\link{fit_c4_aci}}.
}

\value{
  A function with one input argument \code{rc_exdf}, which should be an
  \code{exdf} object representing one C4 CO2 response curve. The return value of
  this function will be a numeric vector with eight elements, representing the
  values of \code{alpha_psii}, \code{gbs}, \code{J_at_25}, \code{RL_at_25},
  \code{rm_frac}, \code{Vcmax_at_25}, \code{Vpmax_at_25}, and \code{Vpr} (in
  that order).
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

# Create the guessing function, using typical values for the alpha_psii, gbs,
# gmc_at_25, and Rm_frac: 0, 0.003, 1, and 0.5
guessing_func <- initial_guess_c4_aci(0, 0.003, 1, 0.5)

# Apply it and see the initial guesses for each curve
print(by(licor_file, licor_file[, 'species_plot'], guessing_func))

# A simple way to visualize the guesses is to "fit" the curves using the null
# optimizer, which simply returns the initial guess
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c4_aci,
  optim_fun = optimizer_null()
))

plot_c4_aci_fit(aci_results, 'species_plot', 'Ci', ylim = c(-10, 100))
}

\concept{exdf}
