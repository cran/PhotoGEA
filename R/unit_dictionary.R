# Specify standard units for some types of quantities
conductance       <- 'mol m^(-2) s^(-1)'
conductance_bar   <- 'mol m^(-2) s^(-1) bar^(-1)'
dimensionless     <- 'dimensionless'
micromol_flux     <- 'micromol m^(-2) s^(-1)'
micromol_fraction <- 'micromol mol^(-1)'
millimol_fraction <- 'mmol mol^(-1)'
pressure          <- 'bar'
temperature       <- 'degrees C'

# Helping function for "normalized" quantities
normalized_units <- function(qname) {
    paste('normalized to', qname, 'at 25 degrees C')
}

# Specify units for some important parameters
unit_dictionary_list <- list(
    A                   = micromol_flux,
    ao                  = dimensionless,
    Ainitial            = micromol_flux,
    alpha_g             = dimensionless,
    alpha_j_at_25       = dimensionless,
    alpha_j_norm        = normalized_units('alpha_j'),
    alpha_old           = dimensionless,
    alpha_psii          = dimensionless,
    alpha_s             = dimensionless,
    alpha_t             = dimensionless,
    bb_index            = conductance,
    c4_curvature        = dimensionless,
    c4_slope            = 'mol m^(-2) s^(-1)',
    Ca                  = micromol_fraction,
    Cc                  = micromol_fraction,
    Ci                  = micromol_fraction,
    CO2_r               = micromol_fraction,
    CO2_s               = micromol_fraction,
    CorrFact            = NA,
    Csurface            = micromol_fraction,
    ETR                 = micromol_flux,
    Flow                = 'micromol s^(-1)',
    gamma_star          = dimensionless,
    Gamma_star_at_25    = micromol_fraction,
    Gamma_star_norm     = normalized_units('Gamma_star'),
    gbs                 = conductance_bar,
    gmc_at_25           = conductance_bar,
    gmc_norm            = normalized_units('gmc'),
    gsc                 = conductance,
    gsw                 = conductance,
    H2O_r               = millimol_fraction,
    H2O_s               = millimol_fraction,
    I2                  = micromol_flux,
    I2_at_25            = micromol_flux,
    I2_tl               = micromol_flux,
    J                   = micromol_flux,
    J_at_25             = micromol_flux,
    J_at_25_lower       = micromol_flux,
    J_at_25_upper       = micromol_flux,
    J_norm              = normalized_units('J'),
    J_tl_avg            = micromol_flux,
    J_tl_avg_lower      = micromol_flux,
    J_tl_avg_upper      = micromol_flux,
    Jmax_at_25          = micromol_flux,
    Jmax_norm           = normalized_units('Jmax'),
    Jmax_tl             = micromol_flux,
    Kc_at_25            = micromol_fraction,
    Kc_norm             = normalized_units('Kc'),
    Ko_at_25            = millimol_fraction,
    Ko_norm             = normalized_units('Ko'),
    oxygen              = 'percent',
    PCm                 = 'microbar',
    PhiPS2              = dimensionless,
    Qin                 = micromol_flux,
    rL                  = micromol_flux,
    RL_at_25            = micromol_flux,
    RL_norm             = normalized_units('RL'),
    Rm_frac             = dimensionless,
    rubisco_specificity = 'M / M',
    S                   = 'cm^2',
    tau                 = dimensionless,
    theta_j_at_25       = dimensionless,
    theta_j_norm        = normalized_units('theta_j'),
    Tleaf_avg           = temperature,
    TleafCnd            = temperature,
    total_pressure      = pressure,
    Tp_at_25            = micromol_flux,
    Tp_norm             = normalized_units('Tp'),
    Vcmax_at_25         = micromol_flux,
    Vcmax_norm          = normalized_units('Vcmax'),
    Vmax                = micromol_flux,
    VPDleaf             = 'kPa',
    Vpmax_at_25         = micromol_flux,
    Vpmax_norm          = normalized_units('Vpmax'),
    Vpr                 = micromol_flux
)

unit_dictionary <- function(quantity_name) {
    if (!quantity_name %in% names(unit_dictionary_list)) {
        msg <- paste0(
            'Units were requested for a quantity named `', quantity_name,
            '`, but it is not included in the unit dictionary'
        )

        stop(msg)
    }

    unit_dictionary_list[[quantity_name]]
}
