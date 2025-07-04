# Specify default fit settings
c4_aci_lower       <- list(alpha_psii = -1,  gbs = -1,    gmc_at_25 = -1, J_at_25 = -50,   RL_at_25 = -10,   Rm_frac = -10, Vcmax_at_25 = -50,   Vpmax_at_25 = -50,   Vpr = -50)
c4_aci_upper       <- list(alpha_psii = 10,  gbs = 10,    gmc_at_25 = 10, J_at_25 = 1000,  RL_at_25 = 100,   Rm_frac = 10,  Vcmax_at_25 = 1000,  Vpmax_at_25 = 1000,  Vpr = 1000)
c4_aci_fit_options <- list(alpha_psii = 0,   gbs = 0.003, gmc_at_25 = 1,  J_at_25 = 1000,  RL_at_25 = 'fit', Rm_frac = 0.5, Vcmax_at_25 = 'fit', Vpmax_at_25 = 'fit', Vpr = 1000)

c4_aci_param <- c('alpha_psii', 'gbs', 'gmc_at_25', 'J_at_25', 'RL_at_25', 'Rm_frac', 'Vcmax_at_25', 'Vpmax_at_25', 'Vpr')

# Fitting function
fit_c4_aci <- function(
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
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_c4_aci requires an exdf object')
    }

    if (sd_A != 'RMSE') {
        stop('At this time, the only supported option for sd_A is `RMSE`')
    }

    # Get the replicate identifier columns
    replicate_identifiers <- identifier_columns(replicate_exdf)

    if (debug_mode) {
        debug_msg('fit_c4_aci curve identifiers:')
        cat('\n')
        utils::str(replicate_identifiers$main_data)
    }

    # Define the total error function; units will also be checked by this
    # function
    total_error_fcn <- error_function_c4_aci(
        replicate_exdf,
        fit_options,
        1, # sd_A
        x_etr,
        a_column_name,
        ao_column_name,
        ci_column_name,
        gamma_star_column_name,
        gmc_norm_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        kp_column_name,
        oxygen_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        vpmax_norm_column_name,
        hard_constraints,
        debug_mode
    )

    # Make sure the required variables are defined and have the correct units;
    # most units have already been checked by error_function_c4_aci
    required_variables <- list()
    required_variables[[ca_column_name]] <- unit_dictionary('Ca')

    check_required_variables(replicate_exdf, required_variables, check_NA = FALSE)

    # Assemble lower, upper, and fit_options
    luf <- assemble_luf(
        c4_aci_param,
        c4_aci_lower, c4_aci_upper, c4_aci_fit_options,
        lower, upper, fit_options
    )

    lower_complete <- luf$lower
    upper_complete <- luf$upper
    fit_options <- luf$fit_options
    fit_options_vec <- luf$fit_options_vec
    param_to_fit <- luf$param_to_fit

    # Get an initial guess for all the parameter values
    alpha_psii_guess <- if (fit_options$alpha_psii == 'fit') {0.1}   else {fit_options$alpha_psii}
    gbs_guess        <- if (fit_options$gbs == 'fit')        {0.003} else {fit_options$gbs}
    gmc_guess        <- if (fit_options$gmc_at_25 == 'fit')  {1.0}   else {fit_options$gmc_at_25}
    rm_frac_guess    <- if (fit_options$Rm_frac == 'fit')    {0.5}   else {fit_options$Rm_frac}

    initial_guess_fun <- initial_guess_c4_aci(
        alpha_psii_guess,
        gbs_guess,
        gmc_guess,
        rm_frac_guess,
        40, # pcm_threshold_rlm
        x_etr,
        a_column_name,
        ci_column_name,
        gmc_norm_column_name,
        j_norm_column_name,
        kp_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        vpmax_norm_column_name,
        debug_mode
    )

    initial_guess <- initial_guess_fun(replicate_exdf)

    if (debug_mode) {
        debug_msg(
            'fit_c4_aci initial_guess:',
            paste(initial_guess, collapse = ', ')
        )
    }

    # Find the best values for the parameters that should be varied
    optim_result <- optim_fun(
        initial_guess[param_to_fit],
        total_error_fcn,
        lower = lower_complete[param_to_fit],
        upper = upper_complete[param_to_fit]
    )

    check_optim_result(optim_result)

    # Get the values of all parameters following the optimization
    best_X <- fit_options_vec
    best_X[param_to_fit] <- optim_result[['par']]

    if (debug_mode) {
        debug_msg(
            'fit_c4_aci best_X:',
            paste(best_X, collapse = ', ')
        )
    }

    # Get the corresponding values of Cc at the best guess
    replicate_exdf <- apply_gm(
        replicate_exdf,
        best_X[3], # gmc_at_25
        'C4',
        TRUE,
        a_column_name,
        ca_column_name,
        ci_column_name,
        gmc_norm_column_name,
        total_pressure_column_name,
        perform_checks = FALSE
    )

    pcm_column_name <- 'PCm'

    # Get the corresponding values of An at the best guess
    aci <- calculate_c4_assimilation(
        replicate_exdf,
        best_X[1], # alpha_psii
        best_X[2], # gbs
        best_X[4], # J_at_25
        best_X[5], # RL_at_25
        best_X[6], # Rm_frac
        best_X[7], # Vcmax_at_25
        best_X[8], # Vpmax_at_25
        best_X[9], # Vpr
        x_etr,
        ao_column_name,
        gamma_star_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        kp_column_name,
        oxygen_column_name,
        pcm_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        vpmax_norm_column_name,
        hard_constraints,
        perform_checks = FALSE
    )

    # Remove any columns in replicate_exdf that are also included in the
    # output from calculate_c4_assimilation
    replicate_exdf <- remove_repeated_colnames(replicate_exdf, aci)

    # Set all categories to `fit_c4_aci` and rename the `An` variable to
    # indicate that it contains fitted values of `a_column_name`
    aci$categories[1,] <- 'fit_c4_aci'
    colnames(aci)[colnames(aci) == 'An'] <- paste0(a_column_name, '_fit')

    # Append the fitting results to the original exdf object
    replicate_exdf <- cbind(replicate_exdf, aci)

    # Get operating point information
    operating_point_info <- estimate_operating_point(
        replicate_exdf,
        Ca_atmospheric,
        type = 'c4',
        a_column_name,
        ca_column_name,
        cc_column_name = NULL,
        ci_column_name,
        pcm_column_name,
        return_list = TRUE
    )

    # Estimate An at the operating point
    operating_An_model <- calculate_c4_assimilation(
        operating_point_info$operating_exdf,
        best_X[1], # alpha_psii
        best_X[2], # gbs
        best_X[4], # J_at_25
        best_X[5], # RL_at_25
        best_X[6], # Rm_frac
        best_X[7], # Vcmax_at_25
        best_X[8], # Vpmax_at_25
        best_X[9], # Vpr
        x_etr,
        ao_column_name,
        gamma_star_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        kp_column_name,
        oxygen_column_name,
        pcm_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        vpmax_norm_column_name,
        hard_constraints,
        perform_checks = FALSE
    )[, 'An']

    # Interpolate onto a finer Cc spacing and recalculate fitted rates
    replicate_exdf_interpolated <- interpolate_assimilation_inputs(
        replicate_exdf,
        c(
            'alpha_psii',
            'gbs',
            'J_at_25',
            'RL_at_25',
            'Rm_frac',
            'Vcmax_at_25',
            'Vpmax_at_25',
            'Vpr',
            ao_column_name,
            ci_column_name,
            gamma_star_column_name,
            j_norm_column_name,
            kc_column_name,
            ko_column_name,
            kp_column_name,
            oxygen_column_name,
            pcm_column_name,
            rl_norm_column_name,
            total_pressure_column_name,
            vcmax_norm_column_name,
            vpmax_norm_column_name
        ),
        ci_column_name,
        c_step = 1
    )

    assim_interpolated <- calculate_c4_assimilation(
        replicate_exdf_interpolated,
        '', # alpha_psii
        '', # gbs
        '', # J_at_25
        '', # RL_at_25
        '', # Rm_frac
        '', # Vcmax_at_25
        '', # Vpmax_at_25
        '', # Vpr
        x_etr,
        ao_column_name,
        gamma_star_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        kp_column_name,
        oxygen_column_name,
        pcm_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        vpmax_norm_column_name,
        hard_constraints,
        perform_checks = FALSE
    )

    # Remove any columns in replicate_exdf_interpolated that are also included
    # in the output from calculate_c4_assimilation
    replicate_exdf_interpolated <- remove_repeated_colnames(
        replicate_exdf_interpolated,
        assim_interpolated
    )

    fits_interpolated <- cbind(
        replicate_exdf_interpolated,
        assim_interpolated
    )

    # If there was a problem, set all the fit results to NA
    fit_failure <- aci[1, 'c4_assimilation_msg'] != ''

    if (fit_failure) {
        best_X[param_to_fit] <- NA
        operating_An_model <- NA

        for (cn in colnames(aci)) {
            if (cn != 'c4_assimilation_msg') {
                replicate_exdf[, cn] <- NA
            }
        }

        for (cn in colnames(assim_interpolated)) {
            if (cn != 'c4_assimilation_msg') {
                fits_interpolated[, cn] <- NA
            }
        }
    }

    if (debug_mode) {
        debug_msg(
            'fit_c4_aci outcome:',
            if (fit_failure) {'failure'} else {'success'}
        )
    }

    # Include the atmospheric CO2 concentration
    replicate_exdf[, 'Ca_atmospheric'] <- Ca_atmospheric

    # Document the new columns that were added
    replicate_exdf <- document_variables(
        replicate_exdf,
        c('fit_c4_aci', 'Ca_atmospheric', 'micromol mol^(-1)')
    )

    # Add a column for the residuals
    replicate_exdf <- calculate_residuals(replicate_exdf, a_column_name)

    # Attach identifiers to interpolated rates, making sure to avoid duplicating
    # any columns
    identifiers_to_keep <-
        colnames(replicate_identifiers)[!colnames(replicate_identifiers) %in% colnames(fits_interpolated)]

    fits_interpolated <- cbind(
        fits_interpolated,
        replicate_identifiers[, identifiers_to_keep, TRUE]
    )

    # Attach the residual stats to the identifiers
    replicate_identifiers <- cbind(
        replicate_identifiers,
        residual_stats(
            replicate_exdf[, paste0(a_column_name, '_residuals')],
            replicate_exdf$units[[a_column_name]],
            length(which(param_to_fit))
        )
    )

    # Attach the best-fit parameters to the identifiers
    replicate_identifiers[, 'alpha_psii']  <- best_X[1]
    replicate_identifiers[, 'gbs']         <- best_X[2]
    replicate_identifiers[, 'gmc_at_25']   <- best_X[3]
    replicate_identifiers[, 'J_at_25']     <- best_X[4]
    replicate_identifiers[, 'RL_at_25']    <- best_X[5]
    replicate_identifiers[, 'Rm_frac']     <- best_X[6]
    replicate_identifiers[, 'Vcmax_at_25'] <- best_X[7]
    replicate_identifiers[, 'Vpmax_at_25'] <- best_X[8]
    replicate_identifiers[, 'Vpr']         <- best_X[9]

    # Attach the average leaf-temperature values of fitting parameters
    replicate_identifiers[, 'gmc_tl_avg']   <- mean(replicate_exdf[, 'gmc_tl'])
    replicate_identifiers[, 'J_tl_avg']     <- mean(replicate_exdf[, 'J_tl'])
    replicate_identifiers[, 'RL_tl_avg']    <- mean(replicate_exdf[, 'RL_tl'])
    replicate_identifiers[, 'Vcmax_tl_avg'] <- mean(replicate_exdf[, 'Vcmax_tl'])
    replicate_identifiers[, 'Vpmax_tl_avg'] <- mean(replicate_exdf[, 'Vpmax_tl'])

    # Attach fitting details
    replicate_identifiers[, 'convergence']         <- optim_result[['convergence']]
    replicate_identifiers[, 'convergence_msg']     <- optim_result[['convergence_msg']]
    replicate_identifiers[, 'feval']               <- optim_result[['feval']]
    replicate_identifiers[, 'optimizer']           <- optim_result[['optimizer']]
    replicate_identifiers[, 'c4_assimilation_msg'] <- replicate_exdf[1, 'c4_assimilation_msg']

    # Attach operating point information
    replicate_identifiers[, 'operating_Ci']        <- operating_point_info$operating_Ci
    replicate_identifiers[, 'operating_PCm']       <- operating_point_info$operating_PCm
    replicate_identifiers[, 'operating_An']        <- operating_point_info$operating_An
    replicate_identifiers[, 'operating_An_model']  <- operating_An_model
    replicate_identifiers[, 'operating_point_msg'] <- operating_point_info$operating_point_msg

    # Get an updated likelihood value using the RMSE
    replicate_identifiers[, 'optimum_val'] <- if (fit_failure) {
        NA
    } else {
        error_function_c4_aci(
            replicate_exdf,
            fit_options,
            replicate_identifiers[, 'RMSE'], # sd_A
            x_etr,
            a_column_name,
            ao_column_name,
            ci_column_name,
            gamma_star_column_name,
            gmc_norm_column_name,
            j_norm_column_name,
            kc_column_name,
            ko_column_name,
            kp_column_name,
            oxygen_column_name,
            rl_norm_column_name,
            total_pressure_column_name,
            vcmax_norm_column_name,
            vpmax_norm_column_name,
            hard_constraints
        )(best_X[param_to_fit])
    }

    # Document the new columns that were added
    replicate_identifiers <- document_variables(
        replicate_identifiers,
        c('fit_c4_aci',               'alpha_psii',          unit_dictionary('alpha_psii')),
        c('fit_c4_aci',               'gbs',                 unit_dictionary('gbs')),
        c('fit_c4_aci',               'gmc_at_25',           unit_dictionary('gmc_at_25')),
        c('fit_c4_aci',               'J_at_25',             'micromol m^(-2) s^(-1)'),
        c('fit_c4_aci',               'RL_at_25',            'micromol m^(-2) s^(-1)'),
        c('fit_c4_aci',               'Rm_frac',             unit_dictionary('Rm_frac')),
        c('fit_c4_aci',               'Vcmax_at_25',         'micromol m^(-2) s^(-1)'),
        c('fit_c4_aci',               'Vpmax_at_25',         'micromol m^(-2) s^(-1)'),
        c('fit_c4_aci',               'Vpr',                 'micromol m^(-2) s^(-1)'),
        c('fit_c4_aci',               'gmc_tl_avg',          unit_dictionary('gmc_at_25')),
        c('fit_c4_aci',               'J_tl_avg',            'micromol m^(-2) s^(-1)'),
        c('fit_c4_aci',               'RL_tl_avg',           'micromol m^(-2) s^(-1)'),
        c('fit_c4_aci',               'Vcmax_tl_avg',        'micromol m^(-2) s^(-1)'),
        c('fit_c4_aci',               'Vpmax_tl_avg',        'micromol m^(-2) s^(-1)'),
        c('estimate_operating_point', 'operating_Ci',        replicate_exdf$units[[ci_column_name]]),
        c('estimate_operating_point', 'operating_PCm',       replicate_exdf$units[[pcm_column_name]]),
        c('estimate_operating_point', 'operating_An',        replicate_exdf$units[[a_column_name]]),
        c('estimate_operating_point', 'operating_point_msg', ''),
        c('fit_c4_aci',               'operating_An_model',  replicate_exdf$units[[a_column_name]]),
        c('fit_c4_aci',               'convergence',         ''),
        c('fit_c4_aci',               'convergence_msg',     ''),
        c('fit_c4_aci',               'feval',               ''),
        c('fit_c4_aci',               'optimum_val',         ''),
        c('fit_c4_aci',               'c4_assimilation_msg', '')
    )

    # Calculate confidence intervals, if necessary
    if (calculate_confidence_intervals) {
        replicate_identifiers <- confidence_intervals_c4_aci(
            replicate_exdf,
            replicate_identifiers,
            lower,
            upper,
            fit_options,
            if (fit_failure) {0} else {replicate_identifiers[, 'RMSE']}, # sd_A
            relative_likelihood_threshold,
            x_etr,
            a_column_name,
            ao_column_name,
            ci_column_name,
            gamma_star_column_name,
            gmc_norm_column_name,
            j_norm_column_name,
            kc_column_name,
            ko_column_name,
            kp_column_name,
            oxygen_column_name,
            rl_norm_column_name,
            total_pressure_column_name,
            vcmax_norm_column_name,
            vpmax_norm_column_name,
            hard_constraints
        )

        # Attach limits for the average leaf-temperature values of fitting
        # parameters
        replicate_identifiers <- confidence_intervals_leaf_temperature(
            replicate_identifiers,
            c('gmc', 'J', 'RL', 'Vcmax', 'Vpmax'),
            'fit_c4_aci'
        )
    }

    # Return the results, including indicators of unreliable parameter estimates
    identify_c4_unreliable_points(
        replicate_identifiers,
        replicate_exdf,
        fits_interpolated,
        remove_unreliable_param
    )
}
