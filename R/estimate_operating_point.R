estimate_operating_point <- function(
    aci_exdf,
    Ca_atmospheric,
    type = 'c3',
    a_column_name = 'A',
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    pcm_column_name = 'PCm',
    return_list = FALSE
)
{
    if (!is.exdf(aci_exdf)) {
        stop('estimate_operating_point requires an exdf object')
    }

    type <- tolower(type)

    if (!type %in% c('c3', 'c4')) {
        stop('`type` must be "c3" or "c4"')
    }

    # Only use points designated for fitting
    aci_exdf <- aci_exdf[points_for_fitting(aci_exdf), , TRUE]

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]  <- unit_dictionary('A')
    required_variables[[ci_column_name]] <- unit_dictionary('Ci')

    if (type == 'c3') {
        required_variables[[cc_column_name]] <- unit_dictionary('Cc')
    } else {
        required_variables[[pcm_column_name]] <- unit_dictionary('PCm')
    }

    check_required_variables(aci_exdf, required_variables)

    # Allow Ca to be NA
    required_variables <- list()
    required_variables[[ca_column_name]] <- unit_dictionary('Ca')
    check_required_variables(aci_exdf, required_variables, check_NA = FALSE)

    # Check to see if we should bypass the calculations
    bypass <- is.na(Ca_atmospheric)

    # Make sure the atmospheric Ca is included in the Ca range
    ca_vals   <- aci_exdf[, ca_column_name]
    all_ca_NA <- all(is.na(ca_vals))

    min_ca <- if (all_ca_NA) {NA} else {min(ca_vals, na.rm = TRUE)}
    max_ca <- if (all_ca_NA) {NA} else {max(ca_vals, na.rm = TRUE)}

    unreliable <- !bypass &&
        (all_ca_NA || Ca_atmospheric < min_ca || Ca_atmospheric > max_ca)

    msg <- if (unreliable) {
        if (all_ca_NA) {
            paste('All values of', ca_column_name, 'are NA')
            paste0(
                'All values of the atmospheric CO2 concentration column (',
                ca_column_name, ') are NA'
            )
        } else {
            paste0(
                'The atmospheric CO2 concentration (', Ca_atmospheric,
                ') is outside the measured range (', min_ca, ' to ', max_ca, ')'
            )
        }
    } else {
        ''
    }

    # Use linear interpolation to estimate the operating point
    operating_An <- if (!bypass && !unreliable) {
        stats::approx(
            aci_exdf[, ca_column_name],
            aci_exdf[, a_column_name],
            Ca_atmospheric,
            ties = list('ordered', mean)
        )[['y']]
    } else {
        NA
    }

    operating_Ci <- if (!bypass && !unreliable) {
        stats::approx(
            aci_exdf[, ca_column_name],
            aci_exdf[, ci_column_name],
            Ca_atmospheric,
            ties = list('ordered', mean)
        )[['y']]
    } else {
        NA
    }

    operating_Cc <- if (!bypass && !unreliable && type == 'c3') {
        stats::approx(
            aci_exdf[, ca_column_name],
            aci_exdf[, cc_column_name],
            Ca_atmospheric,
            ties = list('ordered', mean)
        )[['y']]
    } else {
        NA
    }

    operating_PCm <- if (!bypass && !unreliable && type == 'c4') {
        stats::approx(
            aci_exdf[, ca_column_name],
            aci_exdf[, pcm_column_name],
            Ca_atmospheric,
            ties = list('ordered', mean)
        )[['y']]
    } else {
        NA
    }

    if (return_list) {
        # Prepare an exdf that can be used for calculating the operating An
        # using fit parameters
        if (type == 'c3') {
            cc_seq <- aci_exdf[, cc_column_name]

            operating_row <- if (!bypass && !unreliable) {
                which(abs(cc_seq - operating_Cc) == min(abs(cc_seq - operating_Cc)))
            } else {
                1
            }

            operating_exdf <- aci_exdf[operating_row, , TRUE]

            operating_exdf[, cc_column_name] <- operating_Cc
            operating_exdf[, ci_column_name] <- operating_Ci

            list(
                Ca_atmospheric = Ca_atmospheric,
                operating_An = operating_An,
                operating_Cc = operating_Cc,
                operating_Ci = operating_Ci,
                operating_exdf = operating_exdf,
                operating_point_msg = msg
            )
        } else {
            pcm_seq <- aci_exdf[, pcm_column_name]

            operating_row <- if (!bypass && !unreliable) {
                which(abs(pcm_seq - operating_PCm) == min(abs(pcm_seq - operating_PCm)))
            } else {
                1
            }

            operating_exdf <- aci_exdf[operating_row, , TRUE]

            operating_exdf[, pcm_column_name] <- operating_PCm
            operating_exdf[, ci_column_name] <- operating_Ci

            list(
                Ca_atmospheric = Ca_atmospheric,
                operating_An = operating_An,
                operating_PCm = operating_PCm,
                operating_Ci = operating_Ci,
                operating_exdf = operating_exdf,
                operating_point_msg = msg
            )
        }
    } else {
        # Get the replicate identifier columns
        aci_identifiers <- identifier_columns(aci_exdf)

        # Store the results
        aci_identifiers[, 'Ca_atmospheric'] <- Ca_atmospheric
        aci_identifiers[, 'operating_An']   <- operating_An
        aci_identifiers[, 'operating_Ci']   <- operating_Ci

        if (type == 'c3') {
            aci_identifiers[, 'operating_Cc'] <- operating_Cc
        } else {
            aci_identifiers[, 'operating_PCm'] <- operating_PCm
        }

        aci_identifiers[, 'operating_point_msg'] <- msg

        # Document the new columns that were added and return the exdf
        aci_identifiers <- document_variables(
            aci_identifiers,
            c('estimate_operating_point', 'Ca_atmospheric',      'micromol mol^(-1)'),
            c('estimate_operating_point', 'operating_An',        aci_exdf$units[[a_column_name]]),
            c('estimate_operating_point', 'operating_Ci',        aci_exdf$units[[ci_column_name]]),
            c('estimate_operating_point', 'operating_point_msg', '')
        )

        if (type == 'c3') {
            document_variables(
                aci_identifiers,
                c('estimate_operating_point', 'operating_Cc', aci_exdf$units[[cc_column_name]])
            )
        } else {
            document_variables(
                aci_identifiers,
                c('estimate_operating_point', 'operating_PCm', aci_exdf$units[[pcm_column_name]])
            )
        }
    }
}
