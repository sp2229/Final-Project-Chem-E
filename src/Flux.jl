"""
    calculate_optimal_flux_distribution(S,[Lv,Uv],[Lx,Ux],c; min_flag=true)
    Computes the optimal metabolic flux distribution given the constraints.
    Inputs:
    `S` - stoichiometric_matrix (M x R)
    `[Lv,Uv]` - (R x 2) array of the flux lower (L) and upper (U) bounds
    `[Lx,Ux]` - (M x 2) array of species lower (L) and upper (U) bounds (if at steady-state, L = U = 0)
    `c` - R x 1 vector holding indexes for objective vector
    Outputs:
    `objective_value` - value of the objective function at the optimum
    `calculated_flux_array` - R x 1 flux array at the optimum
    `dual_value_array` - R x 1 dual values
    `uptake_array` - M x 1 array of S*v
    `exit_flag` = 0 if optimal
    `status_flag` = 5 if optimal
"""
function _calculate_optimal_flux_distribution(stoichiometric_matrix::Array{Float64,2}, default_bounds_array::Array{Float64,2},
    species_bounds_array::Array{Float64,2}, objective_coefficient_array::Array{Float64,1}; min_flag::Bool = true)

    try

        # set some constants -
        TIME_RESTART_LIM = 60

        # Get the stoichiometric_matrix from data_dictionary -
        (number_of_species, number_of_fluxes) = size(stoichiometric_matrix)

        # # Setup the GLPK problem -
        lp_problem = GLPK.glp_create_prob()
        GLPK.glp_set_prob_name(lp_problem, "sample")
        GLPK.glp_set_obj_name(lp_problem, "objective")

        # Set solver parameters
        solver_parameters = GLPK.glp_smcp()
        GLPK.glp_init_smcp(solver_parameters)
        solver_parameters.msg_lev = GLPK.GLP_MSG_OFF

        # Are we doing min -or- max?
        if min_flag == true
            GLPK.glp_set_obj_dir(lp_problem, GLPK.GLP_MIN)
        else
            GLPK.glp_set_obj_dir(lp_problem, GLPK.GLP_MAX)
        end

        # Set the number of constraints and fluxes -
        GLPK.glp_add_rows(lp_problem, number_of_species)
        GLPK.glp_add_cols(lp_problem, number_of_fluxes)

        # Setup flux bounds, and objective function -
        (number_of_fluxes, number_of_bounds) = size(default_bounds_array)
        for flux_index = 1:number_of_fluxes

            flux_lower_bound = default_bounds_array[flux_index, 1]
            flux_upper_bound = default_bounds_array[flux_index, 2]

            # Check bounds type ... default is DB -
            if (flux_upper_bound == flux_lower_bound)
                flux_constraint_type = GLPK.GLP_FX
            else
                flux_constraint_type = GLPK.GLP_DB
            end

            # flux symbol? (later use name - for now, fake it)
            flux_symbol = "R_" * string(flux_index)

            # Set the bounds in GLPK -
            GLPK.glp_set_col_name(lp_problem, flux_index, flux_symbol)
            GLPK.glp_set_col_bnds(lp_problem, flux_index, flux_constraint_type, flux_lower_bound, flux_upper_bound)
        end

        # Setup objective function -
        for (flux_index, obj_coeff) in enumerate(objective_coefficient_array)

            # Set the objective function value in GLPK -
            GLPK.glp_set_obj_coef(lp_problem, flux_index, obj_coeff)
        end

        # Setup problem constraints for the metabolites -
        for species_index = 1:number_of_species

            species_lower_bound = species_bounds_array[species_index, 1]
            species_upper_bound = species_bounds_array[species_index, 2]

            # defualt
            species_constraint_type = GLPK.GLP_FX
            if (species_lower_bound != species_upper_bound)
                species_constraint_type = GLPK.GLP_DB
            end

            # set the symbol -
            species_symbol = "x_" * string(species_index)

            # Set the species bounds in GLPK -
            GLPK.glp_set_row_name(lp_problem, species_index, species_symbol)
            GLPK.glp_set_row_bnds(lp_problem, species_index, species_constraint_type, species_lower_bound, species_upper_bound)

        end

        # Setup the stoichiometric array -
        counter = 1
        row_index_array = zeros(Cint, number_of_species * number_of_fluxes)
        col_index_array = zeros(Cint, number_of_species * number_of_fluxes)
        species_index_vector = collect(1:number_of_species)
        flux_index_vector = collect(1:number_of_fluxes)
        flat_stoichiometric_array = zeros(Float64, number_of_species * number_of_fluxes)
        for species_index in species_index_vector
            for flux_index in flux_index_vector
                row_index_array[counter] = species_index
                col_index_array[counter] = flux_index
                flat_stoichiometric_array[counter] = stoichiometric_matrix[species_index, flux_index]
                counter += 1
            end
        end
        GLPK.glp_load_matrix(lp_problem, number_of_species * number_of_fluxes, GLPK.offset(row_index_array), GLPK.offset(col_index_array), GLPK.offset(flat_stoichiometric_array))

        # Call the solver -
        exit_flag = GLPK.glp_simplex(lp_problem, solver_parameters)

        # Get the objective function value -
        objective_value = GLPK.glp_get_obj_val(lp_problem)

        # Get the calculated flux values from GLPK -
        calculated_flux_array = zeros(Float64, number_of_fluxes)
        for flux_index in flux_index_vector
            calculated_flux_array[flux_index] = GLPK.glp_get_col_prim(lp_problem, flux_index)
        end

        # Get the dual values -
        dual_value_array = zeros(Float64, number_of_fluxes)
        for flux_index in flux_index_vector
            dual_value_array[flux_index] = GLPK.glp_get_col_dual(lp_problem, flux_index)
        end

        # is this solution optimal?
        status_flag = GLPK.glp_get_status(lp_problem)

        # Calculate the uptake array -
        uptake_array = stoichiometric_matrix * calculated_flux_array

        # return results_tuple -
        results_tuple = (objective_value = objective_value, calculated_flux_array = calculated_flux_array,
            dual_value_array = dual_value_array, uptake_array = uptake_array, exit_flag = exit_flag,
            status_flag = status_flag)

        # Formulate the return tuple -
        return results_tuple

    catch error

        # what is our error policy? => for now, just print the message
        error_message = sprint(showerror, error, catch_backtrace())
        println(error_message)
    end
end

function compute_optimal_extent(stoichiometric_matrix::Array{Float64,2}, default_bounds_array::Array{Float64,2},
    species_bounds_array::Array{Float64,2}, objective_coefficient_array::Array{Float64,1}; min_flag::Bool = true)

    # compute first pass -
    result = _calculate_optimal_flux_distribution(stoichiometric_matrix, default_bounds_array, species_bounds_array, 
        objective_coefficient_array; min_flag = min_flag);

    # return -
    return result;
end