using GLPK
include("DataDictionary.jl")

    # Inputs:
    # `S` - stoichiometric_matrix (M x R)
    # `[Lv,Uv]` - (R x 2) array of the flux lower (L) and upper (U) bounds
    # `[Lx,Ux]` - (M x 2) array of species lower (L) and upper (U) bounds (if at steady-state, L = U = 0)
    # `c` - R x 1 vector holding indexes for objective vector
    # Outputs:
    # `objective_value` - value of the objective function at the optimum
    # `calculated_flux_array` - R x 1 flux array at the optimum
    # `dual_value_array` - R x 1 dual values
    # `uptake_array` - M x 1 array of S*v
    # `exit_flag` = 0 if optimal
    # `status_flag` = 5 if optimal

function calculate_optimal_flux_distribution(data_dictionary)

    stoichiometric_matrix = data_dictionary["stoichiometric_matrix"];
    metabolic_flux_bounds_array = data_dictionary["metabolic_flux_bounds_array"];
    species_bounds_array = data_dictionary["species_bounds_array"];
    objective_coefficient_array = data_dictionary["objective_coefficient_array"];
    min_flag = data_dictionary["min_flag"];

    # set some constants -
    TIME_RESTART_LIM = 60
    # Get the stoichiometric_matrix from data_dictionary -
    (number_of_species,number_of_fluxes) = size(stoichiometric_matrix);
    # # Setup the GLPK problem -
    lp_problem = GLPK.Prob();
    GLPK.set_prob_name(lp_problem, "sample");
    GLPK.set_obj_name(lp_problem, "objective")
    # Set solver parameters
    solver_parameters = GLPK.SimplexParam();
    GLPK.init_smcp(solver_parameters);
    solver_parameters.msg_lev = GLPK.MSG_OFF;
    # Are we doing min -or- max?
    if min_flag == true
        GLPK.set_obj_dir(lp_problem, GLPK.MIN);
    else
        GLPK.set_obj_dir(lp_problem, GLPK.MAX);
    end
    # Set the number of constraints and fluxes -
    GLPK.add_rows(lp_problem, number_of_species);
    GLPK.add_cols(lp_problem, number_of_fluxes);
    # Setup flux bounds, and objective function -
    (number_of_fluxes,number_of_bounds) = size(metabolic_flux_bounds_array)
    for flux_index = 1:number_of_fluxes
        flux_lower_bound = metabolic_flux_bounds_array[flux_index,1]
        flux_upper_bound = metabolic_flux_bounds_array[flux_index,2]
        # Check bounds type ... default is DB -
        if (flux_upper_bound == flux_lower_bound)
            flux_constraint_type = GLPK.FX
        else
            flux_constraint_type = GLPK.DB
        end
        # flux symbol? (later use name - for now, fake it)
        flux_symbol = "R_"*string(flux_index)
        # Set the bounds in GLPK -
        GLPK.set_col_name(lp_problem, flux_index, flux_symbol);
        GLPK.set_col_bnds(lp_problem, flux_index, flux_constraint_type, flux_lower_bound, flux_upper_bound);
    end
    # Setup objective function -
    for (flux_index,obj_coeff) in enumerate(objective_coefficient_array)
        # Set the objective function value in GLPK -
        GLPK.set_obj_coef(lp_problem, flux_index, obj_coeff);
    end
    # Setup problem constraints for the metabolites -
    for species_index = 1:number_of_species
        species_lower_bound = species_bounds_array[species_index,1]
        species_upper_bound = species_bounds_array[species_index,2]
        # defualt
        species_constraint_type = GLPK.FX
        if (species_lower_bound != species_upper_bound)
            species_constraint_type = GLPK.DB
        end
        # set the symbol -
        species_symbol = "x_"*string(species_index)
        # Set the species bounds in GLPK -
        GLPK.set_row_name(lp_problem, species_index, species_symbol);
        GLPK.set_row_bnds(lp_problem, species_index, species_constraint_type, species_lower_bound, species_upper_bound);
    end
    # Setup the stoichiometric array -
    counter = 1;
    row_index_array = zeros(Int,number_of_species*number_of_fluxes);
    col_index_array = zeros(Int,number_of_species*number_of_fluxes);
    species_index_vector = collect(1:number_of_species);
    flux_index_vector = collect(1:number_of_fluxes);
    flat_stoichiometric_array = zeros(Float64,number_of_species*number_of_fluxes);
    for species_index in species_index_vector
        for flux_index in flux_index_vector
            row_index_array[counter] = species_index;
            col_index_array[counter] = flux_index;
            flat_stoichiometric_array[counter] = stoichiometric_matrix[species_index,flux_index];
            counter+=1;
        end
    end
    GLPK.load_matrix(lp_problem, number_of_species*number_of_fluxes, row_index_array, col_index_array, flat_stoichiometric_array);
    # Call the solver -
    exit_flag = GLPK.simplex(lp_problem, solver_parameters);
    # Get the objective function value -
    objective_value = GLPK.get_obj_val(lp_problem);
    # Get the calculated flux values from GLPK -
    calculated_flux_array = zeros(Float64,number_of_fluxes);
    for flux_index in flux_index_vector
        calculated_flux_array[flux_index] = GLPK.get_col_prim(lp_problem, flux_index);
    end
    # Get the dual values -
    dual_value_array = zeros(Float64,number_of_fluxes);
    for flux_index in flux_index_vector
        dual_value_array[flux_index] = GLPK.get_col_dual(lp_problem, flux_index);
    end
    # is this solution optimal?
    status_flag = GLPK.get_status(lp_problem)
    # Calculate the uptake array -
    uptake_array = stoichiometric_matrix*calculated_flux_array;
    # Formulate the return tuple -
    return (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag);

    # --- SOLVE THE FVA PROBLEM -------------------------------------------- #
    target_value = 0.0
    if (GLPK.get_obj_dir(lp_problem) == GLPK.MAX)
        target_value = floor(objective_value/TOL)*(TOL);
    else
        target_value = ceil(objective_value/TOL)*(TOL);
    end

    # do we have any additional flux constraints?
    number_of_extra_flux_constraints = length(additional_flux_constraints_array)

    # Add an extra row for the objective -
    m = GLPK.add_rows(lp_problem,number_of_extra_flux_constraints+1)
    if (GLPK.get_obj_dir(lp_problem) == GLPK.MAX)
        GLPK.set_row_bnds(lp_problem, m, GLPK.LO,target_value,0.0);
    else
        GLPK.set_row_bnds(lp_problem, m, GLPK.UP,0.0,target_value);
    end

    # populate new row, add the object value -
    index_array = Int[]
    value_array = Float64[]
    for flux_index = 1:number_of_fluxes
        push!(index_array,flux_index)
        push!(value_array,GLPK.get_obj_coef(lp_problem, flux_index))
    end
    GLPK.set_mat_row(lp_problem,m,number_of_fluxes,index_array,value_array)

    # Add additional rows for the extra flux constraints -
    for flux_constraint_index = 1:number_of_extra_flux_constraints

        # what rate do I want to constrain?
        target_flux_index = additional_flux_constraints_array[flux_constraint_index]

        # What is my current row index?
        current_row_index = m + flux_constraint_index

        # What is my target value?
        target_value = calculated_flux_array[target_flux_index]

        # add the row bound -
        if (GLPK.get_obj_dir(lp_problem) == GLPK.MAX)
            GLPK.set_row_bnds(lp_problem, current_row_index, GLPK.LO,target_value,0.0);
        else
            GLPK.set_row_bnds(lp_problem, current_row_index, GLPK.UP,0.0,target_value);
        end

        # Add the extra row to the st array -
        index_array = Int[]
        value_array = Float64[]
        for flux_index = 1:number_of_fluxes
            push!(index_array,flux_index)

            if (flux_index == target_flux_index)
                push!(value_array,1.0)
            else
                push!(value_array,0.0)
            end
        end
        GLPK.set_mat_row(lp_problem,current_row_index,number_of_fluxes,index_array,value_array)
    end

    # zero out the objective coefficients -
    for flux_index = 1:number_of_fluxes
        GLPK.set_obj_coef(lp_problem, flux_index, 0.0);
    end

    solver_parameters.presolve = GLPK.OFF;
    solver_parameters.msg_lev = GLPK.MSG_OFF;
    solver_parameters.tm_lim = 1000*TIME_RESTART_LIM;

    # main FVA loop -
    for round_index = 1:2

        # are we min are maximizing the flux?
        if round_index == 1
            GLPK.set_obj_dir(lp_problem, GLPK.MIN);
        else
            GLPK.set_obj_dir(lp_problem, GLPK.MAX);
        end

        for flux_index = 1:number_of_fluxes

            # run the calc -
            GLPK.set_obj_coef(lp_problem, flux_index, 1.0);
            exit_flag = GLPK.simplex(lp_problem, solver_parameters);

            # reset the objective coefficient -
            GLPK.set_obj_coef(lp_problem,flux_index,0.0)

            if (exit_flag != EXIT_FLAG_SUCCESS)

                solver_parameters.tm_lim = 10000*TIME_RESTART_LIM;
                solver_parameters.presolve = GLPK.ON
                GLPK.adv_basis(lp_problem, 0);
                exit_flag = GLPK.simplex(lp_problem, solver_parameters);
                if (exit_flag != EXIT_FLAG_SUCCESS)
                    # ok, something muy malo is happening. throw an error ...
                    error("ERROR: Unable to solve modified problem. Exiting ... ")
                end
            end

            # reset the problem options -
            solver_parameters.presolve = GLPK.OFF;
            solver_parameters.tm_lim = 1000*TIME_RESTART_LIM;

            if (GLPK.get_obj_dir(lp_problem) == GLPK.MIN)
                calculated_flux_array[flux_index,2] = GLPK.get_obj_val(lp_problem)
            else
                calculated_flux_array[flux_index,3] = GLPK.get_obj_val(lp_problem)
            end
        end
    end
    # ------------------------------------------------------------------------------- #

    # return the array -
    return (calculated_flux_array,dual_value_array)
end
