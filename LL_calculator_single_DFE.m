function [combined_LL, unscaled_gradient_vector, grad_parameter_names] = ...
    LL_calculator_single_DFE(param_vals, input_value_dict, ...
        pre_MLE_output_dict)
    
    % Takes in parameters that apply across all test strains, excluding
        % DFE parameters
    % Estimate the optimal DFE parameters given the current iteration of
        % general parameter values, using a subfunction that optimizes
        % strain parameter values given general parameters and DFE
        % parameters
    % This function writes an output csv file which contains the MLE parameters
        % corresponding to each DFE parameter if the overall LL is higher than
        % previous iterations of this function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    random_effect_names = input_value_dict('random_effect_names');
    mle_parameter_names = input_value_dict('mle_parameter_names');
    gradient_specification = input_value_dict('gradient_specification');
    DFE_parameters = input_value_dict('DFE_parameters');
    current_model = input_value_dict('model');

    strain_list = pre_MLE_output_dict('strain_list');

    parameter_dict = containers.Map(mle_parameter_names, param_vals);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify general parameter values
    
    petite_colony_sigma = parameter_dict('petite_colony_sigma');
        % s.d. of colony GRs of petite distribution
    nonpetite_colony_sigma = parameter_dict('nonpetite_colony_sigma');
        % s.d. of colony GRs of non-petite reference and test strain distribution
    petite_mean = parameter_dict('petite_mean');
        % mean growth rate of petite colonies, regardless of genotype
    ref_mean = parameter_dict('ref_mean');
        % mean growth rate of non-petite ref colonies
    ref_petite_prop = parameter_dict('ref_petite_prop');
        % proportion of petites in ref strain
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create new pre_MLE_output_dict for strain looper LL calculator, which will
        % include general parameter values identified above
    pre_MLE_output_dict_strain_looper = ...
        containers.Map(keys(pre_MLE_output_dict), values(pre_MLE_output_dict));
    pre_MLE_output_dict_strain_looper('petite_colony_sigma') = ...
        petite_colony_sigma;
    pre_MLE_output_dict_strain_looper('nonpetite_colony_sigma') = ...
        nonpetite_colony_sigma;
    pre_MLE_output_dict_strain_looper('petite_mean') = ...
        petite_mean;
    pre_MLE_output_dict_strain_looper('ref_mean') = ...
        ref_mean;
    pre_MLE_output_dict_strain_looper('ref_petite_prop') = ...
        ref_petite_prop;

    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate likelihood of observing petite data given current global parameters

    % replace parameter names in parameter list supplied to LL_mixef_calculator
    input_value_dict_mixef = ...
        containers.Map(keys(input_value_dict),values(input_value_dict));
    parameter_list_mixef = mle_parameter_names;

    ranef_names_from_param_list = strcat(random_effect_names,'_sigma');
    params_to_replace = [{'petite_mean'}, ranef_names_from_param_list];
    replacement_params = [{'petite'}, random_effect_names];
    [~, replacement_idx] = ismember(params_to_replace, parameter_list_mixef);
    parameter_list_mixef(replacement_idx(replacement_idx>0)) = replacement_params;
    input_value_dict_mixef('mle_parameter_names') = parameter_list_mixef;

    [LL_petite, unscaled_gradient_vector_petite, grad_parameter_names_petite] = ...
        LL_mixef_calc(param_vals, input_value_dict_mixef, pre_MLE_output_dict);
    petite_param_gradient_dict = ...
        containers.Map(grad_parameter_names_petite, ...
        unscaled_gradient_vector_petite);
    
    clear input_value_dict_mixef

    % add output of calculation above to pre_MLE_output_dict_strain_looper
    pre_MLE_output_dict_strain_looper('LL_petite') = LL_petite;
    pre_MLE_output_dict_strain_looper('petite_param_gradient_dict') = ...
        petite_param_gradient_dict;
    pre_MLE_output_dict_strain_looper('ranef_names_from_param_list') = ...
        ranef_names_from_param_list;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Identify ML DFE parameters

    % Create a new input_value_dict for strain looper LL calculator,
        % with DFE parameters as top-level parameters
    input_value_dict_strain_looper = ...
        containers.Map(keys(input_value_dict), values(input_value_dict));

    input_value_dict_strain_looper('top_level_parameters') = DFE_parameters;
    input_value_dict_strain_looper('pause_at_end') = false;
    input_value_dict_strain_looper('nonlinear_constraint_function') = ...
        strcat('nonlin_constraint_', current_model);
    input_value_dict_strain_looper('LL_calculator') = ...
        'LL_calculator_strain_looper_pairwise_single_dfe';
    input_value_dict_strain_looper('pre_MLE_function') = 'pre_MLE_nested';
    input_value_dict_strain_looper('post_MLE_function_name') = ...
        'post_MLE_GR_diff_strainwise_with_mixef.m';

    input_value_dict_strain_looper('pre_MLE_output_dict') = ...
        pre_MLE_output_dict_strain_looper;

    % set write_checkpoint to false to avoid overwriting global
        % checkpoint file
    input_value_dict_strain_looper('write_checkpoint') = false;

    % Run MLE_search_executer to identify ML DFE parameters given current
        % general parameters
    [DFE_MLE_param_vals, combined_LL_from_MLE, ~, ~, ~, ...
            DFE_parameters_MLE_output, ~, ~] = ...
        MLE_search_executer(input_value_dict_strain_looper);

    % reorder DFE parameters to match DFE_parameters if necessary
    [~, DFE_parameter_order] = ...
        ismember(DFE_parameters, DFE_parameters_MLE_output);
    DFE_MLE_param_vals = DFE_MLE_param_vals(DFE_parameter_order);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % At ML DFE (and strain) parameters, recalculate LL to get gradient

    if gradient_specification
        % retrieve strain parameter info
        test_strain_ML_file = pre_MLE_output_dict('test_strain_ML_file');
        test_strain_ML_table = readtable(test_strain_ML_file);

        strain_param_names = strcat(strain_list,'_pp');
        strain_param_MLE_vals = test_strain_ML_table{1, strain_param_names};

        % calculate global likelihood gradient given current global
            % parameters and current ML estimates for strain-specific
            % parameters
        % The latter must be provided via fixing their start values in
            % input_value_dict_grad_calc
        parameter_values_grad_calc = DFE_MLE_param_vals;

        input_value_dict_grad_calc = ...
            containers.Map(keys(input_value_dict_strain_looper), ...
                values(input_value_dict_strain_looper));
        pre_MLE_output_dict_grad_calc = ...
            containers.Map(keys(pre_MLE_output_dict_strain_looper), ...
                values(pre_MLE_output_dict_strain_looper));

        fitted_parameter_list_grad_calc = ...
            {'petite_colony_sigma', 'nonpetite_colony_sigma', 'petite_mean', ...
                'ref_mean', 'ref_petite_prop'};
        pre_MLE_output_dict_grad_calc('fitted_parameters') = ...
            fitted_parameter_list_grad_calc;

        % fix all strain parameters (the only ones optimized over within
            % LL_calculator_strain_looper_pairwise_single_dfe)
        parameter_list = input_value_dict_grad_calc('parameter_list');
        strain_param_indices = ...
            parameter_identifier(parameter_list, strain_param_names);
        [~, strain_param_positions] = ...
            ismember(strain_param_names, parameter_list);

        combined_fixed_parameter_array_grad_calc = ...
            input_value_dict_grad_calc('tempfixed_parameter_bool');
        combined_profile_ub_array_unscaled_grad_calc = ...
            input_value_dict_grad_calc('profile_upper_limits');
        combined_profile_lb_array_unscaled_grad_calc = ...
            input_value_dict_grad_calc('profile_lower_limits');
        combined_start_values_array_unscaled_grad_calc = ...
            input_value_dict_grad_calc('starting_parameter_vals');
        combined_length_array_unscaled_grad_calc = ...
            input_value_dict_grad_calc('profile_point_num_list');

        combined_fixed_parameter_array_grad_calc(strain_param_indices) = true;
        combined_profile_ub_array_unscaled_grad_calc(strain_param_positions) = ...
            strain_param_MLE_vals;
        combined_profile_lb_array_unscaled_grad_calc(strain_param_positions) = ...
            strain_param_MLE_vals;
        combined_start_values_array_unscaled_grad_calc(strain_param_positions) = ...
            strain_param_MLE_vals;
        combined_length_array_unscaled_grad_calc(strain_param_positions) = 1;


        input_value_dict_grad_calc('combined_position_array') = {'1'};
        input_value_dict_grad_calc('tempfixed_parameter_bool') = ...
            combined_fixed_parameter_array_grad_calc;
        input_value_dict_grad_calc('profile_upper_limits') = ...
            combined_profile_ub_array_unscaled_grad_calc;
        input_value_dict_grad_calc('profile_lower_limits') = ...
            combined_profile_lb_array_unscaled_grad_calc;
        input_value_dict_grad_calc('starting_parameter_vals') = ...
            combined_start_values_array_unscaled_grad_calc;
        input_value_dict_grad_calc('profile_point_num_list') = ...
            combined_length_array_unscaled_grad_calc;

        % Get final log likelihood and unscaled gradient vector
        [combined_LL, unscaled_gradient_vector, grad_parameter_names] = ...
            LL_calculator_strain_looper_pairwise_single_dfe(...
                parameter_values_grad_calc, input_value_dict_grad_calc, ...
                pre_MLE_output_dict_grad_calc);
    else
        combined_LL = combined_LL_from_MLE;
        unscaled_gradient_vector = [];
        grad_parameter_names = {};
    end

end