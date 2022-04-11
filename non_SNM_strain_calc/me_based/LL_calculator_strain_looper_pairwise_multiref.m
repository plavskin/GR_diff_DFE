function [combined_LL, unscaled_gradient_vector, grad_parameter_names] = ...
    LL_calculator_strain_looper_pairwise_multiref(param_vals,...
    input_value_dict, pre_MLE_output_dict)
	% EP 17-11-07

    % calculates the log likelihood of observing a list of differences
        % between pairs of test and reference strain colony GRs
    % can work with 1+ test strains, 1+ reference strains
    % assumes different GR variance for petite and non-petite colonies,
        % but same GR variance for ref and all test strains

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gradient_specification = input_value_dict('gradient_specification');
    mle_parameter_names = input_value_dict('mle_parameter_names');
    
    strain_list = pre_MLE_output_dict('strain_list');
    fitted_parameters = pre_MLE_output_dict('fitted_parameters');
    unique_test_strains = pre_MLE_output_dict('unique_test_strains');
    unique_ref_strains = pre_MLE_output_dict('unique_ref_strains');
    test_strain_ID_mat = pre_MLE_output_dict('test_strain_ID_mat');
    ref_strain_ID_mat = pre_MLE_output_dict('ref_strain_ID_mat');
    GR_diff_list = pre_MLE_output_dict('GR_diff_list');

    
    parameter_dict = containers.Map(mle_parameter_names, param_vals);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify parameter values

    petite_colony_sigma = parameter_dict('petite_colony_sigma');
        % s.d. of colony GRs of petite distribution
    nonpetite_colony_sigma = parameter_dict('nonpetite_colony_sigma');
        % s.d. of colony GRs of non-petite reference and test strain
            % distribution
    petite_mean = parameter_dict('petite_mean');
        % mean growth rate of petite colonies, regardless of genotype
    ref_mean = parameter_dict('ref_mean');
        % mean growth rate of non-petite ref colonies

    ref_sigma = nonpetite_colony_sigma;
    test_sigma = nonpetite_colony_sigma;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set up gradient dictionary
    unscaled_gradient_dict_combined = containers.Map('KeyType', 'char', ...
        'ValueType', 'any');

    grad_parameter_names = ...
        [{'petite_mean', 'petite_colony_sigma', ...
            'nonpetite_colony_sigma', 'ref_mean'}, ...
        strcat(strain_list, '_pp'), strcat(strain_list, '_me')];

    % Pre-set all gradients to NaN, and then only calculate the actual
        % values if those gradients are for fitted parameters
    for current_param_idx = 1:length(grad_parameter_names)
        current_param_name = grad_parameter_names{current_param_idx};
        unscaled_gradient_dict_combined(current_param_name) = NaN;
    end
        
    test_strain_number = length(unique_test_strains);
    ref_strain_number = length(unique_ref_strains);

    sigma_pp = petite_colony_sigma * sqrt(2);  
    sigma_pr = sqrt(petite_colony_sigma^2 + ref_sigma^2);
    sigma_tp = sqrt(petite_colony_sigma^2 + test_sigma^2);
    sigma_tr = sqrt(ref_sigma^2 + test_sigma^2);

    mean_pp = 0;

    combined_LL = 0;

    % make list of names to be used to determine which gradients to calc
    general_fitted_parameters = {};
    if any(strcmp('petite_colony_sigma', fitted_parameters))
        general_fitted_parameters = ...
            [general_fitted_parameters, 'petite_sigma'];
    end
    if any(strcmp('petite_mean', fitted_parameters))
        general_fitted_parameters = ...
            [general_fitted_parameters, 'petite_mean'];
    end
    if any(strcmp('nonpetite_colony_sigma', fitted_parameters))
        general_fitted_parameters = ...
            [general_fitted_parameters, 'ref_sigma', 'test_sigma'];
    end
    possible_fitted_parameter_general_names = ...
        {'test_petite_prop', 'test_mean', 'ref_petite_prop', 'ref_mean'};

    % loop through every unique combination of reference and test strains,
        % calculate likelihoods
    for test_strain_idx = 1:test_strain_number
        current_test_strain = unique_test_strains{test_strain_idx};
        test_pp_name = strcat(current_test_strain, '_pp');
        test_me_name = strcat(current_test_strain, '_me');
        test_petite_prop = parameter_dict(test_pp_name);
        test_mut_effect = parameter_dict(test_me_name);
        test_mean = ref_mean * exp(test_mut_effect);

        for ref_strain_idx = 1:ref_strain_number

            current_ref_strain = unique_ref_strains{ref_strain_idx};
            ref_pp_name = strcat(current_ref_strain, '_pp');
            ref_me_name = strcat(current_ref_strain, '_me');
            ref_petite_prop = parameter_dict(ref_pp_name);
            current_ref_mut_effect = parameter_dict(ref_me_name);
            current_ref_mean = ref_mean * exp(current_ref_mut_effect);

            % pick out growth rate differences corresponding to current
                % combination of reference and test strain
            current_data_indices = ...
                test_strain_ID_mat(:, test_strain_idx) & ...
                ref_strain_ID_mat(:, ref_strain_idx);

            if sum(current_data_indices, 1) > 0

                % determine which parameters currently need to be fitted
                possible_fitted_parameters = ...
                    {test_pp_name, test_me_name, ref_pp_name, ref_me_name};
                [~, current_fitted_parameter_indices] = ...
                    intersect(possible_fitted_parameters, ...
                        fitted_parameters);

                current_fitted_parameters = ...
                    [general_fitted_parameters, ...
                    possible_fitted_parameter_general_names{current_fitted_parameter_indices}];

                % determine which growth rate differences pertain to
                    % current ref and test strain
                current_GR_diff_list = GR_diff_list(current_data_indices);

                % pre-set up distribution parameters
                lambda_pp = test_petite_prop*ref_petite_prop;
                lambda_pr = test_petite_prop*(1-ref_petite_prop);
                lambda_tp = (1-test_petite_prop)*ref_petite_prop;
                lambda_tr = (1-test_petite_prop)*(1-ref_petite_prop);
                
                mean_pr = petite_mean-current_ref_mean;
                mean_tp = test_mean-petite_mean;
                mean_tr = test_mean-current_ref_mean;

                [current_LL, current_gradient_dict] = ...
                    LL_calc_within_pair_fast(test_mean, current_ref_mean, ...
                        test_petite_prop, ref_petite_prop, petite_mean, ...
                        test_sigma, ref_sigma, petite_colony_sigma, ...
                        current_GR_diff_list, gradient_specification, ...
                        current_fitted_parameters, lambda_pp, ...
                        lambda_pr, lambda_tp, lambda_tr, mean_pp, ...
                        mean_pr, mean_tp, mean_tr, sigma_pp, ...
                        sigma_pr, sigma_tp, sigma_tr);

                combined_LL = combined_LL + current_LL;

                if gradient_specification
                    current_gradient_dict('test_me') = ...
                        current_gradient_dict('test_mean') * ...
                        ref_mean * exp(test_mut_effect);
                    current_gradient_dict('current_ref_me') = ...
                        current_gradient_dict('ref_mean') * ...
                        ref_mean * exp(current_ref_mut_effect);
                    current_ref_mean_grad = ...
                        current_gradient_dict('test_mean') * exp(test_mut_effect) + ...
                        current_gradient_dict('ref_mean') * exp(current_ref_mut_effect);
                    unscaled_gradient_dict_combined('petite_mean') = ...
                        nansum([unscaled_gradient_dict_combined('petite_mean'), ...
                            current_gradient_dict('petite_mean')]);
                    unscaled_gradient_dict_combined('petite_colony_sigma') = ...
                        nansum([unscaled_gradient_dict_combined('petite_colony_sigma'), ...
                            current_gradient_dict('petite_sigma')]);
                    unscaled_gradient_dict_combined('nonpetite_colony_sigma') = ...
                        nansum([unscaled_gradient_dict_combined('nonpetite_colony_sigma'), ...
                            current_gradient_dict('test_sigma'), ...
                            current_gradient_dict('ref_sigma')]);
                    unscaled_gradient_dict_combined(ref_pp_name) = ...
                        nansum([unscaled_gradient_dict_combined(ref_pp_name), ...
                            current_gradient_dict('ref_petite_prop')]);
                    unscaled_gradient_dict_combined(test_pp_name) = ...
                        nansum([unscaled_gradient_dict_combined(test_pp_name), ...
                            current_gradient_dict('test_petite_prop')]);
                    unscaled_gradient_dict_combined(ref_me_name) = ...
                        nansum([unscaled_gradient_dict_combined(ref_me_name), ...
                            current_gradient_dict('current_ref_me')]);
                    unscaled_gradient_dict_combined(test_me_name) = ...
                        nansum([unscaled_gradient_dict_combined(test_me_name), ...
                            current_gradient_dict('test_me')]);
                    unscaled_gradient_dict_combined('ref_mean') = ...
                        nansum([unscaled_gradient_dict_combined('ref_mean'), ...
                            current_ref_mean_grad]);
                end
            end

        end
    end

    unscaled_gradient_vector = cell2mat(values(unscaled_gradient_dict_combined));
    grad_parameter_names = keys(unscaled_gradient_dict_combined);

end