function [combined_LL, unscaled_gradient_vector, grad_parameter_names] = ...
    LL_calculator_strains_pairwise_from_me_dist(param_vals,...
    input_value_dict, pre_MLE_output_dict)
	% EP 17-11-07

    % calculates the log likelihood of observing a list of differences between
        % pairs of test and reference strain colony GRs
    % can work with 1+ test strains, 1+ pairs/strain
    % assumes different GR variance for petite and non-petite colonies,
        % but same GR variance for ref and all test strains

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gradient_specification = input_value_dict('gradient_specification');
    mle_parameter_names = input_value_dict('mle_parameter_names');
%    max_neg_LL_val = input_value_dict('max_neg_LL_val');
    DFE_parameters = input_value_dict('DFE_parameters');
    
    strain_list = pre_MLE_output_dict('strain_list');
    fitted_parameters = pre_MLE_output_dict('fitted_parameters');
    test_strain_list_by_pair = pre_MLE_output_dict('test_strain_list_by_pair');
    GR_diff_list = pre_MLE_output_dict('GR_diff_list');
    me_pdf = pre_MLE_output_dict('me_pdf');
    me_pdf_xvals = pre_MLE_output_dict('me_pdf_xvals');
    strain_pdf_xvals = pre_MLE_output_dict('strain_pdf_xvals');
    me_dist_param_grad_vector_dict = pre_MLE_output_dict('me_dist_param_grad_vector_dict');
%    if gradient_specification
%        d_me_LL_d_me = pre_MLE_output_dict('d_me_LL_d_me');
%    end
    
    parameter_dict = containers.Map(mle_parameter_names, param_vals);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify parameter values
    
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

    test_strain_number = length(strain_list);
        % # of strains besides reference strain whose likelihood being fitted

    ref_sigma = nonpetite_colony_sigma;
    test_sigma = nonpetite_colony_sigma;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop through test strains; within every test strain, loop through fields,
        % and calculate their log likelihood

    % initialize the log likelihood
    combined_LL = 0;
    
    gradient_dict = containers.Map('KeyType', 'char', ...
            'ValueType', 'any');
    gradient_keys = [mle_parameter_names, DFE_parameters];

    % Pre-set all gradients to 0, and then only calculate the actual
        % values if those gradients are for fitted parameters
    for current_param_idx = 1:length(gradient_keys)
        current_param_name = gradient_keys{current_param_idx};
        gradient_dict(current_param_name) = 0;
    end

    fitted_parameters = strrep(fitted_parameters, '_colony_sigma','_sigma');
    if any(strcmp('nonpetite_sigma',fitted_parameters))
        fitted_parameters = [fitted_parameters, {'ref_sigma', 'test_sigma'}];
    end
    
    for strain_idx = 1:test_strain_number
        current_strain = strain_list{strain_idx};
        current_indices = find(strcmp(test_strain_list_by_pair, ...
            current_strain));
        current_sample_size = size(current_indices, 1);
        % only run mle/ll calc if there are samples from current_strain
            % (this is an issue when working with a random subset of data)
        if current_sample_size > 0
            
            current_strain_pp_name = strcat(current_strain,'_pp');
%            current_strain_me_name = strcat(current_strain,'_me');
            test_petite_prop = parameter_dict(current_strain_pp_name);
%            test_mut_effect = parameter_dict(current_strain_me_name);
%            test_mean = ref_mean * (1 + test_mut_effect);

            strain_GR_diff_list = GR_diff_list(current_indices);

%            test_mut_effect_likelihood = ...
%                interp1(me_pdf_xvals, me_pdf, test_mut_effect);
%            if test_mut_effect_likelihood > 0
%                test_mut_effect_LL = ...
%                    log(test_mut_effect_likelihood);
%            else
%                test_mut_effect_LL = - max_neg_LL_val;
%            end

            if gradient_specification
                
%                d_LL_d_test_mut_effect_dist = ...
%                    interp1(me_pdf_xvals, d_me_LL_d_me, test_mut_effect);

                temp_fitted_parameters = {};
                if any(strcmp(current_strain_pp_name, fitted_parameters))
                    temp_fitted_parameters = [temp_fitted_parameters, ...
                        {'test_petite_prop'}];
                end
%                if any(strcmp(current_strain_me_name, fitted_parameters))
%                    temp_fitted_parameters = [temp_fitted_parameters, ...
%                        {'test_mean'}];
%                end
                current_fitted_parameters = [fitted_parameters, ...
                    temp_fitted_parameters];

                [current_LL, current_grad_dict] = ...
                    LL_calculator_observed_diffs_given_me_dist(ref_mean, ...
                        test_petite_prop, ref_petite_prop, petite_mean, ...
                        test_sigma, ref_sigma, petite_colony_sigma, me_pdf, ...
                        me_pdf_xvals, strain_pdf_xvals, ...
                        me_dist_param_grad_vector_dict, ...
                        DFE_parameters, strain_GR_diff_list, ...
                        gradient_specification, current_fitted_parameters);
%                % Likelihood = L(data | mut_effect) * L(mut_effect | me_pdf)
%                current_LL = data_LL + test_mut_effect_LL;

                % add test strain sigma gradient to ref strain sigma to get
                    % nonpetite sigma gradient
                current_grad_dict('nonpetite_sigma') = ...
                    current_grad_dict('ref_sigma') + ...
                    current_grad_dict('test_sigma');

%                % convert derivative of test strain mean to derivative of
%                    % mutational effect, taking into account gradient of
%                    % test_me relative to its source distribution
%                current_grad_dict('test_me') = ...
%                    current_grad_dict('test_mean') * ref_mean + ...
%                    d_LL_d_test_mut_effect_dist;

                gradient_dict('petite_colony_sigma') = ...
                    gradient_dict('petite_colony_sigma') + ...
                    current_grad_dict('petite_sigma');
                gradient_dict('nonpetite_colony_sigma') = ...
                    gradient_dict('nonpetite_colony_sigma') + ...
                    current_grad_dict('nonpetite_sigma');
                gradient_dict('petite_mean') = ...
                    gradient_dict('petite_mean') + ...
                    current_grad_dict('petite_mean');
                gradient_dict('ref_mean') = ...
                    gradient_dict('ref_mean') + ...
                    current_grad_dict('ref_mean');
                gradient_dict('ref_petite_prop') = ...
                    gradient_dict('ref_petite_prop') + ...
                    current_grad_dict('ref_petite_prop');
                gradient_dict(current_strain_pp_name) = ...
                    gradient_dict(current_strain_pp_name) + ...
                    current_grad_dict('test_petite_prop');
%                gradient_dict(current_strain_me_name) = ...
%                    gradient_dict(current_strain_me_name) + ...
%                    current_grad_dict('test_me');

                % Add DFE parameter gradients to gradient_dict
                for DFE_param_idx = 1:length(DFE_parameters)
                    current_param = DFE_parameters{DFE_param_idx};
                    if any(strcmp(fitted_parameters, current_param))
                        gradient_dict(current_param) = ...
                            gradient_dict(current_param) + ...
                            current_grad_dict(current_param);
                    end
                end

            else
                current_fitted_parameters = {};
                current_LL = ...
                    LL_calculator_observed_diffs_given_me_dist(ref_mean, ...
                        test_petite_prop, ref_petite_prop, petite_mean, ...
                        test_sigma, ref_sigma, petite_colony_sigma, me_pdf, ...
                        me_pdf_xvals, strain_pdf_xvals, ...
                        me_dist_param_grad_vector_dict, ...
                        DFE_parameters, strain_GR_diff_list, ...
                        gradient_specification, current_fitted_parameters);
%                % Likelihood = L(data | mut_effect) * L(mut_effect | me_pdf)
%                current_LL = data_LL + test_mut_effect_LL;
            end
            combined_LL = combined_LL + current_LL;
        end
    end
    
    unscaled_gradient_vector = cell2mat(values(gradient_dict));
    grad_parameter_names = keys(gradient_dict);

end