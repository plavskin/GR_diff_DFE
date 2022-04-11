function [combined_LL, unscaled_gradient_vector, grad_parameter_names] = ...
    LL_calculator_me_estimate_single_dfe(param_vals,...
    input_value_dict, pre_MLE_output_dict)
    
    % Loops through test strains, convolving current version of the
    % distribution of fitness effects with the standard error of that
    % strain's mutation effect estimate, and then estimating the likelihood
    % of observing that strain's mutational effect given the convolved
    % distribution
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    max_neg_LL_val = input_value_dict('max_neg_LL_val');
    gradient_specification = input_value_dict('gradient_specification');
    L = input_value_dict('L');
    current_model = input_value_dict('model');

    strain_list = pre_MLE_output_dict('strain_list');
    strain_me_estimate_list = ...
        pre_MLE_output_dict('strain_me_estimate_list');
    strain_me_ste_list = pre_MLE_output_dict('strain_me_ste_list');
    lambda_mult_list = pre_MLE_output_dict('lambda_mult_list');
    lambda_mult_strain_idx_list = pre_MLE_output_dict('lambda_mult_strain_idx_list');
    N = pre_MLE_output_dict('N');
    me_pdf_frequency_vals = pre_MLE_output_dict('me_pdf_frequency_vals');
    me_pdf_xvals = pre_MLE_output_dict('me_pdf_xvals');
    
    test_strain_number = length(strain_list);
        % # of strains besides whose mutational effect has been estimated

    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate likelihood of observing mutation effects given current
        % global parameters (distribution of fitness effects, DFE)
    % First set up function that will calculate fourier transform of DFE
    fourier_dfe_fun_name = strcat('fourier_dfe_', current_model);
    fourier_dfe_function = str2func(fourier_dfe_fun_name);
    % Calculate fourier transform of DFE
    % lambdas may be different, need to account for this!
    fourier_output_dict = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
    for lambda_mult_idx = 1:length(lambda_mult_list)
        curr_lambda_mult = lambda_mult_list(lambda_mult_idx);
        curr_param_vals = param_vals;
        lambda_idx = find(strcmp(input_value_dict('mle_parameter_names'),'lambda_SNM'));
        curr_param_vals(lambda_idx) = param_vals(lambda_idx)*curr_lambda_mult;
        [curr_Fa, curr_Fa_gradient_dict, curr_complete_me_parameter_list, ...
            curr_fitted_parameters_me_fullname] = ...
            fourier_dfe_function(curr_param_vals, input_value_dict, pre_MLE_output_dict);
        fourier_output_dict(lambda_mult_idx) = ...
            {curr_Fa, curr_Fa_gradient_dict, curr_complete_me_parameter_list, ...
            curr_fitted_parameters_me_fullname};
    end
%    [Fa, Fa_gradient_dict, complete_me_parameter_list, ...
%        fitted_parameters_me_fullname] = fourier_dfe_function(...
%            param_vals, input_value_dict, pre_MLE_output_dict);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % initialize the log likelihood
    combined_LL = 0;
    me_dist_param_grad_dict = containers.Map('KeyType', 'char', ...
        'ValueType', 'any');
    
    for strain_idx = 1:test_strain_number

        current_mut_effect = strain_me_estimate_list(strain_idx);
        current_me_ste = strain_me_ste_list(strain_idx);
        current_lambda_idx = lambda_mult_strain_idx_list(strain_idx);

        curr_fourier_output_vals = fourier_output_dict(current_lambda_idx);
        [Fa, Fa_gradient_dict, complete_me_parameter_list, ...
            fitted_parameters_me_fullname] = curr_fourier_output_vals{:};
        
        current_F_kernel = fourier_domain_gauss(me_pdf_frequency_vals, ...
            0, current_me_ste, {}, false);
        % Calculate Fs, which is Fa smoothed with current_F_kernel
        empty_F_kernel_grad_dict = containers.Map('KeyType', 'char', ...
            'ValueType', 'any');
        [Fs, Fs_gradient_dict] = derivative_multiplier(Fa, ...
            current_F_kernel, Fa_gradient_dict, empty_F_kernel_grad_dict);
        % Calculate discrete fourier transform of Fs to get a smoothed
            % version of the pdf of mutational effects per strain (and, if
            % applicable, the gradients of the log of the pdf with respect
            % to each parameter in gradient_param_list)
        [real_me_pdf_smooth, temp_me_pdf_xvals, ...
            current_me_dist_param_grad_vector_dict] = ...
            likelihood_from_fourier(Fs, me_pdf_xvals, N, L, ...
                Fs_gradient_dict, gradient_specification);
            
        current_mut_effect_likelihood = interp1(temp_me_pdf_xvals, ...
            real_me_pdf_smooth, current_mut_effect);
        if current_mut_effect_likelihood > 0
            current_mut_effect_LL = ...
                log(current_mut_effect_likelihood);
        else
            current_mut_effect_LL = - max_neg_LL_val;
        end
        combined_LL = combined_LL + current_mut_effect_LL;
        
        % calculate me_dist_param_grad_dict by interpolating strain me vals
            % in the gradients relative to each me distribution parameter
        if gradient_specification
            % If calculating gradients of LL, need to divide gradient relative
                % to parameter by real_me_pdf_smooth, so remove any values where
                % real_me_pdf_smooth is 0
            non_zero_positions = (real_me_pdf_smooth ~= 0);
            need_to_remove_zero_vals = ...
                length(real_me_pdf_smooth) > sum(non_zero_positions);
            if need_to_remove_zero_vals
                real_me_pdf_smooth = real_me_pdf_smooth(non_zero_positions);
                temp_me_pdf_xvals = temp_me_pdf_xvals(non_zero_positions);
            end
            for current_fitted_param_idx = ...
                    1:length(complete_me_parameter_list)
                current_fitted_param = ...
                    complete_me_parameter_list{current_fitted_param_idx};
                
                if any(strcmp(fitted_parameters_me_fullname, ...
                        current_fitted_param))
                    d_L_d_current_fitted_param = ...
                        current_me_dist_param_grad_vector_dict(...
                            current_fitted_param);

                    if need_to_remove_zero_vals
                        % remove values where real_pdf_vect is 0
                        d_L_d_current_fitted_param = ...
                            d_L_d_current_fitted_param(non_zero_positions);
                    end
                    d_LL_d_current_fitted_param = ...
                        real(d_L_d_current_fitted_param) ./ ...
                            real_me_pdf_smooth;
%                   figure; plot(temp_me_pdf_xvals, d_LL_d_current_fitted_param, movmean(temp_me_pdf_xvals,2), movmean(d_LL_d_current_fitted_param,2)); title(current_fitted_param);               
                    d_LL_d_current_fitted_param_current_data = ...
                        nansum(interp1(temp_me_pdf_xvals, ...
                            d_LL_d_current_fitted_param, ...
                            current_mut_effect));
                else
                    d_LL_d_current_fitted_param_current_data = NaN;
                end
                
                if isKey(me_dist_param_grad_dict, current_fitted_param)
                    me_dist_param_grad_dict(current_fitted_param) = ...
                        nansum(...
                            [me_dist_param_grad_dict(current_fitted_param),...
                            d_LL_d_current_fitted_param_current_data]);
                else
                    me_dist_param_grad_dict(current_fitted_param) = ...
                        d_LL_d_current_fitted_param_current_data;
                end
            end            
        end

    end

    unscaled_gradient_vector = cell2mat(values(me_dist_param_grad_dict));
    grad_parameter_names = keys(me_dist_param_grad_dict);

end