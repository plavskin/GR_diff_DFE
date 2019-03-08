function [LL, gradient_dict] = ...
    LL_calculator_observed_diffs_given_me_dist(ref_mean, test_petite_prop, ...
    ref_petite_prop, petite_mean, test_sigma, ref_sigma, petite_sigma, ...
    me_pdf, me_pdf_xvals, me_dist_param_grad_vector_dict, DFE_parameters, ...
    strain_GR_diff_list, gradient_specification, fitted_parameters)
    % Calculates the likelihood of observed differences in growth rate
        % between colonies of test_strain and ref_strain, given that the
        % distribution of mutational effects on growth me_pdf_xvals is
        % described by real_me_pdf_smooth
    % 'me' is the mutational effect
    % test_strain_mean = ref_mean * exp(me)
    % L (observations ) = 
    %    integral across region of interest(
    %        P(mutational effect = x | theta{DME}) *
    %        product across strain_GR_diff_list(
    %            P(GR_diff | mutational effect = x, theta{general})
    %            ) * dx
    %        )
    
    general_param_names = {'ref_petite_prop', 'test_petite_prop', ...
        'ref_mean', 'petite_mean', 'test_mean', 'ref_sigma', 'petite_sigma', ...
        'test_sigma'};
    fitted_general_param_names = intersect(general_param_names, fitted_parameters);

    test_strain_mean_grid = ref_mean .* exp(me_pdf_xvals);

    % initialize LL_observed_diffs with ones at every position
    LL_observed_diffs = ones(size(test_strain_mean_grid));

    % Although LL_calculator_within_pair_different_sigmas_simple can
        % accept a list of strain_GR_diff_list values, this would result
        % in the creation of a large matrix, so pass strain_GR_diff_list
        % values one at a time
    for current_GR_diff_counter = 1:length(strain_GR_diff_list)
        current_GR_diff = strain_GR_diff_list(current_GR_diff_counter);
        [current_LL_observed_diffs, ...
            current_LL_observed_diffs_grad_dict] = ...
                LL_calc_within_pair(test_strain_mean_grid, ref_mean, ...
                    test_petite_prop, ref_petite_prop, petite_mean, test_sigma, ...
                    ref_sigma, petite_sigma, current_GR_diff, ...
                    gradient_specification, fitted_parameters);
        LL_observed_diffs = ...
            LL_observed_diffs + current_LL_observed_diffs;
        if gradient_specification
            if current_GR_diff_counter == 1 | ...
            	~exist('LL_observed_diffs_grad_dict', 'var')
                LL_observed_diffs_grad_dict = ...
                    current_LL_observed_diffs_grad_dict;
            else
                for param_idx = 1:length(fitted_general_param_names)
                    current_param = fitted_general_param_names{param_idx};
                    LL_observed_diffs_grad_dict(current_param) = ...
                        LL_observed_diffs_grad_dict(current_param) + ...
                        current_LL_observed_diffs_grad_dict(current_param);
                end
            end
        end
    end

    % solve the numerical issue of having huge likelihood values to
        % integrate over by dividing both the integral and the log
        % likelihood by a large correction term, k, which can be done
        % in logspace

    % log( L (observations ) ) - log(k) = 
    %    integral across region of interest(
    %        exp(
    %           log (P(mutational effect = x | theta{DME})) +
    %           log (sum across strain_GR_diff_list(
    %               P(GR_diff | mutational effect = x, theta{general})
    %               )
    %           )
    %           - log( k )
    %        ) * dx
    %    )

    % Set log_k to be rel_log_diff less than the max of the log of the likelihood
        % function to integrate over, so that the numbers the
        % integration is actually performed on are always reasonable
    % This can be a very large number, so no point in calculating k directly
    rel_log_diff = 100;
        % exp(100) can still be evaluated, and using a high number here
        % minimizes the number of points whose value in log space is so low
        % that it exponentiates to 0
    log_likelihood_to_integrate_over = log(me_pdf) + LL_observed_diffs;
    log_k = max(log_likelihood_to_integrate_over) - rel_log_diff;


    likelihood_to_integrate_over_uncorrected = ...
        exp(log_likelihood_to_integrate_over - log_k);

    %figure; plot(me_pdf_xvals, likelihood_to_integrate_over_uncorrected, '-b', me_pdf_xvals, exp(LL_observed_diffs - log_k), '-r')
    
    likelihood_uncorrected = trapz(me_pdf_xvals, likelihood_to_integrate_over_uncorrected);
    LL = log(likelihood_uncorrected) + log_k;
    
    % calculate gradient_dict by interpolating strain me vals
        % in the gradients relative to each me distribution parameter
    gradient_dict = containers.Map('KeyType', 'char', ...
            'ValueType', 'any');

    if gradient_specification
        % Calculate gradient of log likelihood wrt DFE_parameters
        % d ( L (observations)) / d ( theta{DME} ) ) = 
        %    integral across region of interest(
        %        (     d ( P(mutational effect = x | theta{DME}) ) /
        %            d ( theta{DME} ) ) *
        %        product across strain_GR_diff_list(
        %            P(GR_diff | mutational effect = x, theta{general})
        %            ) *
        %        dx
        %        )
        for current_fitted_DFE_param_idx = 1:length(DFE_parameters)
            current_fitted_DFE_param = DFE_parameters{current_fitted_DFE_param_idx};
            if any(strcmp(fitted_parameters, current_fitted_DFE_param))
                d_me_pdf_d_current_fitted_DFE_param = ...
                    me_dist_param_grad_vector_dict(current_fitted_DFE_param);
                d_L_to_int_over_uncorr_d_current_fitted_DFE_param = ...
                    d_me_pdf_d_current_fitted_DFE_param .* ...
                    exp(LL_observed_diffs - log_k);
                d_likelihood_uncorrected_d_current_fitted_DFE_param = ...
                    trapz(me_pdf_xvals, ...
                         d_L_to_int_over_uncorr_d_current_fitted_DFE_param);
                d_LL_d_current_fitted_DFE_param = ...
                    d_likelihood_uncorrected_d_current_fitted_DFE_param / ...
                    likelihood_uncorrected;
                gradient_dict(current_fitted_DFE_param) = ...
                    d_LL_d_current_fitted_DFE_param;
            else
                gradient_dict(current_fitted_DFE_param) = NaN;
            end
        end

        % Calculate gradient of log likelihood wrt 'general' parameters
            % (non-DFE)
        % d ( L (observations)) / d( theta{general} ) ) = 
        %    integral across region of interest(
        %        P(mutational effect = x | theta{DME}) *
        %        product across strain_GR_diff_list(
        %            P(GR_diff | mutational effect = x, theta{general})
        %            ) *
        %        sum across strain_GR_diff_list(
        %            d ( ln ( P(GR_diff | mutational effect = x, theta{general}) ) ) /
        %                d ( theta{general} )
        %            ) * 
        %        dx
        %        )
        for param_idx = 1:length(general_param_names)
            current_param = general_param_names{param_idx};
            d_LL_observed_diffs_d_current_param = ...
                LL_observed_diffs_grad_dict(current_param);
            d_L_to_integrate_over_uncorrected_d_current_param = ...
                likelihood_to_integrate_over_uncorrected .* ...
                d_LL_observed_diffs_d_current_param;
            if strcmp(current_param, 'ref_mean')
                % account for test strain mean's contribution to ref mean
                    % gradient
                    % d(LL)/d(ref_mean) = d(LL)/d(test_mean) *
                        % d(test_mean)/d(ref_mean);
                    % d(test_mean)/d(ref_mean) = exp(test_mut_effect)
                d_LL_observed_diffs_d_test_mean = ...
                	LL_observed_diffs_grad_dict('test_mean');
                d_L_to_integrate_over_uncorrected_d_test_mean = ...
	                likelihood_to_integrate_over_uncorrected .* ...
	                d_LL_observed_diffs_d_test_mean;
                test_mean_contrib_to_ref_mean_grad = ...
                    d_L_to_integrate_over_uncorrected_d_test_mean .* ...
                    exp(me_pdf_xvals);
                d_L_to_integrate_over_uncorrected_d_current_param = ...
                	d_L_to_integrate_over_uncorrected_d_current_param .* ...
                	test_mean_contrib_to_ref_mean_grad;
                % d_L_to_integrate_over_uncorrected_d_current_param still
                % needs to be multiplied by exp(log_k), but doing so at
                % this step could result in numerical issues, so rescale
                % after integration
            end
            %figure; plot(me_pdf_xvals, d_L_to_integrate_over_uncorrected_d_current_param)
            d_likelihood_uncorrected_d_current_param = ...
                trapz(me_pdf_xvals, ...
                    d_L_to_integrate_over_uncorrected_d_current_param);
            d_LL_d_current_param = ...
                d_likelihood_uncorrected_d_current_param / ...
                likelihood_uncorrected;
            if strcmp(current_param, 'ref_mean')
                % rescale from above
                d_LL_d_current_param = exp(log_k) * d_LL_d_current_param;
            end
            gradient_dict(current_param) = d_LL_d_current_param;
        end
    end

end
