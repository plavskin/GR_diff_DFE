function [LL, gradient_dict] = ...
    LL_calculator_observed_diffs_given_me_dist(ref_mean, test_petite_prop, ...
    ref_petite_prop, petite_mean, test_sigma, ref_sigma, petite_sigma, ...
    me_pdf, me_pdf_xvals, strain_pdf_xvals, me_dist_param_grad_vector_dict, DFE_parameters, ...
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set maximum amount of memory that matrices can occupy
    max_matrix_mem_in_bytes = 0.002*1024^3; % 2 Mb
    float_memory = 8;
    single_matrix_mem = max_matrix_mem_in_bytes/length(fitted_parameters);
    single_matrix_floats = single_matrix_mem/float_memory;
    col_num = length(me_pdf_xvals);
    allowed_row_num = max(floor(single_matrix_floats / col_num), 1);
    data_point_num = length(strain_GR_diff_list);
    data_section_num = ceil(data_point_num / allowed_row_num);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    general_param_names = {'ref_petite_prop', 'test_petite_prop', ...
        'ref_mean', 'petite_mean', 'test_mean', 'ref_sigma', 'petite_sigma', ...
        'test_sigma'};

    if any(strcmp('ref_mean', fitted_parameters))
        if ~any(strcmp('test_mean', fitted_parameters))
            fitted_parameters = [fitted_parameters, 'test_mean'];
        end
    end

    fitted_general_param_names = intersect(general_param_names, fitted_parameters);

    test_strain_mean_grid = ref_mean .* exp(strain_pdf_xvals);
    
    % pre-calculate some vectors used in LL_calc_within_pair_fast
    mean_number = length(test_strain_mean_grid);
    if mean_number > 1
        test_strain_mean_grid = ...
            reshape(test_strain_mean_grid, [1 mean_number]);
    end
    
    lambda_pp = test_petite_prop*ref_petite_prop;
	lambda_pr = test_petite_prop*(1-ref_petite_prop);
	lambda_tp = (1-test_petite_prop)*ref_petite_prop;
	lambda_tr = (1-test_petite_prop)*(1-ref_petite_prop);

	mean_pp = 0;
	mean_pr = petite_mean-ref_mean;
	mean_tp = test_strain_mean_grid-petite_mean;
	mean_tr = test_strain_mean_grid-ref_mean;

	sigma_pp = petite_sigma * sqrt(2);	
	sigma_pr = sqrt(petite_sigma^2 + ref_sigma^2);
	sigma_tp = sqrt(petite_sigma^2 + test_sigma^2);
	sigma_tr = sqrt(ref_sigma^2 + test_sigma^2);
    
    % initialize LL_observed_diffs_sparse with ones at every position
    LL_observed_diffs_sparse = ones(size(test_strain_mean_grid));
    
    for current_GR_diff_counter = 1:data_section_num
        start_idx = 1 + allowed_row_num * (current_GR_diff_counter - 1);
        end_idx = min(allowed_row_num * current_GR_diff_counter, ...
            data_point_num);
        current_indices = start_idx:end_idx;
        current_GR_diff = strain_GR_diff_list(current_indices);
        
        [current_LL_observed_diffs, ...
            current_LL_observed_diffs_grad_dict] = ...
                LL_calc_within_pair_fast(test_strain_mean_grid, ref_mean, ...
                    test_petite_prop, ref_petite_prop, petite_mean, test_sigma, ...
                    ref_sigma, petite_sigma, current_GR_diff, ...
                    gradient_specification, fitted_parameters, ...
                    lambda_pp, lambda_pr, lambda_tp, lambda_tr, ...
                    mean_pp, mean_pr, mean_tp, mean_tr, ...
                    sigma_pp, sigma_pr, sigma_tp, sigma_tr);
        LL_observed_diffs_sparse = ...
            LL_observed_diffs_sparse + current_LL_observed_diffs_sparse;
        if gradient_specification
            if current_GR_diff_counter == 1 | ...
            	~exist('LL_observed_diffs_sparse_grad_dict', 'var')
                LL_observed_diffs_sparse_grad_dict = ...
                    current_LL_observed_diffs_sparse_grad_dict;
            else
                for param_idx = 1:length(fitted_general_param_names)
                    current_param = fitted_general_param_names{param_idx};
                    LL_observed_diffs_sparse_grad_dict(current_param) = ...
                        LL_observed_diffs_sparse_grad_dict(current_param) + ...
                        current_LL_observed_diffs_sparse_grad_dict(current_param);
                end
            end
        end
    end
    
    % interpolate LL_observed_diffs from LL_observed_sparse_diffs
    % (and dp the same for LL_observed_diffs_sparse_grad_dict)
    LL_observed_diffs = ...
        interp1(strain_pdf_xvals, LL_observed_diffs_sparse, me_pdf_xvals);
    if gradient_specification
        LL_observed_diffs_grad_dict = containers.Map('KeyType', 'char', ...
            'ValueType', 'any');
        grad_dict_keys = keys(LL_observed_diffs_sparse_grad_dict);
        for param_idx = 1:length(grad_dict_keys)
            current_param = grad_dict_keys{param_idx};
            current_sparse_grad =...
                LL_observed_diffs_sparse_grad_dict(current_param);
            if all(isnan(current_sparse_grad))
                LL_observed_diffs_grad_dict(current_param) = ...
                    current_sparse_grad;
            else
                LL_observed_diffs_grad_dict(current_param) = ...
                    interp1(strain_pdf_xvals, ...
                        LL_observed_diffs_sparse_grad_dict(current_param), ...
                        me_pdf_xvals);
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
    % correct for <0 values in me_pdf, which cause imaginary likelihoods
        % unless something is horribly wrong, the magnitude of these should
        % be tiny
    me_pdf(me_pdf < 0) = 0;
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
            d_LL_d_current_param = ...
                general_param_LL_grad_calculator(...
                    LL_observed_diffs_grad_dict, current_param, ...
                    likelihood_to_integrate_over_uncorrected, me_pdf_xvals, ...
                    likelihood_uncorrected);
            if strcmp(current_param, 'ref_mean')
                % account for test strain mean's contribution to ref mean
                    % gradient
                % For component of d(LL)/d(ref_mean) that comes from the
                    % test strain mean,
                % d(LL)/d(ref_mean) = d(LL)/d(test_mean) *
                    % d(test_mean)/d(ref_mean);
                % d(test_mean)/d(ref_mean) = exp(test_mut_effect)
                L_to_int_over_uncorr_for_test_mean_contrib_to_ref_grad = ...
                    likelihood_to_integrate_over_uncorrected .* ...
                    exp(me_pdf_xvals);
                d_LL_d_test_mean_contrib_to_ref_mean_grad = ...
                    general_param_LL_grad_calculator(...
                        LL_observed_diffs_grad_dict, 'test_mean', ...
                        L_to_int_over_uncorr_for_test_mean_contrib_to_ref_grad, ...
                        me_pdf_xvals, likelihood_uncorrected);
                d_LL_d_current_param = d_LL_d_current_param + ...
                    d_LL_d_test_mean_contrib_to_ref_mean_grad;
            end
            gradient_dict(current_param) = d_LL_d_current_param;
        end
    end
    
    save('~/Documents/yeast stuff/temp.mat')

end
