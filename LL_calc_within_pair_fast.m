function [LL, unscaled_gradient_dict] = LL_calc_within_pair_fast(test_mean, ...
	ref_mean, test_petite_prop, ref_petite_prop, petite_mean, test_sigma, ...
	ref_sigma, petite_sigma, strain_GR_diff_list, gradient_specification, ...
	fitted_parameters, lambda_pp, lambda_pr, lambda_tp, lambda_tr, ...
    mean_pp, mean_pr, mean_tp, mean_tr, sigma_pp, sigma_pr, sigma_tp, sigma_tr)

	% Calculates the log likelihood of the observed difference between
		% reference and test strain colony pair growth rates, accounting
		% for random proportion of petites

	% Equivalent to LL_calculator_within_field, but assumes
		% a single colony per strain
	% The reason the code is different is because LL_calculator_within_field
		% isn't set up to calculate likelihoods for multiple fields (or, in
		% this case, colony pairs) at a time

	% Takes some important parameters as pre-calculated so speed things
		% up when looping over multiple datapoints whose likelihood
		% needs to be estimated

%	tic;

	datapoint_num = length(strain_GR_diff_list);
	if datapoint_num > 1
		strain_GR_diff_list = reshape(strain_GR_diff_list, [datapoint_num 1]);
	end
%	mean_number = length(test_mean);
%	if mean_number > 1
%		test_mean = reshape(test_mean, [1 mean_number]);
%	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% approach below is about 1.5x faster when number of colonies in each sample = 1
%	tic;
%	lambda_pp = test_petite_prop*ref_petite_prop;
%	lambda_pr = test_petite_prop*(1-ref_petite_prop);
%	lambda_tp = (1-test_petite_prop)*ref_petite_prop;
%	lambda_tr = (1-test_petite_prop)*(1-ref_petite_prop);

%	mean_pp = 0;
%	mean_pr = petite_mean-ref_mean;
%	mean_tp = test_mean-petite_mean;
%	mean_tr = test_mean-ref_mean;

%	sigma_pp = petite_sigma * sqrt(2);	
%	sigma_pr = sqrt(petite_sigma^2 + ref_sigma^2);
%	sigma_tp = sqrt(petite_sigma^2 + test_sigma^2);
%	sigma_tr = sqrt(ref_sigma^2 + test_sigma^2);
	
	norm_pp = normpdf(strain_GR_diff_list, mean_pp, sigma_pp);
	norm_pr = normpdf(strain_GR_diff_list, mean_pr, sigma_pr);
	norm_tp = normpdf(strain_GR_diff_list, mean_tp, sigma_tp);
	norm_tr = normpdf(strain_GR_diff_list, mean_tr, sigma_tr);

	L_pp = lambda_pp * norm_pp;
	L_pr = lambda_pr * norm_pr;
	L_tp = lambda_tp * norm_tp;
	L_tr = lambda_tr * norm_tr;

	GR_diff_likelihood_list = L_pp + L_pr + L_tp + L_tr;
%    runtime_short_version = toc;

	LL = sum(log(GR_diff_likelihood_list), 1);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	unscaled_gradient_dict = containers.Map('KeyType', 'char', ...
		'ValueType', 'any');

	grad_parameter_names = ...
		{'ref_petite_prop', 'test_petite_prop', 'ref_mean', 'petite_mean', ...
			'test_mean', 'ref_sigma', 'petite_sigma', 'test_sigma'};

	% Pre-set all gradients to NaN, and then only calculate the actual
		% values if those gradients are for fitted parameters
	for current_param_idx = 1:length(grad_parameter_names)
		current_param_name = grad_parameter_names{current_param_idx};
		unscaled_gradient_dict(current_param_name) = NaN;
	end

	if gradient_specification

		% 0/0, which can easily occur when dividing the derivative of
			%  the LL list by the LL list, returns NaN; when summing
			% across derivatives of likelihoods, use nansum. In effect,
			% this means values at which both the likelihood and its
			% derivative are 0 won't contribute to the gradient
			% calculation, which is correct

		% d(LL)/d(test_petite_prop)
		if any(strcmp('test_petite_prop',fitted_parameters))
			d_lambda_pp_d_test_petite_prop = ref_petite_prop;
			d_lambda_pr_d_test_petite_prop = 1-ref_petite_prop;
			d_lambda_tp_d_test_petite_prop = -ref_petite_prop;
			d_lambda_tr_d_test_petite_prop = -(1-ref_petite_prop);

			d_L_pp_d_test_petite_prop = ...
				d_lambda_pp_d_test_petite_prop * norm_pp;
			d_L_pr_d_test_petite_prop = ...
				d_lambda_pr_d_test_petite_prop * norm_pr;
			d_L_tp_d_test_petite_prop = ...
				d_lambda_tp_d_test_petite_prop * norm_tp;
			d_L_tr_d_test_petite_prop = ...
				d_lambda_tr_d_test_petite_prop * norm_tr;

			d_GR_diff_likelihood_list_d_test_petite_prop = ...
				d_L_pp_d_test_petite_prop + d_L_pr_d_test_petite_prop + ...
				d_L_tp_d_test_petite_prop + d_L_tr_d_test_petite_prop;

			d_GR_diff_LL_list_d_test_petite_prop = ...
				d_GR_diff_likelihood_list_d_test_petite_prop./...
				GR_diff_likelihood_list;
			d_LL_d_test_petite_prop = ...
				nansum(d_GR_diff_LL_list_d_test_petite_prop, 1);
			unscaled_gradient_dict('test_petite_prop') = ...
				d_LL_d_test_petite_prop;
		end

		% d(LL)/d(ref_petite_prop)
		if any(strcmp('ref_petite_prop',fitted_parameters))
			d_lambda_pp_d_ref_petite_prop = test_petite_prop;
			d_lambda_pr_d_ref_petite_prop = -test_petite_prop;
			d_lambda_tp_d_ref_petite_prop = 1-test_petite_prop;
			d_lambda_tr_d_ref_petite_prop = -(1-test_petite_prop);

			d_L_pp_d_ref_petite_prop = ...
				d_lambda_pp_d_ref_petite_prop * norm_pp;
			d_L_pr_d_ref_petite_prop = ...
				d_lambda_pr_d_ref_petite_prop * norm_pr;
			d_L_tp_d_ref_petite_prop = ...
				d_lambda_tp_d_ref_petite_prop * norm_tp;
			d_L_tr_d_ref_petite_prop = ...
				d_lambda_tr_d_ref_petite_prop * norm_tr;

			d_GR_diff_likelihood_list_d_ref_petite_prop = ...
				d_L_pp_d_ref_petite_prop + d_L_pr_d_ref_petite_prop + ...
				d_L_tp_d_ref_petite_prop + d_L_tr_d_ref_petite_prop;

			d_GR_diff_LL_list_d_ref_petite_prop = ...
				d_GR_diff_likelihood_list_d_ref_petite_prop./...
				GR_diff_likelihood_list;
			d_LL_d_ref_petite_prop = ...
				nansum(d_GR_diff_LL_list_d_ref_petite_prop, 1);
			unscaled_gradient_dict('ref_petite_prop') = ...
				d_LL_d_ref_petite_prop;
		end

		% gradients relatvie to means
		if ~isempty(intersect({'petite_mean','ref_mean','test_mean'}, ...
			fitted_parameters))

			d_L_pr_d_mean_pr = d_LL_d_mu_norm_calc(strain_GR_diff_list, mean_pr, sigma_pr) .* L_pr;
			d_L_tp_d_mean_tp = d_LL_d_mu_norm_calc(strain_GR_diff_list, mean_tp, sigma_tp) .* L_tp;

			% d(LL)/d(petite_mean)
            if any(strcmp('petite_mean',fitted_parameters))
				
            	d_mean_pr_d_petite_mean = 1;
            	d_mean_tp_d_petite_mean = -1;

            	d_L_pr_d_petite_mean = d_L_pr_d_mean_pr .* d_mean_pr_d_petite_mean;
            	d_L_tp_d_petite_mean = d_L_tp_d_mean_tp .* d_mean_tp_d_petite_mean;

				d_GR_diff_likelihood_list_d_petite_mean = ...
					d_L_pr_d_petite_mean + d_L_tp_d_petite_mean;
				d_GR_diff_LL_list_d_petite_mean = ...
					d_GR_diff_likelihood_list_d_petite_mean./...
					GR_diff_likelihood_list;
				d_LL_d_petite_mean = nansum(d_GR_diff_LL_list_d_petite_mean, 1);
				unscaled_gradient_dict('petite_mean') = d_LL_d_petite_mean;

            end
            
            if ~isempty(intersect({'ref_mean','test_mean'}, ...
				fitted_parameters))

            	d_L_tr_d_mean_tr = d_LL_d_mu_norm_calc(strain_GR_diff_list, mean_tr, sigma_tr) .* L_tr;

            	% gradient relative to test_mean contributes to
						% ref_mean gradient if test_mean is specified as an
						% effect on ref_mean, and thus needs to be
						% calculated if ref_mean is fitted, regardless of
						% whether test_mean is fitted
				d_mean_tp_d_test_mean = 1;
				d_mean_tr_d_test_mean = 1;

            	d_L_tp_d_test_mean = d_L_tp_d_mean_tp .* d_mean_tp_d_test_mean;
            	d_L_tr_d_test_mean = d_L_tr_d_mean_tr .* d_mean_tr_d_test_mean;
				
				d_GR_diff_likelihood_list_d_test_mean = ...
					d_L_tp_d_test_mean + d_L_tr_d_test_mean;
				d_GR_diff_LL_list_d_test_mean = ...
					d_GR_diff_likelihood_list_d_test_mean./...
					GR_diff_likelihood_list;
				d_LL_d_test_mean = nansum(d_GR_diff_LL_list_d_test_mean, 1);
				unscaled_gradient_dict('test_mean') = d_LL_d_test_mean;

				% d(LL)/d(ref_mean)
				if any(strcmp('ref_mean',fitted_parameters))
					d_mean_pr_d_ref_mean = -1;
					d_mean_tr_d_ref_mean = -1;

	            	d_L_pr_d_ref_mean = d_L_pr_d_mean_pr * d_mean_pr_d_ref_mean;
	            	d_L_tr_d_ref_mean = d_L_tr_d_mean_tr * d_mean_tr_d_ref_mean;
					
					d_GR_diff_likelihood_list_d_ref_mean = ...
						d_L_pr_d_ref_mean + d_L_tr_d_ref_mean;
					d_GR_diff_LL_list_d_ref_mean = ...
						d_GR_diff_likelihood_list_d_ref_mean./...
						GR_diff_likelihood_list;
					d_LL_d_ref_mean = nansum(d_GR_diff_LL_list_d_ref_mean, 1);
					unscaled_gradient_dict('ref_mean') = d_LL_d_ref_mean;
				end
            end

		end

		% gradients relatvie to sigmas
        if ~isempty(intersect({'petite_sigma','ref_sigma','test_sigma'}, ...
			fitted_parameters))

			d_L_pr_d_sigma_pr = ...
				d_LL_d_sigma_norm_calc(strain_GR_diff_list, mean_pr, ...
					sigma_pr) .* L_pr;
			d_L_tp_d_sigma_tp = ...
				d_LL_d_sigma_norm_calc(strain_GR_diff_list, mean_tp, ...
					sigma_tp) .* L_tp;	

			% d(LL)/d(petite_sigma)
            if any(strcmp('petite_sigma',fitted_parameters))

            	d_L_pp_d_sigma_pp = ...
            		d_LL_d_sigma_norm_calc(strain_GR_diff_list, mean_pp, ...
            			sigma_pp) .* L_pp;

            	d_sigma_pp_d_petite_sigma = sqrt(2);
            	d_sigma_pr_d_petite_sigma = petite_sigma/sigma_pr;
            	d_sigma_tp_d_petite_sigma = petite_sigma/sigma_tp;

            	d_L_pp_d_petite_sigma = ...
            		d_L_pp_d_sigma_pp * d_sigma_pp_d_petite_sigma;
            	d_L_pr_d_petite_sigma = ...
            		d_L_pr_d_sigma_pr * d_sigma_pr_d_petite_sigma;
            	d_L_tp_d_petite_sigma = ...
            		d_L_tp_d_sigma_tp * d_sigma_tp_d_petite_sigma;

				d_GR_diff_likelihood_list_d_petite_sigma = ...
					d_L_pp_d_petite_sigma + ...
					d_L_pr_d_petite_sigma + ...
					d_L_tp_d_petite_sigma;
				d_GR_diff_LL_list_d_petite_sigma = ...
					d_GR_diff_likelihood_list_d_petite_sigma./...
					GR_diff_likelihood_list;
				d_LL_d_petite_sigma = nansum(d_GR_diff_LL_list_d_petite_sigma, 1);
				unscaled_gradient_dict('petite_sigma') = d_LL_d_petite_sigma;
            end
            
            if ~isempty(intersect({'ref_sigma','test_sigma'}, ...
				fitted_parameters))

            	d_L_tr_d_sigma_tr = ...
					d_LL_d_sigma_norm_calc(strain_GR_diff_list, mean_tr, ...
						sigma_tr) .* L_tr;

            	% gradient relative to test_sigma contributes to
					% ref_sigma gradient if test_sigma is specified as
					% an effect on ref_sigma, or if both are actually
					% equal to a general nonpetite sigma, and thus needs
					% to be calculated if ref_mean is fitted, regardless
					% of whether test_mean is fitted
				% d(LL)/d(test_sigma)
				d_sigma_tr_d_test_sigma = test_sigma/sigma_tr;
            	d_sigma_tp_d_test_sigma = test_sigma/sigma_tp;

            	d_L_tr_d_test_sigma = ...
            		d_L_tr_d_sigma_tr * d_sigma_tr_d_test_sigma;
            	d_L_tp_d_test_sigma = ...
            		d_L_tp_d_sigma_tp * d_sigma_tp_d_test_sigma;

				d_GR_diff_likelihood_list_d_test_sigma = ...
					d_L_tr_d_test_sigma + ...
					d_L_tp_d_test_sigma;
				d_GR_diff_LL_list_d_test_sigma = ...
					d_GR_diff_likelihood_list_d_test_sigma./...
					GR_diff_likelihood_list;
				d_LL_d_test_sigma = nansum(d_GR_diff_LL_list_d_test_sigma, 1);
				unscaled_gradient_dict('test_sigma') = d_LL_d_test_sigma;

				% d(LL)/d(ref_sigma)
				if any(strcmp('ref_sigma',fitted_parameters))
					d_sigma_tr_d_ref_sigma = ref_sigma/sigma_tr;
	            	d_sigma_pr_d_ref_sigma = ref_sigma/sigma_pr;

	            	d_L_tr_d_ref_sigma = ...
	            		d_L_tr_d_sigma_tr * d_sigma_tr_d_ref_sigma;
	            	d_L_pr_d_ref_sigma = ...
	            		d_L_pr_d_sigma_pr * d_sigma_pr_d_ref_sigma;

					d_GR_diff_likelihood_list_d_ref_sigma = ...
						d_L_tr_d_ref_sigma + ...
						d_L_pr_d_ref_sigma;
					d_GR_diff_LL_list_d_ref_sigma = ...
						d_GR_diff_likelihood_list_d_ref_sigma./...
						GR_diff_likelihood_list;
					d_LL_d_ref_sigma = nansum(d_GR_diff_LL_list_d_ref_sigma, 1);
					unscaled_gradient_dict('ref_sigma') = d_LL_d_ref_sigma;
				end

            end
        end

	end

%	disp(runtime_long_version);
%	disp(runtime_short_version);


end