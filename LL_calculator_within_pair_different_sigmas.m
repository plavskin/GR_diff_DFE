function [LL,gradient_vector_strain,gradient_vector_global]=LL_calculator_within_pair_diff_sigmas(test_mean,...
	ref_mean,test_petite_prop,ref_petite_prop,petite_mean,test_sigma,ref_sigma,petite_sigma,...
	strain_GR_diff_list,return_all_grads)
	% EP 17-08-21

	% Calculates the log likelihood of the observed difference between
		% reference and test strain colony pair growth rates, accounting
		% for random proportion of petites

	% Equivalent to LL_calculator_within_field, but assumes
		% a single colony per strain
	% The reason the code is different is because LL_calculator_within_field
		% isn't set up to calculate likelihoods for multiple fields (or, in
		% this case, colony pairs) at a time

%	tic;

	strain_GR_diff_list = reshape(strain_GR_diff_list,[1 length(strain_GR_diff_list)]);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% As in LL_calculator_within_field, but the number of reference
		% and test colonies is assumed to be 1 each
	ref_colony_num = 1;
	test_colony_num = 1;

	% Create a vector of possible number of petite and non-petite colony
		% numbers for reference and test strain
	ref_petite_colony_nums = 0:ref_colony_num;
	test_petite_colony_nums = 0:test_colony_num;
	ref_nonpetite_colony_nums = ref_colony_num-ref_petite_colony_nums;
	test_nonpetite_colony_nums = test_colony_num-test_petite_colony_nums;

	% Create matrices of the propoprtion of petite to non-petite
		% colonies for the reference and test strains within each
		% subfield (for a pair, this is [1 0; 1 0] or [1 1; 0 0]
		% for the ref and test strains, respectively)
	[ref_petite_colony_num_matrix,test_petite_colony_num_matrix] = ...
		meshgrid(ref_petite_colony_nums,test_petite_colony_nums);

	ref_nonpetite_colony_num_matrix = ref_colony_num - ref_petite_colony_num_matrix;
	test_nonpetite_colony_num_matrix = test_colony_num - test_petite_colony_num_matrix;

	ref_petite_colony_prop_matrix = ...
		ref_petite_colony_num_matrix/ref_colony_num;
	test_petite_colony_prop_matrix = ...
		test_petite_colony_num_matrix/test_colony_num;

	ref_nonpetite_colony_prop_matrix = 1 - ref_petite_colony_prop_matrix;
	test_nonpetite_colony_prop_matrix = 1 - test_petite_colony_prop_matrix;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Create a matrix with the expected difference in mean GRs between
		% the test and reference colonies within each pair, for each
		% combination of possible petite colony numbers for those two
		% strains

	% For each possible number of petites, calculate the mean growth
		% rate of the reference and test colonies in the pair
	ref_pair_mean_matrix = ...
		petite_mean*ref_petite_colony_prop_matrix + ...
		ref_mean*ref_nonpetite_colony_prop_matrix;
	test_pair_mean_matrix = ...
		petite_mean*test_petite_colony_prop_matrix + ...
		test_mean*test_nonpetite_colony_prop_matrix;

	% calculate the expected difference in growth rate for each
		% combination of petite colony number within test and ref strains
	pair_mean_diff_matrix = test_pair_mean_matrix-ref_pair_mean_matrix;

	% Calculate variance of mean growth rate estimates within the pair
		% for each strain, for each possible number of petites

	ref_pair_var_matrix = ...
		petite_sigma^2*ref_petite_colony_prop_matrix + ...
		ref_sigma^2*ref_nonpetite_colony_prop_matrix;
	test_pair_var_matrix = ...
		petite_sigma^2*test_petite_colony_prop_matrix + ...
		ref_sigma^2*test_nonpetite_colony_prop_matrix;

	% calculate the expected sd in growth rate differences for each
		% combination of petite colony number within test and ref strains
	pair_sd_diff_matrix = sqrt(test_pair_var_matrix+ref_pair_var_matrix);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Create a matrix with the probability of each possible petite
		% colony number combination for the mean and reference strain

	ref_petite_num_probabilities = binopdf(ref_petite_colony_nums,...
		ref_colony_num,ref_petite_prop);
	test_petite_num_probabilities = binopdf(test_petite_colony_nums,...
		test_colony_num,test_petite_prop);

	petite_num_probability_matrix = ...
		test_petite_num_probabilities'*ref_petite_num_probabilities;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Create a matrix in which every column is a list of likelihoods of
		% every observed difference contingent on a combination of number
		% of colonies coming from each of the two distributions (petite or 
		% non-petite) for each strain
	% For a single colony from each strain, these combinations are simply
		% petite reference and test colonies, non-petite reference and test
		% colonies, and either petite ref and non-petite test or vice versa
	% The likelihoods in the columns of this matrix are contingent on the 

	pair_mean_diff_vector = pair_mean_diff_matrix(:)';
	pair_sd_diff_vector = pair_sd_diff_matrix(:)';

	contingent_GR_diff_probability_matrix = normpdf(strain_GR_diff_list',...
		pair_mean_diff_vector,pair_sd_diff_vector);

	petite_num_probability_vector = petite_num_probability_matrix(:);
		% probability of each possible petite colony number combination
			% for the mean and reference strain

	GR_diff_likelihood_list = contingent_GR_diff_probability_matrix*...
		petite_num_probability_vector;

	if nargout > 1

		% Create matrices for ref and test colonies corresponding to
			% numbers of non-petites from each strain in each postion
			% in mean matrices

		d_test_petite_num_probabilities_d_test_petite_prop = ...
			1/(test_petite_prop*(1-test_petite_prop))*...
			test_petite_num_probabilities.*...
			((1-test_petite_prop)*test_petite_colony_nums-...
				test_petite_prop*test_nonpetite_colony_nums);

		d_petite_num_probability_matrix_d_test_petite_prop = ...
			d_test_petite_num_probabilities_d_test_petite_prop'*...
			ref_petite_num_probabilities;

		[mean_diff_across_subfields_matrix,subfield_GR_diff_matrix] = ...
			meshgrid(pair_mean_diff_vector,strain_GR_diff_list);
		variance_matrix = repmat(pair_sd_diff_vector.^2,[length(strain_GR_diff_list),1]);
		distance_matrix = subfield_GR_diff_matrix-mean_diff_across_subfields_matrix;
        var_scaled_distance_matrix = ...
			distance_matrix./variance_matrix;
		% d_contingent_GR_diff_probability_matrix_d_mean_diff is not the
			% true derivative of d_contingent_GR_diff_probability_matrix,
			% it still needs to be scaled by petite_num_probability_vector
			% (which happens 2 steps below)
        
        % d(GR diff contingent on # of petite colonies of ref and test
        	% strain for each possible number of petite colonies) / 
			% d(diff between means of ref and test colonies in subfield)
		d_contingent_GR_diff_probability_matrix_d_mean_diff = ...
			contingent_GR_diff_probability_matrix.*...
			var_scaled_distance_matrix;
		% d(diff between means of ref and test colonies in subfield) / 
			% d(test strain mean GR)	
		d_mean_diff_d_test_mean = test_nonpetite_colony_prop_matrix(:);
		d_GR_diff_likelihood_list_d_test_mean = ...
			d_contingent_GR_diff_probability_matrix_d_mean_diff*...
			(petite_num_probability_vector.*...
				d_mean_diff_d_test_mean);
		d_GR_diff_LL_list_d_test_mean = ...
			d_GR_diff_likelihood_list_d_test_mean./...
			GR_diff_likelihood_list;

		d_petite_num_probability_vector_d_test_petite_prop = ...
			d_petite_num_probability_matrix_d_test_petite_prop(:);
		d_GR_diff_likelihood_list_d_test_petite_prop = ...
			(contingent_GR_diff_probability_matrix*...
			d_petite_num_probability_vector_d_test_petite_prop);
		d_GR_diff_LL_list_d_test_petite_prop = ...
			d_GR_diff_likelihood_list_d_test_petite_prop./...
			GR_diff_likelihood_list;

		% calculate gradient relative to test sigma
		d_pair_sd_diff_matrix_d_test_sigma = ...
				test_sigma*(test_nonpetite_colony_prop_matrix)./...
				pair_sd_diff_matrix;

		sd_matrix = sqrt(variance_matrix);
		sd_modified_distance_matrix = ...
			(distance_matrix-sd_matrix).*...
			(distance_matrix+sd_matrix)./sd_matrix.^3;
		% d_contingent_GR_diff_probability_matrix_d_pair_diff_sd is not the
			% true derivative of d_contingent_GR_diff_probability_matrix,
			% it still needs to be scaled by petite_num_probability_vector
			% (which happens 2 steps below)
		d_contingent_GR_diff_probability_matrix_d_pair_diff_sd = ...
			contingent_GR_diff_probability_matrix.*...
			sd_modified_distance_matrix;
		d_GR_diff_likelihood_list_d_test_sigma = ...
			d_contingent_GR_diff_probability_matrix_d_pair_diff_sd*...
			(petite_num_probability_vector.*...
				d_pair_sd_diff_matrix_d_test_sigma(:));
		d_GR_diff_LL_list_d_test_sigma = ...
			d_GR_diff_likelihood_list_d_test_sigma./...
			GR_diff_likelihood_list;

		if return_all_grads
			% Calculate and return gradients for every parameter, not
				% just mutation effect and test strain petite proportion

			% gradient relative to ref_petite_prop
			d_ref_petite_num_probabilities_d_ref_petite_prop = ...
				1/(ref_petite_prop*(1-ref_petite_prop))*...
				ref_petite_num_probabilities.*...
				((1-ref_petite_prop)*ref_petite_colony_nums-...
					ref_petite_prop*ref_nonpetite_colony_nums);

			d_petite_num_probability_matrix_d_ref_petite_prop = ...
				test_petite_num_probabilities'*...
				d_ref_petite_num_probabilities_d_ref_petite_prop;

			d_petite_num_probability_vector_d_ref_petite_prop = ...
				d_petite_num_probability_matrix_d_ref_petite_prop(:);
			d_GR_diff_likelihood_list_d_ref_petite_prop = ...
				(contingent_GR_diff_probability_matrix*...
				d_petite_num_probability_vector_d_ref_petite_prop);
			d_GR_diff_LL_list_d_ref_petite_prop = ...
				d_GR_diff_likelihood_list_d_ref_petite_prop./...
				GR_diff_likelihood_list;

			% gradient relative to petite_mean
			d_mean_diff_d_petite_mean = ...
				(test_petite_colony_prop_matrix(:)-...
						ref_petite_colony_prop_matrix(:));
			d_GR_diff_likelihood_list_d_petite_mean = ...
				d_contingent_GR_diff_probability_matrix_d_mean_diff*...
				(petite_num_probability_vector.*...
					d_mean_diff_d_petite_mean);
			d_GR_diff_LL_list_d_petite_mean = ...
				d_GR_diff_likelihood_list_d_petite_mean./...
				GR_diff_likelihood_list;

			% gradient relative to ref_mean
%			d_mean_diff_d_ref_mean = ...
%				((test_mean/ref_mean)*test_nonpetite_colony_prop_matrix(:)-...
%						ref_nonpetite_colony_prop_matrix(:));
			d_mean_diff_d_ref_mean = ...
				-ref_nonpetite_colony_prop_matrix(:);
			d_GR_diff_likelihood_list_d_ref_mean = ...
				d_contingent_GR_diff_probability_matrix_d_mean_diff*...
				(petite_num_probability_vector.*...
					d_mean_diff_d_ref_mean);
			d_GR_diff_LL_list_d_ref_mean = ...
				d_GR_diff_likelihood_list_d_ref_mean./...
				GR_diff_likelihood_list;

			% gradient relative to ref and petite sigmas
			d_pair_sd_diff_matrix_d_petite_sigma = ...
				petite_sigma*...
				(ref_petite_colony_prop_matrix+test_petite_colony_prop_matrix)./...
				pair_sd_diff_matrix;
			d_pair_sd_diff_matrix_d_ref_sigma = ...
				ref_sigma*(ref_nonpetite_colony_prop_matrix)./...
				pair_sd_diff_matrix;

			d_GR_diff_likelihood_list_d_ref_sigma = ...
				d_contingent_GR_diff_probability_matrix_d_pair_diff_sd*...
				(petite_num_probability_vector.*...
					d_pair_sd_diff_matrix_d_ref_sigma(:));
			d_GR_diff_LL_list_d_ref_sigma = ...
				d_GR_diff_likelihood_list_d_ref_sigma./...
				GR_diff_likelihood_list;
			d_GR_diff_likelihood_list_d_petite_sigma = ...
				d_contingent_GR_diff_probability_matrix_d_pair_diff_sd*...
				(petite_num_probability_vector.*...
					d_pair_sd_diff_matrix_d_petite_sigma(:));
			d_GR_diff_LL_list_d_petite_sigma = ...
				d_GR_diff_likelihood_list_d_petite_sigma./...
				GR_diff_likelihood_list;

		end

	end

%	runtime_long_version = toc;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% approach below is about 1.5x faster when number of colonies in each sample = 1
%	tic;
%	likelihood_list = test_petite_prop*ref_petite_prop*(normpdf(strain_GR_diff_list,0,sigma_colony*sqrt(2)))+...
%    	test_petite_prop*(1-ref_petite_prop)*(normpdf(strain_GR_diff_list,(petite_mean-ref_mean),sigma_colony*sqrt(2)))+...
%    	(1-test_petite_prop)*ref_petite_prop*(normpdf(strain_GR_diff_list,(test_mean-petite_mean),sigma_colony*sqrt(2)))+...
%    	(1-test_petite_prop)*(1-ref_petite_prop)*(normpdf(strain_GR_diff_list,(test_mean-ref_mean),sigma_colony*sqrt(2)));
%    runtime_short_version = toc;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	LL = sum(log(GR_diff_likelihood_list));
%	LL = sum(log(likelihood_list))

	if nargout > 1

		% 0/0, which can easily occur when dividing the derivative of
			%  the LL list by the LL list, returns NaN; when summing
			% across derivatives of likelihoods, use nansum. In effect,
			% this means values at which both the likelihood and its
			% derivative are 0 won't contribute to the gradient
			% calculation, which is correct
		d_LL_d_test_petite_prop = nansum(d_GR_diff_LL_list_d_test_petite_prop);
		d_LL_d_test_mean = nansum(d_GR_diff_LL_list_d_test_mean);
		d_LL_d_test_sigma = nansum(d_GR_diff_LL_list_d_test_sigma);

		gradient_vector_strain = [d_LL_d_test_mean, d_LL_d_test_petite_prop, d_LL_d_test_sigma];

		if return_all_grads
			d_LL_d_ref_petite_prop = nansum(d_GR_diff_LL_list_d_ref_petite_prop);
			d_LL_d_petite_mean = nansum(d_GR_diff_LL_list_d_petite_mean);
			d_LL_d_ref_mean = nansum(d_GR_diff_LL_list_d_ref_mean);
			d_LL_d_ref_sigma = nansum(d_GR_diff_LL_list_d_ref_sigma);
			d_LL_d_petite_sigma = nansum(d_GR_diff_LL_list_d_petite_sigma);

			gradient_vector_global = [d_LL_d_petite_sigma,d_LL_d_ref_sigma,d_LL_d_petite_mean,d_LL_d_ref_mean,d_LL_d_ref_petite_prop];
			
		end

	end

%	disp(runtime_long_version);
%	disp(runtime_short_version);


end