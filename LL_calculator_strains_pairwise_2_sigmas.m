function [neg_LL_val,neg_gradient_vector_partial_strain,neg_gradient_vector_global]=LL_calculator_strains_pairwise_2_sigmas(...
    strainwise_parameter_vals_partial,fixed_parameter_indices,fixed_parameter_values,...
    strain_current_list,test_strain_current_list,GR_diff_current_list,return_all_grads,max_neg_LL_val,...
    current_scaling_values,current_logspace_bool)
	% EP 17-11-07

    % calculates the log likelihood of observing a list of differences between
        % pairs of test and reference strain colony GRs
    % can work with 1+ test strains, 1+ pairs/strain
    % assumes different GR variance for petite and non-petite colonies,
        % but same GR variance for ref and all test strains

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % assign current data to cell arrays

%    well_current_list = current_data_table{1};
%    plate_current_list = current_data_table{2};
%    ref_strain_current_list = current_data_table{3};
%    ref_count_current_list = current_data_table{4};
%    test_strain_current_list = current_data_table{5};
%    test_count_current_list = current_data_table{6};
%    GR_diff_current_list = current_data_table{7};
%    field_count_var_current_list = current_data_table{8};

%    unique_test_strains = unique(strain_current_list);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify parameter values

    % some parameters may be 'fixed' (i.e. not fitted by current iteration of
        % MLE); these are provided to the function in a separate list
    % compile a list of all parameters, fixed or fitted, and identify values
        % belonging to each individual parameter using that list
    param_vals = NaN(size(fixed_parameter_indices));
    param_vals(fixed_parameter_indices) = fixed_parameter_values(fixed_parameter_indices);
    param_vals(~fixed_parameter_indices) = strainwise_parameter_vals_partial;
    param_vals = reverse_value_scaler(param_vals,current_logspace_bool,current_scaling_values);

    petite_sigma = param_vals(1);
        % s.d. of colony GRs of petite distribution
    nonpetite_sigma = param_vals(2);
        % s.d. of colony GRs of non-petite reference and test strain distribution
    petite_mean = param_vals(3);
        % mean growth rate of petite colonies, regardless of genotype
    ref_mean = param_vals(4);
        % mean growth rate of non-petite ref colonies
    ref_petite_prop = param_vals(5);
        % proportion of petites in ref strain

    ref_sigma = nonpetite_sigma;
    test_sigma = nonpetite_sigma;

    global_param_number = 5;
    test_strain_number = length(strain_current_list);
        % number of strains besides reference strain whose likelihood being fitted
    global_param_indices = 1:global_param_number;
    mut_effect_indices = (global_param_number+1):(global_param_number+test_strain_number);
    petite_prop_indices = (global_param_number+test_strain_number+1):(global_param_number+2*test_strain_number);
    non_global_indices = [mut_effect_indices,petite_prop_indices];

    test_strain_mut_effects = param_vals(mut_effect_indices);
    test_strain_means = ref_mean*exp(test_strain_mut_effects);
        % mean growth rates of non-petite colonies for each test strain
    test_petite_prop_list = param_vals(petite_prop_indices);
        % proportion of petites in each test strain

%    NaNs_in_pp = sum(isnan(test_petite_prop_list));
%    if NaNs_in_pp > 0 
%        disp('warning')
%        NaNs_in_pp
%        param_vals
%        return
%    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop through test strains; within every test strain, loop through fields,
        % and calculate their log likelihood

    % initialize the log likelihood
    LL_val = 0;
    if nargout > 1
        gradient_vector_strain = zeros([1,test_strain_number*2]);

        if return_all_grads
            gradient_vector_global = [0,0,0,0,0];
        end

    end
    
    for strain_idx = 1:test_strain_number
        current_strain = strain_current_list{strain_idx};
        current_indices = find(strcmp(test_strain_current_list,current_strain));
        % only run mle/ll calc if there are samples from current_strain
            % (this is an issue when working with a random subset of data)
        if size(current_indices, 1) > 0

            test_petite_prop = test_petite_prop_list(strain_idx);
            test_mean = test_strain_means(strain_idx);

            strain_GR_diff_list = GR_diff_current_list(current_indices);

            if nargout > 1
    %            test_mean
    %            ref_mean
    %            test_petite_prop
    %            ref_petite_prop
    %            petite_mean
    %            sigma_colony
    %            %strain_GR_diff_list

                if return_all_grads
                    [current_LL,current_strain_gradient,current_global_gradient] = ...
                        LL_calculator_within_pair_different_sigmas(test_mean,ref_mean,...
                            test_petite_prop,ref_petite_prop,petite_mean,test_sigma,...
                            ref_sigma,petite_sigma,strain_GR_diff_list,return_all_grads);
                    % account for the fact that test_mean is a function of ref_mean
                    current_contrib_to_ref_mean_grad = current_strain_gradient(1)*test_mean/ref_mean;
                    current_global_gradient(4) = current_global_gradient(4)+current_contrib_to_ref_mean_grad;
                    % add test strain sigma gradient to nonpetite sigma gradient
                    current_global_gradient(2) = current_global_gradient(2)+current_strain_gradient(3);

                    gradient_vector_global = gradient_vector_global+current_global_gradient;

                else
                    [current_LL,current_strain_gradient] = ...
                        LL_calculator_within_pair_different_sigmas(test_mean,ref_mean,test_petite_prop,...
                            ref_petite_prop,petite_mean,test_sigma,ref_sigma,...
                            petite_sigma,strain_GR_diff_list,return_all_grads);
                end

                % convert derivative of test strain mean to derivative of
                    % mutational effect
                current_strain_gradient(1) = current_strain_gradient(1)*ref_mean*exp(test_strain_mut_effects(strain_idx));
                gradient_vector_strain(strain_idx) = current_strain_gradient(1);
                gradient_vector_strain(test_strain_number+strain_idx) = current_strain_gradient(2);
            else
                current_LL = LL_calculator_within_pair_different_sigmas(test_mean,ref_mean,...
                    test_petite_prop,ref_petite_prop,petite_mean,test_sigma,...
                    ref_sigma,petite_sigma,strain_GR_diff_list,return_all_grads);
            end

            LL_val = LL_val+current_LL;
        end

    end

    neg_LL_val = -LL_val;

    if neg_LL_val > max_neg_LL_val
        neg_LL_val = max_neg_LL_val;
    end

    if nargout > 1
        % rescale gradient vector back to space being used by MLE
        gradient_vector_strain_scaled = ...
            gradient_value_rescaler(gradient_vector_strain,...
            param_vals(non_global_indices),...
            current_logspace_bool(non_global_indices),...
            current_scaling_values(non_global_indices));

        if return_all_grads
            gradient_vector_partial_strain = gradient_vector_strain_scaled;
            % rescale gradient vector back to space being used by MLE
            gradient_vector_global_scaled = ...
                gradient_value_rescaler(gradient_vector_global,...
                param_vals(global_param_indices),...
                current_logspace_bool(global_param_indices),...
                current_scaling_values(global_param_indices));
            neg_gradient_vector_global = -gradient_vector_global_scaled;
                % account for fixed parameters outside of this function
            neg_gradient_vector_global(neg_gradient_vector_global>max_neg_LL_val) = max_neg_LL_val;
            neg_gradient_vector_global(neg_gradient_vector_global<-max_neg_LL_val) = -max_neg_LL_val;
        else
            gradient_vector_partial_strain = ...
            gradient_vector_strain_scaled(~fixed_parameter_indices(non_global_indices));
        end
        neg_gradient_vector_partial_strain = -gradient_vector_partial_strain;

        neg_gradient_vector_partial_strain(neg_gradient_vector_partial_strain>max_neg_LL_val) = max_neg_LL_val;
        neg_gradient_vector_partial_strain(neg_gradient_vector_partial_strain<-max_neg_LL_val) = -max_neg_LL_val;

    end
end