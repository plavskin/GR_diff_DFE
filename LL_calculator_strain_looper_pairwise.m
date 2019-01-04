function [combined_LL, unscaled_gradient_vector, grad_parameter_names] = ...
    LL_calculator_strain_looper_pairwise(param_vals,...
    input_value_dict, pre_MLE_output_dict)
    
    % EP 17-11-07

    % Takes in parameters that apply across all test strains
    % Loops through test strains to estimate the optimal strain parameters
        % given the current iteration of general parameter values
    % For each test strain, calculates the log likelihood of observing a list
        % of differences between test and reference strain GRs within
        % pairs of colonies
    % test strain GR = (ref strain GR)*exp(test strain mut_effect)
    % The input into this function can hold one of the test strain mut_effects or
        % petite proportions fixed
    % This function writes an output csv file which contains the MLE parameters
        % corresponding to each strain if the overall LL is higher than previous
        % iterations of this function

    % V2 - updated how the function deals with fixed parameters; it now takes a
        % list of indices and values for fixed mutation effects and petite
        % proportions, the same size as the list of test strains; each position
        % corresponds to a position in strain_list, and NaN indicates that that
        % parameter is not fixed
    % V3 - uses different sigma values for petite and non-petite colonies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    strainwise_search_type = input_value_dict('strainwise_search_type');
    tolx_val = input_value_dict('x_tolerance');
    tolfun_val = input_value_dict('fun_tolerance');
    max_neg_LL_val = input_value_dict('max_neg_LL_val');
    random_effect_names = input_value_dict('random_effect_names');
    global_mle_parameter_names = input_value_dict('mle_parameter_names');
    
    strain_list = pre_MLE_output_dict('strain_list');
    test_strain_list_by_pair = pre_MLE_output_dict('test_strain_list_by_pair');
    GR_diff_list = pre_MLE_output_dict('GR_diff_list');
    test_strain_ML_file = pre_MLE_output_dict('test_strain_ML_file');
    strain_LL_table_file = pre_MLE_output_dict('strain_LL_table_file');

    parameter_dict = containers.Map(global_mle_parameter_names,param_vals);
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

    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate likelihood of observing petite data given current global parameters

    % replace parameter names in parameter list supplied to LL_mixef_calculator
    input_value_dict_mixef = ...
        containers.Map(keys(input_value_dict),values(input_value_dict));
    parameter_list_mixef = global_mle_parameter_names;

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop through test strains and identify ML test strain parameters for each
    % The sum of LL across all strains is the combined log likelihood for the
        % data, given the current parameter values
    % Because test strain-specific parameters are independent of each other in
        % terms of their effect on likelihood, we can calculate likelihood for 
        % one test strain at a time (given reference strain and general
        % parameters), greatly reducing the dimensionality of the test strain
        % parameter search

    % don't calculate global parameter gradients until ML parameters identified
        % for each strain
    return_all_grads = false;

    strain_tolx_val = tolx_val/(test_strain_number*1.1);
    strain_tolfun_val = tolfun_val/(test_strain_number*1.1);

    if strcmp(strainwise_search_type,'global')
    % set up structure for global solution search

        gs = GlobalSearch;
        gs.TolFun = strain_tolfun_val;
        gs.TolX = strain_tolx_val;
        gs.NumStageOnePoints = 2^2; % 9
        gs.NumTrialPoints = 3^2; % 90
        gs.StartPointsToRun='bounds';
        gs.Display = 'off';
    end
    
    %fmincon_opts = optimset('TolX',10^-7,'TolFun',10^-7,'Algorithm','interior-point','MaxIter',5000,'MaxFunEvals',12000,'UseParallel',true,'GradObj','on','GradConstr','on');
    % 2016 and beyond version:
    fmincon_opts = optimoptions('fmincon','TolX',strain_tolx_val,'TolFun',strain_tolfun_val,...
        'Algorithm','interior-point','MaxIter',5000,'MaxFunEvals',12000,...
        'SpecifyObjectiveGradient',true,'CheckGradients',false,'Display','off');

    current_iter_parameter_values = [petite_colony_sigma,nonpetite_colony_sigma,petite_mean,...
        ref_mean,ref_petite_prop];
        % parameter values fixed over this iteration

    % Initialize a matrix for storing test strain mut_effects and petite proportions
    test_strain_ML_params = NaN(test_strain_number,2);

    % global parameters have already been descaled and de-logspaced, so
        % use these arrays when checking whether to descale them or
        % de-logspace them in later LL functions
    global_scaling_for_strain_mle = ones(size(current_iter_parameter_values));
    global_logspace_for_strain_mle = false(size(current_iter_parameter_values));

    % initialize the negative log likelihood, which is minimized in
        % iterations of this function
    % the log likelihood here is the max LL across all the MA strains
        % plus the LL of the petite data above
        
    combined_LL = 0 + LL_petite;

    for strain_idx = 1:test_strain_number

        current_strain = strain_list{strain_idx};
        % set up list of data to pass to likelihood estimation function
        % find which data in list corresponds to current_strain
        current_indices = find(strcmp(test_strain_list_by_pair,current_strain));
        % only run mle/ll calc if there are samples from current_strain
            % (this is an issue when working with a random subset of data)
            
        if size(current_indices, 1) > 0
            
            strain_current_list = {current_strain};
            test_strain_current_list = test_strain_list_by_pair(current_indices);
            GR_diff_current_list = GR_diff_list(current_indices);
            
            strain_param_current_list = [strcat(strain_current_list,'_me'), ...
                strcat(strain_current_list,'_pp')];
            
            [~, ~, current_strain_param_logspace_bool, ...
                current_strain_param_scaling_values, current_strain_fixed_values, current_strain_fixed_indices, ...
                current_lb_values, current_ub_values, current_start_values] = ...
                parameter_array_subsetter(strain_param_current_list, input_value_dict);

            current_scaling_values = [global_scaling_for_strain_mle,...
                current_strain_param_scaling_values];
            current_logspace_bool = [global_logspace_for_strain_mle,...
                current_strain_param_logspace_bool];

            % 'global' parameters provided to this function are also fixed
            fixed_parameter_indices = [true(size(current_iter_parameter_values)),...
                current_strain_fixed_indices];
            fixed_parameter_values = [current_iter_parameter_values,...
                current_strain_fixed_values];

            % Estimate likelihood for all pairs of current strain
            % If there are no parameters to fit, just run the likelihood estimation

            if sum(~fixed_parameter_indices) > 0
                current_min_problem = createOptimProblem('fmincon','objective',...
                    @(strainwise_parameter_vals_partial) ...
                    LL_calculator_strains_pairwise_2_sigmas(strainwise_parameter_vals_partial,...
                        fixed_parameter_indices, fixed_parameter_values,...
                        strain_current_list, test_strain_current_list, GR_diff_current_list,...
                        return_all_grads, max_neg_LL_val, current_scaling_values,...
                        current_logspace_bool),'x0',current_start_values,...
                    'lb',current_lb_values,'ub',current_ub_values,'options',fmincon_opts);

                if strcmp(strainwise_search_type,'global')
                    [current_out_param_vals, neg_LL_current]=...
                        run(gs,current_min_problem);
                else
                    [current_out_param_vals, neg_LL_current]=...
                        fmincon(current_min_problem);
                end
            else
                neg_LL_current = LL_calculator_strains_pairwise_2_sigmas(NaN,...
                    fixed_parameter_indices,fixed_parameter_values,...
                    strain_current_list,test_strain_current_list,GR_diff_current_list,...
                    return_all_grads,max_neg_LL_val,current_scaling_values,...
                        current_logspace_bool);
                current_out_param_vals = NaN;
            end

            combined_LL = combined_LL - neg_LL_current;

            corrected_current_param_vals = NaN(1,2);
            corrected_current_param_vals(current_strain_fixed_indices) = ...
                current_strain_fixed_values(current_strain_fixed_indices);
            corrected_current_param_vals(~current_strain_fixed_indices) = ...
                current_out_param_vals;
            % record output parameter values before unscaling them
            unscaled_current_param_vals = ...
                reverse_value_scaler(corrected_current_param_vals, ...
                current_strain_param_logspace_bool, ...
                current_strain_param_scaling_values);
            
            test_strain_ML_params(strain_idx,:) = unscaled_current_param_vals;
        end

    end

    % calculate global likelihood gradient given current global
        % parameters and current ML estimates for strain-specific
        % parameters
    return_all_grads = true;
    fixed_parameter_values_grad_calc = [petite_colony_sigma,nonpetite_colony_sigma,petite_mean,...
        ref_mean,ref_petite_prop,test_strain_ML_params(:,1)',...
        test_strain_ML_params(:,2)'];
    scaling_array_grad_calc = ones(size(fixed_parameter_values_grad_calc));
    logspace_array_grad_calc = false(size(fixed_parameter_values_grad_calc));
    fixed_parameter_indices_grad_calc = true(size(fixed_parameter_values_grad_calc));
    [~,neg_strain_param_gradient_vector,neg_global_gradient_vector_no_petite_params] = ...
        LL_calculator_strains_pairwise_2_sigmas(NaN,...
                fixed_parameter_indices_grad_calc,fixed_parameter_values_grad_calc,...
                strain_list,test_strain_list_by_pair,GR_diff_list,return_all_grads,...
                max_neg_LL_val,scaling_array_grad_calc,logspace_array_grad_calc);
    global_gradient_vector_no_petite_params = ...
        -[neg_global_gradient_vector_no_petite_params, ...
        neg_strain_param_gradient_vector];
        % In this case, the strain_gradient_vector output is blank,
            % which happens when all strain parameters are fixed in
            % LL_calculator_strains_pairwise
    
    strain_param_names = [strcat(strain_list,'_me'), ...
        strcat(strain_list,'_pp')];
    strain_param_gradient = global_gradient_vector_no_petite_params(6:end);
    unscaled_gradient_dict = containers.Map(strain_param_names, strain_param_gradient);
    % d_LL_d_nonpetite_colony_sigma
    unscaled_gradient_dict('nonpetite_colony_sigma') = global_gradient_vector_no_petite_params(2);
    % d_LL_d_petite_mean
    unscaled_gradient_dict('petite_mean') = ...
        petite_param_gradient_dict('petite') + ...
        global_gradient_vector_no_petite_params(3);
    % d_LL_d_ref_mean
    unscaled_gradient_dict('ref_mean') = global_gradient_vector_no_petite_params(4);
    % d_LL_d_ref_petite_prop
    unscaled_gradient_dict('ref_petite_prop') = global_gradient_vector_no_petite_params(5);
    % d_LL_d_ranef
    ranef_gradient_vals = values(petite_param_gradient_dict, ...
        random_effect_names);
    ranef_gradient_dict = containers.Map(ranef_names_from_param_list, ...
        ranef_gradient_vals);
    
    unscaled_gradient_dict_complete = ...
        [unscaled_gradient_dict; ranef_gradient_dict];

    % d_LL_d_petite_colony_sigma
    unscaled_gradient_dict_complete('petite_colony_sigma') = ...
        petite_param_gradient_dict('petite_colony') + ...
        global_gradient_vector_no_petite_params(1);
    
    unscaled_gradient_vector = cell2mat(values(unscaled_gradient_dict_complete));
    grad_parameter_names = keys(unscaled_gradient_dict_complete);
    
    runtime = toc;

    % need to first read previous file and see whether current combined_LL is higher than one in prev file

    % default is to write output to a file, but an suppress this if previous
        % run of MLE already produced a better log likelihood

    write_output = true;

    if exist(strain_LL_table_file,'file') == 2
        previous_LL_table = readtable(strain_LL_table_file);

        if combined_LL < previous_LL_table.Current_Best_LL(1)
            write_output = false;
        end
    end

    if write_output
        Strain_Names = {strain_list{:}}';
        Mutation_Effects = test_strain_ML_params(:,1);
        Petite_Proportions = test_strain_ML_params(:,2);
        export_table = table(Strain_Names,Mutation_Effects,Petite_Proportions);

        Current_Best_LL = combined_LL;
        Best_LL_Runtime = runtime;
        strain_LL_table = table(Current_Best_LL, Best_LL_Runtime);

        writetable(export_table,test_strain_ML_file);
        writetable(strain_LL_table,strain_LL_table_file);

    end

end
    
