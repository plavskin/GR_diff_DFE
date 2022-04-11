function [combined_LL, unscaled_gradient_vector, grad_parameter_names] = ...
    LL_calculator_strain_looper_pairwise_multiref_with_mixef(param_vals,...
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    strainwise_search_type = input_value_dict('strainwise_search_type');
    tolx_val = input_value_dict('x_tolerance');
    tolfun_val = input_value_dict('fun_tolerance');
    max_neg_LL_val = input_value_dict('max_neg_LL_val');
    random_effect_names = input_value_dict('random_effect_names');
    mle_parameter_names = input_value_dict('mle_parameter_names');
    gradient_specification = input_value_dict('gradient_specification');

    strain_list = pre_MLE_output_dict('strain_list');
    test_strain_ML_file = pre_MLE_output_dict('test_strain_ML_file');
    strain_LL_table_file = pre_MLE_output_dict('strain_LL_table_file');
    
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify ML test strain parameters for each strain

    %fmincon_opts = optimset('TolX',10^-7,'TolFun',10^-7,'Algorithm','interior-point','MaxIter',5000,'MaxFunEvals',12000,'UseParallel',true,'GradObj','on','GradConstr','on');
    % 2016 and beyond version:
    fmincon_opts = optimoptions('fmincon','TolX', tolx_val,'TolFun', tolfun_val,...
        'Algorithm','interior-point','MaxIter',5000,'MaxFunEvals',12000,...
        'SpecifyObjectiveGradient',gradient_specification,'CheckGradients',false,'Display','off');

    current_iter_parameter_values = ...
        [petite_colony_sigma, nonpetite_colony_sigma, petite_mean, ref_mean];
        % parameter values fixed over this iteration

    % global parameters have already been descaled and de-logspaced, so
        % use these arrays when checking whether to descale them or
        % de-logspace them in later LL functions
    global_scaling_for_strain_mle = ones(size(current_iter_parameter_values));
    global_logspace_for_strain_mle = false(size(current_iter_parameter_values));
    global_mle_parameter_names_for_strain_mle = {'petite_colony_sigma', ...
        'nonpetite_colony_sigma', 'petite_mean', 'ref_mean'};

    % initialize the log likelihood, which is maximized in
        % iterations of this function
    % the log likelihood here is the max LL across all the MA strains
        % plus the LL of the petite data above
        
%    combined_LL = 0 + LL_petite;
            
    strain_param_list = [strcat(strain_list,'_me'), ...
        strcat(strain_list,'_pp')];
            
    [~, mle_parameter_names_strain_only, strain_param_logspace_bool, strain_param_scaling_values, ...
        strain_fixed_values, strain_fixed_indices, lb_values, ub_values, ...
        start_values] = parameter_array_subsetter(strain_param_list, ...
            input_value_dict);

    scaling_values = [global_scaling_for_strain_mle,...
        strain_param_scaling_values];
    logspace_bool = [global_logspace_for_strain_mle,...
        strain_param_logspace_bool];

    % 'global' parameters provided to this function are also fixed
    fixed_parameter_indices = [true(size(current_iter_parameter_values)),...
        strain_fixed_indices];
    fixed_parameter_values = [current_iter_parameter_values,...
        strain_fixed_values];

    mle_parameter_names_strain_mle = ...
        [global_mle_parameter_names_for_strain_mle, mle_parameter_names_strain_only];

    pre_MLE_output_dict_strain = ...
        containers.Map(keys(pre_MLE_output_dict), values(pre_MLE_output_dict));

    input_value_dict_strain = ...
        containers.Map({'LL_calculator', 'gradient_specification', ...
                'mle_parameter_names', 'write_checkpoint'}, ...
            {'LL_calculator_strain_looper_pairwise_multiref', ...
                gradient_specification, mle_parameter_names_strain_mle, ...
                false});

    % Estimate likelihood for all pairs of current strain
    % If there are no parameters to fit, just run the likelihood estimation

    if sum(~fixed_parameter_indices) > 0

        min_problem = createOptimProblem('fmincon','objective',...
            @(strainwise_parameter_vals_partial) ...
            LL_calculator(strainwise_parameter_vals_partial,...
                mle_parameter_names_strain_mle, fixed_parameter_indices, ...
                fixed_parameter_values, logspace_bool, ...
                scaling_values, max_neg_LL_val, ...
                input_value_dict_strain, ...
                pre_MLE_output_dict_strain), ...
            'x0', start_values,...
            'lb', lb_values, 'ub', ub_values, ...
            'options', fmincon_opts);
        
        gs = GlobalSearch;
        gs.TolFun = tolfun_val;
        gs.TolX = tolx_val;
        gs.StartPointsToRun='bounds';
        gs.Display='off';
        gs.NumStageOnePoints = 30;
        gs.NumTrialPoints = 100;

        if strcmp(strainwise_search_type,'global')
            [out_param_vals]=...
                run(gs,min_problem);
        else
            [out_param_vals]=...
                fmincon(min_problem);
        end
    else
        out_param_vals = start_values;
    end

    % record output parameter values before unscaling them
    corrected_param_vals = NaN(size(strain_fixed_indices));
    corrected_param_vals(strain_fixed_indices) = ...
        strain_fixed_values(strain_fixed_indices);
    corrected_param_vals(~strain_fixed_indices) = ...
        out_param_vals;
    
    unscaled_strain_param_vals = ...
        reverse_value_scaler(corrected_param_vals, ...
        strain_param_logspace_bool, ...
        strain_param_scaling_values);

    % calculate global likelihood gradient given current global
        % parameters and current ML estimates for strain-specific
        % parameters
    parameter_values_final = [petite_colony_sigma, ...
        nonpetite_colony_sigma, petite_mean, ref_mean, unscaled_strain_param_vals];

    pre_MLE_output_dict_strain('fitted_parameters') = mle_parameter_names_strain_mle;

    [LL_GR_diff, unscaled_gradient_vector_gr_diff, grad_parameter_names_gr_diff] = ...
        LL_calculator_strain_looper_pairwise_multiref(parameter_values_final, ...
            input_value_dict_strain, pre_MLE_output_dict_strain);
    
    unscaled_gradient_dict_complete = ...
        containers.Map(grad_parameter_names_gr_diff, ...
        unscaled_gradient_vector_gr_diff);

    % d_LL_d_ranef from petite_param_gradient
    ranef_gradient_vals = values(petite_param_gradient_dict, ...
        random_effect_names);
    ranef_gradient_dict = containers.Map(ranef_names_from_param_list, ...
        ranef_gradient_vals);
    % remove petite_colony_sigma from ranef gradient dictionary to
        % prevent overwriting petite_colony_sigma gradient value in
        % unscaled_gradient_dict_complete
    remove(ranef_gradient_dict, 'petite_colony_sigma');
    unscaled_gradient_dict_complete = ...
        [unscaled_gradient_dict_complete; ranef_gradient_dict];

    % d_LL_d_petite_mean includes a component from petite_param_gradient
    unscaled_gradient_dict_complete('petite_mean') = ...
        unscaled_gradient_dict_complete('petite_mean') + ...
        petite_param_gradient_dict('petite');

    % d_LL_d_petite_colony_sigma includes a component from petite_param_gradient
    unscaled_gradient_dict_complete('petite_colony_sigma') = ...
        unscaled_gradient_dict_complete('petite_colony_sigma') + ...
        petite_param_gradient_dict('petite_colony');
    
    unscaled_gradient_vector = cell2mat(values(unscaled_gradient_dict_complete));
    grad_parameter_names = keys(unscaled_gradient_dict_complete);

    combined_LL = LL_GR_diff + LL_petite;
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
        export_table_data = num2cell(unscaled_strain_param_vals);
        export_table = ...
            table(export_table_data{:},'VariableNames', mle_parameter_names_strain_only);

        Current_Best_LL = combined_LL;
        Best_LL_Runtime = runtime;
        strain_LL_table = table(Current_Best_LL, Best_LL_Runtime);

        writetable(export_table,test_strain_ML_file);
        writetable(strain_LL_table,strain_LL_table_file);

    end

end
    
