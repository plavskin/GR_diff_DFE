function [combined_LL, unscaled_gradient_vector, grad_parameter_names] = ...
    LL_calculator_strain_looper_pairwise_lambdafit(param_vals,...
    input_value_dict, pre_MLE_output_dict)
    
    % EP 17-11-07

    % Takes in parameters that apply across all test strains
    % Loops through test strains to estimate the optimal strain parameters
        % given the current iteration of general parameter values
    % For each test strain, calculates the log likelihood of observing a list
        % of differences between test and reference strain GRs within
        % pairs of colonies
    % test strain GR = (ref strain GR)*exp(test strain mut_effect)
    % test strain mut_effect drawn from poisson-distributed number of
        % reflected gamma distributions
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
    L = input_value_dict('L');

    strain_list = pre_MLE_output_dict('strain_list');
    test_strain_list_by_pair = pre_MLE_output_dict('test_strain_list_by_pair');
    GR_diff_list = pre_MLE_output_dict('GR_diff_list');
    test_strain_ML_file = pre_MLE_output_dict('test_strain_ML_file');
    strain_LL_table_file = pre_MLE_output_dict('strain_LL_table_file');
    N = pre_MLE_output_dict('N');
    me_pdf_frequency_vals = pre_MLE_output_dict('me_pdf_frequency_vals');
    me_pdf_xvals = pre_MLE_output_dict('me_pdf_xvals');
    F_kernel = pre_MLE_output_dict('F_kernel');
    test_strain_occurances = pre_MLE_output_dict('test_strain_occurances');
    
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

    % mutational effect distribution parameters
    mu_SNM = parameter_dict('mu_SNM');
    shape_SNM = parameter_dict('shape_SNM');
    prop_pos_SNM = parameter_dict('prop_pos_SNM');
    lambda_SNM = parameter_dict('lambda_SNM');
    
    test_strain_number = length(strain_list);
        % # of strains besides reference strain whose likelihood being fitted

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
    fitted_parameters_me_idx = contains(mle_parameter_names, '_SNM');
    fitted_parameters_me = ...
        strrep(mle_parameter_names(fitted_parameters_me_idx), '_SNM', '');

    % Calculate likelihood of observing mutation effects given current
        % global parameters
    % Calculate Fz, the fourier transform of the pdf of observing a
        % single mutational effect from a digamma distribution
    % Then calculate Fa, the fourier transform of the pdf of observing
        % a poisson random number of mutational effects from the
        % digamma distribution in Fz
    if gradient_specification
        [Fz, Fz_gradient_dict] = fourier_domain_digamma(...
            me_pdf_frequency_vals, mu_SNM, shape_SNM, prop_pos_SNM, ...
            fitted_parameters_me, gradient_specification);
        [Fa, Fa_gradient_dict] = fourier_domain_poisson_dist_num(...
            Fz, me_pdf_frequency_vals, lambda_SNM, Fz_gradient_dict, ...
            fitted_parameters_me, gradient_specification);
    else
        Fz = fourier_domain_digamma(...
            me_pdf_frequency_vals, mu_SNM, shape_SNM, prop_pos_SNM, ...
            fitted_parameters_me, gradient_specification);
        Fz_gradient_dict = NaN;
        Fa = fourier_domain_poisson_dist_num(...
            Fz, me_pdf_frequency_vals, lambda_SNM, Fz_gradient_dict, ...
            fitted_parameters_me, gradient_specification);
    end

    % Calculate Fs, which is Fa smoothed with F_kernel
    Fs = Fa .* F_kernel;
    % Calculate discrete fourier transform of Fs to get a smoothed
        % version of the pdf of mutational effects per strain
    me_pdf_smooth = fourier_inverter(Fs, N, L);
    % Need to remove stray imaginary and negative values from
        % me_pdf_smooth
    real_me_pdf_smooth = real(me_pdf_smooth);
    
%   figure; plot(me_pdf_xvals, fourier_inverter(Fa, N, L)); hold on; plot(me_pdf_xvals, me_pdf_smooth, '-r'); hold off;

    if gradient_specification
        % If calculating gradients of LL, need to divide gradient relative
            % to parameter by me_pdf_smooth, so remove any values where
            % me_pdf_smooth is 0
        non_zero_positions = (me_pdf_smooth ~= 0);
        need_to_remove_zero_vals = ...
            length(me_pdf_smooth) > sum(non_zero_positions);
        if need_to_remove_zero_vals
            me_pdf_smooth = me_pdf_smooth(non_zero_positions);
            Fs = Fs(non_zero_positions);
            me_pdf_frequency_vals = me_pdf_frequency_vals(non_zero_positions);
            me_pdf_xvals = me_pdf_xvals(non_zero_positions);
        end
        % Calculate gradient of log likelihood of mutational effect
            % distribution relative to mutational effect
        d_Fs_d_x = - 1i * me_pdf_frequency_vals .* Fs;
        d_me_pdf_smooth_d_me = fourier_inverter(d_Fs_d_x, N, L);
        d_me_LL_smooth_d_me = d_me_pdf_smooth_d_me ./ me_pdf_smooth;
        % Need to remove stray imaginary values from d_me_LL_smooth_d_me
        real_d_me_LL_smooth_d_me = real(d_me_LL_smooth_d_me);
        % Create a dictionary of gradients of LL relative to each me
            % distribution parameter
%        me_dist_param_grad_vector_dict = containers.Map(keys(Fa_gradient_dict), ...
%            repmat({NaN},size(keys(Fa_gradient_dict))));
        me_dist_param_grad_vector_dict = containers.Map();
        for current_fitted_param_idx = 1:length(fitted_parameters_me)
            current_fitted_param = ...
                fitted_parameters_me{current_fitted_param_idx};
            d_Fa_d_current_fitted_param = ...
                Fa_gradient_dict(current_fitted_param);
            % Adjust gradients to account for smoothing
            d_Fs_d_current_fitted_param = ...
                d_Fa_d_current_fitted_param .* F_kernel;
            if need_to_remove_zero_vals
                % remove values where me_pdf_smooth is 0
                d_Fs_d_current_fitted_param = ...
                    d_Fs_d_current_fitted_param(non_zero_positions);
            end
            d_me_pdf_smooth_d_current_fitted_param = ...
                fourier_inverter(d_Fs_d_current_fitted_param, N, L);
            %d_LL_d_current_fitted_param = ...
            %    real(d_me_pdf_smooth_d_current_fitted_param ./ ...
            %        me_pdf_smooth);
            d_LL_d_current_fitted_param = ...
                real(d_me_pdf_smooth_d_current_fitted_param) ./ ...
                    real_me_pdf_smooth;
%           figure; plot(me_pdf_xvals, d_LL_d_current_fitted_param, movmean(me_pdf_xvals,2), movmean(d_LL_d_current_fitted_param,2)); title(current_fitted_param);
           me_dist_param_grad_vector_dict(current_fitted_param) = ...
                d_LL_d_current_fitted_param;
        end
        corrected_me_dist_param_grad_vector_dict_keys = ...
            strcat(keys(me_dist_param_grad_vector_dict), '_SNM');
        me_dist_param_grad_vector_dict = ...
            containers.Map(corrected_me_dist_param_grad_vector_dict_keys, ...
                values(me_dist_param_grad_vector_dict));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop through test strains and identify ML test strain parameters for each
    % The sum of LL across all strains is the combined log likelihood for the
        % data, given the current parameter values
    % Because test strain-specific parameters are independent of each other in
        % terms of their effect on likelihood, we can calculate likelihood for 
        % one test strain at a time (given reference strain and general
        % parameters), greatly reducing the dimensionality of the test strain
        % parameter search

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
    global_mle_parameter_names_for_strain_mle = {'petite_colony_sigma', ...
        'nonpetite_colony_sigma', 'petite_mean', 'ref_mean', 'ref_petite_prop'};

    % initialize the log likelihood, which is maximized in
        % iterations of this function
    % the log likelihood here is the max LL across all the MA strains
        % plus the LL of the petite data above
        
%    combined_LL = 0 + LL_petite;

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
            if sum(~fixed_parameter_indices) > 0
                current_mle_parameter_names = ...
                    [global_mle_parameter_names_for_strain_mle, ...
                        strain_param_current_list];

                pre_MLE_output_dict_current_strain = ...
                    containers.Map({'strain_list', 'test_strain_list_by_pair', ...
                            'GR_diff_list', 'me_pdf_xvals', 'me_pdf'}, ...
                        {strain_current_list, test_strain_current_list, ...
                            GR_diff_current_list, me_pdf_xvals, real_me_pdf_smooth});
                if gradient_specification
                    pre_MLE_output_dict_current_strain('d_me_LL_d_me') = ...
                        real_d_me_LL_smooth_d_me;
                end

                input_value_dict_curent_strain = ...
                    containers.Map({'LL_calculator', 'gradient_specification', ...
                            'mle_parameter_names'}, ...
                        {'LL_calculator_strains_pairwise_from_me_dist', ...
                            gradient_specification, current_mle_parameter_names});

                current_min_problem = createOptimProblem('fmincon','objective',...
                    @(strainwise_parameter_vals_partial) ...
                    LL_calculator(strainwise_parameter_vals_partial,...
                        current_mle_parameter_names, fixed_parameter_indices, ...
                        fixed_parameter_values, current_logspace_bool, ...
                        current_scaling_values, max_neg_LL_val, ...
                        input_value_dict_curent_strain, ...
                        pre_MLE_output_dict_current_strain), ...
                    'x0', current_start_values,...
                    'lb', current_lb_values, 'ub', current_ub_values, ...
                    'options', fmincon_opts);

                if strcmp(strainwise_search_type,'global')
                    [current_out_param_vals]=...
                        run(gs,current_min_problem);
                else
                    [current_out_param_vals]=...
                        fmincon(current_min_problem);
                end
            else
                current_out_param_vals = current_start_values;
            end

%            combined_LL = combined_LL - neg_LL_current;

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
    parameter_values_grad_calc = [petite_colony_sigma, ...
        nonpetite_colony_sigma, petite_mean, ref_mean, ref_petite_prop, ...
        test_strain_ML_params(:,1)', test_strain_ML_params(:,2)'];

    strain_param_names = [strcat(strain_list,'_me'), ...
        strcat(strain_list,'_pp')];

    mle_parameter_list_grad_calc = ...
        [global_mle_parameter_names_for_strain_mle, ...
            strain_param_names];

    pre_MLE_output_dict_gr_diff_grad = ...
        containers.Map({'strain_list', 'test_strain_list_by_pair', ...
                'GR_diff_list', 'me_pdf_xvals', 'me_pdf','fitted_parameters'}, ...
            {strain_list, test_strain_list_by_pair, GR_diff_list, ...
                me_pdf_xvals, real_me_pdf_smooth, mle_parameter_list_grad_calc});
    if gradient_specification
        pre_MLE_output_dict_gr_diff_grad('d_me_LL_d_me') = ...
            real_d_me_LL_smooth_d_me;
    end

    input_value_dict_gr_diff_grad = ...
        containers.Map({'gradient_specification', 'mle_parameter_names'}, ...
            {gradient_specification, mle_parameter_list_grad_calc});

    [LL_GR_diff, unscaled_gradient_vector_gr_diff, grad_parameter_names_gr_diff] = ...
        LL_calculator_strains_pairwise_from_me_dist(parameter_values_grad_calc, ...
            input_value_dict_gr_diff_grad, pre_MLE_output_dict_gr_diff_grad);
    
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

    % calculate me_dist_param_grad_dict by interpolating strain me vals
        % in the gradients relative to each me distribution parameter
    me_dist_param_grad_dict = containers.Map();
    if gradient_specification
        for current_fitted_param_idx = ...
                1:length(corrected_me_dist_param_grad_vector_dict_keys)
            current_fitted_param = ...
                corrected_me_dist_param_grad_vector_dict_keys{...
                    current_fitted_param_idx};
            d_LL_d_current_fitted_param = ...
                me_dist_param_grad_vector_dict(current_fitted_param);
            d_LL_d_current_fitted_param_current_data = ...
                interp1(me_pdf_xvals, d_LL_d_current_fitted_param, ...
                    test_strain_ML_params(:,1)');
            me_dist_param_grad_dict(current_fitted_param) = ...
                nansum(d_LL_d_current_fitted_param_current_data .* ...
                    test_strain_occurances);
%            me_dist_param_grad_dict(current_fitted_param) = ...
%                nansum(d_LL_d_current_fitted_param_current_data);
        end            
    else
        me_dist_param_grad_dict = {};
    end
    % Combined unscaled_gradient_dict_complete with me_dist_param_grad_dict
    unscaled_gradient_dict_complete = ...
        [unscaled_gradient_dict_complete; me_dist_param_grad_dict];
    
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
    
