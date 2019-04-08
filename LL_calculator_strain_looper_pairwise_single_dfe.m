function [combined_LL, unscaled_gradient_vector, grad_parameter_names] = ...
    LL_calculator_strain_looper_pairwise_single_dfe(param_vals,...
    input_value_dict, pre_MLE_output_dict)
    
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
    DFE_parameters = input_value_dict('DFE_parameters');
    L = input_value_dict('L');
    current_model = input_value_dict('model');

    strain_list = pre_MLE_output_dict('strain_list');
    test_strain_list_by_pair = pre_MLE_output_dict('test_strain_list_by_pair');
    GR_diff_list = pre_MLE_output_dict('GR_diff_list');
    test_strain_ML_file = pre_MLE_output_dict('test_strain_ML_file');
    strain_LL_table_file = pre_MLE_output_dict('strain_LL_table_file');
    N = pre_MLE_output_dict('N');
    me_pdf_frequency_vals = pre_MLE_output_dict('me_pdf_frequency_vals');
    me_pdf_xvals = pre_MLE_output_dict('me_pdf_xvals');
    F_kernel = pre_MLE_output_dict('F_kernel');
    ranef_names_from_param_list = ...
        pre_MLE_output_dict('ranef_names_from_param_list');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify general parameter values
    
    petite_colony_sigma = pre_MLE_output_dict('petite_colony_sigma');
        % s.d. of colony GRs of petite distribution
    nonpetite_colony_sigma = pre_MLE_output_dict('nonpetite_colony_sigma');
        % s.d. of colony GRs of non-petite reference and test strain distribution
    petite_mean = pre_MLE_output_dict('petite_mean');
        % mean growth rate of petite colonies, regardless of genotype
    ref_mean = pre_MLE_output_dict('ref_mean');
        % mean growth rate of non-petite ref colonies
    ref_petite_prop = pre_MLE_output_dict('ref_petite_prop');
        % proportion of petites in ref strain
    
    test_strain_number = length(strain_list);
        % # of strains besides reference strain whose likelihood being fitted

    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get previously calculated likelihood of observing petite data given
        % current global parameters, as well as its gradient

    LL_petite = pre_MLE_output_dict('LL_petite');
    petite_param_gradient_dict = ...
        pre_MLE_output_dict('petite_param_gradient_dict');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate likelihood of observing mutation effects given current
        % global parameters (distribution of fitness effects, DFE)
    % First set up function that will calculate fourier transform of DFE
    fourier_dfe_fun_name = strcat('fourier_dfe_', current_model);
    fourier_dfe_function = str2func(fourier_dfe_fun_name);
    % Calculate fourier transform of DFE
    [Fa, Fa_gradient_dict, ~, ...
        fitted_parameters_me_fullname] = fourier_dfe_function(...
            param_vals, input_value_dict, pre_MLE_output_dict);
    % Calculate Fs, which is Fa smoothed with F_kernel
    [Fs, Fs_gradient_dict] = fourier_space_smoother(Fa, F_kernel, ...
        me_pdf_frequency_vals, Fa_gradient_dict, ...
        gradient_specification);
    % Calculate discrete fourier transform of Fs to get a smoothed
        % version of the pdf of mutational effects per strain (and, if
        % applicable, the gradients of the log of the pdf with respect
        % to each parameter in gradient_param_list)
    [real_me_pdf_smooth, me_pdf_xvals, me_dist_param_grad_vector_dict] = ...
        likelihood_from_fourier(Fs, me_pdf_xvals, N, L, ...
            Fs_gradient_dict, gradient_specification);

%    if gradient_specification
%        % Pull out gradient of log likelihood of mutational effect
%            % distribution with respect to to mutational effect
%        real_d_me_LL_smooth_d_me = ...
%            me_dist_param_grad_vector_dict('random_variable');
%    end

%   figure; plot(me_pdf_xvals, fourier_inverter(Fa, N, L)); hold on; plot(me_pdf_xvals, real_me_pdf_smooth, '-r'); hold off;

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
    test_strain_ML_params = NaN(test_strain_number, 1);

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
            
            strain_param_current_list = strcat(strain_current_list,'_pp');
            
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
                            'GR_diff_list', 'me_pdf_xvals', 'me_pdf', ...
                            'me_dist_param_grad_vector_dict'}, ...
                        {strain_current_list, test_strain_current_list, ...
                            GR_diff_current_list, me_pdf_xvals, real_me_pdf_smooth, ...
                            me_dist_param_grad_vector_dict});
%                if gradient_specification
%                    pre_MLE_output_dict_current_strain('d_me_LL_d_me') = ...
%                        real_d_me_LL_smooth_d_me;
%                end

                input_value_dict_curent_strain = ...
                    containers.Map({'LL_calculator', ...
                        'gradient_specification', ...
                            'mle_parameter_names', 'DFE_parameters', ...
                            'max_neg_LL_val', 'write_checkpoint'}, ...
                        {'LL_calculator_strains_pairwise_from_me_dist', ...
                            gradient_specification, ...
                            current_mle_parameter_names, ...
                            DFE_parameters, max_neg_LL_val, false});

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

            corrected_current_param_val = NaN;
            corrected_current_param_val(current_strain_fixed_indices) = ...
                current_strain_fixed_values(current_strain_fixed_indices);
            corrected_current_param_val(~current_strain_fixed_indices) = ...
                current_out_param_vals;
            % record output parameter values before unscaling them
            unscaled_current_param_val = ...
                reverse_value_scaler(corrected_current_param_val, ...
                current_strain_param_logspace_bool, ...
                current_strain_param_scaling_values);
            
            test_strain_ML_params(strain_idx) = unscaled_current_param_val;
        end

    end

    % calculate global likelihood gradient given current global
        % parameters and current ML estimates for strain-specific
        % parameters
    parameter_values_grad_calc = [petite_colony_sigma, ...
        nonpetite_colony_sigma, petite_mean, ref_mean, ref_petite_prop, ...
        test_strain_ML_params'];

    strain_param_names = strcat(strain_list,'_pp');

    mle_parameter_list_grad_calc = ...
        [global_mle_parameter_names_for_strain_mle, ...
            fitted_parameters_me_fullname, strain_param_names];
    mle_parameter_names_for_grad_calc = ...
        [global_mle_parameter_names_for_strain_mle, strain_param_names];

    pre_MLE_output_dict_gr_diff_grad = ...
        containers.Map({'strain_list', 'test_strain_list_by_pair', ...
                'GR_diff_list', 'me_pdf_xvals', 'me_pdf', ...
                'me_dist_param_grad_vector_dict', 'fitted_parameters'}, ...
            {strain_list, test_strain_list_by_pair, GR_diff_list, ...
                me_pdf_xvals, real_me_pdf_smooth, ...
                me_dist_param_grad_vector_dict, mle_parameter_list_grad_calc});
%    if gradient_specification
%        pre_MLE_output_dict_gr_diff_grad('d_me_LL_d_me') = ...
%            real_d_me_LL_smooth_d_me;
%    end

    input_value_dict_gr_diff_grad = ...
        containers.Map({'gradient_specification', 'mle_parameter_names', ...
            'DFE_parameters', 'max_neg_LL_val'}, ...
            {gradient_specification, mle_parameter_names_for_grad_calc, ...
                DFE_parameters, max_neg_LL_val});

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

%    % calculate me_dist_param_grad_dict by interpolating strain me vals
%        % in the gradients relative to each me distribution parameter
%    me_dist_param_grad_dict = containers.Map('KeyType', 'char', ...
%            'ValueType', 'any');
%    if gradient_specification
%        for current_fitted_param_idx = ...
%                1:length(complete_me_parameter_list)
%            current_fitted_param = ...
%                complete_me_parameter_list{current_fitted_param_idx};
%            if any(strcmp(fitted_parameters_me_fullname, ...
%                    current_fitted_param))
%                d_LL_d_current_fitted_param = ...
%                    me_dist_param_grad_vector_dict(current_fitted_param);
%                d_LL_d_current_fitted_param_current_data = ...
%                    interp1(me_pdf_xvals, d_LL_d_current_fitted_param, ...
%                        test_strain_ML_params(:,1)');
%                me_dist_param_grad_dict(current_fitted_param) = ...
%                    nansum(d_LL_d_current_fitted_param_current_data);
%            else
%                me_dist_param_grad_dict(current_fitted_param) = NaN;
%            end
%        end            
%    end
%    % Combined unscaled_gradient_dict_complete with me_dist_param_grad_dict
%    unscaled_gradient_dict_complete = ...
%        [unscaled_gradient_dict_complete; me_dist_param_grad_dict];
    
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
        strain_names = {strain_list{:}}';
        pp_mle_vals = test_strain_ML_params(:);

        strain_param_name_list = strcat(strain_names,'_pp')';
        strain_table_data = num2cell(pp_mle_vals);
        strain_table = table(strain_table_data{:},'VariableNames', strain_param_name_list);

        % retrieve DFE parameter info
        [~, DFE_param_positions] = ...
            ismember(DFE_parameters, mle_parameter_names);
        DFE_table_data = num2cell(param_vals(DFE_param_positions)');
        DFE_table = table(DFE_table_data{:}, 'VariableNames', DFE_parameters);

        export_table = [DFE_table strain_table];

        Current_Best_LL = combined_LL;
        Best_LL_Runtime = runtime;
        strain_LL_table = table(Current_Best_LL, Best_LL_Runtime);

        writetable(export_table,test_strain_ML_file);
        writetable(strain_LL_table,strain_LL_table_file);

    end

end