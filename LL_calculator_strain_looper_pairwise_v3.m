function [neg_combined_LL,global_gradient_vector_partial] = LL_calculator_strain_looper_pairwise_v3(param_vals_partial,...
    global_fixed_parameter_indices,global_fixed_parameter_values,...
    mut_effect_fixed_parameter_values,petite_prop_fixed_parameter_values,...
    test_strain_ML_file,strain_LL_table_file,...
    strain_mut_effect_lower_bounds,strain_mut_effect_upper_bounds,...
    strain_petite_proportion_lower_bounds,strain_petite_proportion_upper_bounds,...
    strain_mut_effect_start_vals,strain_petite_proportion_start_vals,strain_list,...
    test_strain_list_by_pair,GR_diff_list,strainwise_search_type,max_neg_LL_val,...
    response_vector,fixef_ID_mat,ranef_corr_struct,unique_fixefs,...
    unique_ranef_categories,order_vector,block_start_positions,tolx_val, tolfun_val)
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
    % Set up csv filename where results of test strain parameters will be saved
%    csv_folder = strcat(output_folder,'/csv_output');
%    csv_file_prename = strcat(csv_folder,'/test_strain_params-',...
%        csvname_string_addition);
%    csv_filename = strcat(csv_folder,'/test_strain_params-',...
%        csvname_string_addition,'.csv');
%    strain_LL_table_file = strcat(csv_folder,'/current_best_strainlooper_LL-',...
%        csvname_string_addition,'.csv');

    % create csv_folder in output_directory if it doesn't currently exist
%    if exist(csv_folder,'dir') ~= 7
%        mkdir(csv_folder);
%    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify parameter values

    % some parameters may be 'fixed' (i.e. not fitted by current iteration of
        % MLE); these are provided to the function in a separate list
    % compile a list of all parameters, fixed or fitted, and identify values
        % belonging to each individual parameter using that list
    param_vals = NaN(size(global_fixed_parameter_indices));
    param_vals(global_fixed_parameter_indices) = global_fixed_parameter_values(~isnan(global_fixed_parameter_values));
    param_vals(~global_fixed_parameter_indices) = param_vals_partial;
    
    petite_sigma = param_vals(1);
        % s.d. of colony GRs of petite distribution
    nonpetite_sigma = param_vals(2);
        % s.d. of colony GRs of non-petite reference and test strain distribution
    petite_GR = param_vals(3);
        % mean growth rate of petite colonies, regardless of genotype
    ref_GR = param_vals(4);
        % mean growth rate of non-petite ref colonies
    ref_petite_proportion = param_vals(5);
        % proportion of petites in ref strain
    plate_sigma = param_vals(6);
        % s.d. of plate effect on petite GR
    well_sigma = param_vals(7);
        % s.d. of well effect on petite GR
    field_sigma = param_vals(8);
        % s.d. of field effect on petite GR

    test_strain_number = length(strain_list);
        % # of strains besides reference strain whose likelihood being fitted

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate likelihood of observing petite data given current global parameters
    tic;
    calculate_gradient = true;
    nonref_indices = [false,true,true,true,true,true];
    petite_parameters_total = [plate_sigma,well_sigma,field_sigma,...
        petite_GR*(~nonref_indices),...
        petite_sigma*(~nonref_indices)];
    petite_fixed_parameter_indices = [false,false,false,...
        nonref_indices,nonref_indices];
    petite_parameters = petite_parameters_total(~petite_fixed_parameter_indices);
    petite_fixed_parameter_values = ...
        petite_parameters_total(petite_fixed_parameter_indices);
        % all parameters that we want to get back gradient for must be
            % labeled as unfixed
    [neg_LL_petite,gradient_vector_partial_petite] = ...
        LL_mixef_calc_v5b(petite_parameters,petite_fixed_parameter_indices,...
            petite_fixed_parameter_values,response_vector,ranef_corr_struct,...
            fixef_ID_mat,nonref_indices,unique_ranef_categories,...
            calculate_gradient,order_vector,block_start_positions);

%    toc;
%    disp('done with petite LL calc')

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

    if strcmp(strainwise_search_type,'global')
    % set up structure for global solution search

        gs = GlobalSearch;
        gs.TolFun = tolfun_val;
        gs.TolX = tolx_val;
        gs.NumStageOnePoints = 2^2; % 9
        gs.NumTrialPoints = 3^2; % 90
        gs.StartPointsToRun='bounds';
        gs.Display = 'off';
    end
    
    %fmincon_opts = optimset('TolX',10^-7,'TolFun',10^-7,'Algorithm','interior-point','MaxIter',5000,'MaxFunEvals',12000,'UseParallel',true,'GradObj','on','GradConstr','on');
    % 2016 and beyond version:
    fmincon_opts = optimoptions('fmincon','TolX',10^-5,'TolFun',10^-4,...
        'Algorithm','interior-point','MaxIter',5000,'MaxFunEvals',12000,...
        'SpecifyObjectiveGradient',true,'CheckGradients',false,'Display','off');

    current_iter_parameter_values = [petite_sigma,nonpetite_sigma,petite_GR,...
        ref_GR,ref_petite_proportion];
        % parameter values fixed over this iteration

    % Initialize a matrix for storing test strain mut_effects and petite proportions
    test_strain_ML_params = NaN(test_strain_number,2);

    tic;

    % initialize the negative log likelihood, which is minimized in
        % iterations of this function
    % the log likelihood here is the max LL across all the MA strains
        % plus the LL of the petite data above
        
    neg_combined_LL = 0+neg_LL_petite;

    for strain_idx = 1:test_strain_number

        current_strain = strain_list{strain_idx};
        
        % set up starting vals and lower and upper bounds for current strain
        current_strain_mut_effect_start_val = strain_mut_effect_start_vals(strain_idx);
        current_strain_petite_proportion_start_val = ...
            strain_petite_proportion_start_vals(strain_idx);
        current_strain_mut_effect_lower_bound = strain_mut_effect_lower_bounds(strain_idx);
        current_strain_petite_proportion_lower_bound = ...
            strain_petite_proportion_lower_bounds(strain_idx);
        current_strain_mut_effect_upper_bound = strain_mut_effect_upper_bounds(strain_idx);
        current_strain_petite_proportion_upper_bound = ...
            strain_petite_proportion_upper_bounds(strain_idx);
        
        % set up list of data to pass to likelihood estimation function
        % find which data in list corresponds to current_strain
        current_indices = find(strcmp(test_strain_list_by_pair,current_strain));
        % only run mle/ll calc if there are samples from current_strain
            % (this is an issue when working with a random subset of data)
        if size(current_indices, 1) > 0
            strain_current_list = {current_strain};
            test_strain_current_list = test_strain_list_by_pair(current_indices);
            GR_diff_current_list = GR_diff_list(current_indices);

            % set up starting parameters
            current_start_values = [current_strain_mut_effect_start_val,...
                current_strain_petite_proportion_start_val];
            current_lb_values = [current_strain_mut_effect_lower_bound,...
                current_strain_petite_proportion_lower_bound];
            current_ub_values = [current_strain_mut_effect_upper_bound,...
                current_strain_petite_proportion_upper_bound];

            % figure out whether the current strain has a fixed petite 
                % proportion or mutation effect
            % mut_effect_fixed_parameter_values and 
                % petite_prop_fixed_parameter_values contain NaN at
                % non-fixed positions
            strain_param_values = [mut_effect_fixed_parameter_values(strain_idx),...
                petite_prop_fixed_parameter_values(strain_idx)];
            fixed_test_strain_params = ~isnan(strain_param_values);
            current_start_values = current_start_values(~fixed_test_strain_params);
            current_lb_values = current_lb_values(~fixed_test_strain_params);
            current_ub_values = current_ub_values(~fixed_test_strain_params);

            % 'global' parameters provided to this function are also fixed
            fixed_parameter_indices = [true(size(current_iter_parameter_values)),...
                fixed_test_strain_params];
            fixed_parameter_values = [current_iter_parameter_values,...
                strain_param_values(fixed_test_strain_params)];

            % Estimate likelihood for all pairs of current strain
            % If there are no parameters to fit, just run the likelihood estimation

            if sum(~fixed_parameter_indices) > 0
                current_min_problem = createOptimProblem('fmincon','objective',...
                    @(strainwise_parameter_vals_partial) LL_calculator_strains_pairwise_2_sigmas(strainwise_parameter_vals_partial,...
                        fixed_parameter_indices,fixed_parameter_values,...
                        strain_current_list,test_strain_current_list,GR_diff_current_list,return_all_grads,max_neg_LL_val),'x0',current_start_values,...
                    'lb',current_lb_values,'ub',current_ub_values,'options',fmincon_opts);

                if strcmp(strainwise_search_type,'global')
                    [current_out_param_vals,neg_LL_current,exitflag_current,...
                        output_current,solutions_current]=...
                        run(gs,current_min_problem);
                else
                    [current_out_param_vals,neg_LL_current,exitflag_current]=...
                        fmincon(current_min_problem);
                end
            else
                neg_LL_current = LL_calculator_strains_pairwise_2_sigmas(NaN,...
                    fixed_parameter_indices,fixed_parameter_values,...
                    strain_current_list,test_strain_current_list,GR_diff_current_list,...
                    return_all_grads,max_neg_LL_val);
                current_out_param_vals = NaN;
            end

            neg_combined_LL = neg_combined_LL+neg_LL_current;

            corrected_current_param_vals = NaN(1,2);
            corrected_current_param_vals(fixed_test_strain_params) = ...
                strain_param_values(fixed_test_strain_params);
            corrected_current_param_vals(~fixed_test_strain_params) = ...
                current_out_param_vals;

            test_strain_ML_params(strain_idx,:) = corrected_current_param_vals;
        end
%    toc
%    disp(current_strain)

    end

    % calculate global likelihood gradient given current global
        % parameters and current ML estimates for strain-specific
        % parameters
    return_all_grads = true;
    fixed_parameter_values_grad_calc = [petite_sigma,nonpetite_sigma,petite_GR,...
        ref_GR,ref_petite_proportion,test_strain_ML_params(:,1)',...
        test_strain_ML_params(:,2)'];
    fixed_parameter_indices_grad_calc = true(size(fixed_parameter_values_grad_calc));
    [neg_LL_current,strain_gradient_vector,global_gradient_vector_no_petites] = LL_calculator_strains_pairwise_2_sigmas(NaN,...
                fixed_parameter_indices_grad_calc,fixed_parameter_values_grad_calc,...
                strain_list,test_strain_list_by_pair,GR_diff_list,return_all_grads,max_neg_LL_val);
    % In this case, the strain_gradient_vector output is blank,
            % which happens when all strain parameters are fixed in
            % LL_calculator_strains_pairwise
    global_gradient_vector = zeros(size(param_vals));
    % d_LL_d_petite_sigma
    global_gradient_vector(1) = gradient_vector_partial_petite(5) + ...
        global_gradient_vector_no_petites(1);
    % d_LL_d_nonpetite_sigma
    global_gradient_vector(2) = global_gradient_vector_no_petites(2);
    % d_LL_d_petite_GR
    global_gradient_vector(3) = gradient_vector_partial_petite(4) + ...
        global_gradient_vector_no_petites(3);
    % d_LL_d_ref_GR
    global_gradient_vector(4) = global_gradient_vector_no_petites(4);
    % d_LL_d_ref_petite_proportion
    global_gradient_vector(5) = global_gradient_vector_no_petites(5);
    % d_LL_d_ranef
    global_gradient_vector(6:8) = gradient_vector_partial_petite(1:3);

    global_gradient_vector_partial = global_gradient_vector(~global_fixed_parameter_indices);

    global_gradient_vector_partial(global_gradient_vector_partial>max_neg_LL_val) = max_neg_LL_val;
    global_gradient_vector_partial(global_gradient_vector_partial<-max_neg_LL_val) = -max_neg_LL_val;

    if neg_combined_LL > max_neg_LL_val
        neg_combined_LL = max_neg_LL_val;
    end
    
    runtime = toc;

    % need to first read previous file and see whether current combined_LL is higher than one in prev file

    % default is to write output to a file, but an suppress this if previous
        % run of MLE already produced a better log likelihood

    write_output = true;

    if exist(strain_LL_table_file,'file') == 2
        previous_LL_table = readtable(strain_LL_table_file);

        if -neg_combined_LL < previous_LL_table.Current_Best_LL(1)
            write_output = false;
        end
    end

    if write_output
        Strain_Names = {strain_list{:}}';
        Mutation_Effects = test_strain_ML_params(:,1);
        Petite_Proportions = test_strain_ML_params(:,2);
        export_table = table(Strain_Names,Mutation_Effects,Petite_Proportions);

        Current_Best_LL = -neg_combined_LL;
        Best_LL_Runtime = runtime;
        strain_LL_table = table(Current_Best_LL,Best_LL_Runtime);

        writetable(export_table,test_strain_ML_file);
        writetable(strain_LL_table,strain_LL_table_file);

    end

 %   if runtime < 120
 %       pausetime=120-runtime;
 %       pause(pausetime)
 %   end
end
    
