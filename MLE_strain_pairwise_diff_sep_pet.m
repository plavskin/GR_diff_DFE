function MLE_strain_pairwise_diff_sep_pet(key_list, value_list)

    % EP 18-07-27

    % This function writes an output csv file which contains the MLE values for
        % each global parameter
    % In addition, an internal function writes a csv file which contains the MLE
        % values of each strain's parameters, given the global parameters that
        % result in the highest LL

    % fixed_global_parameter_array contains NaN in the position of
        % experiment-wide parameters that are not fixed, and 

    %%%%%%%%%%%%%% ??? %%%%%%%%%%%%%%
    % FOR NOW, ASSUMES A SINGLE REFERENCE STRAIN
    % SIGMA VALUE MUST BE ABOVE 0!
    %%%%%%%%%%%%%% ??? %%%%%%%%%%%%%%

    % sep_pet version: optimizes different gr sigma values for petite
        % vs non-petite colonies
    % uses linear mixed effect model to separately calculate likelihood
        % of observing petite distribution given their mean and sigma
        % growth rate parameters; thus includes parameters for field,
        % well, and plate effects

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % To avoid overflow issues, set a max value that log likelihood can
        % be equal to throughout the MLE
    %max_neg_LL_val = realmax/(strain_number*2);
    max_neg_LL_val = 10^50;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    parameter_dict = containers.Map(key_list,value_list);

    external_counter = str2double(parameter_dict('external_counter'));
    combined_fixed_parameter_array = parameter_dict('combined_fixed_parameter_array');
    combined_min_array_unscaled = parameter_dict('combined_min_array');
    combined_max_array_unscaled = parameter_dict('combined_max_array');
    combined_length_array = parameter_dict('combined_length_array');
    combined_position_array = cellfun(@str2num,parameter_dict('combined_position_array'));
    combined_start_values_array_unscaled = parameter_dict('combined_start_values_array');
    combined_scaling_array = parameter_dict('combined_scaling_array');
    parameter_list = parameter_dict('parameter_list');
    output_file = parameter_dict('output_file');
    parallel_processors = parameter_dict('parallel_processors');
    ms_positions = parameter_dict('ms_positions');
    combined_profile_ub_array_unscaled = parameter_dict('combined_profile_ub_array');
    combined_profile_lb_array_unscaled = parameter_dict('combined_profile_lb_array');
    ms_grid_parameter_array = parameter_dict('ms_grid_parameter_array');
    combined_logspace_parameters = parameter_dict('combined_logspace_parameters');
    datafile_path = parameter_dict('datafile_path');
    output_id_parameter = parameter_dict('output_id_parameter');
    % optional parameters
    phenotype_file = parameter_dict('phenotype_file');
    petite_file = parameter_dict('petite_file');
    % process parameter name arrays into bool arrays
    combined_logspace_array = parameter_identifier(parameter_list,combined_logspace_parameters);
    indices_to_multistart = parameter_identifier(parameter_list,ms_grid_parameter_array);
    test_strain_ML_file = fullfile(datafile_path,strcat('test_strain_params-',...
        output_id_parameter,'_',int2str(external_counter),'.csv'));
    strain_LL_table_file = fullfile(datafile_path,strcat('current_best_strainlooper_LL-',...
        output_id_parameter,'_',int2str(external_counter),'.csv'));
    % set default parameters for the following if they are not provided
    if isKey(parameter_dict,'ms_starting_point_file')
        ms_starting_point_file = parameter_dict('ms_starting_point_file');
        if exist(ms_starting_point_file, 'file')
            ms_starting_point_mat = csvread(ms_starting_point_file);
        else
            ms_starting_point_mat = NaN;
        end
    else
        ms_starting_point_file = NaN;
        ms_starting_point_mat = NaN;
    end
    if isKey(parameter_dict,'tolx_val')
        tolx_val = parameter_dict('tolx_val');
    else
        tolx_val = 10^-3;
    end
    if isKey(parameter_dict,'tolfun_val')
        tolfun_val = parameter_dict('tolfun_val');
    else
        tolfun_val = 10^-2;
    end
    if isKey(parameter_dict,'initial_data_fraction')
        initial_data_fraction = parameter_dict('initial_data_fraction');
    else
        initial_data_fraction = 1;
    end
    if isKey(parameter_dict,'write_solutions_output')
        write_solutions_output = parameter_dict('write_solutions_output');
    else
        write_solutions_output = false;
    end

    % rescale parameters and convert to logspace as needed
    combined_min_array = value_rescaler(combined_min_array_unscaled,combined_logspace_array,combined_scaling_array);
    combined_max_array = value_rescaler(combined_max_array_unscaled,combined_logspace_array,combined_scaling_array);
    combined_start_values_array = value_rescaler(combined_start_values_array_unscaled,combined_logspace_array,combined_scaling_array);
    combined_profile_ub_array = value_rescaler(combined_profile_ub_array_unscaled,combined_logspace_array,combined_scaling_array);
    combined_profile_lb_array = value_rescaler(combined_profile_lb_array_unscaled,combined_logspace_array,combined_scaling_array);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%    csv_output_prename,output_folder,


    % Read phenotype data for petites
    response_name = 'GR';
    fixef_name = 'Strain';
    random_effect_names = {'plate','well','field'};
    block_effect_name = 'plate';

    [response_vector,fixef_ID_mat,ranef_corr_struct,unique_fixefs,...
        unique_ranef_categories,order_vector,block_start_positions] = mixef_data_reader(petite_file,response_name,...
        fixef_name,random_effect_names,block_effect_name,initial_data_fraction);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read GR difference phenotype data

    data_table = readtable(phenotype_file);
    indices_to_use = logical(binornd(1,initial_data_fraction,[1 size(data_table,1)]));
    data_table = data_table(indices_to_use,:);

    % assign current data to cell arrays

%    well_list = data_table.well;
%    plate_list = data_table.plate;
%    ref_strain_list = data_table.ref_Strain;
    test_strain_list_by_pair = data_table.MA_Strain;
    GR_diff_list = data_table.GR_diff;

    % parameter_list (and combined parameter arrays) follow format:
        % global parameters, strainwise mutation effects,
        % strainwise petite effects

    global_param_number = 8;
    strain_params = parameter_list((global_param_number+1):end);
    strain_number = length(strain_params)/2;
    mut_effect_parameter_list = strain_params(1:strain_number);
    %strain_list = regexprep(mut_effect_parameter_list,'-.*$','');
    strain_list = regexprep(mut_effect_parameter_list,'_.\w$','');
        % removes everything after final underscore followed only by word
            % characters (including the underscore)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if length(combined_position_array)==1
        combined_position_array = repmat(combined_position_array(1),size(parameter_list));
    end

    global_index_list = 1:global_param_number;
    mut_effect_index_list = (global_param_number+1):...
        (global_param_number+strain_number);
    petite_prop_index_list = (global_param_number+strain_number+1):...
        (global_param_number+2*strain_number);

    global_fixed_parameter_array = ...
        combined_fixed_parameter_array(global_index_list);
    global_min_array = ...
        combined_min_array(global_index_list);
    global_max_array = ...
        combined_max_array(global_index_list);
    global_length_array = ...
        combined_length_array(global_index_list);
    global_position_array = ...
        combined_position_array(global_index_list);
    global_start_values = ...
        combined_start_values_array(global_index_list);
    global_logspace_array = ...
        combined_logspace_array(global_index_list);
    global_scaling_array = ...
        combined_scaling_array(global_index_list);

    mut_effect_fixed_parameter_array = ...
        combined_fixed_parameter_array(mut_effect_index_list);
    mut_effect_min_array = ...
        combined_min_array(mut_effect_index_list);
    mut_effect_max_array = ...
        combined_max_array(mut_effect_index_list);
    mut_effect_length_array = ...
        combined_length_array(mut_effect_index_list);
    mut_effect_position_array = ...
        combined_position_array(mut_effect_index_list);
    mut_effect_start_values = ...
        combined_start_values_array(mut_effect_index_list);
    mut_effect_logspace_array = ...
        combined_logspace_array(mut_effect_index_list);
    mut_effect_scaling_array = ...
        combined_scaling_array(mut_effect_index_list);

    petite_prop_fixed_parameter_array = ...
        combined_fixed_parameter_array(petite_prop_index_list);
    petite_prop_min_array = ...
        combined_min_array(petite_prop_index_list);
    petite_prop_max_array = ...
        combined_max_array(petite_prop_index_list);
    petite_prop_length_array = ...
        combined_length_array(petite_prop_index_list);
    petite_prop_position_array = ...
        combined_position_array(petite_prop_index_list);
    petite_prop_start_values = ...
        combined_start_values_array(petite_prop_index_list);
    petite_prop_logspace_array = ...
        combined_logspace_array(petite_prop_index_list);
    petite_prop_scaling_array = ...
        combined_scaling_array(petite_prop_index_list);

    global_profile_lb_array = ...
        combined_profile_lb_array(global_index_list);
    global_profile_ub_array = ...
        combined_profile_ub_array(global_index_list);
    mut_effect_profile_lb_array = ...
        combined_profile_lb_array(mut_effect_index_list);
    mut_effect_profile_ub_array = ...
        combined_profile_ub_array(mut_effect_index_list);
    petite_prop_profile_lb_array = ...
        combined_profile_lb_array(petite_prop_index_list);
    petite_prop_profile_ub_array = ...
        combined_profile_ub_array(petite_prop_index_list);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set up a list of fixed effect values for 'global' parameters (these
        % are: sigma_colony , petite_GR, ref_GR, ref_petite_proportion), as well
        % as strain_mut_effects and strain_petite_props

    % Create an array indicating which fixed parameters need to be
        % created on a log scale, rather than a linear scale
%    global_logspace_array = false(size(global_fixed_parameter_array));
%    mut_effect_logspace_array = false(size(mut_effect_fixed_parameter_array));
%    petite_prop_logspace_array = false(size(petite_prop_fixed_parameter_array));


    global_fixed_parameter_values = ...
        fixed_parameter_processor(global_fixed_parameter_array,global_profile_lb_array,...
            global_profile_ub_array,global_length_array,global_position_array);
    mut_effect_fixed_parameter_values = ...
        fixed_parameter_processor(mut_effect_fixed_parameter_array,mut_effect_profile_lb_array,...
            mut_effect_profile_ub_array,mut_effect_length_array,mut_effect_position_array);
    petite_prop_fixed_parameter_values = ...
        fixed_parameter_processor(petite_prop_fixed_parameter_array,petite_prop_profile_lb_array,...
            petite_prop_profile_ub_array,petite_prop_length_array,petite_prop_position_array);

    global_fixed_parameter_indices = logical(global_fixed_parameter_array);
%    mut_effect_fixed_parameter_indices = logical(mut_effect_fixed_parameter_array);
%    petite_prop_fixed_parameter_indices = logical(petite_prop_fixed_parameter_array);

%    % use linear space if fitted parameter is q or p0; otherwise, use log space
%    if (strcmp('q_1',current_fixed_parameter) || strcmp('p0_1',current_fixed_parameter))
%        fixed_parameter_values = fixed_parameter_values*linear_multiplier;
%    else
%        fixed_parameter_values = log(fixed_parameter_values);
%    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up lower and upper bounds for parameters

    global_lower_bounds_fitted = global_min_array(~global_fixed_parameter_indices);
    global_upper_bounds_fitted = global_max_array(~global_fixed_parameter_indices);
    global_start_vals_fitted = global_start_values(~global_fixed_parameter_indices);

%    strain_mut_effect_lower_bounds = mut_effect_min_array;
%    strain_mut_effect_upper_bounds = mut_effect_max_array;
%    strain_mut_effect_start_vals = mut_effect_start_values;

%    strain_petite_proportion_lower_bounds = petite_prop_min_array;
%    strain_petite_proportion_upper_bounds = petite_prop_max_array;
%    strain_petite_proportion_start_vals = petite_prop_start_values;

%    global_param_number_fitted = length(global_lower_bounds_fitted);

%    global_profile_lb_fitted = global_profile_lb_array(~global_fixed_parameter_indices);
%    global_profile_ub_fitted = global_profile_ub_array(~global_fixed_parameter_indices);
%    strain_mut_effect_profile_lb = mut_effect_profile_lb_array;
%    strain_mut_effect_profile_ub = mut_effect_profile_ub_array;
%    strain_petite_proportion_profile_lb = petite_prop_profile_lb_array;
%    strain_petite_proportion_profile_ub = petite_prop_profile_ub_array;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up global search structure
    gs = GlobalSearch;
    gs.TolFun = tolfun_val;
    gs.TolX = tolx_val;
    gs.NumStageOnePoints = 2^global_param_number; % default = 200; sufficient is 50
    gs.NumTrialPoints = 3^global_param_number; % default = 1000; sufficient is 200
    gs.StartPointsToRun='bounds';
    gs.Display='final';

    ms = MultiStart(gs);
        % assign same Tol, StartPointsToRun, and Display properties to ms as to gs
    if parallel_processors > 1
        ms.UseParallel = 1;
    else
        ms.UseParallel = 0;
    end

    fmincon_opts = optimoptions('fmincon','TolX',tolx_val,'TolFun',tolfun_val,...
        'Algorithm','interior-point','MaxIter',5000,'MaxFunEvals',12000,...
        'SpecifyObjectiveGradient',true,'CheckGradients',false,'Display','off');

    % To find parameters corresponding to each strain, use simple
        % space search with one start point or global search with
        % multiple starts?
    %strainwise_search_type = 'local';   
    strainwise_search_type = 'global';    

    global_start_time = tic;

    min_problem_fixed_params = createOptimProblem('fmincon','objective',...
        @(v) LL_calculator_strain_looper_pairwise_v3(v,...
            global_fixed_parameter_indices,global_fixed_parameter_values,...
            mut_effect_fixed_parameter_values,...
            petite_prop_fixed_parameter_values,...
            test_strain_ML_file,strain_LL_table_file,...
            mut_effect_min_array,mut_effect_max_array,...
            petite_prop_min_array,petite_prop_max_array,...
            mut_effect_start_values,petite_prop_start_values,strain_list,...
            test_strain_list_by_pair,GR_diff_list,strainwise_search_type,max_neg_LL_val,...
            response_vector,fixef_ID_mat,ranef_corr_struct,unique_fixefs,...
            unique_ranef_categories,order_vector,block_start_positions,...
            tolx_val,tolfun_val,...
            global_logspace_array,mut_effect_logspace_array,petite_prop_logspace_array,...
            global_scaling_array,mut_effect_scaling_array,petite_prop_scaling_array),...
        'x0',global_start_vals_fitted,'lb',global_lower_bounds_fitted,'ub',global_upper_bounds_fitted,...
        'options',fmincon_opts);
    
    if ms_positions == 0
        % use globalsearch algorithm
        [vout_current,fval_current,~,~,solutions_current]=...
            run(gs,min_problem_fixed_params);
    elseif ms_positions == 1
        % only use single startpoint, at global_start_values
        [vout_current,fval_current,~]=...
            fmincon(min_problem_fixed_params);
    else

        if parallel_processors > 1
            parpool('local',parallel_processors);
        end

        %ms_positions_across_dimensions = ms_positions^global_param_number+1;
%        indices_to_multistart = [false true true true true false false false];
        indices_to_multistart_fitted = indices_to_multistart(~global_fixed_parameter_indices);

        startpoints = MS_Startposition_Generator_v2(indices_to_multistart_fitted,ms_positions,...
            global_start_vals_fitted,global_lower_bounds_fitted,global_upper_bounds_fitted,ms_starting_point_mat);
        [vout_current,fval_current,~,~,solutions_current]=...
            run(ms,min_problem_fixed_params,startpoints);
    end
    
    vout_current_corrected = zeros(size(global_fixed_parameter_indices));
    vout_current_corrected(global_fixed_parameter_indices) = global_fixed_parameter_values(~isnan(global_fixed_parameter_values));
    vout_current_corrected(~global_fixed_parameter_indices) = vout_current;
    vout_current_corrected = reverse_value_scaler(vout_current_corrected,global_logspace_array,global_scaling_array);

    % return log-space parameters back to their real values
%    for parameter_name_counter = 1:length(parameter_names)
%        if (strcmp('q_1',parameter_names{parameter_name_counter}) || strcmp('p0_1',parameter_names{parameter_name_counter}))
%            vout_current_corrected(parameter_name_counter)=vout_current_corrected(parameter_name_counter)/linear_multiplier;
%        elseif (strcmp('combined_mean_2',parameter_names{parameter_name_counter}) || strcmp('combined_sd_2',parameter_names{parameter_name_counter}))
%            vout_current_corrected(parameter_name_counter)=vout_current_corrected(parameter_name_counter); % do nothing
%        else
%            vout_current_corrected(parameter_name_counter)=exp(vout_current_corrected(parameter_name_counter));
%        end
%    end
    
    runtime = toc(global_start_time);

    % close current parallel pool
    delete(gcp('nocreate'));

    % get data for strain parameters
    test_strain_ML_param_table = readtable(test_strain_ML_file);
    
    table_data = num2cell([-fval_current,runtime,vout_current_corrected,...
        test_strain_ML_param_table.Mutation_Effects',...
        test_strain_ML_param_table.Petite_Proportions']);

    % if doing mutlistart or global search and one solution didn't converge, replace export_data with NA
    if exist('solutions_current','var') == 1
        for solution_counter = 1:length(solutions_current)
            temp_exit = solutions_current(solution_counter).Exitflag;
            if temp_exit < 1
                table_data = NaN(size(table_data));
            end
        end
        % write content of solutions_current to ms_starting_point_file
        if ~isnan(ms_starting_point_file) & write_solutions_output
            number_fitted_parameters = size(solutions_current(1).X,2);
            num_solutions = size(solutions_current,2);
            point_mat_unsorted = reshape([solutions_current(:).X],[number_fitted_parameters,num_solutions])';
%            [~,sorted_by_LL_indices] = sort(-[solutions_current(:).Fval],'descend');
%            point_mat = point_mat_unsorted(sorted_by_LL_indices);
            point_mat = point_mat_unsorted;
            compressed_point_mat = Combine_Nearby_Points(point_mat, tolx_val*1.01);
            csvwrite(ms_starting_point_file,compressed_point_mat);
        end
    end



    table_headers = [{'LL','runtime_in_secs'} ...
        regexprep(parameter_list,'[\.,-]','_')];
        % need to replace periods and dashes in parameter_list with
            % underscores to make them valid variable names
    T = table(table_data{:},'VariableNames',table_headers);
    writetable(T,output_file);
    
end

