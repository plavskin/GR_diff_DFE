function output_dict = pre_MLE_GR_diff(input_value_dict)
    % requires mixef_data_reader, effect_position_finder

    % Read phenotype data
    external_counter = str2num(input_value_dict('external_counter'));
    parameter_list = input_value_dict('parameter_list');
    output_file_path = input_value_dict('output_path');
    output_file_label = input_value_dict('output_file_label');

    phenotype_file = input_value_dict('phenotype_file');

    if isKey(input_value_dict,'initial_data_fraction')
        initial_data_fraction = input_value_dict('initial_data_fraction');
    else
        initial_data_fraction = 1;
    end

    data_table = readtable(phenotype_file);
    indices_to_use = logical(binornd(1,initial_data_fraction,[1 size(data_table,1)]));
    data_table = data_table(indices_to_use,:);

    test_strain_list_by_pair = data_table.MA_Strain;
    GR_diff_list = data_table.GR_diff;
    strain_list_from_data = unique(test_strain_list_by_pair);

    % Check whether parameters in parameter_list with _pp endings and _me endings match
%    me_parameters_in_list = intersect(parameter_list, mut_effect_mle_parameters);
%    pp_parameters_in_list = intersect(parameter_list, petite_prop_mle_parameters);
    me_parameter_indices = ~cellfun(@isempty,regexp(parameter_list,'_me$'));
    pp_parameter_indices = ~cellfun(@isempty,regexp(parameter_list,'_pp$'));
    me_parameters_in_list = parameter_list(me_parameter_indices);
    pp_parameters_in_list = parameter_list(pp_parameter_indices);

    strain_list_from_me_param_list = regexprep(me_parameters_in_list,'(_me)$','');
    strain_list_from_pp_param_list = regexprep(pp_parameters_in_list,'(_pp)$','');

    set_diff_me_pp = [setdiff(strain_list_from_me_param_list, strain_list_from_pp_param_list), ...
        setdiff(strain_list_from_pp_param_list, strain_list_from_me_param_list)];
    if ~isempty(set_diff_me_pp)
        disp('Error:')
        disp(set_diff_me_pp)
        error('Error: lists of strains with *_pp and *_me parameter names dont match');
    end

    % Check whether there are any [strain]_me parameters in parameter_list that don't exist in strain_list_from_data
    set_diff_param_vs_data_strains = setdiff(strain_list_from_me_param_list, strain_list_from_data);
    if ~isempty(set_diff_param_vs_data_strains)
        disp('Error:')
        disp(set_diff_param_vs_data_strains)
        error('The subset of strains with *_pp and *_me parameter names in parameter_list above cant be found in phenotype data');
    end

    test_strain_ML_file = fullfile(output_file_path,strcat('test_strain_params-',...
        output_file_label,'_',int2str(external_counter),'.csv'));
    strain_LL_table_file = fullfile(output_file_path,strcat('current_best_strainlooper_LL-',...
        output_file_label,'_',int2str(external_counter),'.csv'));

    output_dict = containers.Map({'strain_list', 'test_strain_list_by_pair', 'GR_diff_list', ...
        'test_strain_ML_file', 'strain_LL_table_file'}, ...
        {strain_list_from_me_param_list, test_strain_list_by_pair, GR_diff_list, ...
        test_strain_ML_file, strain_LL_table_file});

end

