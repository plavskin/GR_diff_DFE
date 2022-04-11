function output_dict = pre_MLE_GR_diff_non_SNM_strainwise_direct(input_value_dict)
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
    ref_strain_list_by_pair = data_table.ref_Strain;
    GR_diff_list = data_table.GR_diff;

    % make matrix of indices of every unique test strain and reference strain name
    [test_strain_ID_mat, unique_test_strains] = ...
        effect_position_finder(test_strain_list_by_pair);
    [ref_strain_ID_mat, unique_ref_strains] = ...
        effect_position_finder(ref_strain_list_by_pair);

    strain_list_from_data = ...
        [reshape(unique_test_strains, [1 length(unique_test_strains)]), ...
        reshape(unique_ref_strains, [1 length(unique_ref_strains)])];

    % Check whether parameters in parameter_list with _pp endings and _gr endings match
    gr_parameter_indices = ~cellfun(@isempty,regexp(parameter_list,'_gr$'));
    pp_parameter_indices = ~cellfun(@isempty,regexp(parameter_list,'_pp$'));
    gr_parameters_in_list = parameter_list(gr_parameter_indices);
    pp_parameters_in_list = parameter_list(pp_parameter_indices);

    strain_list_from_gr_param_list = regexprep(gr_parameters_in_list,'(_gr)$','');
    strain_list_from_pp_param_list = regexprep(pp_parameters_in_list,'(_pp)$','');

    set_diff_gr_pp = [setdiff(strain_list_from_gr_param_list, strain_list_from_pp_param_list), ...
        setdiff(strain_list_from_pp_param_list, strain_list_from_gr_param_list)];
    if ~isempty(set_diff_gr_pp)
        disp('Error:')
        disp(set_diff_gr_pp)
        error('Error: lists of strains with *_pp and *_gr parameter names dont match');
    end

    % Check whether there are any [strain]_gr parameters in parameter_list that don't exist in strain_list_from_data
    set_diff_param_vs_data_strains = setdiff(strain_list_from_gr_param_list, strain_list_from_data);
    if ~isempty(set_diff_param_vs_data_strains)
        disp('Error:')
        disp(set_diff_param_vs_data_strains)
        error('The subset of strains with *_pp and *_gr parameter names in parameter_list above cant be found in phenotype data');
    end

    test_strain_ML_file = fullfile(output_file_path,strcat('test_strain_params-',...
        output_file_label,'_',int2str(external_counter),'.csv'));
    strain_LL_table_file = fullfile(output_file_path,strcat('current_best_strainlooper_LL-',...
        output_file_label,'_',int2str(external_counter),'.csv'));

    output_dict = containers.Map({'strain_list', 'unique_test_strains', ...
            'unique_ref_strains', 'test_strain_ID_mat', 'ref_strain_ID_mat', ...
            'GR_diff_list', 'test_strain_ML_file', 'strain_LL_table_file'}, ...
        {strain_list_from_gr_param_list, unique_test_strains, ...
            unique_ref_strains, test_strain_ID_mat, ref_strain_ID_mat, ...
            GR_diff_list, test_strain_ML_file, strain_LL_table_file});

end

