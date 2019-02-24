function output_dict = post_MLE_GR_diff_strainwise_with_mixef(input_value_dict, pre_MLE_output_dict, output_table)

    output_file = input_value_dict('output_file');
    global_output_table = readtable(output_file);

    lower_level_param_ML_file = pre_MLE_output_dict('test_strain_ML_file');
    lower_level_param_ML_table = readtable(lower_level_param_ML_file);
    lower_level_param_ML_table_names = ...
        lower_level_param_ML_table.Properties.VariableNames;

    % retrieve DFE parameter info
    DFE_parameters = input_value_dict('DFE_parameters');

    final_table = [global_output_table lower_level_param_ML_table];
    final_table = ...
        movevars(final_table,'point_num', 'After', ...
            lower_level_param_ML_table_names{end});
    final_table = movevars(final_table, DFE_parameters, 'After', 'runtime_in_secs');

    writetable(final_table,output_file);

end