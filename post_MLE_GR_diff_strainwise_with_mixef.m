function output_dict = post_MLE_GR_diff_strainwise_with_mixef(input_value_dict, pre_MLE_output_dict, output_table)

    output_file = input_value_dict('output_file');
    test_strain_ML_file = pre_MLE_output_dict('test_strain_ML_file');
    test_strain_ML_table = readtable(test_strain_ML_file);

    strain_names = test_strain_ML_table.Strain_Names';
    me_mle_vals = test_strain_ML_table.Mutation_Effects;
    pp_mle_vals = test_strain_ML_table.Petite_Proportions;

    strain_param_name_list = [strcat(strain_names,'_me'), ...
        strcat(strain_names,'_pp')];
    strain_table_data = num2cell([me_mle_vals', pp_mle_vals']);
    strain_table = table(strain_table_data{:},'VariableNames', strain_param_name_list);

    final_table = [global_output_table strain_table];
    final_table = movevars(final_table,'point_num','After',strain_param_name_list{end});

    writetable(final_table,output_file);

end