function output_dict = pre_MLE_me_estimate(input_value_dict)
    % requires mixef_data_reader, effect_position_finder

    % Read phenotype data
    phenotype_file = input_value_dict('phenotype_file');

    data_table = readtable(phenotype_file);

    strain_list = data_table.MA_Strain;
    strain_me_estimate_list = data_table.param_MLE;
    strain_me_ste_list = data_table.param_ste;

    output_dict = containers.Map({'strain_list', 'strain_me_estimate_list', ...
            'strain_me_ste_list'}, ...
        {strain_list, strain_me_estimate_list, strain_me_ste_list});

end

