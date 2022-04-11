function output_dict = pre_MLE_me_estimate(input_value_dict)
    % requires mixef_data_reader, effect_position_finder

    % Read phenotype data
    phenotype_file = input_value_dict('phenotype_file');

    data_table = readtable(phenotype_file);

    strain_list = data_table.MA_Strain;
    strain_me_estimate_list = data_table.param_MLE;
    strain_me_ste_list = data_table.param_ste;
    % Add 'lambda multiplier', e.g. for cases when MA across multiple generation times
    if sum(ismember(data_table.Properties.VariableNames,'lambda_multiplier'))
        [lambda_mult_list,~,lambda_mult_strain_idx_list] = ...
            unique(data_table.lambda_multiplier);
    else
        [lambda_mult_list,~,lambda_mult_strain_idx_list] = unique(ones(length(strain_me_estimate_list)));
    end

    output_dict = containers.Map({'strain_list', 'strain_me_estimate_list', ...
            'strain_me_ste_list','lambda_mult_list','lambda_mult_strain_idx_list'}, ...
        {strain_list, strain_me_estimate_list, strain_me_ste_list,lambda_mult_list,lambda_mult_strain_idx_list});

end

