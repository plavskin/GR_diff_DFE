function output_dict = pre_MLE_mixef(input_value_dict)
    % requires mixef_data_reader, effect_position_finder

    % Read phenotype data for petites
    response_name = input_value_dict('response_name');
    fixef_category = input_value_dict('fixef_category');
    random_effect_names = input_value_dict('random_effect_names');
    block_effect_name = input_value_dict('block_effect_name');
    data_file = input_value_dict('data_file');

    if isKey(input_value_dict,'initial_data_fraction')
        initial_data_fraction = input_value_dict('initial_data_fraction');
    else
        initial_data_fraction = 1;
    end

    [response_vector, fixef_ID_mat, ranef_corr_struct, unique_fixefs,...
        unique_ranef_categories, order_vector, block_start_positions] = ...
        mixef_data_reader(data_file, response_name, fixef_category, ...
        random_effect_names, block_effect_name, initial_data_fraction);

    output_dict = containers.Map({'response_vector', 'fixef_ID_mat', ...
        'ranef_corr_struct', 'unique_fixefs','unique_ranef_categories', ...
        'order_vector', 'block_start_positions'}, ...
        {response_vector, fixef_ID_mat, ranef_corr_struct, unique_fixefs,...
        unique_ranef_categories, order_vector, block_start_positions});