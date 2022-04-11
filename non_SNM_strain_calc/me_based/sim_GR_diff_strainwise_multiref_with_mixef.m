function sim_GR_diff_strainwise_multiref_with_mixef(key_list, value_list)

	input_value_dict_mixef = ...
        containers.Map(key_list, value_list);
    parameter_list_mixef = input_value_dict_mixef('parameter_list');
    random_effect_names = input_value_dict_mixef('random_effect_names');
    ranef_names_from_param_list = strcat(random_effect_names,'_sigma');
    params_to_replace = [{'petite_mean'}, ranef_names_from_param_list];
    replacement_params = [{'petite'}, random_effect_names];
    [~, replacement_idx] = ismember(params_to_replace, parameter_list_mixef);
    parameter_list_mixef(replacement_idx(replacement_idx>0)) = replacement_params;
    input_value_dict_mixef('parameter_list') = parameter_list_mixef;

    mixef_keys = keys(input_value_dict_mixef);
    mixef_vals = values(input_value_dict_mixef);
	sim_mixef(mixef_keys, mixef_vals);

	GR_diff_keys = [key_list, {'output_path', 'output_file_label'}];
	GR_diff_vals = [value_list, {'', ''}];
	sim_pairwise_GR_diff_multiref(GR_diff_keys, GR_diff_vals);

end