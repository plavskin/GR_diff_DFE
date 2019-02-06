function sim_GR_diff_strainwise_with_mixef_and_mut_eff(key_list, value_list)

    % simulate linear mixed effect modeled portion
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

    % simulate mutation effects using appropriate distribution of
        % fitness effects
    input_value_dict_GR = containers.Map(key_list, value_list);
    if strcmp(input_value_dict_GR('model'), 'lambdafit')
        new_param_values = sim_mut_effect_lambdafit(key_list, value_list);
    elseif strcmp(input_value_dict_GR('model'), 'poissonmut')
        new_param_values = sim_mut_effect_poissonmut(key_list, value_list);
    end
    input_value_dict_GR('starting_parameter_vals') = new_param_values;

    % simulate individual growth rates based on global parameters and
        % using the mutation effects simulated above
	GR_diff_keys = [keys(input_value_dict_GR), {'output_path', 'output_file_label'}];
	GR_diff_vals = [values(input_value_dict_GR), {'', ''}];
	sim_pairwise_GR_diff(GR_diff_keys, GR_diff_vals);

end