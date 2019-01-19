function output_dict = pre_MLE_GR_diff_strainwise_mixef_me_dist(input_value_dict)

    output_dict_GR_diff = pre_MLE_GR_diff(input_value_dict);
    output_dict_mixef = pre_MLE_mixef(input_value_dict);
    output_dict_me_dist = pre_MLE_me_dist(input_value_dict);

    output_dict = [output_dict_GR_diff; output_dict_mixef; ...
    	output_dict_me_dist];

end