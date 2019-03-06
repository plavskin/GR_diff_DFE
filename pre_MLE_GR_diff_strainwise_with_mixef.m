function output_dict = pre_MLE_GR_diff_strainwise_with_mixef(input_value_dict)

    output_dict_GR_diff = pre_MLE_GR_diff_strainwise(input_value_dict);
    output_dict_mixef = pre_MLE_mixef(input_value_dict);

    output_dict = [output_dict_GR_diff; output_dict_mixef];

end