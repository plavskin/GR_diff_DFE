function output_dict = pre_MLE_me_estimate_and_dist(input_value_dict)

    output_dict_me_estimate = pre_MLE_me_estimate(input_value_dict);
    output_dict_me_dist = pre_MLE_me_dist(input_value_dict);

    output_dict = [output_dict_me_estimate; output_dict_me_dist];

end