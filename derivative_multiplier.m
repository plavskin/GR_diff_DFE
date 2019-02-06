function [multiplied_vector, d_multiplied_v_dict] = derivative_multiplier(v1, v2, d_v1_dict, ...
    d_v2_dict)
    % calculates pointwise product of v1 and v1, and for derivatives of
    % v1 and v2 with respect to any parameters (listed in d_v1_dict and
    % d_v2_dict), calculates the derivative of the product of v1 and v2
    % with respect to that parameter

    multiplied_vector = v1 .* v2;
    
    d_v1_keys = keys(d_v1_dict);
    d_v2_keys = keys(d_v2_dict);
    combined_keys = union(d_v1_keys, d_v2_keys);

    d_multiplied_v_dict = containers.Map();

    for current_key_idx = 1:length(combined_keys)
        current_key = combined_keys{current_key_idx};
        if any(strcmp(current_key, d_v1_keys))
            d_v1_component = d_v1_dict(current_key) .* v2;
        else
            d_v1_component = 0;
        end
        if any(strcmp(current_key, d_v2_keys))
            d_v2_component = d_v2_dict(current_key) .* v1;
        else
            d_v2_component = 0;
        end
        d_multiplied_v_dict(current_key) = d_v1_component + d_v2_component;
    end

end