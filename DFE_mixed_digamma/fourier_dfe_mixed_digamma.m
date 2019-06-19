function [Fa, Fa_gradient_dict, complete_me_parameter_list, ...
        fitted_parameters_me_fullname] = ...
    fourier_dfe_mixed_digamma(param_vals, input_value_dict, pre_MLE_output_dict)
    
    % Creates Fa, the fourier transform of a distribution of a poisson
    % number of mutational effects each drawn from a reflected gamma
    % distribution in which a proportion of the variants are neutral, and
    % Fa_gradient_dict, a Map object containing the fourier transform of
    % the gradient of the likelihood of the distribution of fitness effects
    % with respect to each one of the distribution's parameters
    % Also returns lists of dfe-relevant parameters and fitted parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mle_parameter_names = input_value_dict('mle_parameter_names');
    gradient_specification = input_value_dict('gradient_specification');

    me_pdf_frequency_vals = pre_MLE_output_dict('me_pdf_frequency_vals');
    fitted_parameters = pre_MLE_output_dict('fitted_parameters');
    
    parameter_dict = containers.Map(mle_parameter_names, param_vals);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify parameter values
    
    % mutational effect distribution parameters
    mu_SNM = parameter_dict('mu_SNM');
    shape_SNM = parameter_dict('shape_SNM');
    prop_pos_SNM = parameter_dict('prop_pos_SNM');
    lambda_SNM = parameter_dict('lambda_SNM');
    prop_neut_SNM = parameter_dict('prop_neut_SNM');

    mu_unseq_muts = parameter_dict('mu_unseq_muts');
    shape_unseq_muts = parameter_dict('shape_unseq_muts');
    prop_pos_unseq_muts = parameter_dict('prop_pos_unseq_muts');
    lambda_unseq_muts = parameter_dict('lambda_unseq_muts');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    complete_me_parameter_list = ...
        [mle_parameter_names(contains(mle_parameter_names, '_SNM')), ...
        mle_parameter_names(contains(mle_parameter_names, '_unseq_muts'))];

    fitted_parameters_SNM_fullname = ...
        fitted_parameters(contains(fitted_parameters, '_SNM'));
    fitted_parameters_SNM = ...
        strrep(fitted_parameters_SNM_fullname, '_SNM', '');

    fitted_parameters_unseq_muts_fullname = ...
        fitted_parameters(contains(fitted_parameters, '_unseq_muts'));
    fitted_parameters_unseq_muts = ...
        strrep(fitted_parameters_unseq_muts_fullname, '_unseq_muts', '');

    fitted_parameters_me_fullname = intersect(complete_me_parameter_list, ...
        fitted_parameters);

    % Calculate likelihood of observing mutation effects given current
        % global parameters

    % Calculate Fz_SNM, the fourier transform of the pdf of observing a
        % single mutational effect from a digamma distribution
    % Then calculate Fa_SNM, the fourier transform of the pdf of observing
        % a poisson random number of mutational effects from the
        % digamma distribution in Fz_SNM
    [Fz_SNM, Fz_SNM_gradient_dict_temp_names] = fourier_domain_digamma(...
        me_pdf_frequency_vals, mu_SNM, shape_SNM, prop_pos_SNM, ...
        fitted_parameters_SNM, gradient_specification);
    [Fa_SNM, Fa_SNM_gradient_dict_temp_names] = ...
        fourier_domain_poisson_dist_num_with_neutrals(Fz_SNM, ...
            me_pdf_frequency_vals, lambda_SNM, prop_neut_SNM, ...
                Fz_SNM_gradient_dict_temp_names, fitted_parameters_SNM, ...
                gradient_specification);


    % Calculate Fz_unseq_muts, the fourier transform of the pdf of observing a
        % single mutational effect from a digamma distribution with unsequenced
        % mutation parameters
    % Then calculate Fa_unseq_muts, the fourier transform of the pdf of observing
        % a poisson random number of mutational effects from the
        % digamma distribution in Fz_unseq_muts
    [Fz_unseq_muts, Fz_unseq_muts_gradient_dict_temp_names] = ...
        fourier_domain_digamma(me_pdf_frequency_vals, mu_unseq_muts, ...
            shape_unseq_muts, prop_pos_unseq_muts, ...
        fitted_parameters_unseq_muts, gradient_specification);
    [Fa_unseq_muts, Fa_unseq_muts_gradient_dict_temp_names] = ...
        fourier_domain_poisson_dist_num(Fz_unseq_muts, ...
            me_pdf_frequency_vals, lambda__unseq_muts, ...
            Fz_unseq_muts_gradient_dict_temp_names, ...
            fitted_parameters_unseq_muts, gradient_specification);
    
    % Change name of dictionary keys to original names of parameters
    if gradient_specification
        corrected_SNM_gradient_dict_keys =  ...
            strcat(keys(Fa_SNM_gradient_dict_temp_names), '_SNM');
        Fa_SNM_gradient_dict = containers.Map(corrected_SNM_gradient_dict_keys, ...
            values(Fa_SNM_gradient_dict_temp_names));
        corrected_unseq_muts_gradient_dict_keys =  ...
            strcat(keys(Fa_unseq_muts_gradient_dict_temp_names), '_unseq_muts');
        Fa_unseq_muts_gradient_dict = ...
            containers.Map(corrected_unseq_muts_gradient_dict_keys, ...
                values(Fa_unseq_muts_gradient_dict_temp_names));
    else
        Fa_SNM_gradient_dict = containers.Map('KeyType', 'char', ...
            'ValueType', 'any');
        Fa_unseq_muts_gradient_dict = containers.Map('KeyType', 'char', ...
            'ValueType', 'any');
    end

    % Create fourier transform of distribution of mutational effects in
        % a strain containing both SNMs and unsequenced mutations
    [Fa, Fa_gradient_dict] = derivative_multiplier(Fa_SNM, Fa_unseq_muts, ...
        Fa_SNM_gradient_dict, Fa_unseq_muts_gradient_dict);

end
    
