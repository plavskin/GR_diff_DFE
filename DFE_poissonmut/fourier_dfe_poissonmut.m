function [Fa, Fa_gradient_dict, complete_me_parameter_list, ...
        fitted_parameters_me_fullname] = ...
    fourier_dfe_poissonmut(param_vals, input_value_dict, pre_MLE_output_dict)
    
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    complete_me_parameter_list = ...
        mle_parameter_names(contains(mle_parameter_names, '_SNM'));
    fitted_parameters_me_fullname = ...
        fitted_parameters(contains(fitted_parameters, '_SNM'));
    fitted_parameters_me = ...
        strrep(fitted_parameters_me_fullname, '_SNM', '');

    % Calculate likelihood of observing mutation effects given current
        % global parameters
    % Calculate Fz, the fourier transform of the pdf of observing a
        % single mutational effect from a digamma distribution
    % Then calculate Fa, the fourier transform of the pdf of observing
        % a poisson random number of mutational effects from the
        % digamma distribution in Fz
    [Fz, Fz_gradient_dict_temp_names] = fourier_domain_digamma(...
        me_pdf_frequency_vals, mu_SNM, shape_SNM, prop_pos_SNM, ...
        fitted_parameters_me, gradient_specification);
    [Fa, Fa_gradient_dict_temp_names] = ...
        fourier_domain_poisson_dist_num_with_neutrals(Fz, ...
            me_pdf_frequency_vals, lambda_SNM, prop_neut_SNM, ...
                Fz_gradient_dict_temp_names, fitted_parameters_me, ...
                gradient_specification);
    
    % Change name of dictionary keys to original names of parameters
    if gradient_specification
        corrected_gradient_dict_keys =  ...
            strcat(keys(Fa_gradient_dict_temp_names), '_SNM');
        Fa_gradient_dict = containers.Map(corrected_gradient_dict_keys, ...
            values(Fa_gradient_dict_temp_names));
    else
        Fa_gradient_dict = containers.Map('KeyType', 'char', ...
            'ValueType', 'any');
    end
    

end
    
