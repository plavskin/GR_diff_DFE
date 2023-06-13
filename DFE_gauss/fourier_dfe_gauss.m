function [Fa, Fa_gradient_dict, complete_me_parameter_list, ...
        fitted_parameters_me_fullname] = ...
    fourier_dfe_gauss(param_vals, input_value_dict, pre_MLE_output_dict)
    
    % Creates Fa, the fourier transform of a gaussian distribution, and
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
    mu_unseq_muts = parameter_dict('mu_unseq_muts');
    sigma_unseq_muts = parameter_dict('sigma_unseq_muts');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    complete_me_parameter_list = ...
        mle_parameter_names(contains(mle_parameter_names, '_unseq_muts'));

    fitted_parameters_unseq_muts_fullname = ...
        fitted_parameters(contains(fitted_parameters, '_unseq_muts'));
    fitted_parameters_unseq_muts = ...
        strrep(fitted_parameters_unseq_muts_fullname, '_unseq_muts', '');

    fitted_parameters_me_fullname = intersect(complete_me_parameter_list, ...
        fitted_parameters);

    % Calculate likelihood of observing mutation effects given current
        % global parameters

    % Calculate Fa, the fourier transform of the pdf of
        % observing combined unsequenced mutational effects drawn from a
        % normal distribution with a mean of mu_unseq_muts per strain
        % and standard deviation of sigma_unseq_muts across strains
    [Fa, Fa_gradient_dict_temp_names] = ...
        fourier_domain_gauss(me_pdf_frequency_vals, mu_unseq_muts, ...
            sigma_unseq_muts, fitted_parameters_unseq_muts, ...
            gradient_specification);
    
    % Change name of dictionary keys to original names of parameters
    if gradient_specification
        corrected_unseq_muts_gradient_dict_keys =  ...
            strcat(keys(Fa_gradient_dict_temp_names), '_unseq_muts');
        Fa_gradient_dict = ...
            containers.Map(corrected_unseq_muts_gradient_dict_keys, ...
                values(Fa_gradient_dict_temp_names));
    else
        Fa_gradient_dict = containers.Map('KeyType', 'char', ...
            'ValueType', 'any');
    end

end
    
