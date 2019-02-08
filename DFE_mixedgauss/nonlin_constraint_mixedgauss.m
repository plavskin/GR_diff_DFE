function [c, ceq, c_gradient, ceq_gradient, gradient_key] = ...
            nonlin_constraint_mixedgauss(param_vals, input_value_dict, ...
                ~)

    % given gamma distribution parameters and added error value,
        % calculates non-linear constraint
    % also returns gradients along all fitted parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mle_parameter_names = input_value_dict('mle_parameter_names');
    gradient_specification = input_value_dict('gradient_specification');
    L = input_value_dict('L');
    sigma_kernel = input_value_dict('smoothing_kernel_sigma');
    
    parameter_dict = containers.Map(mle_parameter_names, param_vals);

    % mutational effect distribution parameters
    mu_SNM = parameter_dict('mu_SNM');
    shape_SNM = parameter_dict('shape_SNM');
    lambda_SNM = parameter_dict('lambda_SNM');
    prop_neut_SNM = parameter_dict('prop_neut_SNM');

    mu_unseq_muts = parameter_dict('mu_unseq_muts');
    sigma_unseq_muts = parameter_dict('sigma_unseq_muts');

    effective_lambda_SNM = (1 - prop_neut_SNM) * lambda_SNM;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % no equality constraint
    ceq = [];
    ceq_gradient=[];

    % inequality constraint: c must be less than or equal to 0
    % this is a conservative limit in cases when prop_pos_SNM < ~1,
        % since it doesn't account for the decrease in density on a
        % given side of 0

    a = mu_SNM ^ 2 * effective_lambda_SNM / shape_SNM + sigma_kernel ^ 2 + ...
        sigma_unseq_muts ^ 2;
    sqrt_a = sqrt(a);
    c = mu_SNM * effective_lambda_SNM - mu_unseq_muts +  4 * sqrt_a - L / 2;

    if gradient_specification % check that gradient required

        d_c_d_prop_pos_SNM = 0;
        d_c_d_effective_lambda_SNM = mu_SNM + (2 / sqrt_a) * mu_SNM ^ 2 / shape_SNM;
        d_effective_lambda_SNM_d_lambda_SNM = 1 - prop_neut_SNM;
        d_effective_lambda_SNM_d_prop_neut_SNM = -lambda_SNM;
        d_c_d_lambda_SNM = d_c_d_effective_lambda_SNM * d_effective_lambda_SNM_d_lambda_SNM;
        d_c_d_prop_neut_SNM = d_c_d_effective_lambda_SNM * d_effective_lambda_SNM_d_prop_neut_SNM;
        d_c_d_mu_SNM = ...
            effective_lambda_SNM + (2 / sqrt_a) * 2 * mu_SNM * effective_lambda_SNM / shape_SNM;
        d_c_d_shape_SNM = ...
            (2 / sqrt_a) * mu_SNM ^ 2 * effective_lambda_SNM * (-1) / shape_SNM ^ 2;
        d_c_d_mu_unseq_muts = -1;
        d_c_d_sigma_unseq_muts = (2 / sqrt_a) * 2 * sigma_kernel;

        c_gradient = ...
            [d_c_d_prop_pos_SNM, d_c_d_lambda_SNM, d_c_d_mu_SNM, ...
            d_c_d_shape_SNM, d_c_d_prop_neut_SNM, d_c_d_mu_unseq_muts, ...
            d_c_d_sigma_unseq_muts];
        gradient_key = {'prop_pos_SNM', 'lambda_SNM', 'mu_SNM', 'shape_SNM', ...
            'prop_neut_SNM', 'mu_unseq_muts', 'sigma_unseq_muts'};

    end

end