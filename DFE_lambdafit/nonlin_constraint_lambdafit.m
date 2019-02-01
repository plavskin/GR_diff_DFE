function [c, ceq, c_gradient, ceq_gradient, gradient_key] = ...
            nonlin_constraint_lambdafit(param_vals, input_value_dict, ...
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % no equality constraint
    ceq = [];
    ceq_gradient=[];

    % inequality constraint: c must be less than or equal to 0
    % this is a conservative limit in cases when prop_pos_SNM < ~1,
        % since it doesn't account for the decrease in density on a
        % given side of 0

    a = mu_SNM ^ 2 * lambda_SNM / shape_SNM + sigma_kernel ^ 2;
    sqrt_a = sqrt(a);
    c = mu_SNM * lambda_SNM + 4 * sqrt_a - L / 2;

    if gradient_specification % check that gradient required

        d_c_d_prop_pos_SNM = 0;
        d_c_d_lambda_SNM = mu_SNM + (2 / sqrt_a) * mu_SNM ^ 2 / shape_SNM;
        d_c_d_mu_SNM = ...
            lambda_SNM + (2 / sqrt_a) * 2 * mu_SNM * lambda_SNM / shape_SNM;
        d_c_d_shape_SNM = ...
            (2 / sqrt_a) * mu_SNM ^ 2 * lambda_SNM * (-1) / shape_SNM ^ 2;

        c_gradient = ...
            [d_c_d_prop_pos_SNM, d_c_d_lambda_SNM, d_c_d_mu_SNM, d_c_d_shape_SNM];
        gradient_key = {'prop_pos_SNM', 'lambda_SNM', 'mu_SNM', 'shape_SNM'};

    end

end