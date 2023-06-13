function [c, ceq, c_gradient, ceq_gradient, gradient_key] = ...
            nonlin_constraint_gauss(param_vals, input_value_dict, ...
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
    mu_unseq_muts = parameter_dict('mu_unseq_muts');
    sigma_unseq_muts = parameter_dict('sigma_unseq_muts');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % no equality constraint
    ceq = [];
    ceq_gradient=[];

    % inequality constraint: c must be less than or equal to 0

    a = sigma_kernel ^ 2 + sigma_unseq_muts ^ 2;
    sqrt_a = sqrt(a);
    c = mu_unseq_muts +  4 * sqrt_a - L / 2;

    if gradient_specification % check that gradient required

        d_c_d_mu_unseq_muts = 1;
        d_c_d_sigma_unseq_muts = 8 * sigma_unseq_muts;

        c_gradient = [d_c_d_mu_unseq_muts, d_c_d_sigma_unseq_muts];
        gradient_key = {'mu_unseq_muts', 'sigma_unseq_muts'};

    end

end