function [Fn, grad_dict] = fourier_domain_gauss(...
    f, mu, sigma, fitted_parameters, gradient_specification)
    
    % fourier transform of continuous portion of the distribution of
        % 'true' phenotypic effects of each individual mutation,
        % for normally-distributed mutations
    % produces a vector of imaginary values

    variance = sigma ^ 2;

    % In order to speed up runtime for gradients, break up the above
        % equation for Fn into reusable components

    % Fn is a gaussian expressed in terms of a 'mu'
        % and a 'variance' (mean and variance of all mut effects
        % together, respectively)
    Fn = exp(1i * f * mu - variance * f .^ 2 / 2);

    grad_dict = containers.Map('KeyType', 'char', 'ValueType', 'any');

    if gradient_specification % check that gradient required

        % Pre-set all gradients to NaN, and then only calculate the actual
            % values if those gradients are for fitted parameters
        d_Fn_d_mu = NaN;
        d_Fn_d_sigma = NaN;
        d_Fn_d_random_variable = NaN;

        grad_dict = containers.Map('KeyType', 'char', ...
            'ValueType', 'any');

        % Calculate derivative of Fn with respect to each parameter
        if any(strcmp('sigma',fitted_parameters))
            %d_Fn_d_variance = -Fn.*(2*pi*f).^2/2;
            d_Fn_d_sigma = - Fn .* f .^2 * sigma;
        end

        if ~isempty(intersect({'mu', 'random_variable'}, fitted_parameters))

            d_Fn_d_mu = 1i * f .* Fn;

            if any(strcmp('random_variable',fitted_parameters))
                %d_Fn_d_random_variable = Fn.*(2*pi*1i*f);
                d_Fn_d_random_variable = - d_Fn_d_mu;
            end
        end

        grad_dict('mu') = d_Fn_d_mu;
        grad_dict('sigma') = d_Fn_d_sigma;
        grad_dict('random_variable') = d_Fn_d_random_variable;
    end

   

end