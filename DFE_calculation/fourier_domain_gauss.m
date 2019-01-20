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
    Fn = exp(j * f * mu - variance * f .^ 2 / 2);

    if gradient_specification % check that gradient required

        % Pre-set all gradients to NaN, and then only calculate the actual
            % values if those gradients are for fitted parameters
        d_Fn_d_mu = NaN;
        d_Fn_d_sigma = NaN;
        d_Fn_d_x = NaN;

        % Calculate derivative of Fn with respect to each parameter
        if any(strcmp('sigma',fitted_parameters))
            %d_Fn_d_variance = -Fn.*(2*pi*f).^2/2;
            d_Fn_d_sigma = - Fn .* f .^2 * sigma;
        end

        if ~isempty(intersect({'mu', 'x'}, fitted_parameters))

            d_Fn_d_mu = j * f .* Fn;

            if any(strcmp('x',fitted_parameters))
                %d_Fn_d_x = Fn.*(2*pi*1i*f);
                d_Fn_d_x = - d_Fn_d_mu;
            end
        end

        unscaled_gradient_vector = {d_Fn_d_mu, ...
            d_Fn_d_sigma, d_Fn_d_x};
        grad_parameter_names = {'mu', 'sigma', 'x'};
        grad_dict = containers.Map(grad_parameter_names, unscaled_gradient_vector);

    end

   

end