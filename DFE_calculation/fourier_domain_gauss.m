function [continuous_Fz, grad_dict] = fourier_domain_gauss(...
    f, mu, sigma, fitted_parameters, gradient_specification)
    
    % fourier transform of continuous portion of the distribution of
        % 'true' phenotypic effects of each individual mutation,
        % for normally-distributed mutations
    % produces a vector of imaginary values

    variance = sigma^2;

    % In order to speed up runtime for gradients, break up the above
        % equation for continuous_Fz into reusable components

    % continuous_Fz is a gaussian expressed in terms of a 'mu'
        % and a 'variance' (mean and variance of all mut effects
        % together, respectively)
    continuous_Fz = exp(-1i*2*pi*f*mu-variance*(2*pi*f).^2/2);

    if gradient_specification % check that gradient required

        % Pre-set all gradients to NaN, and then only calculate the actual
            % values if those gradients are for fitted parameters
        d_continuous_Fz_d_mu = NaN;
        d_continuous_Fz_d_sigma = NaN;
        d_continuous_Fz_d_x = NaN;

        % Calculate derivative of continuous_Fz with respect to each parameter
        if any(strcmp('sigma',fitted_parameters))
            %d_continuous_Fz_d_variance = -continuous_Fz.*(2*pi*f).^2/2;
            d_continuous_Fz_d_sigma = -continuous_Fz.*(2*pi*f).^2*sigma;
        end

        if ~isempty(intersect({'mu', 'x'}, fitted_parameters))

            d_continuous_Fz_d_mu = -1i*2*pi*continuous_Fz.*f;

            if any(strcmp('x',fitted_parameters))
                %d_continuous_Fz_d_x = continuous_Fz.*(2*pi*1i*f);
                d_continuous_Fz_d_x = -d_continuous_Fz_d_mu;
            end
        end

        unscaled_gradient_vector = {d_continuous_Fz_d_mu, ...
            d_continuous_Fz_d_sigma, d_continuous_Fz_d_x};
        grad_parameter_names = {'mu', 'sigma', 'x'};
        grad_dict = containers.Map(grad_parameter_names, unscaled_gradient_vector);

    end

   

end