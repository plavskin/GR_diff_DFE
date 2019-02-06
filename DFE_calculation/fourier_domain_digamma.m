function [Fz, grad_dict]=fourier_domain_digamma(...
    f, mu, shape, prop_pos, fitted_parameters, gradient_specification)

    %Fz = prop_pos*(1./(1+1i*2*pi*f*theta_pos)).^shape+(1-prop_pos)*(1./(1-1i*2*pi*f*theta_neg)).^shape;
    
    % fourier transform of continuous portion of the distribution of
        % 'true' phenotypic effects of a single individual mutation,
        % drawn from a reflected gamma distribution with unequal sides
    % produces a vector of imaginary values

    % In order to speed up runtime for gradients, break up the above
        % equation for Fz into reusable components

    C = 1i * f;
    A = 1 - C .* mu/shape;
    B = 1 + C .* mu/shape;
    A_shape = A.^(-shape);
    B_shape = B.^(-shape);

    Fz = prop_pos * A_shape + (1 - prop_pos) * B_shape;

    grad_dict = containers.Map('KeyType', 'char', ...
            'ValueType', 'any');
    
    if gradient_specification

        % Calculate derivative of Fz with respect to each
            % parameter

        d_Fz_d_prop_pos = NaN;
        d_Fz_d_mu = NaN;
        d_Fz_d_shape = NaN;
        d_Fz_d_random_variable = NaN;

        if any(strcmp('random_variable',fitted_parameters))
            d_Fz_d_random_variable = -1i * Fz .* f;
        end

        if any(strcmp('prop_pos',fitted_parameters))
            d_Fz_d_prop_pos = A_shape - B_shape;
        end

        if ~isempty(intersect({'mu', 'shape'}, fitted_parameters))
            
            d_Fz_d_mu = prop_pos * C .* A .^ (-(shape + 1)) - ...
                (1 - prop_pos) * C .* B .^ (-(shape + 1));

            if any(strcmp('shape',fitted_parameters))
                d_Fz_d_shape = (-mu / shape) * d_Fz_d_mu + prop_pos * ...
                    (-log(A)) .* A_shape + (1 - prop_pos) * (-log(B)) .* B_shape;
            end
        end

        grad_dict('mu') = d_Fz_d_mu;
        grad_dict('shape') = d_Fz_d_shape;
        grad_dict('prop_pos') = d_Fz_d_prop_pos;
        grad_dict('random_variable') = d_Fz_d_random_variable;

    end


end