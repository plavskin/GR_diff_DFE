function [Fa, grad_dict] = fourier_domain_poisson_dist_num(...
    Fz, f, lambda, initial_grad_dict, fitted_parameters, gradient_specification)

    % takes an input distribution in fourier space and its gradients and
        % returns the fourier transforms of the distribution and
        % gradients of a poisson random number of draws with mean lambda
        % drawn from the initial distribution

    % Fa is the fourier transform of the pdf of the combined
        % mutational effects of all muts in a strain (across all
        % possible mutation numbers)
    Fa = exp(-lambda * (1 - Fz));
    % The discrete inverse fourier transform of Fa results in a point
        % mass at 0, since Fa approaches exp(-lambda), not 0, as f
        % approaches infinity
    % This is problematic because it means the fourier transform of Fa
        % isn't smooth; therefore, need to convolve inverse fourier
        % transform of Fa with a gaussian kernel with a small variance
        % to smooth it out (do this outside of this function)

    grad_dict = containers.Map('KeyType', 'char', ...
        'ValueType', 'any');
    
    if gradient_specification

        % Calculate derivative of Fa with respect to each parameter

        if any(strcmp('lambda',fitted_parameters))
            d_Fa_d_lambda = Fa .* (Fz - 1);
        else
            d_Fa_d_lambda = NaN;
        end

        grad_dict('lambda') = d_Fa_d_lambda;

        gradient_names = keys(initial_grad_dict);
        for current_grad_idx=1:length(gradient_names)
            current_grad_name = gradient_names{current_grad_idx};
            % if current_grad_name in fitted_parameters, use gradient of
                % Fz relative to that gradient to calculate gradient Fa
                % relative to that gradient, and save it in grad_dict
            if any(strcmp(current_grad_name, fitted_parameters))
                if strcmp(current_grad_name, 'random_variable')
                    grad_dict('random_variable') = -1i * Fa .* f;
                else
                    d_Fz_d_current_grad = ...
                        initial_grad_dict(current_grad_name);
                    grad_dict(current_grad_name) = ...
                        lambda * Fa .* d_Fz_d_current_grad;
                end
            else
                grad_dict(current_grad_name) = NaN;
            end
        end
    end

end