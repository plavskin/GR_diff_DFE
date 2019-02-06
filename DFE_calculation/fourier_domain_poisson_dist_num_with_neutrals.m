function [Fa, grad_dict] = fourier_domain_poisson_dist_num_with_neutrals(...
    Fz, f, lambda, prop_neut, initial_grad_dict, fitted_parameters, gradient_specification)

    % takes an input distribution in fourier space and its gradients and
        % returns the fourier transforms of the distribution and
        % gradients of a poisson random number of draws with mean lambda
        % and prop_neut mutations that don't have an effect is drawn
        % from the initial distribution

    % Fa is the fourier transform of the pdf of the combined
        % mutational effects of all muts in a strain (across all
        % possible mutation numbers)
    effective_lambda = (1 - prop_neut) * lambda;
    Fa = exp(-effective_lambda * (1 - Fz));
    % The discrete inverse fourier transform of Fa results in a point
        % mass at 0, since Fa approaches exp(-lambda), not 0, as f
        % approaches infinity
    % This is problematic because it means the fourier transform of Fa
        % isn't smooth; therefore, need to convolve inverse fourier
        % transform of Fa with a gaussian kernel with a small variance
        % to smooth it out (do this outside of this function)

    grad_dict = containers.Map('KeyType', 'char', 'ValueType', 'any');

    if gradient_specification

        % Calculate derivative of Fa with respect to each parameter
            
        d_Fa_d_lambda = NaN;
        d_Fa_d_prop_neut = NaN;
        if ~isempty(intersect({'lambda', 'prop_neut'}, fitted_parameters))

            d_Fa_d_effective_lambda = Fa .* (Fz - 1);

            if any(strcmp('lambda',fitted_parameters))
                d_effective_lambda_d_lambda = 1 - prop_neut;
                d_Fa_d_lambda = d_Fa_d_effective_lambda * ...
                    d_effective_lambda_d_lambda;
            end

            if any(strcmp('prop_neut',fitted_parameters))
                d_effective_lambda_d_prop_neut = -lambda;
                d_Fa_d_prop_neut = d_Fa_d_effective_lambda * ...
                    d_effective_lambda_d_prop_neut;
            end

        end
        
        grad_dict('lambda') = d_Fa_d_lambda;
        grad_dict('prop_neut') = d_Fa_d_prop_neut;

        gradient_names = keys(initial_grad_dict);
        for current_grad_idx=1:length(gradient_names)
            current_grad_name = gradient_names{current_grad_idx};
            % if current_grad_name in fitted_parameters, use gradient of
                % Fz relative to that gradient to calculate gradient Fa
                % relative to that gradient, and save it in
                % grad_dict
            if any(strcmp(current_grad_name, fitted_parameters))
                if strcmp(current_grad_name, 'random_variable')
                    grad_dict('random_variable') = -1i * Fa .* f;
                else
                    d_Fz_d_current_grad = ...
                        initial_grad_dict(current_grad_name);
                    grad_dict(current_grad_name) = ...
                        effective_lambda * Fa .* d_Fz_d_current_grad;
                end
            else
                grad_dict(current_grad_name) = NaN;
            end
        end
    end

end