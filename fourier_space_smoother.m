function [Fs, Fs_gradient_dict] = fourier_space_smoother(F_vector, ...
    F_kernel, freq_vals, F_vector_gradient_dict, gradient_specification)
    
    % Calculate Fs, which is F_vector smoothed with F_kernel
    Fs = F_vector .* F_kernel;
    
    if gradient_specification
        
        % Create a dictionary of fourier transforms of the gradients of
            % each parameter in gradient_parameters relative to the
            % smoothed log likelihood
        % Use keys currently in F_vector_gradient_dict and set every value
            % to NaN, in case gradient_parameters are for some reason just
            % a subset of F_vector_gradient_dict keys
        gradient_dict_keys = keys(F_vector_gradient_dict);
        Fs_gradient_dict = containers.Map();
        
        for current_fitted_param_idx = 1:length(gradient_dict_keys)
            current_fitted_param = ...
                gradient_dict_keys{current_fitted_param_idx};
            d_F_vector_d_current_fitted_param = ...
                F_vector_gradient_dict(current_fitted_param);
            if all(isnan(d_F_vector_d_current_fitted_param))
                Fs_gradient_dict(current_fitted_param) = NaN;
            else
                % Adjust gradients to account for smoothing
                Fs_gradient_dict(current_fitted_param) = ...
                    d_F_vector_d_current_fitted_param .* F_kernel;
            end
        end
        
        % Calculate fourier transform of the gradient of log likelihood of
            % the distribution relative to the random variable
            % (i.e. x-axis of pdf)
        Fs_gradient_dict('random_variable') = - 1i * freq_vals .* Fs;
    end
end