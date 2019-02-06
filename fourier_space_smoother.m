function [Fs, Fs_gradient_dict] = fourier_space_smoother(F_vector, ...
    F_kernel, freq_vals, F_vector_gradient_dict, gradient_specification)
    
    % Calculate Fs, which is F_vector smoothed with F_kernel

    if gradient_specification
        % Create a dictionary of fourier transforms of the gradients of
            % each parameter in gradient_parameters relative to the
            % smoothed log likelihood
        empty_F_kernel_grad_dict = containers.Map();
        [Fs, Fs_gradient_dict] = derivative_multiplier(F_vector, F_kernel, ...
            F_vector_gradient_dict, empty_F_kernel_grad_dict);
        
        % Calculate fourier transform of the gradient of log likelihood of
            % the distribution relative to the random variable
            % (i.e. x-axis of pdf)
        Fs_gradient_dict('random_variable') = - 1i * freq_vals .* Fs;
    else
        Fs = F_vector .* F_kernel;
        Fs_gradient_dict = containers.Map();
    end
end