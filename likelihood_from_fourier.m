function [real_pdf_vect, pdf_xvals, dist_param_grad_vector_dict] = ...
    likelihood_from_fourier(F_vector, pdf_xvals, N, L, ...
    F_vector_gradient_dict, gradient_specification) 
    
    % Calculate discrete fourier transform of F_vector to get pdf
    pdf_vect = fourier_inverter(F_vector, N, L);
    % Need to remove stray imaginary and negative values from
        % pdf_vect
    real_pdf_vect = real(pdf_vect);
    
    dist_param_grad_vector_dict = containers.Map('KeyType', 'char', ...
        'ValueType', 'any');
    
    if gradient_specification
        
        % Create a dictionary of gradients of LL relative to each
            % distribution parameter
        gradient_dict_keys = keys(F_vector_gradient_dict);
        for current_fitted_param_idx = 1:length(gradient_dict_keys)
            current_fitted_param = ...
                gradient_dict_keys{current_fitted_param_idx};
            d_F_vector_d_current_fitted_param = ...
                F_vector_gradient_dict(current_fitted_param);
            if all(isnan(d_F_vector_d_current_fitted_param))
                dist_param_grad_vector_dict(current_fitted_param) = NaN;
            else
                d_pdf_vect_d_current_fitted_param = ...
                    real(fourier_inverter(d_F_vector_d_current_fitted_param, N, L));
    %           figure; plot(pdf_xvals, d_pdf_vect_d_current_fitted_param, movmean(pdf_xvals,2), movmean(d_pdf_vect_d_current_fitted_param,2)); title(current_fitted_param);
                dist_param_grad_vector_dict(current_fitted_param) = ...
                    d_pdf_vect_d_current_fitted_param;
            end
        end
    end