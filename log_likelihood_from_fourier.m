function [real_pdf_vect, pdf_xvals, dist_param_grad_vector_dict] = ...
    log_likelihood_from_fourier(F_vector, pdf_xvals, N, L, ...
    F_vector_gradient_dict, gradient_specification) 
    
    % Calculate discrete fourier transform of F_vector to get pdf
    pdf_vect = fourier_inverter(F_vector, N, L);
    % Need to remove stray imaginary and negative values from
        % pdf_vect
    real_pdf_vect = real(pdf_vect);
    
    if gradient_specification
        % If calculating gradients of LL, need to divide gradient relative
            % to parameter by real_pdf_vect, so remove any values where
            % real_pdf_vect is 0
        non_zero_positions = (real_pdf_vect ~= 0);
        need_to_remove_zero_vals = ...
            length(real_pdf_vect) > sum(non_zero_positions);
        if need_to_remove_zero_vals
            real_pdf_vect = real_pdf_vect(non_zero_positions);
            pdf_xvals = pdf_xvals(non_zero_positions);
        end
        
        % Create a dictionary of gradients of LL relative to each
            % distribution parameter
        gradient_dict_keys = keys(F_vector_gradient_dict);
        dist_param_grad_vector_dict = containers.Map();
        for current_fitted_param_idx = 1:length(gradient_dict_keys)
            current_fitted_param = ...
                gradient_dict_keys{current_fitted_param_idx};
            d_F_vector_d_current_fitted_param = ...
                F_vector_gradient_dict(current_fitted_param);
            if all(isnan(d_F_vector_d_current_fitted_param))
                dist_param_grad_vector_dict(current_fitted_param) = NaN;
            else
                d_pdf_vect_d_current_fitted_param = ...
                    fourier_inverter(d_F_vector_d_current_fitted_param, N, L);
                if need_to_remove_zero_vals
                    % remove values where real_pdf_vect is 0
                    d_pdf_vect_d_current_fitted_param = ...
                        d_pdf_vect_d_current_fitted_param(non_zero_positions);
                end
                d_LL_d_current_fitted_param = ...
                    real(d_pdf_vect_d_current_fitted_param) ./ ...
                        real_pdf_vect;
    %           figure; plot(pdf_xvals, d_LL_d_current_fitted_param, movmean(pdf_xvals,2), movmean(d_LL_d_current_fitted_param,2)); title(current_fitted_param);
                dist_param_grad_vector_dict(current_fitted_param) = ...
                    d_LL_d_current_fitted_param;
            end
        end
    end