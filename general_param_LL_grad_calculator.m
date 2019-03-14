function d_LL_d_current_param = ...
    general_param_LL_grad_calculator(LL_observed_diffs_grad_dict, ...
    current_param, likelihood_to_integrate_over_uncorrected, ...
    me_pdf_xvals, likelihood_uncorrected)

    % Calculate gradient of log likelihood wrt 'general' parameters
        % (non-DFE)
    % d ( L (observations)) / d( theta{general} ) ) = 
    %    integral across region of interest(
    %        P(mutational effect = x | theta{DME}) *
    %        product across strain_GR_diff_list(
    %            P(GR_diff | mutational effect = x, theta{general})
    %            ) *
    %        sum across strain_GR_diff_list(
    %            d ( ln ( P(GR_diff | mutational effect = x, theta{general}) ) ) /
    %                d ( theta{general} )
    %            ) * 
    %        dx
    %        )
    
    d_LL_observed_diffs_d_current_param = ...
        LL_observed_diffs_grad_dict(current_param);
    d_L_to_integrate_over_uncorrected_d_current_param = ...
        likelihood_to_integrate_over_uncorrected .* ...
        d_LL_observed_diffs_d_current_param;
    %figure; plot(me_pdf_xvals, d_L_to_integrate_over_uncorrected_d_current_param)
    d_likelihood_uncorrected_d_current_param = ...
        trapz(me_pdf_xvals, ...
            d_L_to_integrate_over_uncorrected_d_current_param);
    d_LL_d_current_param = ...
        d_likelihood_uncorrected_d_current_param / ...
        likelihood_uncorrected;

end