function output_dict = pre_MLE_me_dist(input_value_dict)

    % Read phenotype data for petites
    gridpower = input_value_dict('gridpower');
    smoothing_kernel_sigma = input_value_dict('smoothing_kernel_sigma');
    L = input_value_dict('L');
    if (isKey(input_value_dict, 'strain_mean_gridpower'))
        strain_mean_gridpower = input_value_dict('strain_mean_gridpower');
    else
        strain_mean_gridpower = smoothing_kernel_sigma;
    end

    N = 2^gridpower;
    [me_pdf_frequency_vals, me_pdf_xvals] = fourier_space_creator(N, L);
    % Calculate Fn, the fourier transform of the smoothing kernel
    F_kernel = fourier_domain_gauss(me_pdf_frequency_vals, 0, ...
        smoothing_kernel_sigma, {}, false);
    
    % min and max values of strain_pdf_xvals must be at or beyond min
        % and max of me_pdf_xvals
    strain_pdf_xvals = linspace(min(me_pdf_xvals), max(me_pdf_xvals), ...
        2^strain_mean_gridpower);

    output_dict = ...
        containers.Map({'N', 'me_pdf_frequency_vals', 'me_pdf_xvals', ...
            'F_kernel', 'strain_pdf_xvals'}, ...
        {N, me_pdf_frequency_vals, me_pdf_xvals, F_kernel, strain_pdf_xvals});

end