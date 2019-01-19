function output_dict = pre_MLE_me_dist(input_value_dict)

    % Read phenotype data for petites
    gridpower = input_value_dict('gridpower');
    smoothing_kernel_sigma = input_value_dict('smoothing_kernel_sigma');
    L = input_value_dict('L');

    N = 2^gridpower;
    [me_pdf_frequency_vals, me_pdf_xvals] = fourier_space_creator(N, L);
    % Calculate Fn, the fourier transform of the smoothing kernel
    Fn = fourier_domain_gauss(me_pdf_frequency_vals, 0, ...
        smoothing_kernel_sigma, {}, false);
    
    output_dict = ...
        containers.Map({'N', 'me_pdf_frequency_vals', 'me_pdf_xvals', 'Fn'}, ...
        {N, me_pdf_frequency_vals, me_pdf_xvals, Fn});

end