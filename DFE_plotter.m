function point_likelihoods = DFE_plotter(gridpower, L, smoothing_kernel_sigma, ...
    current_model, points_to_eval, mle_parameter_names, ...
    mle_parameter_vals)

    gradient_specification = false;
    
    input_value_dict = ...
        containers.Map('KeyType', 'char', 'ValueType', 'any');
    input_value_dict('gridpower') = gridpower;
    input_value_dict('strain_mean_gridpower') = gridpower;
    input_value_dict('smoothing_kernel_sigma') = smoothing_kernel_sigma;
    input_value_dict('mle_parameter_names') = mle_parameter_names;
    input_value_dict('gradient_specification') = gradient_specification;
    input_value_dict('L') = L;

	pre_MLE_output_dict = pre_MLE_me_dist(input_value_dict);
    pre_MLE_output_dict('fitted_parameters') = {};


    N = pre_MLE_output_dict('N');
    me_pdf_frequency_vals = pre_MLE_output_dict('me_pdf_frequency_vals');
    me_pdf_xvals = pre_MLE_output_dict('me_pdf_xvals');
    F_kernel = pre_MLE_output_dict('F_kernel');

	fourier_dfe_fun_name = strcat('fourier_dfe_', current_model);
    fourier_dfe_function = str2func(fourier_dfe_fun_name);
    % Calculate fourier transform of DFE
    [Fa, Fa_gradient_dict] = fourier_dfe_function(...
            mle_parameter_vals, input_value_dict, pre_MLE_output_dict);
    % Calculate Fs, which is Fa smoothed with F_kernel
    [Fs, Fs_gradient_dict] = fourier_space_smoother(Fa, F_kernel, ...
        me_pdf_frequency_vals, Fa_gradient_dict, ...
        gradient_specification);
    % Calculate discrete fourier transform of Fs to get a smoothed
        % version of the pdf of mutational effects per strain (and, if
        % applicable, the gradients of the log of the pdf with respect
        % to each parameter in gradient_param_list)
    [real_me_pdf_smooth, me_pdf_xvals] = ...
        likelihood_from_fourier(Fs, me_pdf_xvals, N, L, ...
            Fs_gradient_dict, gradient_specification);

    point_likelihoods = ...
        interp1(me_pdf_xvals, real_me_pdf_smooth, points_to_eval);
        
%    output_table = table(points_to_eval', point_likelihoods', ...
%    	'VariableNames', {'phenotype', 'density'});
%    writetable(output_table, output_file);






