function [neg_LL,gradient_vector_partial] = LL_mixef_calc_v5b(fitted_parameters,...
    fixed_parameter_indices,fixed_parameter_values,response_vector,ranef_corr_struct,...
    fixef_ID_mat,nonref_indices,unique_ranef_categories,calculate_gradient,...
    order_vector,block_start_positions)

    % first num_ranef elements of fitted_parameters are current guesses for
        % ranef std
    % calculates LL as a function of Lambda, a matrix of sds multiplied
        % by ranef identities for all ranefs
    % allows an intercept and slopees for both the fixef effects and the
        % inter-individual s.d. relative to a 'reference' strain
    % allows separate processing of matrices along block-diagonal

%    disp('starting lme')

    param_vals = NaN(size(fixed_parameter_indices));
    param_vals(fixed_parameter_indices) = fixed_parameter_values(~isnan(fixed_parameter_values));
    param_vals(~fixed_parameter_indices) = fitted_parameters;

    ranef_names = fieldnames(ranef_corr_struct);
    num_ranef = length(ranef_names);
    num_responses = length(response_vector);
    num_fixef = size(fixef_ID_mat,2);
    
    ranef_guess_list = param_vals(1:num_ranef);
    fixef_guess_list = param_vals((num_ranef+1):(num_ranef+num_fixef));
    sigma_guess_list = param_vals((num_ranef+num_fixef+1):(num_ranef+2*num_fixef));

    fixef_intercept_guess = fixef_guess_list(~nonref_indices);
    fixef_slopes_guess = fixef_guess_list(nonref_indices);
    sigma_intercept_guess = sigma_guess_list(~nonref_indices);
    sigma_slopes_guess = sigma_guess_list(nonref_indices);

    % get a vector of population (expected) means and s.d. for growth rates of each colony
    fixef_guess_relative = zeros(size(nonref_indices));
    fixef_guess_relative(nonref_indices) = fixef_slopes_guess;
    fixef_guess = fixef_guess_relative+fixef_intercept_guess;

    sigma_guess_relative = zeros(size(nonref_indices));
    sigma_guess_relative(nonref_indices) = sigma_slopes_guess;
    sigma_guess = sigma_guess_relative+sigma_intercept_guess;

    mean_vector = fixef_ID_mat*fixef_guess';
    sigma_vector = fixef_ID_mat*sigma_guess';

    %tic;
    
    % loop through ranef cov structure matrices, multiply each by expected
        % var value, add structure together to form cov matrix

%    disp('starting ranef calculations')

    unique_ranef_categories_with_sigma = [unique_ranef_categories,num_responses];
    total_unique_ranef_states = sum(unique_ranef_categories_with_sigma);
    total_nonzero_elements = (num_ranef+1)*num_responses;
    ranef_column_start_positions = 1+[0,cumsum(unique_ranef_categories_with_sigma)];
    ranef_column_end_positions = cumsum(unique_ranef_categories_with_sigma);

    Lambda = spalloc(num_responses,total_unique_ranef_states,total_nonzero_elements);
    
%    disp('lambda allocated')

    % include random effect matrices in combined Lambda matrix
    for current_ranef_index = 1:num_ranef
        current_ranef_name = ranef_names(current_ranef_index);
        current_sd_guess = ranef_guess_list(current_ranef_index);
        current_ranef_cov_struct = ranef_corr_struct.(current_ranef_name{:})*current_sd_guess;

        current_first_column = ranef_column_start_positions(current_ranef_index);
        current_last_column = ranef_column_end_positions(current_ranef_index);

        Lambda(:,current_first_column:current_last_column) = current_ranef_cov_struct;

    end

    clear current_ranef_cov_struct
%    clear ranef_corr_struct

    % include residual (strain-specific sigma) effects in Lambda matrix
    current_first_column = ranef_column_start_positions(num_ranef+1);
    current_last_column = ranef_column_end_positions(num_ranef+1);

    Lambda(:,current_first_column:current_last_column) = ...
        spdiags(sigma_vector,0,num_responses,num_responses);

    total_cov = Lambda*Lambda';

    clear Lambda

%    disp('total_cov calculated')

    [order_vector,block_structure] = Block_Finder(total_cov,order_vector,block_start_positions);
    
    neg_LL = 0;
    block_number = size(block_structure,2);

    if calculate_gradient
        gradient_vector = zeros(size(param_vals));
    end

    for block_counter=1:block_number
%        disp(block_counter)
        current_cov = block_structure(block_counter).block;
        current_start = block_structure(block_counter).start;
        current_end = block_structure(block_counter).end;

        current_indices = order_vector(current_start:current_end);
        current_response_vector = response_vector(current_indices);
        current_mean_vector = mean_vector(current_indices);
        current_num_responses = length(current_indices);

%        disp(current_num_responses)

        current_diff_vector = current_response_vector - current_mean_vector;
        
%        disp('calculating log_det')

        current_log_det_cov = full(log_determinant_calculator(current_cov));

%        disp('log_det calculated')

        current_diff_vector_over_cov = current_diff_vector'/current_cov;
        current_LL = -(current_num_responses*log(2*pi) + current_log_det_cov)/2+...
            ((-1/2)*(current_diff_vector_over_cov)*current_diff_vector);

%        disp('LL calculated')
        neg_LL = neg_LL - current_LL;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if calculate_gradient

            current_fixef_ID_mat = fixef_ID_mat(current_indices,:);

            % calculate gradient relative to fixef
            d_LL_d_fixef_guess = current_diff_vector_over_cov*current_fixef_ID_mat;
            d_LL_d_fixef_intercept = sum(d_LL_d_fixef_guess);
            d_LL_d_fixef_slopes = d_LL_d_fixef_guess(nonref_indices);
            d_LL_d_fixef = NaN(size(nonref_indices));
            d_LL_d_fixef(nonref_indices) = d_LL_d_fixef_slopes;
            d_LL_d_fixef(~nonref_indices) = d_LL_d_fixef_intercept;

            % calculate gradient relative to ranef (excluding residuals)
            d_LL_d_ranef = NaN(1,num_ranef);

            for current_ranef_index = 1:num_ranef

                current_first_column = ranef_column_start_positions(current_ranef_index);
                current_last_column = ranef_column_end_positions(current_ranef_index);

                current_ranef_name = ranef_names(current_ranef_index);
                total_ranef_cov_positions = ranef_corr_struct.(current_ranef_name{:});
                current_ranef_cov_positions = total_ranef_cov_positions(current_indices,:);
                current_ranef_val = ranef_guess_list(current_ranef_index);

                current_cov_mat_positions = current_ranef_cov_positions*current_ranef_cov_positions';

            %    d_LL_d_current_ranef = -current_ranef_val*...
            %        (trace(current_ranef_cov_positions*current_ranef_cov_positions'/total_cov)-...
            %        diff_vector'/total_cov*current_ranef_cov_positions*current_ranef_cov_positions'/...
            %            total_cov*diff_vector);
                d_LL_d_current_ranef = d_LL_d_cov_component(current_ranef_val,...
                    current_cov_mat_positions,current_cov,current_diff_vector,...
                    current_diff_vector_over_cov);

                d_LL_d_ranef(current_ranef_index) = d_LL_d_current_ranef;

            end

            % calculate gradient relative to inter-individual sd (residuals)
            d_LL_d_sigma_guess = NaN(1,length(sigma_guess));

            for current_sigma_index = 1:length(sigma_guess)

                current_sigma = sigma_guess(current_sigma_index);
                current_diagonal = current_fixef_ID_mat(:,current_sigma_index);
                current_cov_mat_positions = spdiags(current_diagonal,0,...
                    current_num_responses,current_num_responses);

                d_LL_d_current_sigma = d_LL_d_cov_component(current_sigma,...
                    current_cov_mat_positions,current_cov,current_diff_vector,...
                    current_diff_vector_over_cov);

                d_LL_d_sigma_guess(current_sigma_index) = d_LL_d_current_sigma;

            end

            d_LL_d_sigma_intercept = sum(d_LL_d_sigma_guess);
            d_LL_d_sigma_slopes = d_LL_d_sigma_guess(nonref_indices);
            d_LL_d_sigma = NaN(size(nonref_indices));
            d_LL_d_sigma(nonref_indices) = d_LL_d_sigma_slopes;
            d_LL_d_sigma(~nonref_indices) = d_LL_d_sigma_intercept;

            % combined gradients into single vector
            gradient_vector = gradient_vector-[d_LL_d_ranef,d_LL_d_fixef,d_LL_d_sigma];
        else
            gradient_vector = NaN(size(param_vals));
        end


    end

    clear current_cov

    gradient_vector_partial = gradient_vector(~fixed_parameter_indices);

end