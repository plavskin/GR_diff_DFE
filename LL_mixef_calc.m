function [LL, unscaled_gradient_vector, grad_parameter_names] = LL_mixef_calc(param_vals,...
    input_value_dict, pre_MLE_output_dict)

    % first num_ranef elements of fitted_parameters are current guesses for
        % ranef std
    % calculates LL as a function of Lambda, a matrix of sds multiplied
        % by ranef identities for all ranefs
    % allows an intercept and slopees for both the fixef effects and the
        % inter-individual s.d. relative to a 'reference' strain
    % allows separate processing of matrices along block-diagonal

%    disp('starting lme')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ranef_corr_struct = pre_MLE_output_dict('ranef_corr_struct');
    response_vector = pre_MLE_output_dict('response_vector');
    fixef_ID_mat = pre_MLE_output_dict('fixef_ID_mat');
    unique_ranef_categories = pre_MLE_output_dict('unique_ranef_categories');
    order_vector = pre_MLE_output_dict('order_vector');
    block_start_positions = pre_MLE_output_dict('block_start_positions');
    fixef_names = pre_MLE_output_dict('unique_fixefs');

    calculate_gradient = input_value_dict('gradient_specification');
    parameter_list = input_value_dict('parameter_list');
    ranef_names = input_value_dict('random_effect_names');
    intercept_parameter = input_value_dict('intercept_parameter');
    mle_parameter_names = input_value_dict('mle_parameter_names');

    parameter_dict = containers.Map(mle_parameter_names, param_vals);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    num_ranef = length(ranef_names);
    num_responses = length(response_vector);

    ranef_guess_list = cell2mat(values(parameter_dict, ranef_names));
    fixef_guess_list = cell2mat(values(parameter_dict, fixef_names));
    
    if isempty(intercept_parameter)
        slope_bool = true(size(fixef_names));
        fixef_intercept_guess = 0;
    else
        slope_bool = ~ismember(fixef_names,intercept_parameter);
        fixef_intercept_guess = fixef_guess_list(~slope_bool);
    end
   
%    sigma_guess_list = param_vals((num_ranef+num_fixef+1):(num_ranef+2*num_fixef));

    fixef_slopes_guess = fixef_guess_list(slope_bool);
%    sigma_intercept_guess = sigma_guess_list(~slope_bool);
%    sigma_slopes_guess = sigma_guess_list(slope_bool);

    % get a vector of population (expected) means and s.d. for growth rates of each colony
    fixef_guess_relative = zeros(size(slope_bool));
    fixef_guess_relative(slope_bool) = fixef_slopes_guess;
    fixef_guess = fixef_guess_relative+fixef_intercept_guess;

%    sigma_guess_relative = zeros(size(slope_bool));
%    sigma_guess_relative(slope_bool) = sigma_slopes_guess;
%    sigma_guess = sigma_guess_relative+sigma_intercept_guess;

    mean_vector = fixef_ID_mat*fixef_guess';
%    sigma_vector = fixef_ID_mat*sigma_guess';

    %tic;
    
    % loop through ranef cov structure matrices, multiply each by expected
        % var value, add structure together to form cov matrix

%    disp('starting ranef calculations')

%    unique_ranef_categories_with_sigma = [unique_ranef_categories,num_responses];
%    total_unique_ranef_states = sum(unique_ranef_categories_with_sigma);
%    total_nonzero_elements = (num_ranef)*num_responses;
%    ranef_column_start_positions = 1+[0,cumsum(unique_ranef_categories_with_sigma)];
%    ranef_column_end_positions = cumsum(unique_ranef_categories_with_sigma);
    total_unique_ranef_states = sum(unique_ranef_categories);
    total_nonzero_elements = (num_ranef)*num_responses;
    ranef_column_start_positions = 1+[0,cumsum(unique_ranef_categories)];
    ranef_column_end_positions = cumsum(unique_ranef_categories);


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

%   disp('Lambda calculated')

    clear current_ranef_cov_struct
%    clear ranef_corr_struct

%    % include residual (strain-specific sigma) effects in Lambda matrix
%    current_first_column = ranef_column_start_positions(num_ranef+1);
%    current_last_column = ranef_column_end_positions(num_ranef+1);

%    Lambda(:,current_first_column:current_last_column) = ...
%        spdiags(sigma_vector,0,num_responses,num_responses);

    


%    disp('total_cov calculated')

    if (isnan(any(order_vector)) | isnan(any(block_start_positions)))
        total_cov = Lambda*Lambda';
        [order_vector,block_structure] = ...
            Block_Finder(total_cov,order_vector,block_start_positions);
        block_start_positions = [block_structure(:).start];
        block_end_positions = [block_structure(:).end];
        calculate_each_cov = false;
    else
        block_end_positions = ...
            [(block_start_positions(2:end)-1), num_responses];
        calculate_each_cov = true;
    end

    LL = 0;
    block_number = length(block_start_positions);

    if calculate_gradient
        unscaled_gradient_vector = zeros([1 length(ranef_names)+length(fixef_names)]);
    end

    Lambda_ordered = Lambda(order_vector, :);
    
    grad_parameter_names = [ranef_names, fixef_names];

    for block_counter=1:block_number
        %disp(block_counter)
        current_start = block_start_positions(block_counter);
        current_end = block_end_positions(block_counter);

        if calculate_each_cov
            current_Lambda = Lambda_ordered(current_start:current_end, :);
            current_cov = full(current_Lambda * current_Lambda');
        else
            current_cov = block_structure(block_counter).block;
        end
        current_indices = order_vector(current_start:current_end);
        current_response_vector = response_vector(current_indices);
        current_mean_vector = mean_vector(current_indices);
        current_num_responses = length(current_indices);

%        disp(current_num_responses)

        current_diff_vector = current_response_vector - current_mean_vector;
        
%        disp('calculating log_det')

        current_log_det_cov = log_determinant_calculator(current_cov);

%        disp('log_det calculated')

        current_diff_vector_over_cov = current_diff_vector'/current_cov;
        current_LL = -(current_num_responses*log(2*pi) + current_log_det_cov)/2+...
            ((-1/2)*(current_diff_vector_over_cov)*current_diff_vector);

%        disp('LL calculated')
        LL = LL + current_LL;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if calculate_gradient

            current_fixef_ID_mat = fixef_ID_mat(current_indices,:);

            % calculate gradient relative to fixef
            d_LL_d_fixef_guess = current_diff_vector_over_cov*current_fixef_ID_mat;
            d_LL_d_fixef_intercept = sum(d_LL_d_fixef_guess);
            d_LL_d_fixef_slopes = d_LL_d_fixef_guess(slope_bool);
            d_LL_d_fixef = NaN(size(slope_bool));
            d_LL_d_fixef(slope_bool) = d_LL_d_fixef_slopes;
            d_LL_d_fixef(~slope_bool) = d_LL_d_fixef_intercept;

            % calculate gradient relative to ranef (excluding residuals)
            d_LL_d_ranef = NaN(1,num_ranef);
            
            for current_ranef_index = 1:num_ranef

%                current_first_column = ranef_column_start_positions(current_ranef_index);
%                current_last_column = ranef_column_end_positions(current_ranef_index);

                current_ranef_name = ranef_names(current_ranef_index)
                total_ranef_cov_positions = ranef_corr_struct.(current_ranef_name{:});
                current_ranef_cov_positions = full(total_ranef_cov_positions(current_indices,:));
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

%            % calculate gradient relative to inter-individual sd (residuals)
%            d_LL_d_sigma_guess = NaN(1,length(sigma_guess));
%
%            for current_sigma_index = 1:length(sigma_guess)
%
%                current_sigma = sigma_guess(current_sigma_index);
%                current_diagonal = current_fixef_ID_mat(:,current_sigma_index);
%                current_cov_mat_positions = spdiags(current_diagonal,0,...
%                    current_num_responses,current_num_responses);
%
%                d_LL_d_current_sigma = d_LL_d_cov_component(current_sigma,...
%                    current_cov_mat_positions,current_cov,current_diff_vector,...
%                    current_diff_vector_over_cov);
%
%                d_LL_d_sigma_guess(current_sigma_index) = d_LL_d_current_sigma;
%
%            end
%
%            d_LL_d_sigma_intercept = sum(d_LL_d_sigma_guess);
%            d_LL_d_sigma_slopes = d_LL_d_sigma_guess(slope_bool);
%            d_LL_d_sigma = NaN(size(slope_bool));
%            d_LL_d_sigma(slope_bool) = d_LL_d_sigma_slopes;
%            d_LL_d_sigma(~slope_bool) = d_LL_d_sigma_intercept;

            % combined gradients into single vector
%            gradient_vector = gradient_vector-[d_LL_d_ranef,d_LL_d_fixef,d_LL_d_sigma];
            unscaled_gradient_vector = unscaled_gradient_vector+[d_LL_d_ranef,d_LL_d_fixef];
        else
            unscaled_gradient_vector = NaN(size(grad_parameter_names));
        end
    end
end