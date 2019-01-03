function sim_mixef(key_list, value_list)

    % takes input from a data file, and simulates the same number of datapoints

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    input_value_dict = containers.Map(key_list,value_list);
    
    external_counter = str2num(input_value_dict('external_counter'));
    combined_start_values_array_unscaled = input_value_dict('starting_parameter_vals');
    parameter_list = input_value_dict('parameter_list');
    original_data_file = input_value_dict('original_data_file');
    data_file = input_value_dict('data_file');
    pause_at_end = input_value_dict('pause_at_end');

    ranef_names = input_value_dict('random_effect_names');
    intercept_parameter = input_value_dict('intercept_parameter');
    response_name = input_value_dict('response_name');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_value_dict_for_pre_MLE_function = input_value_dict;
    input_value_dict_for_pre_MLE_function('data_file') = original_data_file;
    pre_MLE_output_dict = pre_MLE_mixef(input_value_dict_for_pre_MLE_function);

    ranef_corr_struct = pre_MLE_output_dict('ranef_corr_struct');
    fixef_ID_mat = pre_MLE_output_dict('fixef_ID_mat');
    unique_ranef_categories = pre_MLE_output_dict('unique_ranef_categories');
    fixef_names = pre_MLE_output_dict('unique_fixefs');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use a random seed that's the sum of the current time in seconds and
        % external_counter, so that mutliple jobs starting at the same time have
        % different seeds
    rng('shuffle')
    rng_shuffle = rng;
    random_seed = rng_shuffle.Seed + external_counter;
        % note that if random_seed exceeds 2^32, it maxes out
    if random_seed > 2^32
        random_seed = rng_shuffle.Seed - external_counter;
    end
    rng(random_seed);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get data
    tic;
    data_table = readtable(original_data_file);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parameter_dict = containers.Map(parameter_list,...
        combined_start_values_array_unscaled);
    
    num_ranef = length(ranef_names);

    ranef_guess_list = cell2mat(values(parameter_dict, ranef_names));
    fixef_guess_list = cell2mat(values(parameter_dict, fixef_names));

    if isempty(intercept_parameter)
        slope_bool = true(size(fixef_names));
        fixef_intercept_guess = 0;
    else
        slope_bool = ~ismember(fixef_names,intercept_parameter);
        fixef_intercept_guess = fixef_guess_list(~slope_bool);
    end
   
    fixef_slopes_guess = fixef_guess_list(slope_bool);

    % get a vector of population (expected) means and s.d. for growth rates of each colony
    fixef_guess_relative = zeros(size(slope_bool));
    fixef_guess_relative(slope_bool) = fixef_slopes_guess;
    fixef_guess = fixef_guess_relative+fixef_intercept_guess;

    mean_vector = fixef_ID_mat*fixef_guess';

    sim_ranef_vector = zeros(size(mean_vector));

    for current_ranef_index = 1:num_ranef
        current_ranef_name = ranef_names(current_ranef_index);
        current_sd_guess = ranef_guess_list(current_ranef_index);
        current_ranef_ID_mat = ranef_corr_struct.(current_ranef_name{:});

        current_unique_ranef_categories = unique_ranef_categories(current_ranef_index);
        current_ranef_values = random('normal', 0, current_sd_guess, [1 current_unique_ranef_categories]);

        current_ranef_vector = current_ranef_ID_mat*current_ranef_values';

        sim_ranef_vector = sim_ranef_vector + current_ranef_vector;
    end

    sim_response_vector = mean_vector + sim_ranef_vector;

    data_table{:,response_name} = sim_response_vector;
    writetable(data_table, data_file);

    runtime = toc;

    if pause_at_end & runtime < 120
        pausetime=120-runtime;
        pause(pausetime)
    end

end

