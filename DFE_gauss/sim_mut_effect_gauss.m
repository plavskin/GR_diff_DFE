function new_param_values = sim_mut_effect_gauss(key_list, value_list)

    % takes input from a data file, and simulates the same number of datapoints

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    input_value_dict = containers.Map(key_list,value_list);
    
    external_counter = str2num(input_value_dict('external_counter'));
    combined_start_values_array_unscaled = input_value_dict('starting_parameter_vals');
    parameter_list = input_value_dict('parameter_list');
    original_phenotype_file = input_value_dict('original_phenotype_file');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_value_dict_for_pre_MLE_function = input_value_dict;
    input_value_dict_for_pre_MLE_function('phenotype_file') = original_phenotype_file;
    pre_MLE_output_dict = pre_MLE_GR_diff(input_value_dict_for_pre_MLE_function);

    strain_list = pre_MLE_output_dict('strain_list');

    strain_me_names = strcat(strain_list,'_me');
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
    % identify parameter values
    parameter_dict = containers.Map(parameter_list,...
        combined_start_values_array_unscaled);
    
    % mutational effect distribution parameters
    mu_unseq_muts = parameter_dict('mu_unseq_muts');
    sigma_unseq_muts = parameter_dict('sigma_unseq_muts');

    test_strain_number = length(strain_list);
        % # of strains besides reference strain whose likelihood being fitted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate mutational effects per strain coming from unsequenced
        % mutations, whose strainwise effects are drawn from a gaussian
    total_mut_effects = ...
        normrnd(mu_unseq_muts, sigma_unseq_muts, [test_strain_number, 1]);

    strain_me_indices = ...
        parameter_identifier(parameter_list, strain_me_names);
    new_param_values = combined_start_values_array_unscaled;
    new_param_values(strain_me_indices) = total_mut_effects;

end

