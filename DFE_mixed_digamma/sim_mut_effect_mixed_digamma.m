function new_param_values = sim_mut_effect_mixed_digamma(key_list, value_list)

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
    mu_SNM = parameter_dict('mu_SNM');
    shape_SNM = parameter_dict('shape_SNM');
    prop_pos_SNM = parameter_dict('prop_pos_SNM');
    lambda_SNM = parameter_dict('lambda_SNM');
    prop_neut_SNM = parameter_dict('prop_neut_SNM');
    theta_SNM = mu_SNM/shape_SNM;

    mu_unseq_muts = parameter_dict('mu_unseq_muts');
    shape_unseq_muts = parameter_dict('shape_unseq_muts');
    prop_pos_unseq_muts = parameter_dict('prop_pos_unseq_muts');
    lambda_unseq_muts = parameter_dict('lambda_unseq_muts');
    theta_unseq_muts = mu_unseq_muts/shape_unseq_muts;

    test_strain_number = length(strain_list);
        % # of strains besides reference strain whose likelihood being fitted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop through test strains; within every test strain, loop through fields,
        % and calculate their log likelihood

    % Generate a list containing the total mutation number in each strain,
        % based on the known mutation number and lambda in the 'unsequenced'
        % portion of the genome
    total_mut_numbers_SNM = poissrnd(lambda_SNM, [test_strain_number, 1]);
    total_mut_numbers_unseq_muts = ...
        poissrnd(lambda_unseq_muts, [test_strain_number, 1]);

    % Make a vector of random non-neutral mutation numbers, dependent on p0, in each strain
    non_neutral_muts_SNM = binornd(total_mut_numbers_SNM,(1-prop_neut_SNM));

    % Make a vector of mutation numbers causing a positive effect on the phenotype
    positive_muts_SNM = binornd(non_neutral_muts_SNM, prop_pos_SNM);
    negative_muts_SNM = non_neutral_muts_SNM - positive_muts_SNM;

    positive_muts_unseq_muts = ...
        binornd(total_mut_numbers_unseq_muts, prop_pos_unseq_muts);
    negative_muts_unseq_muts = ...
        total_mut_numbers_unseq_muts - positive_muts_unseq_muts;

    % Make a vector of positive and negative mutational effects
    % Use the fact that the joint distribution of multiple
        % gamma-distributed variables with the same scale parameter is
        % a gamma distribution with that scale parameter and the sum of
        % each distribution's shape parameters
    neg_gamma_shapes_SNM = shape_SNM * negative_muts_SNM;
    pos_gamma_shapes_SNM = shape_SNM * positive_muts_SNM;

    neg_gamma_shapes_unseq_muts = ...
        shape_unseq_muts * negative_muts_unseq_muts;
    pos_gamma_shapes_unseq_muts = ...
        shape_unseq_muts * positive_muts_unseq_muts;

    neg_mut_effects_SNM = -gamrnd(neg_gamma_shapes_SNM, theta_SNM);
    pos_mut_effects_SNM = gamrnd(pos_gamma_shapes_SNM, theta_SNM);
        % for cases with 0 mutations, shape_SNM will be 0, which should
            % produce 0 mut effects
    neg_mut_effects_unseq_muts = ...
        -gamrnd(neg_gamma_shapes_unseq_muts, theta_unseq_muts);
    pos_mut_effects_unseq_muts = ...
        gamrnd(pos_gamma_shapes_unseq_muts, theta_unseq_muts);

    % Combined positive and negative effects in each strain
    total_SNM_mut_effects = neg_mut_effects_SNM + pos_mut_effects_SNM;
    total_unseq_muts_mut_effects = ...
        neg_mut_effects_unseq_muts + pos_mut_effects_unseq_muts;

    % combined mutational effects for each strain are the total SNM
        % effects for that strain plus the total unsequenced mutation
        % effects for that strain
    total_mut_effects = total_SNM_mut_effects + total_unseq_mut_effects;

    strain_me_indices = ...
        parameter_identifier(parameter_list, strain_me_names);
    new_param_values = combined_start_values_array_unscaled;
    new_param_values(strain_me_indices) = total_mut_effects;

end

