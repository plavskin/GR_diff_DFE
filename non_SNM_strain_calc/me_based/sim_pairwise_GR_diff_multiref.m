function sim_pairwise_GR_diff_multiref(key_list, value_list)

    % takes input from a data file, and simulates the same number of datapoints

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    input_value_dict = containers.Map(key_list,value_list);
    
    external_counter = str2num(input_value_dict('external_counter'));
    combined_start_values_array_unscaled = input_value_dict('starting_parameter_vals');
    parameter_list = input_value_dict('parameter_list');
    original_phenotype_file = input_value_dict('original_phenotype_file');
    phenotype_file = input_value_dict('phenotype_file');
    pause_at_end = input_value_dict('pause_at_end');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_value_dict_for_pre_MLE_function = input_value_dict;
    input_value_dict_for_pre_MLE_function('phenotype_file') = original_phenotype_file;
    pre_MLE_output_dict = pre_MLE_GR_diff_non_SNM_strainwise(input_value_dict_for_pre_MLE_function);

    strain_list = pre_MLE_output_dict('strain_list');
    unique_test_strains = pre_MLE_output_dict('unique_test_strains');
    unique_ref_strains = pre_MLE_output_dict('unique_ref_strains');
    test_strain_ID_mat = pre_MLE_output_dict('test_strain_ID_mat');
    ref_strain_ID_mat = pre_MLE_output_dict('ref_strain_ID_mat');
    GR_diff_list = pre_MLE_output_dict('GR_diff_list');
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
    data_table = readtable(original_phenotype_file);
    GR_diff_list_new = NaN(size(GR_diff_list));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify parameter values
    parameter_dict = containers.Map(parameter_list,...
        combined_start_values_array_unscaled);
    
    petite_colony_sigma = parameter_dict('petite_colony_sigma');
        % s.d. of colony GRs of petite distribution
    nonpetite_colony_sigma = parameter_dict('nonpetite_colony_sigma');
        % s.d. of colony GRs of non-petite reference and test strain distribution
    petite_mean = parameter_dict('petite_mean');
        % mean growth rate of petite colonies, regardless of genotype
    ref_mean = parameter_dict('ref_mean');
        % mean growth rate of non-petite ref colonies

    test_strain_number = length(strain_list);
        % # of strains besides reference strain whose likelihood being fitted

    ref_sigma = nonpetite_colony_sigma;
    test_sigma = nonpetite_colony_sigma;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop through every unique combination of reference and test strains,
        % simulate growth rate differences
    for test_strain_idx = 1:test_strain_number
        current_test_strain = unique_test_strains(test_strain_idx);
        test_pp_name = strcat(current_test_strain, '_pp');
        test_me_name = strcat(current_test_strain, '_me');
        test_petite_prop = parameter_dict(test_pp_name);
        test_mut_effect = parameter_dict(test_me_name);
        test_mean = ref_mean * exp(test_mut_effect);

        for ref_strain_idx = 1:ref_strain_number

            current_ref_strain = unique_ref_strains{ref_strain_idx};
            ref_pp_name = strcat(current_ref_strain, '_pp');
            ref_me_name = strcat(current_ref_strain, '_me');
            ref_petite_prop = parameter_dict(ref_pp_name);
            ref_mut_effect = parameter_dict(ref_me_name);
            current_ref_mean = ref_mean * exp(ref_mut_effect);

            % pick out growth rate differences corresponding to current
                % combination of reference and test strain
            current_data_indices = ...
                test_strain_ID_mat(:, test_strain_idx) .* ...
                ref_strain_ID_mat(:, ref_strain_idx);

            if size(current_data_indices, 1) > 0

                current_GR_diff_list = ...
                    sim_within_pair_different_sigmas(test_mean,...
                        current_ref_mean, test_petite_prop, ref_petite_prop, ...
                        petite_mean, test_sigma, ref_sigma, ...
                        petite_colony_sigma, current_datapoint_number);

                GR_diff_list_new(current_data_indices) = ...
                    current_GR_diff_list;

            end
        end
    end
    
    data_table.GR_diff = GR_diff_list_new;
    writetable(data_table, phenotype_file);

    runtime = toc;

    if pause_at_end & runtime < 120
        pausetime=120-runtime;
        pause(pausetime)
    end

end

