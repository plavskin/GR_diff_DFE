function GR_diff_list = sim_within_pair_different_sigmas(test_mean, ...
	ref_mean, test_petite_prop, ref_petite_prop, petite_mean, test_sigma, ...
	ref_sigma, petite_sigma, datapoint_number)
	
	% simulates datapoint_number datapoints of differences in growth
		% rate between a test and reference strain, where growth rates
		% in both strains come from two distinct distributions: the
		% reference/test strain distribution, and the petite strain
		% distribution

    ref_GR_list = NaN([1 datapoint_number]);
    test_GR_list = NaN([1 datapoint_number]);

	ref_petite_colony_indices = ...
        logical(binornd(ones(size(ref_GR_list)), ref_petite_prop));
	test_petite_colony_indices = ...
        logical(binornd(ones(size(test_GR_list)), test_petite_prop));
	
    ref_petite_colony_num = sum(ref_petite_colony_indices);
    test_petite_colony_num = sum(test_petite_colony_indices);
    
    ref_petite_GR_list = ...
        normrnd(petite_mean, petite_sigma, [1 ref_petite_colony_num]);
    test_petite_GR_list = ...
        normrnd(petite_mean, petite_sigma, [1 test_petite_colony_num]);
    
    ref_nonpetite_GR_list = ...
        normrnd(ref_mean, ref_sigma, ...
        	[1 (datapoint_number-ref_petite_colony_num)]);
    test_nonpetite_GR_list = ...
        normrnd(test_mean, test_sigma, ...
        	[1 (datapoint_number-test_petite_colony_num)]);
    
    ref_GR_list(ref_petite_colony_indices) = ref_petite_GR_list;
    ref_GR_list(~ref_petite_colony_indices) = ref_nonpetite_GR_list;
    
    test_GR_list(test_petite_colony_indices) = test_petite_GR_list;
    test_GR_list(~test_petite_colony_indices) = test_nonpetite_GR_list;
    
    GR_diff_list = test_GR_list - ref_GR_list;

end