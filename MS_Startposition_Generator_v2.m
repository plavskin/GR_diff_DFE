function startpoints = MS_Startposition_Generator_v2(indices_to_multistart,...
	ms_positions,global_start_values,global_min_array,global_max_array)

	% Generates a custom start point set object for multistart that
		% includes the assigned start position but is otherwise spread
		% evenly across allowed parameter space

	global_param_number = length(indices_to_multistart);

	if ms_positions > 1
		N = cell([1,global_param_number]);
		[N{indices_to_multistart}] = ndgrid(linspace(0,1,ms_positions));
		[N(~indices_to_multistart)] = num2cell(global_start_values(~indices_to_multistart));

		number_multistart_parameters = sum(indices_to_multistart);
		pt_matrix = zeros([(ms_positions^number_multistart_parameters+1),global_param_number]);
		pt_matrix(1,:) = global_start_values;

		for counter = 1:global_param_number
			if indices_to_multistart(counter)
				current_lb = global_min_array(counter);
				current_ub = global_max_array(counter);
				current_interval = current_ub - current_lb;
				N{counter} = N{counter}*current_interval+current_lb;
			end
			current_mat = N{counter};
			pt_matrix(2:end,counter) = current_mat(:);
		end
	else
		pt_matrix = global_start_values;
	end

	startpoints = CustomStartPointSet(pt_matrix);
end