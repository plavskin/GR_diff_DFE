function compressed_point_mat = Combine_Nearby_Points(point_mat, tolerance)
	% Takes a matrix point_mat where each row is a point in 
		% column-dimnesional space and combines points that are within
		% tolerance of each other in each dimension

	space_for_improvement = true;
	while space_for_improvement
		% initialize adjacency_matrix with every point interconnected, and
			% update iteratively for each dimension
		% if points aren't interconnected in any dimension, set
			% interconnection to 0
		point_num = size(point_mat,1);
		adjacency_matrix = ones(point_num);

		for column_counter = 1:size(point_mat,2)
			current_column = point_mat(:,column_counter);
			current_distance_mat = ...
				abs(repmat(current_column,[1,point_num])-...
					repmat(current_column',[point_num,1]));
			within_tolerance = current_distance_mat <= abs(tolerance);
			adjacency_matrix = adjacency_matrix.*within_tolerance;
		end

		if sum(adjacency_matrix(:)) > point_num
			% Get list of largest group of points within tolerance of each other
			interconnected_points = Block_Finder_Interconnected(adjacency_matrix);
			% Replace interconnected points in point_mat with a single
				% point that is the average location of the
				% interconnected points
			interconnected_mat = point_mat(interconnected_points,:);
			average_point = mean(interconnected_mat,1)
			point_mat = removerows(point_mat,'ind',interconnected_points);
			point_mat = [point_mat; average_point];
		else
			space_for_improvement = false;
		end
	end

	compressed_point_mat = point_mat;
end

