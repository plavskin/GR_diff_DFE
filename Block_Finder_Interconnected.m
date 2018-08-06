function biggest_interconnected_block_indices = Block_Finder_Interconnected(M)

	% Finds largest group of mutually inerconnected vertices in a graph using
		% its adjacency matrix
	% determines whether symmetric matrix M is block-diagonal and, if
		% so, creates a structure containing each contiguous (i.e. no
		% zeros) block and its start and end positions
	% M(order_vector,order_vector) contains clear blocks that start at
		% position M(block_start_positions(i),block_start_positions(i))
		% and end at M(block_start_positions(i+1),block_start_positions(i+1))
	% Based on Block_Finder, which was initially inspired by
		% L-RCM: Matrix block detection algorithm by Miguel Rebollo

	M_with_diag = M - diag(diag(M)) + diag(ones(size(M,1),1));
		% M_with_diag is just M, but with ones along the diagonal
	order_vector = symamd(M_with_diag);

	M_ordered = M(order_vector,order_vector);
	M_ordered_no_diag = M_ordered-diag(diag(M_ordered));
	M_upper_triag = triu(M_ordered_no_diag);

	% get rid of non-block connections in M_upper_triag
	M_filled_bottom = M_upper_triag + tril(ones(size(M_upper_triag)));
	% shift M_filled_bottom to the right and fill the missing column with zeros
	left_filter = zeros(size(M_filled_bottom));
	left_filter(:,2:end) = M_filled_bottom(:,1:end-1);
	% shift M_filled_bottom up and fill the missing row with zeros
	bottom_filter = zeros(size(M_filled_bottom));
	bottom_filter(1:end-1,:) = M_filled_bottom(2:end,:);
	M_upper_triag_cleared = M_upper_triag.*left_filter.*bottom_filter;

	% find a subset of all the interconnected blocks, which should include the largest block
	block_sizes_unfiltered = sum(M_upper_triag_cleared,2);
	block_start_positions = find(block_sizes_unfiltered>0);
	block_sizes = block_sizes_unfiltered(block_start_positions) + 1;
	block_end_positions = block_start_positions + block_sizes - 1;

	[~,biggest_block_idx] = max(block_sizes);
	biggest_block_start = block_start_positions(biggest_block_idx);
	biggest_block_end = block_end_positions(biggest_block_idx);

	biggest_interconnected_block_indices = order_vector(biggest_block_start:biggest_block_end);

%	is_block_diagonal = true;
%
%	block_structure = struct;
%	for current_block_num = 1:length(block_start_positions)
%		current_start = ...
%			block_start_positions(current_block_num);
%		current_end = ...
%			block_end_positions(current_block_num);
%		current_size = block_sizes(current_block_num);
%
%		block_structure(current_block_num).start = current_start;
%		block_structure(current_block_num).end = current_end;
%		block_structure(current_block_num).size = current_size;
%		% isolate current block on the diagonal
%		% since block-diagonal blocks unlikely to be sparse, convert to full
%		block_structure(current_block_num).block = full(M_ordered(current_start:current_end,...
%			current_start:current_end));
%
%	end

end