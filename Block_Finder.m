%function [is_block_diagonal,order_vector,block_structure] = Block_Finder(M)
function [order_vector,block_structure] = Block_Finder(M,order_vector,block_start_positions)

	% determines whether symmetric matrix M is block-diagonal and, if
		% so, creates a structure containing each block and its start
		% and end positions
	% M(order_vector,order_vector) contains clear blocks that start at
		% position M(block_start_positions(i),block_start_positions(i))
		% and end at M(block_start_positions(i+1),block_start_positions(i+1))
	% inspired by L-RCM: Matrix block detection algorithm by Miguel Rebollo

	if isnan(order_vector)
		order_vector = symamd(M);
	end
%	disp('order_vector created')
	M_ordered = M(order_vector,order_vector);
	if isnan(block_start_positions)
%		disp('M_ordered created')
		M_ordered_no_diag = M_ordered-diag(diag(M_ordered));
%		disp('M_ordered_no_diag created')
		M_upper_triag = triu(M_ordered_no_diag);
%		disp('M_upper_triag created')
		clear M_ordered_no_diag
		block_start_positions = find(sum(M_upper_triag)==0);
		clear M_upper_triag
	end

	block_end_positions = [(block_start_positions(2:end)-1),size(M,2)];

%	is_block_diagonal = true;

	block_structure = struct;
	for current_block_num = 1:length(block_start_positions)
		current_start = ...
			block_start_positions(current_block_num);
		current_end = ...
			block_end_positions(current_block_num);

		block_structure(current_block_num).start = current_start;
		block_structure(current_block_num).end = current_end;
		% isolate current block on the diagonal
		% since block-diagonal blocks unlikely to be sparse, convert to full
		block_structure(current_block_num).block = full(M_ordered(current_start:current_end,...
			current_start:current_end));

		% check whether there are any off-diagonal elements; if so,
			% return blank structure
%		if current_block_num ~= length(block_start_positions)
%			current_off_diag_horizontal_block = ...
%				M_ordered(current_start:current_end,(current_end+1):end);
%			current_off_diag_horizontal_elements = ...
%				sum(abs(current_off_diag_horizontal_block(:)));
%			if current_off_diag_horizontal_elements > 0
%				block_structure = struct;
%				is_block_diagonal = false;
%				break
%			end
%		end
		
	end

end