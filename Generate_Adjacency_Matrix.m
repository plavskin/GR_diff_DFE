function M = Generate_Adjacency_Matrix(p_connection, n)

	% Generates a symmetrical adjacency matrix based on a probability
		% of connection between vertices is p_connection, and which has
		% n vertices total

	initial_mat = binornd(1,p_connection,[n,n]);
	ones_matrix = ones(n);

	% reflect matrix about diagonal
	M = triu(initial_mat) + triu(initial_mat)'.*(ones_matrix-triu(ones_matrix));
	
