function d_LL_d_sigma = d_LL_d_cov_component(sigma,cov_mat_positions,...
    total_cov,diff_vector,diff_vector_over_cov)

	% Calculates the gradient of one of the random effects that makes
		% up the covariance matrix with respect to the log of the
		% likelihood (for normally-distributed values)
		
	% sigma is the s.d. of the current random effect
	% cov_mat_positions are the positions in the covariance matrix that
		% sigma^2 contributes to
	% total_cov is the covariance matrix
	% diff_vector is the difference between each measurement and its
		%expected value (in a vector equal to the number of measurements)

	if isnan(diff_vector_over_cov)
		diff_vector_over_cov = diff_vector'/total_cov;
	end

    d_LL_d_sigma = -sigma*...
        (trace(cov_mat_positions/total_cov)-...
        diff_vector_over_cov*cov_mat_positions/total_cov*diff_vector);
    
end