function log_det = log_determinant_calculator(current_matrix)

    log_det = 2*sum(log(diag(chol(current_matrix))));

end