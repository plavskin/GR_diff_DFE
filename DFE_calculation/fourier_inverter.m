function output_vector = fourier_inverter(fourier_space_vector, N, L)

	fourier_space_vector( (N/2+1+1):N ) = ...
		conj(fourier_space_vector( (N/2+1-1):-1:2) );
	output_vector = (N/L) * ifft(fourier_space_vector);
	output_vector = output_vector([ (N/2+1):N 1:(N/2) ]);

end