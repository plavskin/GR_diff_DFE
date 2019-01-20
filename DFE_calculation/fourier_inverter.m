function output_vector = fourier_inverter(fourier_space_vector, N, L)
	
	% Transforms frequency space (as set up by fourier_space_creator)
		% to time space
	output_vector = fftshift(ifft(fftshift(N/L * fourier_space_vector)));

end