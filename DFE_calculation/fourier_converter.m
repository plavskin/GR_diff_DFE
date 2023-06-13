function output_vector = fourier_converter(time_space_vector, N, L)
	
	% Transforms frequency space (as set up by fourier_space_creator)
		% to time space
	output_vector = fftshift(fft(fftshift(L/N * time_space_vector)));

end