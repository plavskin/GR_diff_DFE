function [f, g] = fourier_space_creator(N, L)

	% Sets up frequncy space to be used to calculate discrete characteristic functions

	% discretized increment between frquencies in the discrete
        % Fourier transform

	df=1/L;
	% discretized increment between values of the argument for the
		% density fucntion
    dg=L/N;
	    
	N_space = (-N/2):(N/2-1);

	% array of argument values for the density function
	g = N_space * dg;

	% discrete frequencies for Fourier transform
	f = -2 * pi * df * N_space;

end