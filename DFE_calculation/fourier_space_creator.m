function [f, g] = fourier_space_creator(N, L)

	df=1/L;
    % discretized increment between frquencies in the discrete
        % Fourier transform
	f=(0:(N-1))*df;
	    % discrete frequencies for Fourier transform
	dg=L/N;
	    % discretized increment between values of the argument for the
	        % density fucntion
	g=( (-N/2):(N/2-1) )*dg;
	    % array of argument values for the density function

end