function [ Phi ] = fourier_dict( N, L )
%FOURIER_DICT An l-fold overdetermined fourier dictionary with waveforms
%(atoms) of length n
%   See Chen, Donoho and Saunders (1998)
P = N*L;
n = [0:(N-1)];
Pc = (0:(P/2));
Ps = (1:(P/2 - 1));
Phi = [cos(2*pi/P*n'*Pc) sin(2*pi/P*n'*Ps)];
end

