function [a0,a,b] = fourier_coeff(x,M)
% Calculate Fourier coefficients using FFT

coeff = fft(x)/length(x);
coeff = coeff(1:M+1);
a0 =  coeff(1);
a  =  2*real(coeff(2:end));
b  = -2*imag(coeff(2:end));

end
