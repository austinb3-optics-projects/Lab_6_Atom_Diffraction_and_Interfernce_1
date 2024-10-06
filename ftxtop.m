function [fp] = ftxtop(fx,dx,hbar)

% function ftptox computes the fourier transform fp 
% (a wavfunction in the momentum variable p)  of the position 
% wavefunction fx.  Constants are chosen such that normalization
% is preserved.  The definitions of the fourier transform and inverse 
% transform assumed here have a coefficient of sqrt(1/(2*pi*hbar))
% for both the transform and inverse transform integrals.
%
% dx is the grid spacing of the position array.
% hbar is treated as an input, so that by setting hbar=1, a
% function over wavenumber (k) can be obtained instead of the momentum
% wavefunction.  To work in standard (SI) units, set hbar=1.044e-34.

N  = length(fx);
fp = dx*(1/2/pi/hbar)^(1/2)*fftshift(fft(ifftshift(fx)));

end

