function [fx] = iftptox(fp,dp,hbar)

% function iftptox computes the inverse fourier transform fx 
% (a wavfunction in the position variable x)  of the momentum 
% wavefunction fp.  Constants are chosen such that normalization
% is preserved.  The definitions of the fourier transform and inverse 
% transform assumed here have a coefficient of sqrt(1/(2*pi*hbar))
% for both the transform and inverse transform integrals.
%
% dp is the grid spacing of the momentum array.
% hbar is treated as an input, so that by setting hbar=1, a
% function over wavenumber (k) can be used instead of the momentum
% wavefunction.  To work in standard (SI) units, set hbar=1.044e-34.

N  = length(fp);
fx = (N)*dp*(1/2/pi/hbar)^(1/2)*ifftshift(ifft(fftshift(fp)));

end
