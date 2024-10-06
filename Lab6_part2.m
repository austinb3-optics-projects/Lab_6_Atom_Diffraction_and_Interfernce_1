%% inputs 
clear; close all;

% defin "hbar". The first value is the SI unit value, the second is if we
% are FTing from space to k-space and not momentum space.
hbar = [1.044e-34,1]; 

%% grid setup

% Number of grid points 
N = 2^14;
% maximum spatial grid extent
xmax = 8e-3;        % size of grid is 100 times expected maximum.
% spatial grid increment
dx = xmax/N;
% unit step vector to define grid.
n = 0:1:(N-1); p = n;
% spatial grid vector
xn = -xmax/2 + n*dx;

% momentum space grid.
% from DFFT.pdf dK = 2pi/xmax. To translate to momentum, use p = hbar*K
% This is Nyquist sampling
% define maximum momentum on grid
pmax = 2*pi*hbar(1)/dx;
% momentum step of grid.
dp = 2*pi*hbar(1)/xmax;
% define the p-space grid
pn = -pmax/2 + p.*dp;

%% important parameters from Ref1
z_sa = 0.96;                    % [m]
z_ad = 1;                       % [m]
source_slit = 2.5e-9;           % [nm]
diffract_slit = 25.4e-6;        % [um] nominal value (will vary)
detector_size = 80e-6;          % half distance is 40um.
lam_dB = 0.175e-10;             % [Angstrom]deBroigle wavelength
m_n = 1.675e-27;                % [kg] Neutron mass

%% Setting up wavefunctions for II: Diffraction of an atomic beam due to a narrow slit
% In your simulation, you will calculate probability density distributions
% that correspond to those of Figs. 2 and 3 of Ref1.
m_k = 39*m_n;
% time evolution operator setup
% define a propagation vector for z-direction
Nz = 1000;
zend = 1;
% define a step dz = dt*v where v is particle velocity.
% let dz = dx
dz = zend/Nz;
z = 0:dz:zend;
% define velocity as h/m/lam_dB where h is hbar/2/pi. Also assume t0 = 0
v = 2*pi*hbar(1)/m_k/lam_dB;        % [m/s]
% define time step
dt = dz/v;
% define differential time evolution operator dU
% dU = exp(-1i*pn.^2.*dt./2./m_n./hbar(1));
U = dU(dt,pn,m_k,hbar(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that the wavefunction psi0 is the wavefunction immediately         %  
% following the diffracting slit.                                         %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% wavefunctions




%% Functions 
function [U] = dU(dt,p,m,hbar)
    U = exp(-1i.*(p.^2).*dt./2./m./hbar);
end

function y = rect(t,bound)
    y = abs(t) <= bound/2;
    y = double(y);  % make sure the output datatype is double and not logical
end

function [fp] = ftxtop(fx,dx,hbar)
    N  = length(fx);
    fp = dx*(1/2/pi/hbar)^(1/2)*fftshift(fft(ifftshift(fx)));
end

function [fx] = iftptox(fp,dp,hbar)
    N  = length(fp);
    fx = (N)*dp*(1/2/pi/hbar)^(1/2)*ifftshift(ifft(fftshift(fp)));
end

