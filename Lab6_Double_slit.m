%% inputs 
clear; close all;

% defin "hbar". The first value is the SI unit value, the second is if we
% are FTing from space to k-space and not momentum space.
hbar = [1.044e-34,1]; 

%% grid setup

% Number of grid points need 2^16 or 2^17 to resolve the grating structure
N = 2^17;
% maximum spatial grid extent
xmax = 10e-3;        % size of grid is 100 times expected maximum.
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

%% important parameters from Ref2
z_sa = 64e-2;                    % [m]
z_ad = 64e-2;                       % [m]
source_slit = 2e-6;             % [um]
%%%%%%%
% the detector size might be 10periods * 8um/period = 80um
detector_size = 80e-6;          % half distance is 40um.
lam_dB = 1.03e-10;             % [Angstrom]deBroigle wavelength
m_n = 1.675e-27;                % [kg] Neutron mass
d = 8e-6;                     % periodicity
diffract_slit = 1e-6;        % [um] corresponds to a duty cycle of 10%
slit_separation = 8e-6;     % center to center separation of the two slits.
fill_factor = 50;

%% Propagation vectors, grating setup
% In your simulation, you will calculate probability density distributions
% He atoms
m_He = 2*m_n;
% time evolution operator setup
% define a propagation vector for z-direction
Nz = 5000;
zend = z_ad;
% define a step dz = dt*v where v is particle velocity.
% let dz = dx
dz = zend/Nz;
z = 0:dz:zend;
% define velocity as h/m/lam_dB where h is hbar/2/pi. Also assume t0 = 0
v = 2*pi*hbar(1)/m_He/lam_dB;        % [m/s]
% define time step
dt = dz/v;
% define differential time evolution operator dU
% dU = exp(-1i*pn.^2.*dt./2./m_n./hbar(1));
U = dU(dt,pn,m_He,hbar(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that the wavefunction psi0 is the wavefunction immediately         %  
% following the diffracting slit.                                         %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grating setup
rect1 = @(x) double(abs(x) <= diffract_slit/2);  % This creates a rect function from -0.5 to 0.5
dslit = rect1(xn - slit_separation/2) + rect1(xn + slit_separation/2);
% figure;
% plot(xn,dslit,'LineWidth',1.5);
x_grating = (2*pi/d).*xn;
% 
% % grating vector. The amplitude scaling and shifting by 0.5 is to ensure
% % that the wavefrom is between 0 and 1 with an amplitude of 1.
gv = 0.5.*square(x_grating,fill_factor) + 0.5;


%% Spherical wave calculations

% diffracting aperture function
slit = dslit;
% time step for propagation from source slit to diffracting slit
dt = z_sa/v;
% source wavefunction
psi0source = rect(xn,source_slit);
psi0spherical = slit.*iftptox(dU(dt,pn,m_He,hbar(1)).*ftxtop(psi0source,dx,hbar(1)),dp,hbar(1));

psi0spherical_p = ftxtop(psi0spherical,dx,hbar(1));
Psi0spherical_p = zeros(length(z),length(psi0spherical_p));
Psi0spherical_x = zeros(length(z),length(psi0spherical_p));
for zi = 1:length(z)
    dt = z(zi)./v;
    Psi0spherical_p(zi,:) = dU(dt,pn,m_He,hbar(1)).*psi0spherical_p;
    Psi0spherical_x(zi,:) = iftptox(Psi0spherical_p(zi,:),dp,hbar(1));
end

Psi0 = gv.*Psi0spherical_x(end,:);
ap = (xn >= -detector_size/2 & xn <= detector_size/2);
x_limited = xn(ap);
figure;
imagesc([0,z_ad],[-detector_size/2, detector_size/2],abs(Psi0spherical_x(:,ap)').^2);
set(gca,'FontSize',15);
xlabel('z');
ylabel('x');
title('$|\Psi(z,t)|^2$ : spherical wave',Interpreter='latex');
colormap turbo;
colorbar;

% wavefunction at detector plane
% THIS ISN'T QUITE CORRECT YET. ASK QUESTIONS.
figure;
plot(xn,abs(Psi0).^2,'LineWidth',1);
xlabel('x_n');
ylabel('$|\Psi(1,t)|^2$',Interpreter='latex');
title('Wavefunction at detector plane (spherical)');
xlim([-detector_size/2,detector_size/2]);
set(gca,'FontSize',15);


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

