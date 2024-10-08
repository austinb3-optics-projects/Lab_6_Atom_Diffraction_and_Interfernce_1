%% inputs 
clear; close all;

% defin "hbar". The first value is the SI unit value, the second is if we
% are FTing from space to k-space and not momentum space.
hbar = [1.044e-34,1]; 

%% grid setup

% Number of grid points need 2^16 or 2^17 to resolve the grating structure
N = 2^18;
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
z_12 = 0.25;                    % [m]
z_23 = 0.25;                       % [m]
source_slit = 2.5e-6;             % [um]
%%%%%%%
% the detector size might be 10periods * 8um/period = 80um
detector_size = 100e-6;          % half distance is 40um.
lam_dB = 0.12e-10;             % [Angstrom]deBroigle wavelength
m_n = 1.675e-27;                % [kg] Neutron mass
d = 8e-6;                     % periodicity
diffract_slit = 1e-6;        % [um] corresponds to a duty cycle of 10%
slit_separation = 8e-6;     % center to center separation of the two slits.
fill_factor = 50;
lam = 811e-9;               % wavelength of light for the phase gratings
k_vec = 2*pi/lam;           % photon k-vector
%% Propagation vectors, grating setup
% In your simulation, you will calculate probability density distributions
% Argon atoms
m_Ar = 40*m_n;
% time evolution operator setup
% define a propagation vector for z-direction
Nz = 999;
zend1 = z_12;
zend2 = z_23;
% define a step dz = dt*v where v is particle velocity.
% let dz = dx
dz = zend1/Nz;
z1 = 0:dz:zend1;
z2 = 0:dz:zend2;
% define velocity as h/m/lam_dB where h is hbar/2/pi. Also assume t0 = 0
v = 850;        % [m/s]
% define time step
dt = dz/v;
% define differential time evolution operator dU
% dU = exp(-1i*pn.^2.*dt./2./m_n./hbar(1));
U = dU(dt,pn,m_Ar,hbar(1));
phi01 = 2.56;
phi02 = 4.34;       % maximum phase shifts of the phase gratings
dphi1 = phi01.*(cos(k_vec.*xn)).^2;
dphi2 = phi02.*(cos(k_vec.*xn)).^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that the wavefunction psi0 is the wavefunction immediately         %  
% following the diffracting slit.                                         %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grating setup
rect1 = @(x) double(abs(x) <= source_slit/2);  % This creates a rect function from -0.5 to 0.5
source = rect1(xn);
phase_grating1 = exp(1i.*dphi1);
% figure;
% plot(xn,phase_grating1);
% xlim([-0.5e-4,0.5e-4]);
phase_grating2 = exp(1i.*dphi2);
% figure;
% plot(xn,phase_grating2);
% phi0*sin(4*pi*x/lam)
% xlim([-0.5e-4,0.5e-4]);
%% Facts from the paper

% - period of the standing wave phase gratings is 405nm (lam/2)
% - collimation slits have a size of 5um. This will correspond the the width
% of the plane wave input.
% - distance to first grating is 85cm.
% - distance from grating 1 to 2 = 25cm
% - distance from grating 2 to 3 = 25cm
% - phase shift experienced dphi = phi0*cos^2(kx)
% - The modulation is akin to a sinusoidal phase grating with period lam/2
% from the light.
% - phi0 for 1 and 3 is 2.56
% - phi0 for 2 is 4.34



psi0_pw = source;       % spatial wavefunction immediately after the aperture
% figure; 
% plot(xn,psi0_pw);
psi1_x = phase_grating1.*psi0_pw;  % wavefunction immediately after grating 1.
% figure;
% plot(xn,psi1_x);
psi2_p = ftxtop(psi1_x,dx,hbar(1));   % fourier transform after grating 1.

% wavefunctions to hold z dependence up to the second grating.
Psi0_p1 = zeros(length(z1),length(psi2_p));
Psi0_x1 = zeros(length(z1),length(psi2_p));

for zi = 1:length(z1)
    dt = z1(zi)./v;
    Psi0_p1(zi,:) = dU(dt,pn,m_Ar,hbar(1)).*psi2_p;
    Psi0_x1(zi,:) = iftptox(Psi0_p1(zi,:),dp,hbar(1));
end

ap = (xn >= -detector_size/2 & xn <= detector_size/2);
x_limited = xn(ap);
% figure;
% % imagesc([0,z_ad],[-detector_size/2, detector_size/2],abs(Psi0_pw_x').^2);
% imagesc([0,z_12],[-detector_size/2, detector_size/2],abs(Psi0_x1(:,ap)').^2);
% set(gca,'FontSize',15);
% xlabel('z');
% ylabel('x');
% title('$|\Psi(z,t)|^2$ : plane wave',Interpreter='latex');
% colormap turbo;
% colorbar;


% Second grating propagation
psi0_2 = Psi0_x1(end,:);        % take the final row of the propagation to grating 2 as the initial wavefunction before grating 2.
psi1_2 = phase_grating2.*psi0_2;        % wavefunction immediately following grating 2.
psi1_2p = ftxtop(psi1_2,dx,hbar(1));

% wavefunctions to hold z dependence up to the second grating.
Psi1_p1 = zeros(length(z2),length(psi1_2p));
Psi1_x1 = zeros(length(z2),length(psi1_2p));

for zi = 1:length(z2)
    dt = z2(zi)./v;
    Psi1_p1(zi,:) = dU(dt,pn,m_Ar,hbar(1)).*psi1_2p;
    Psi1_x1(zi,:) = iftptox(Psi1_p1(zi,:),dp,hbar(1));
end

% figure;
% % imagesc([0,z_ad],[-detector_size/2, detector_size/2],abs(Psi0_pw_x').^2);
% imagesc([0,z_23],[-detector_size/2, detector_size/2],abs(Psi1_x1(:,ap)').^2);
% set(gca,'FontSize',15);
% xlabel('z');
% ylabel('x');
% title('$|\Psi(z,t)|^2$ : plane wave',Interpreter='latex');
% colormap turbo;
% colorbar;

Psi = vertcat(Psi0_x1,Psi1_x1);
% ap = (xn >= -detector_size/2 & xn <= detector_size/2);
% x_limited = xn(ap);
% 
% % wavefunction evolution over propagation distance
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This takes awhile to run                                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
imagesc([0,(z_12 + z_23)],[-detector_size/2, detector_size/2],abs(Psi(:,ap)').^2);
hold on;
xline(z_12,'LineWidth',3);
set(gca,'FontSize',15);
xlabel('z');
ylabel('x');
title('$|\Psi(z,t)|^2$ : plane wave',Interpreter='latex');
colormap turbo;
colorbar;

figure;
plot(xn,abs(Psi(end,:)).^2,'LineWidth',1.5);
xlabel('x_n');
ylabel('$|\Psi(1,t)|^2$',Interpreter='latex');
title('Wavefunction at talbot distance (plane)');
set(gca,'FontSize',15);
xlim([-detector_size/2,detector_size/2]);


%% Questions

% - What happens when you reduce the width of the source by a factor of
% two. What do the arms look like.
% Answer: The arms of the interferometer are more obvious.
% though there is reduced contrast of the arms. Maybe it's better to
% describe it as decreased diffraction efficiency.

%% Functions 
function [U] = dU(dt,p,m,hbar)
    U = exp(-1i.*(p.^2).*dt./2./m./hbar);
end

function [fp] = ftxtop(fx,dx,hbar)
    N  = length(fx);
    fp = dx*(1/2/pi/hbar)^(1/2)*fftshift(fft(ifftshift(fx)));
end

function [fx] = iftptox(fp,dp,hbar)
    N  = length(fp);
    fx = (N)*dp*(1/2/pi/hbar)^(1/2)*ifftshift(ifft(fftshift(fp)));
end

