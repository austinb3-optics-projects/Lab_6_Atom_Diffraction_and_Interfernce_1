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

%% Plane wave wavefunction
% dim is the dimension over which to perform the fft
dim = 2;        % performed on the rows.

psi0_pw = rect(xn,diffract_slit);       % spatial wavefunction immediately after the aperture
psi0_pw_p = ftxtop(psi0_pw,dx,hbar(1));
Psi0_pw_p = zeros(length(z),length(psi0_pw_p));
Psi0_pw_x = zeros(length(z),length(psi0_pw_p));
for zi = 1:length(z)
    dt = z(zi)./v;
    Psi0_pw_p(zi,:) = dU(dt,pn,m_k,hbar(1)).*psi0_pw_p;
    Psi0_pw_x(zi,:) = iftptox(Psi0_pw_p(zi,:),dp,hbar(1));
end
ap = (xn >= -detector_size/2 & xn <= detector_size/2);
x_limited = xn(ap);
% wavefunction evolution over propagation distance
% figure;
% imagesc([0,z_ad],[-detector_size/2, detector_size/2],abs(Psi0_pw_x(:,ap)').^2);
% set(gca,'FontSize',15);
% xlabel('z');
% ylabel('x');
% title('$|\Psi(z,t)|^2$ : plane wave',Interpreter='latex');
% colormap turbo;
% colorbar;

% wavefunction at detector plane
% figure;
% plot(xn,abs(Psi0_pw_x(end,:)).^2,'LineWidth',1.5);
% xlabel('x_n');
% ylabel('$|\Psi(1,t)|^2$',Interpreter='latex');
% title('Wavefunction at detector plane (plane)');
% set(gca,'FontSize',15);
% xlim([-detector_size/2,detector_size/2]);

%% Spherical wave implementation
% diffracting aperture function
slit = rect(xn,diffract_slit);
% time step for propagation from source slit to diffracting slit
dt = z_sa/v;
% source wavefunction
psi0source = rect(xn,source_slit);
psi0spherical = slit.*iftptox(dU(dt,pn,m_k,hbar(1)).*ftxtop(psi0source,dx,hbar(1)),dp,hbar(1));

psi0spherical_p = ftxtop(psi0spherical,dx,hbar(1));
Psi0spherical_p = zeros(length(z),length(psi0spherical_p));
Psi0spherical_x = zeros(length(z),length(psi0spherical_p));
for zi = 1:length(z)
    dt = z(zi)./v;
    Psi0spherical_p(zi,:) = dU(dt,pn,m_k,hbar(1)).*psi0spherical_p;
    Psi0spherical_x(zi,:) = iftptox(Psi0spherical_p(zi,:),dp,hbar(1));
end

% Psi0spherical_x_limited = Psi0spherical_x(:,((-detector_size/2 <= xn) && (xn <= detector_size/2)));
% wavefunction evolution over propagation distance
% figure;
% imagesc([0,z_ad],[-detector_size/2, detector_size/2],abs(Psi0spherical_x(:,ap)').^2);
% set(gca,'FontSize',15);
% xlabel('z');
% ylabel('x');
% title('$|\Psi(z,t)|^2$ : spherical wave',Interpreter='latex');
% colormap turbo;
% colorbar;

% wavefunction at detector plane
% figure;
% plot(xn,abs(Psi0spherical_x(end,:)).^2,'LineWidth',1);
% xlabel('x_n');
% ylabel('$|\Psi(1,t)|^2$',Interpreter='latex');
% title('Wavefunction at detector plane (spherical)');
% xlim([-detector_size/2,detector_size/2]);
% set(gca,'FontSize',15);

%% Questions from lab handout:

% 1. I notice how much less resolved the outputs are when I decrease the
% number of points by a factor of two, and on the other hand, how much more
% resolution I get when I increase the number of points by a factor two.
% 2. Diffraction of a spherical wave from an aperture is different than
% diffraction by an incident plane wave. Something to do with the fact that
% the wavefronts of the plane wave are flat at the aperture while the
% spherical wave has curved wavefronts.
% 3. changing the source size and aperture size
source_slit = 2.5e-6;
diffract_slit = 6.35e-6;
% diffracting aperture function
slit = rect(xn,diffract_slit);
% time step for propagation from source slit to diffracting slit
dt = z_sa/v;
% source wavefunction
psi0source = rect(xn,source_slit);
psi0spherical = slit.*iftptox(dU(dt,pn,m_k,hbar(1)).*ftxtop(psi0source,dx,hbar(1)),dp,hbar(1));

psi0spherical_p = ftxtop(psi0spherical,dx,hbar(1));
Psi0spherical_p = zeros(length(z),length(psi0spherical_p));
Psi0spherical_x = zeros(length(z),length(psi0spherical_p));
for zi = 1:length(z)
    dt = z(zi)./v;
    Psi0spherical_p(zi,:) = dU(dt,pn,m_k,hbar(1)).*psi0spherical_p;
    Psi0spherical_x(zi,:) = iftptox(Psi0spherical_p(zi,:),dp,hbar(1));
end

% wavefunction at detector plane
figure;
plot(xn,abs(Psi0spherical_x(end,:)).^2,'LineWidth',1);
xlabel('x_n');
ylabel('$|\Psi(1,t)|^2$',Interpreter='latex');
title('Wavefunction at detector plane (spherical)');
xlim([-detector_size/2,detector_size/2]);
set(gca,'FontSize',15);

% diffracting aperture size = 2.54um
diffract_slit = 2.54e-6;
slit = rect(xn,diffract_slit);
psi0spherical = slit.*iftptox(dU(dt,pn,m_k,hbar(1)).*ftxtop(psi0source,dx,hbar(1)),dp,hbar(1));

psi0spherical_p = ftxtop(psi0spherical,dx,hbar(1));
Psi0spherical_p = zeros(length(z),length(psi0spherical_p));
Psi0spherical_x = zeros(length(z),length(psi0spherical_p));
for zi = 1:length(z)
    dt = z(zi)./v;
    Psi0spherical_p(zi,:) = dU(dt,pn,m_k,hbar(1)).*psi0spherical_p;
    Psi0spherical_x(zi,:) = iftptox(Psi0spherical_p(zi,:),dp,hbar(1));
end

% wavefunction at detector plane
figure;
plot(xn,abs(Psi0spherical_x(end,:)).^2,'LineWidth',1);
xlabel('x_n');
ylabel('$|\Psi(1,t)|^2$',Interpreter='latex');
title('Wavefunction at detector plane (spherical)');
xlim([-detector_size/2,detector_size/2]);
set(gca,'FontSize',15);

%% movie section

% for i = 1:50
% plot(xn,abs(Psi0spherical_x(i,:)).^2,'LineWidth',1.5);
% xlabel('x_n');
% ylabel('$|\Psi(1,t)|^2$',Interpreter='latex');
% title("(spherical)"+ num2str(i));
% xlim([-detector_size/2,detector_size/2]);
% set(gca,'FontSize',15);
% 
% drawnow;
% M(i) = getframe;
% end
% movie(M,2,1)

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

