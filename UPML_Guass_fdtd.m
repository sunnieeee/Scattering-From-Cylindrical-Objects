% -*- coding: utf-8 -*-

%Created on Wed Dec 30 15:40:12 2020
%
%@author: David Lyu

%2D-FDTD With PML, excited by Gauss source
close all;
clear all;

% Define Units
meters = 1;
cm     = 1e-2 * meters;
mm     = 1e-3 * meters;
seconds = 1;
hertz   = 1 / seconds;
khz     = 1e3 * hertz;
mhz     = 1e6 * hertz;
ghz     = 1e9 * hertz;

% Parameters For Free Space
e0 = 8.85418782e-12;
u0 = 1.25663706e-6;
N0 = sqrt(u0/e0);
c0 = 299792458 * meters/seconds;
% Parameters For Cylindrical Objects
rad = 20 * meters;
er = 6.25;
ur = 1;
% Parameters for Mesh Grids
nmax = sqrt(er);
NRES = 10;
NPML = [40 40 40 40];

% Parameters For Excitation Pulse
fmax = 15 * ghz;
NFREQS = 500;
FREQS  = linspace(0,fmax,NFREQS);
f0 = 5 * ghz;
lam0 = c0/f0;

lam0_min = c0/fmax/nmax;
dx       = lam0_min/NRES;
dy       = lam0_min/NRES;
Nx = 2*NPML(1) + 100;
Ny = 2*NPML(3) + 100;
Sx = Nx*dx;
Sy = Ny*dy;
% Initialize Cylindrical Area
URxx = ones(Nx,Ny);
URyy = ones(Nx,Ny);
ERzz = ones(Nx,Ny);
sig  = zeros(Nx,Ny);
for ny = 1:Ny
    for nx = 1:Nx
        dist = sqrt((nx-100).^2+(ny-100).^2);
        if dist <= rad
            ERzz(nx,ny) = er;
            sig(nx,ny)  = 0.02;
        end
    end
end

% Calculate stepsize for stability
dmin = min([dx,dy]);
dt = dmin/(2*c0);

period = 1/f0;
Nt     = ceil(period/dt);
dt     = period/Nt;

% Set Source origin
nx1_src = NPML(1) + 1;
nx2_src = Nx - NPML(2);
ny_src  = NPML(3) + 2;

tau  = 0.5/fmax;
t0   = 6*tau;
A    = -sqrt(ERzz(1,ny_src)/URyy(1,ny_src));
delt = 0.5*dy/c0 + dt/2;

d = sqrt(Sx^2 + Sy^2);
proptime = nmax*d/c0;
simtime   = 2*t0 + 10*proptime;
STEPS    = ceil(simtime/dt)-9500;

% Set Guass Source
ta = (0:STEPS-1)*dt;
Ez_src = exp(-((ta-t0)/tau).^2);
Hx_src = A *exp(-((ta-t0+delt)/tau).^2);

% Initial field 
K  = exp(-1i*2*pi*dt*FREQS);
K0 = exp(-1i*2*pi*dt*f0);
Eref  = zeros(Nx,NFREQS);
Etrn  = zeros(Nx,NFREQS);
SRC   = zeros(1,NFREQS);
Eref0 = zeros(Nx,1);
Etrn0 = zeros(Nx,1);
ssSRC = 0;

% Define boundary for calculating
ny_ref = NPML(1) +1;
ny_trn = Ny - NPML(1);
% PML Thickness
Nx2 = 2*Nx;
Ny2 = 2*Ny;

% Calculate sigx
sigx = zeros(Nx2,Ny2);
for nx = 1:2*NPML(1)
    nx1 = 2*NPML(1) - nx + 1;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(1))^3;
end
for nx = 1:2*NPML(2)
    nx1 = Nx2 - 2*NPML(2) + nx;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(2))^3;
end

% Update Parameters for Hy
% Calculate sigy
sigy = zeros(Nx2,Ny2);
for ny = 1:2*NPML(3)
    ny1 = 2*NPML(3) - ny + 1;
    sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(3))^3;
end
for ny = 1:2*NPML(4)
    ny1 = Ny2 - 2*NPML(4) + ny;
    sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(4))^3;
end

% Calculate Parameters for Hx,Hy,Dz,Ez
% Calculate Parameters for Hx
sigHx = sigx(1:2:Nx2,2:2:Ny2);
sigHy = sigy(1:2:Nx2,2:2:Ny2);
mHx0  = (1/dt) + sigHy/(2*e0);
mHx1  = ((1/dt) - sigHy/(2*e0))./mHx0;
mHx2  = -c0./URxx./mHx0;
mHx3  = -(c0*dt/e0) * sigHx./URxx./mHx0;
sigHx = sigx(2:2:Nx2,1:2:Ny2);
sigHy = sigy(2:2:Ny2,1:2:Ny2);
mHy0  = (1/dt) + sigHx/(2*e0);
mHy1  = ((1/dt) - sigHx/(2*e0))./mHx0;
mHy2  = - c0./URyy./mHx0;
mHy3  = - (c0*dt/e0) * sigHy./URyy./mHx0;
% Calculate Parameters for Dz
sigDx = sigx(1:2:Nx2,1:2:Ny2);
sigDy = sigy(1:2:Nx2,1:2:Ny2);
mDz0  = (1/dt) + (sigDx+ sigDy)/(2*e0) + (sigDx.*sigDy*(dt/4/e0^2));
mDz1  = ((1/dt) - (sigDx+ sigDy)/(2*e0) - (sigDx.*sigDy*(dt/4/e0^2)))./mDz0;
mDz2  = c0./mDz0;
mDz3  = - (dt/e0^2)*sigDx.*sigDy./mDz0;
% Calculate Parameters for Ez
mEz1  = sig*dt/e0;

% Initialize Hx,Hy,Dz,Ez
Hx = zeros(Nx,Ny);
Hy = zeros(Nx,Ny);
Dz = zeros(Nx,Ny);
Ez = zeros(Nx,Ny);
% Initialize curl field
CEx = zeros(Nx,Ny);
CEy = zeros(Nx,Ny);
CHz = zeros(Nx,Ny);
% Initialize integral field
ICEx = zeros(Nx,Ny);
ICEy = zeros(Nx,Ny);
IDz  = zeros(Nx,Ny);
IEz  = zeros(Nx,Ny);

%gifname = 'CylinScat_UPML_Guass.gif';

for T = 1:STEPS
    % Update CEx
    for ny = 1:Ny-1
        for nx = 1:Nx
            CEx(nx,ny) = (Ez(nx,ny+1) - Ez(nx,ny))/dy;
        end
    end
    for nx = 1:Nx
        CEx(nx,Ny) =(Ez(nx,1) - Ez(nx,Ny))/dy;
    end
    % Update CEy
    for nx = 1:Nx-1
        for ny= 1:Ny
            CEy(nx,ny) = -(Ez(nx+1,ny) -Ez(nx,ny))/dx;
        end
    end
    for ny = 1:Ny
        CEy(Nx,ny) = -(Ez(1,ny) - Ez(Nx,ny))/dx;
    end
    % Import Source Ez
    CEx(nx1_src:nx2_src,ny_src-1) = CEx(nx1_src:nx2_src,ny_src-1) - Ez_src(T)/dy;
    % Update ICEx,ICEy
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;
    % Update Hx,Hy
    Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
    Hy = mHy1.*Hy + mHy2.*CEy + mHy3.*ICEy;
    % Calculate CHz
    CHz(1,1) = (Hy(1,1) - Hy(Nx,1))/dx - (Hx(1,1) - Hx(1,Ny))/dy;
    for nx = 2:Nx
        CHz(nx,1) = (Hy(nx,1) - Hy(nx-1,1))/dx - (Hy(nx,1) - Hx(nx,Ny))/dy;
    end
    for ny = 2:Ny
        CHz(1,ny) = (Hy(1,ny) - Hy(Nx,ny))/dx - (Hx(1,ny) - Hx(1,ny-1))/dy;
        for nx = 2:Nx
            CHz(nx,ny) = (Hy(nx,ny) - Hy(nx-1,ny))/dx - (Hx(nx,ny) - Hx(nx,ny-1))/dy;
        end
    end
    
    % Add Source Hx
    CHz(nx1_src:nx2_src,ny_src) = CHz(nx1_src:nx2_src,ny_src) - Hx_src(T)/dy;
    % Update IDz,IEz
    IDz = IDz + Dz;
    IEz = IEz + Ez;
    % Update Dz
    Dz = mDz1.*Dz + mDz2.*CHz + mDz3.*IDz;
    % Update Ez
    Ez = (Dz - mEz1.*IEz)./(ERzz + mEz1);
    % Update equations
    for nfreq = 1:NFREQS
        Eref(:,nfreq) = Eref(:,nfreq) + (K(nfreq)^T) * Ez(:,ny_ref) * dt;
        Etrn(:,nfreq) = Etrn(:,nfreq) + (K(nfreq)^T) * Ez(:,ny_trn) * dt;
        SRC(nfreq) = SRC(nfreq) + (K(nfreq)^T) *Ez_src(T) * dt;
    end
    % Update equations at f0
    Eref0 = Eref0 + (K0^T) *Ez(:,ny_ref) * dt;
    Etrn0 = Etrn0 + (K0^T) *Ez(:,ny_trn) * dt;
    ssSRC = ssSRC + (K0^T) *Ez_src(T) * dt;
    
    % Background Color For Cylindrical Objects and Surroundings
    ERkk = Ez;
    for ny = 1:Ny
        for nx = 1:Nx
            dist = sqrt((nx-100).^2+(ny-100).^2);
            if dist <= rad
                ERkk(nx,ny)=-0.05;
            end
        end
    end
    % Visualization
    % Plot at each 20 steps
    if mod(T,20) == 0
        timestep = int2str(T);
        hold on;
        s1 = pcolor(Ez);
        hold on;
        s2 = pcolor(ERkk);
        hold off;
        alpha(s2,0.1);
        shading flat;
        axis([1 Nx 1 Ny]);
        c = colorbar;
        caxis([-0.05,0.05]);
        set(get(c,'title'),'string','V/m');
        axis image;
        axis equal tight;
        title(['Ez at time step = ',timestep]);
        pause(0.01);
% Plot gif
%{        
        frame = getframe;
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256);
        if T == 20
            imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.2);
        else
            imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.2);
        end
%}
    end
end





