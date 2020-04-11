clear
close all
%% Test Setup
SourceType = 1; % 1 = monopole, 2 = vertical dipole, 3 = down-range horizontal dipole
f = 100;        % Source frequency, Hz
zs = 30;        % Source depth, m
rzfact = 2;     % Ratio of range step to depth step
R0 = 1;         % Reference range, m
Rmax = 4000;    % Maximum calculation range, m
Hmax = 400;     % Maximum calculation depth, m 
c0 = 1500;      % Reference sound speed, m/s
rho0 = 1000;    % Reference density, kg/m^3
%% Environment Setup
c = 1500;       % Water sound speed, m/s
rho = 1000;     % Water density, m/s
c1 = 1700;      % Bottom sound speed, m/s
rho1 = 1500;    % Bottom density, kg/m^3
alpha1 = 0.0;   % Bottom attenuation, dB/wavelength
% Create ocean and bottom profile structures
ocean.ssp_zr = scatteredInterpolant([0;Hmax;0;Hmax],[0;0;Rmax;Rmax;],c*ones(4,1),'natural','nearest');
ocean.rho = rho;
bottom.bathy_r = griddedInterpolant([0;1e3;2e3;3e3],[100;100;200;200],'linear','nearest');
bottom.c = c1;
bottom.rho = rho1;
bottom.alpha = alpha1;
bottom.cs = 0;
bottom.alphas = 0;
% Create rough surface realization
surface = RoughSurface1D(2,4e3,8192,100,2,1,0,1,10);
% Clean up
clear c rho c1 rho1 alpha1
%% Run M3PE2D
[p,vz,vr,z,r,surf,bathy] = MNPE2D(SourceType,f,zs,rzfact,R0,Rmax,c0,rho0,ocean,bottom,surface);
%% Calculate Transmission Loss
TLp = -20*log10(abs(p));
TLr = -20*log10(abs(vr*rho0*c0));
TLz = -20*log10(abs(vz*rho0*c0));
%% Plot Transmission Loss Field
figure
imagesc(r/1e3,z,TLp)
hold on
plot(r/1e3,surf,'w','LineWidth',2)
plot(r/1e3,bathy,'w','LineWidth',2)
hold off
colorbar
colormap(flipud(colormap));
caxis([0 120])
xlim([0 Rmax/1000])
ylim([min(z) 250])
xlabel('Range (km)')
ylabel('Depth (m)')
title('Pressure Transmission Loss (dB//R0)')
%% Plot Transmission Loss at 50 m
zout = 50;
nout = find(z>=zout,1);
if (zout-z(nout-1))<(z(nout)-zout)
    nout = nout-1;
end
figure
plot(r/1000,TLp(nout,:),'r');
axis ij;
grid on
xlim([0 Rmax/1000])
ylim([30 110])
xlabel('Range (km)')
ylabel('Transmission Loss (dB//R0)')
title(['Transmission Loss at ' num2str(zout) ' m'])
%% Plot Velocity Transmission Loss Fields
figure
imagesc(r/1e3,z,TLr)
hold on
plot(r/1e3,surf,'w','LineWidth',2)
plot(r/1e3,bathy,'w','LineWidth',2)
hold off
colorbar
colormap(flipud(colormap));
caxis([0 120])
ylim([min(z) 250])
xlim([0 Rmax/1000])
xlabel('Range (km)')
ylabel('Depth (m)')
title('Down-Range Horizontal Velocity Transmission Loss (dB//R0)')
figure
imagesc(r/1e3,z,TLz)
hold on
plot(r/1e3,surf,'w','LineWidth',2)
plot(r/1e3,bathy,'w','LineWidth',2)
hold off
colorbar
colormap(flipud(colormap));
caxis([0 120])
ylim([min(z) 250])
xlim([0 Rmax/1000])
xlabel('Range (km)')
ylabel('Depth (m)')
title('Vertical Velocity Transmission Loss (dB//R0)')