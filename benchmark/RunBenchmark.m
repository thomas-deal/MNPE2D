clear
close all
%% Path Setup
path(pathdef)
BENCHPath = pwd;
cd('..')
MNPEPath = pwd;
BINPath = fullfile(MNPEPath,'bin');
MDLPath = fullfile(MNPEPath,'mdl');
cd(BENCHPath)
addpath(BINPath)
addpath(MDLPath)
%% Plot Setup
cax = [30 90];
yax = [0 250];
%% Test Setup
SourceType = 1; % 1 = monopole, 2 = vertical dipole, 3 = down-range horizontal dipole
f = 200;        % Source frequency, Hz
zs = 100;       % Source depth, m
rzfact = 2;     % Ratio of range step to depth step
R0 = 1;         % Reference range, m
Rmax = 4000;    % Maximum calculation range, m
Hmax = 400;     % Maximum calculation depth, m 
c0 = 1485;      % Reference sound speed, m/s
rho0 = 1000;    % Reference density, kg/m^3
%% Environment Setup
H = 200;        % Water depth, m
c = 1500;       % Water sound speed, m/s
rho = 1000;     % Water density, kg/m^3
c1 = 1700;      % Bottom sound speed, m/s
rho1 = 1500;    % Bottom density, kg/m^3
alphaf = 0.333; % Bottom attenuation, dB/m/kHz
alpha1 = ConvertAttenuationUnits(alphaf,'dB/m/kHz','dB/lambda',f,c1);
                % Bottom attenuation, dB/wavelength
%% Run Benchmark
[B.TLp, B.TLvz, B.r, B.z] = ModeTLwithLoss(f,zs,H,c,rho,c1,rho1,alphaf,Rmax);
B.surf = zeros(1,length(B.r));
B.bathy = H*ones(1,length(B.r));
%% Run Matlab MNPE2D
test_mnpe2d
M.TLp = TLp;
M.TLvz = TLvz;
M.r = r;
M.z = z;
M.surf = surf;
M.bathy = bathy;
%% Run MNPE2D
% Write temporary file to answer MNPE command line prompt for accuracy/efficiency
PromptFile='tmp.txt';fid=fopen(PromptFile,'w');fprintf(fid,'1\n');fclose(fid);
system([fullfile(BINPath,'MNPE2D') ' < ' PromptFile]);
delete(PromptFile);
doPlots = false;
scriptedinput = true;
pefile = 'press.bin'; %#ok<*NASGU>
opt1 = 2;
savedat = 'n';
peout1
peout2
F.TLp = tlpress;
pefile = 'apvz.bin';
peout1
peout2
F.TLvz = tlpress;
F.r = rng*1e3;
F.z = dep';
F.surf = surfout';
F.bathy = bathout';
%% Plot Benchmark Pressure Transmission Loss Field
figure
imagesc(B.r/1e3,B.z,B.TLp)
hold on
plot(B.r/1e3,B.surf,'w','LineWidth',2)
plot(B.r/1e3,B.bathy,'w','LineWidth',2)
hold off
colorbar('direction','reverse')
colormap(gray);
caxis(cax)
xlim([0 Rmax/1000])
ylim(yax)
xlabel('Range (km)')
ylabel('Depth (m)')
title('Benchmark Pressure Transmission Loss (dB//R0)')
%% Plot Benchmark Vertical Velocity Transmission Loss Field
figure
imagesc(B.r/1e3,B.z,B.TLvz)
hold on
plot(B.r/1e3,B.surf,'w','LineWidth',2)
plot(B.r/1e3,B.bathy,'w','LineWidth',2)
hold off
colorbar('direction','reverse')
colormap(gray);
caxis(cax)
xlim([0 Rmax/1000])
ylim(yax)
xlabel('Range (km)')
ylabel('Depth (m)')
title('Benchmark Vertical Velocity Transmission Loss (dB//R0)')
%% Plot Matlab MNPE2D Pressure Transmission Loss Field
figure
imagesc(M.r/1e3,M.z,M.TLp)
hold on
plot(M.r/1e3,M.surf,'w','LineWidth',2)
plot(M.r/1e3,M.bathy,'w','LineWidth',2)
hold off
colorbar('direction','reverse')
colormap(gray);
caxis(cax)
xlim([0 Rmax/1000])
ylim(yax)
xlabel('Range (km)')
ylabel('Depth (m)')
title('Matlab MNPE2D Pressure Transmission Loss (dB//R0)')
%% Plot Matlab MNPE2D Vertical Velocity Transmission Loss Field
figure
imagesc(M.r/1e3,M.z,M.TLvz)
hold on
plot(M.r/1e3,M.surf,'w','LineWidth',2)
plot(M.r/1e3,M.bathy,'w','LineWidth',2)
hold off
colorbar('direction','reverse')
colormap(gray);
caxis(cax)
xlim([0 Rmax/1000])
ylim(yax)
xlabel('Range (km)')
ylabel('Depth (m)')
title('Matlab MNPE2D Vertical Velocity Transmission Loss (dB//R0)')
%% Plot MNPE2D Pressure Transmission Loss Field
figure
imagesc(F.r/1e3,F.z,F.TLp)
hold on
plot(F.r/1e3,F.surf,'w','LineWidth',2)
plot(F.r/1e3,F.bathy,'w','LineWidth',2)
hold off
colorbar('direction','reverse')
colormap(gray);
caxis(cax)
xlim([0 Rmax/1000])
ylim(yax)
xlabel('Range (km)')
ylabel('Depth (m)')
title('MNPE2D Pressure Transmission Loss (dB//R0)')
%% Plot MNPE2D Vertical Velocity Transmission Loss Field
figure
imagesc(F.r/1e3,F.z,F.TLvz)
hold on
plot(F.r/1e3,F.surf,'w','LineWidth',2)
plot(F.r/1e3,F.bathy,'w','LineWidth',2)
hold off
colorbar('direction','reverse')
colormap(gray);
caxis(cax)
xlim([0 Rmax/1000])
ylim(yax)
xlabel('Range (km)')
ylabel('Depth (m)')
title('MNPE2D Vertical Velocity Transmission Loss (dB//R0)')
%% Plot Pressure Transmission Loss at 30 m
zout = 30;
B.nout = find(B.z>=zout,1);
if (zout-B.z(B.nout-1))<(B.z(B.nout)-zout)
    B.nout = B.nout-1;
end
M.nout = find(M.z>=zout,1);
if (zout-M.z(M.nout-1))<(M.z(M.nout)-zout)
    M.nout = M.nout-1;
end
F.nout = find(F.z>=zout,1);
if (zout-F.z(F.nout-1))<(F.z(F.nout)-zout)
    F.nout = F.nout-1;
end
figure
plot(B.r/1000,B.TLp(B.nout,:),'k');
hold on
plot(M.r/1000,M.TLp(M.nout,:),'b');
plot(F.r/1000,F.TLp(F.nout,:),'r');
hold off
axis ij;
grid on
xlim([0 Rmax/1000])
ylim(cax)
legend('Benchmark','Matlab MNPE2D','MNPE2D','location','northeast')
xlabel('Range (km)')
ylabel('Transmission Loss (dB//R0)')
title(['Pressure Transmission Loss at ' num2str(zout) ' m'])
%% Plot Vertical Velocity Transmission Loss at 30 m
figure
plot(B.r/1000,B.TLvz(B.nout,:),'k');
hold on
plot(M.r/1000,M.TLvz(M.nout,:),'b');
plot(F.r/1000,F.TLvz(F.nout,:),'r');
hold off
axis ij;
grid on
xlim([0 Rmax/1000])
ylim(cax)
legend('Benchmark','Matlab MNPE2D','MNPE2D','location','northeast')
xlabel('Range (km)')
ylabel('Transmission Loss (dB//R0)')
title(['Vertical Velocity Transmission Loss at ' num2str(zout) ' m'])