function [k0,dz,dr,Nz,Nr,z,kz,r,zf,kf,D] = SetupCalculation(f,rzfact,R0,Rmax,c0,bottom)
%% function [k0,dz,dr,Nz,Nr,z,kz,r,zf,kf,D] = SetupCalculation(f,rzfact,R0,Rmax,c0,bottom)
%
% Determines computational grid dimensions and sampling for M3PE 2D
% calculations.
%

%% Determine maximum depth from environment
H = max(bottom.bathy_r.Values);
%% Calculate size of computational grid
k0 = 2*pi*f/c0;                 % Reference wavenumber, 1/m
D = 4*H;                        % Initial maximum calculation depth, m
lambda = c0/f;                  % Reference wavelength, m
dz = lambda/10;                 % Depth step, m
dr = dz*rzfact;                 % Range step, m
Nz = 2^nextpow2(round(2*D/dz)); % Number of depth cells 
Nr = floor(Rmax/dr);            % Number of range cells
z = (-Nz/2:Nz/2-1)'*dz;         % Depth vector
r = R0 + (0:Nr-1)*dr;           % Range vector
D = Nz/2*dz;                    % Actual maximum calculation depth, m
dkz = pi/D;                     % Vertical wavenumber step, 1/m
kz = (-Nz/2:(Nz/2-1))'*dkz;     % Vertical wavenumber vector
kz = fftshift(kz);              % Convert to FFT order
zf = D/2;                       % Depth to begin filtering, m
kf = 3/4*k0;                    % Wavenumber to begin filtering, 1/m