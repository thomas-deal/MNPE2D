function [p,vz,vr,z,r,surf,bathy] = MNPE2D(SourceType,f,zs,rzfact,R0,Rmax,c0,rho0,ocean,bottom,surface)
%% function [p,vz,vr,z,r,surf,bathy] = MNPE2D(SourceType,f,zs,rzfact,R0,Rmax,c0,rho0,ocean,bottom,surface)
%
% Matlab Monterey-Newport Parabolic Equation 2 Dimensional (MNPE2D) 
% split-step Fourier parabolic equation propagation model. Outputs pressure
% and velocity components for monopole, vertical dipole, and horizontal 
% dipole sources. Currently supports a single bottom interface (no deep 
% bottom).
%

%% Calculation Setup
[k0,dz,dr,Nz,Nr,z,kz,r,zf,kf,D] = SetupCalculation(f,rzfact,R0,Rmax,c0,bottom);
p = zeros(Nz,Nr);
vz = zeros(Nz,Nr);
vr = zeros(Nz,Nr);
surf = zeros(1,Nr);
bathy = zeros(1,Nr);
%% Wavenumber Domain Filter
FK = WavenumberFilter(k0,kz,kf);
%% Depth Filter
FZ = DepthFilter(D,dz,zf,z);
%% Wavenumber Propagator
[PROPK,Top] = WavenumberPropagator(k0,kz,dr);
%% Lens Propagator
ir = 1;
[PROPZ,n2,~,bathy(ir),surf(ir),detadr] = LensPropagator(ir,f,c0,dz,dr,z,r,ocean,bottom,surface);
%% Initialize psi
psi = SourceFunction(SourceType,R0,k0,kz,zs);
psi = psi/sqrt(rho0);     % Field transformation for density term
psi = FK.*psi;
psir = psi;
psi = ifftshift(ifft(psi));
psiz = 1i*(kz/k0).*psir;
psiz = FK.*psiz;
psiz = ifftshift(ifft(psiz));
psir = Top.*psir;
psir = FK.*psir;
psir = ifftshift(ifft(psir));
psir = (sqrt(n2)+1i./(2*k0*r(ir))).*psi - psir;
p(:,ir) = sqrt(R0/r(ir))*exp(1i*k0*r(ir))*psi;               % Pressure field
vz(:,ir) = 1/(c0*rho0)*sqrt(R0/r(ir))*exp(1i*k0*r(ir))*psiz; % Vertical velocity field
vr(:,ir) = 1/(c0*rho0)*sqrt(R0/r(ir))*exp(1i*k0*r(ir))*psir; % Radial velocity field
%% Propagate
while ir<Nr
    % Step Solution Forward One Range Step
    psi = PROPZ.*psi;                       % Apply first half-step lens propagator
    psi = FZ.*psi;                          % Apply depth filter  
    psi = fft(psi);                         % Transform to wavenumber domain
    [PROPKIO,~] = WavenumberPropagator(k0,kz,dr,detadr);
    psiio = PROPKIO.*psi;                   % Apply image ocean wavenumber propagator
    psiio = FK.*psiio;                      % Apply wavenumber filter
    psiio = ifft(psiio);                    % Transform to depth domain
    psi = PROPK.*psi;                       % Apply wavenumber propagator
    psi = FK.*psi;                          % Apply wavenumber filter
    psi = ifft(psi);                        % Transform to depth domain
    psi(z<surf(ir)) = psiio(z<surf(ir));    % Copy image ocean to psi
    ir = ir + 1;                            % Increment range step
    [PROPZ,n2,rhoz,bathy(ir),surf(ir),detadr] = LensPropagator(ir,f,c0,dz,dr,z,r,ocean,bottom,surface);                         
                                            % Calculate second half-step lens propagator
    psi = PROPZ.*psi;                       % Apply second half-step lens propagator
    psi = FZ.*psi;                          % Apply depth filter   
    % Calculate Velocity Potentials
    psir = fft(psi);
    psiz = (kz/k0).*psir;
    psiz = FK.*psiz;
    psiz = ifft(psiz);
    psir = Top.*psir;
    psir = FK.*psir;
    psir = ifft(psir);
    psir = (sqrt(n2)+1i./(2*k0*r(ir))).*psi - psir;
    % Calculate Pressure and Velocity
    p(:,ir) = sqrt(R0/r(ir))*exp(1i*k0*r(ir))*(psi.*sqrt(rhoz));               % Pressure field
    vz(:,ir) = 1/(c0*rho0)*sqrt(R0/r(ir))*exp(1i*k0*r(ir))*(psiz.*sqrt(rhoz)); % Vertical velocity field
    vr(:,ir) = 1/(c0*rho0)*sqrt(R0/r(ir))*exp(1i*k0*r(ir))*(psir.*sqrt(rhoz)); % Radial velocity field
end
%% Eliminate Image Ocean
for ir=1:Nr
    p(z<surf(ir),ir) = 0;
    vr(z<surf(ir),ir) = 0;
    vz(z<surf(ir),ir) = 0;
end
zmin = min(surf);
izmin = find(z<zmin,1,'last');
z = z(izmin:end);
p = p(izmin:end,:);
vz = vz(izmin:end,:);
vr = vr(izmin:end,:);