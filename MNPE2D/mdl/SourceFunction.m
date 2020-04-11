function psihat = SourceFunction(SourceType,R0,k0,kz,zs)
%% function psihat = SourceFunction(SourceType,R0,k0,kz,zs)
%
% Outputs the kz-domain version of the PE starter field in fft order for a
% monopole, vertical dipole, or horizontal dipole source.
%
% Inputs:   SourceType - 1 = monopole, 2 = vertical dipole, 3 = down-range horizontal dipole
%           R0         - Reference range, m
%           k0         - Reference wavenumber, 1/m
%           kz         - Vertical wavenumber array, 1/m
%           zs         - Source depth, m
%
% Outputs:  psihat     - Starter field in kz-domain, fft order
%

%% Normalization and phase terms
Nz = length(kz);
dkz = abs(diff(kz(1:2)));
alpha = Nz*dkz*sqrt(R0/(2*pi*k0));  % Source normalization
arg1 = min(0.9999,max(-0.9999,kz/k0));
taper = (1-arg1.^2).^(-1/4);        % Wide-angle PE taper function
phase = 1;
dz = 0.1/k0;                        % Dipole spacing
Az = 1i/(2*k0*dz);                  % Dipole equivalent amplitude
%% Monopole source
if SourceType == 1
    psihat = sin(kz*zs);
end
%% Vertical dipole source
if SourceType == 2
    psihat = Az.*sin(kz*(zs+dz)) - ...
        Az.*sin(kz*(zs-dz));
end
%% Horizontal dipole source
if SourceType == 3
    C1 = -1/8*Az^4;
    C2 = -1/2*Az^2 + 1/2*Az^4;
    C3 = 1 + Az^2 - 3/4*Az^4;
    psihat = C1.*sin(kz*(zs-4*dz)) + ...
        C2.*sin(kz*(zs-2*dz)) + ...
        C3.*sin(kz*(zs)) + ...
        C2.*sin(kz*(zs+2*dz)) + ...
        C1.*sin(kz*(zs+4*dz));
end
psihat = -2i*alpha.*taper.*phase.*psihat;

    