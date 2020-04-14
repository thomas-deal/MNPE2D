function FK = WavenumberFilter(k0,kz,kf)
%% function FK = WavenumberFilter(k0,kz,kf)
%
% Calculates the wavenumber-domain filter to eliminate portions of the
% spectrum propagating at angles steeper than supported by the wide-angle
% PE approximation. Vector is returned in FFT order.
%

%% Generate wavenumber filter
kz = fftshift(kz);  % Convert to natural order
FK = zeros(size(kz));
FK(kz>-k0) = (cos((2*pi/k0)*(kz(kz>-k0)-kf))).^2;
FK(kz>-kf) = 1;
FK(kz>kf) = (cos((2*pi/k0)*(kz(kz>kf)-kf))).^2;
FK(kz>k0) = 0;
FK = fftshift(FK);