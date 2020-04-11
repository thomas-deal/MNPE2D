function [PROPK,Top] = WavenumberPropagator(k0,kz,dr,varargin)
%% function [PROPK,Top] = WavenumberPropagator(k0,kz,dr)
% function [PROPK,Top] = WavenumberPropagator(k0,kz,dr,detadr)
%
% Produces the wavenumber propagator and kinetic energy operator for M3PE
% 2D calculations. When argument detadr is omitted, outputs are for the
% real ocean. When it is included, outputs are for the image ocean.
% Outputs are in FFT order.

%% Generate kinetic energy operator and wavenumber propagator
Top = 1-sqrt(1-(kz/k0).^2);
if nargin > 3
    detadr = varargin{1};
    Top = Top + 2*detadr*kz/k0;
end
PROPK = exp(-1i*k0*dr*Top);