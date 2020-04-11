function [PROPZ,n2,rhoz,H,eta,detadr] = LensPropagator(ir,f,c0,dz,dr,z,r,ocean,bottom,surface)
%% function [PROPZ,n2,rhoz,H,eta,detadr] = LensPropagator(ir,f,c0,dz,dr,z,r,ocean,bottom,surface)
%
% Produces the lens propagator, effective sound speed ratio, and 
% interpolated depth at current range step for M3PE 2D calculations. 
% Vectors are in FFT order.
%

%% Calculate derived parameters
Nz = length(z);
lambda = c0/f;
k0 = 2*pi*f/c0;
%% Interpolate ocean sound speed profile
zreal = z(Nz/2+1:Nz);
cw = ocean.ssp_zr(zreal,r(ir)*ones(Nz/2,1));
rhow = ocean.rho;
% Attenuation 
attnconv=1.094/(1000.0*8.686);
fsq = (f/1e3)^2;
attnw = 0.003+0.1*fsq/(1+fsq)+40*fsq/(4100+fsq)+2.75e-4*fsq;
attnw = attnconv*attnw;
attnw = attnw*(1-6.46e-5*zreal);
%% Interpolate bottom properties and bathymetry
cb = bottom.c*ones(Nz/2,1);
alphab = bottom.alpha*ones(Nz/2,1);
attnb = alphab.*f./cb/(20*log10(exp(1)));
% Interpolate bottom depth at this range
H = bottom.bathy_r(r(ir));
nb = floor(H/dz)+1;
% Convert shear wave effect to equivalent complex density
c1 = ocean.ssp_zr(H,r(ir));
if (bottom.cs>0)&&(bottom.cs<c1)
    c2 = bottom.c;
    rho2 = bottom.rho;
    a2 = bottom.alpha;
    cs = bottom.cs;
    as = bottom.alphas;  
    w = 2*pi*f;
    rhob = rho2*((1 - 2/(c1/cs + 1i*as*c1/w)^2)^2 + ...
        4i*(1 - (c1/c2 + 1i*a2*c1/w)^2)^(1/2) * ...
        ((c1/cs + 1i*as*c1/w)^2-1)^(1/2) / ...
        (c1/cs + 1i*as*c1/w)^4);
else
    rhob = bottom.rho;
end
rhob = rhob*ones(Nz/2,1);
%% Interpolate surface between defined ranges
surface.f_r.Values = surface.eta;
eta = surface.f_r(r(ir));
surface.f_r.Values = surface.detadr;
detadr = surface.f_r(r(ir));
%% Mix ocean and bottom profiles
zeta = z(Nz/2+1:Nz)-H;
[Hc,Hr,Hrp,Hrpp] = MixingFunction(zeta,lambda,Nz,dz);  
cw(nb:end) = cw(nb-1);
cb(1:(nb-1)) = cb(nb);
cz = cw + (cb-cw).*Hc;
attn = attnw + (attnb-attnw).*Hc;
rhoz = rhow + (rhob-rhow).*Hr;
rhozp = (rhob-rhow).*Hrp;
rhozpp = (rhob-rhow).*Hrpp;
%% Create Image Ocean
cz = [flipud(cz); cz];
attn = [flipud(attn); attn];
rhoz = [flipud(rhoz); rhoz];
rhozp = [flipud(rhozp); rhozp];
rhozpp = [flipud(rhozpp); rhozpp];
%% Generate effective sound speed ratio, potential energy operator, and lens propagator
n2 = (c0./cz).^2 + 1/(2*k0^2)*(1./rhoz.*rhozpp - 3/2*(1./rhoz.*rhozp).^2); 
n2 = EnforceSymmetry(z,eta,n2,0,0);
Uop = 1-sqrt(n2);
attn = EnforceSymmetry(z,eta,attn,0,0);
PROPZ = exp(-1i*k0*dr/2*Uop).*exp(-attn*dr/2);