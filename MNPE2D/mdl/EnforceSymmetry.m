function x = EnforceSymmetry(z,eta,x,odd,shift)
%% function psi = EnforceSymmetry(z,eta,psi,odd,shift)
%
% Applies odd or even symmetry to vector accounting for the displacement of 
% the rough surface.
%

%% Apply symmetry
Nz = length(z);
if shift
    x = fftshift(x);
end
if odd
    s = -1;
else
    s = 1;
end
ns = find(z<=eta,1,'last');
iz = ns+1;
x(ns) = 0;
for i=(ns-1):-1:1
    x(i) = s*x(iz);
    iz = min(iz+1,Nz); 
end
if shift
    x = fftshift(x);
end
    