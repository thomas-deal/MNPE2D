function [Hc,Hr,Hrp,Hrpp] = MixingFunction(zeta,lambda,Nz,dz)
% function [Hc,Hr,Hrp,Hrpp] = MixingFunction(zeta,lambda,Nz,dz)
%
% Calculates the mixing functions for sound speed and density
% discontinuities at the water/bottom interface for M3PE 2D calculations.
%

%% Sound speed mixing function
Lc = lambda/10; % Sound speed mixing length, m
Lc = max(Lc,dz);
Hc = 1/2*(1+tanh(zeta/(2*Lc)));
%% Density mixing functions
Lr = 2*lambda;  % Density mixing length, m 
Lr = max(Lr,2*dz);
Hr = zeros(Nz/2,1);
Hrp = zeros(Nz/2,1);
Hrpp = zeros(Nz/2,1);
for iz = 1:Nz/2
    zz = zeta(iz);
    if (zz <= -Lr)
        Hr(iz) = 0;
        Hrp(iz) = 0;
        Hrpp(iz) = 0;
    elseif (zz <= -Lr/2)
        Hr(iz) = 2/3*(1+zz/Lr)^3;
        Hrp(iz) = 2/Lr*(1+zz/Lr)^2;
        Hrpp(iz) = 4/Lr^2*(1+zz/Lr);
    elseif (zz <= Lr/2)
        Hr(iz) = 1/2+zz/Lr-2/3*(zz/Lr)^3;
        Hrp(iz) = 1/Lr*(1-2*(zz/Lr)^2);
        Hrpp(iz) = -4/Lr^2*(zz/Lr);
    elseif (zz <= Lr)
        Hr(iz) = 1-2/3*(1-zz/Lr)^3;
        Hrp(iz) = 2/Lr*(1-zz/Lr)^2;
        Hrpp(iz) = -4/Lr^2*(1-zz/Lr);
    else
        Hr(iz) = 1;
        Hrp(iz) = 0;
        Hrpp(iz) = 0;
    end
end