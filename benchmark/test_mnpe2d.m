%% Setup Environment
% Create ocean and bottom profile structures
ocean.ssp_zr = scatteredInterpolant([0;Hmax;0;Hmax],[0;0;Rmax;Rmax;],c*ones(4,1),'natural','nearest');
ocean.rho = rho;
bottom.bathy_r = griddedInterpolant([0;Rmax],[H;H],'linear','nearest');
bottom.c = c1;
bottom.rho = rho1;
bottom.alpha = alpha1;
bottom.cs = 0;
bottom.alphas = 0;
% Create flat surface realization
surface = RoughSurface1D(0,Rmax,8192,H,0);
%% Run M3PE2D
[p,vz,vr,z,r,surf,bathy] = MNPE2D(SourceType,f,zs,rzfact,R0,Rmax,c0,rho0,ocean,bottom,surface);
%% Calculate Transmission Loss
TLp = -20*log10(abs(p));
TLvz = -20*log10(abs(vz*rho*c));