function FZ = DepthFilter(D,dz,zf,z)
%% function FZ = DepthFilter(D,dz,zf,z)
%
% Calculates the depth-domain filter to eliminate reflections from the
% bottom of the computational grid.
%

%% Generate depth filter
FZ = ones(size(z));
Nz = length(FZ);
FZ(z>zf) = (cos((pi/(2*(D-zf)))*(z(z>zf)-zf))).^2;
FZ(z>D-dz) = 0;
FZ(Nz/2:-1:2)=FZ(Nz/2+2:Nz);
FZ(1)=0;