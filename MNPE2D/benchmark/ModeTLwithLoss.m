function [TLp, TLvz, r, z] = ModeTLwithLoss(f,zs,H,c,rho,c1,rho1,alphaf,Rmax)
%  Compute TL based on modes for Pekeris waveguide
% clear

%  define parameters
ssp(1,1)=0.;
% ssp(1,2)=input('Enter water column sound speed in m/s: ');
% sspb(1,1)=input('Enter bottom depth in m: ');
% sspb(1,2)=input('Enter bottom sound speed in m/s: ');
ssp(1,2) = c;
sspb(1,1) = H;
sspb(1,2) = c1;
ssp(2,1)=sspb(1,1); ssp(2,2)=ssp(1,2);
% rhob=input('Enter bottom density in g/cc: ');
% atten=input('Enter bottom attenuation in dB/m/kHz: ');
% z_s=input('Enter source depth in m: ');
% freq=input('Enter source frequency in Hz: ');
rhob = rho1/1e3;
atten = alphaf;
z_s = zs;
freq = f;
w=2*pi*freq;
% r_max=input('Enter max range for TL plot (km): ')*1000;
r_max = Rmax;
attene=atten/8.686*(freq./1000);

N=1024;

zbot=sspb(1,1);
nssp=size(ssp,1);
if (ssp(nssp,1)<zbot)
   ssp(nssp+1,2)=ssp(nssp,2)+(ssp(nssp,2)-ssp(nssp-1,2))*(zbot-ssp(nssp,1))/(ssp(nssp,1)-ssp(nssp-1,1));
   ssp(nssp+1,1)=zbot;
   nssp=nssp+1;
end
  
cDm=ssp(nssp,2);
cDp=sspb(1,2);
rhow=1;

dzw=zbot/(3*N/4-1);
zw=[0:dzw:zbot];
dzb=dzw*sqrt(rhob/rhow);
z_max=zbot+dzb*N/4;
zb=[zbot+dzb:dzb:z_max];
sspb(2,1)=z_max;sspb(2,2)=sspb(1,2);
  
cw=interp1q(ssp(:,1),ssp(:,2),zw');
cb=interp1q(sspb(:,1),sspb(:,2),zb');

inanw=find(isnan(cw)>0, 1);
inanb=find(isnan(cb)>0, 1);
if ~isempty(inanw) || ~isempty(inanb)
   %   c(inan)=c(inan(1)-1);
   disp('NAN found in interpolated ssp')
end


%  define matrix and compute eigenfunctions, eigenvalues, then normalize
ew=[0, dzw^(-2).*ones(1,3*N/4-2)];
eb=[dzb^(-2).*ones(1,N/4-1), 0];
e=[ew, eb];
dw=-2/dzw^2 + w^2./cw(1:3*N/4-1).^2;
dw(3*N/4)=-rhow*(1/rhow+1/rhob)/dzw^2 + w^2*(1/cDm^2+1/cDp^2)/2;
db=-2/dzb^2 + w^2./cb.^2;
d=[dw', db'];
A=diag(d) + diag(e,1) + diag(e,-1);

[psi, K_squ]=eig(A);

[K_squ, norder]=sort(diag(K_squ)');
K_squ=fliplr(K_squ); norder=fliplr(norder);
K=sqrt(K_squ);

norm=dzw/rhow*sum([psi(1,:).^2./2; psi(2:3*N/4-1,:).^2; psi(3*N/4,:).^2./2]);
norm=norm+dzb/rhob*sum([psi(3*N/4,:).^2./2; psi(3*N/4+1:N-1,:).^2; psi(N,:).^2./2]);
norm=sqrt(1./norm);

nmodes=floor(2*zbot*freq*sqrt(cb(1).^2-cw(1).^2)/(cw(1).*cb(1)) + 0.5);
disp(['Number of propagating modes = ' num2str(nmodes) newline]);

psit=zeros(size(psi,1),nmodes);
for n=1:nmodes
  psit(:,n)=norm(norder(n))*psi(:,norder(n));
  if (psit(2,n)<0)
    psit(:,n)=-psit(:,n);
  end
end
psi=psit;clear psit
psidel=zeros(nmodes,1);
for n=1:nmodes
  psidel(n)=dzb*sum([psi(3*N/4,n).^2./2; psi(3*N/4+1:N-1,n).^2; psi(N,n).^2./2]);
  psidel(n)=psidel(n)*2*pi*freq*attene/cb(1)/K(n);
end

z=[zw, zb];

%  plot first 10 modes
% figure
% for n=1:min(nmodes,5)
%   subplot(1,5,n), plot(psi(:,n),z); axis('ij'); axis([min(psi(:,n)) max(psi(:,n)) 0 z_max]);
%   if n ~= 1
%     set(gca,'ytick',[1000],'ygrid','on')
%   end
%   if n == 1
%     ylabel('Depth (m)')
%   end
%   xlabel(['K=',num2str(K(n)),', del=',num2str(psidel(n))])
%   title(['Mode #',num2str(n)])
% end
% if nmodes>5
%     figure
%     for n=6:min(nmodes,10)
%       subplot(1,5,n-5), plot(psi(:,n),z); axis('ij'); axis([min(psi(:,n)) max(psi(:,n)) 0 z_max]);
%       if n ~= 6
%         set(gca,'ytick',[1000],'ygrid','on')
%       end
%       if n == 6
%         ylabel('Depth (m)')
%       end
%       xlabel(['K=',num2str(K(n)),', del=',num2str(psidel(n))])
%       title(['Mode #',num2str(n)])
%     end
% end

psi_s=zeros(1,nmodes);
for n=1:nmodes
  psi_s(n)=interp1(z,psi(:,n),z_s,'pchip');
end

%  Plot TL field
Nr=1001;
dr=r_max/(Nr-1);
rngm=[0:dr:r_max];

TL=zeros(size(psi,1),Nr);
TLvz = TL;
for m=1:Nr
  press=zeros(N,1);
  for n=1:nmodes
    press=press+(psi_s(1,n).*psi(:,n).*exp(1i*K(n)*rngm(m)).*exp(-psidel(n)*rngm(m))./sqrt(K(n)));
  end
  TL(:,m)=-20.*log10(max(10^(-10),abs(1/rhow.*sqrt(2*pi/max(1,rngm(m))).*press)));
  vz = 1/(1i*2*pi*f*rho)*diff(press)./diff(z')*rho*c;
  vz = [0; vz];
  TLvz(:,m)=-20.*log10(max(10^(-10),abs(1/rhow.*sqrt(2*pi/max(1,rngm(m))).*vz)));
end

TLp = TL;
r = rngm;
z = z';

% figure
% imagesc(rngm/1000.,z,TL);colormap(flipud(jet));set(colorbar,'YDir','reverse');
% caxis([min(min(TL)) min(min(TL))+80])
% hold on;plot([0 r_max],[zbot zbot]);
% %v=axis;v(4)=4000;axis(v);
% xlabel('Range (km)'); ylabel('Depth (m)');
% title('Transmission Loss');
% 
% disp(' ');
% modefile=input('Enter name of mode file to save: ','s');
% psi=psi(:,1:nmodes);K=K(1:nmodes);Nm=nmodes;
% eval(['save ' modefile ' psi K psidel z zw dzw zb dzb cw cb rhow rhob freq Nm'])
