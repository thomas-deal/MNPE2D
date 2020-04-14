if ~exist('scriptedinput','var')
    clear pefile;
end

if ~exist('pefile','var')
    pefile=input('Enter name of binary MMPE output file: ','s');
else
    disp(pefile)
end
disp(' ');
fid=fopen(pefile);

% Read in first record of data
nrecs=fread(fid,1,'int32');
c0=fread(fid,1,'float32');
nf=fread(fid,1,'int32');cfreq=fread(fid,1,'float32');freqbw=fread(fid,1,'float32');
irng=1;
nrout=fread(fid,1,'int32');rngmin=fread(fid,1,'float32')/1000;
rngmax=fread(fid,1,'float32')/1000;dr=fread(fid,1,'float32')/1000;
nzout=fread(fid,1,'int32');depmin=fread(fid,1,'float32');depmax=fread(fid,1,'float32');
nrad=fread(fid,1,'int32');drad=fread(fid,1,'float32');
bdint(irng)=fread(fid,1,'float32');dbdint(irng)=fread(fid,1,'float32'); %#ok<*NASGU>
sd=fread(fid,1,'float32');itype=fread(fid,1,'int32');
% Read in bathymetry
fseek(fid,(4*(2*nzout*(1+(1+nrout*nrad)*nf)-18)),0);
for irad=1:nrad
 bath1(:,irad)=fread(fid,nrout+1,'float32'); %#ok<*SAGROW>
 dbath1(:,irad)=fread(fid,nrout+1,'float32');
 surf1(:,irad)=fread(fid,nrout+1,'float32');
end
fclose(fid);
dep=(depmin:(depmax-depmin)/(nzout-1):depmax);
rng0=0.;
% NOTE: If rngmin = 0, then first output (after starting field) is
% actually at rngmin = dr (assuming ranges between successive outputs > dr).
if nrout > 1
  rng=(rngmin:(rngmax-rngmin)/(nrout-1):rngmax);
else
  rng=rngmax;
end
if rngmin == 0 && nrout > 1
  rng(1)=dr;
end
if nrout > 1
  delr=rng(size(rng,2))-rng(size(rng,2)-1);
else
  delr=rng(1);
end
if nrad > 1
  radmin=-(nrad/2-1)*drad;
  radmax=(nrad/2)*drad;
  rad=radmin:drad:radmax;
else
  rad=0;
end
if nf > 1
  dfreq=freqbw/(nf-1);
  freqmin=cfreq-(nf/2)*dfreq;
  freqmax=cfreq+(nf/2-1)*dfreq;
  freq=freqmin:dfreq:freqmax;
  dtime=1/freqbw;
  timmin=-(nf/2-.5)*dtime;
  timmax=(nf/2-.5)*dtime;
  time=timmin:dtime:timmax;
else
  freq=cfreq;
  time=0;
end
bath0=bath1(1,:); dbath0=dbath1(1,:); surf0=surf1(1,:);
bath=bath1(2:nrout+1,:); dbath=dbath1(2:nrout+1,:); surf=surf1(2:nrout+1,:);
clear bath1 dbath1 surf1;

% Output parameter info for this file.
if itype==1
    disp('This file contains scalar pressure data with the following parameters:');
elseif itype==2
    disp('This file contains radial component particle velocity data with the following parameters:');
elseif itype==3
    disp('This file contains vertical component particle velocity data with the following parameters:');
end
if nf > 1
  disp(['   ' num2str(nf) ' frequencies over ' num2str(freqbw) ' Hz centered at ' num2str(cfreq) ' Hz']);
else
  disp(['   frequency = ' num2str(cfreq) ' Hz for a source at depth ' num2str(sd) ' meters']);
end
if nrad > 1
  disp(['   ' num2str(nrad) ' radials over ' num2str((nrad-1)*drad) ' degrees azimuth']);
else
  disp('   single radial');
end
disp(['   ' num2str(nzout) ' points in depth from ' num2str(depmin) ' m to ' num2str(depmax) ' m']);
disp(['   ' num2str(nrout) ' points in range from ' num2str(rngmin) ' km to ' num2str(rngmax) ' km']);
disp(' ');

header=2*nzout;
recl=2*nzout;

clear ans bdint dbdint irad irng;
