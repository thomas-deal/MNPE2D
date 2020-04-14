clear psi0 psi psid press pressd tlpsi0 tlpress tlpressd gaps tmp r0 k0;
clear i nradout nrngout irad irng bathout dbathout;

if ~exist('scriptedinput','var')
    clear opt1 opt2 optt opt1t optt1 optout savedat freqout 
    clear radout rngout fileout depout depmin depmax elmin elmax
end
if ~exist('doPlots','var')
    doPlots = true;
end

% =============================================================================

disp(' ');
if ~exist('pefile','var')
  disp('Must first initialize with "peout1" command.')
else

% Offer options for calculations.
disp('Choose one of the following options.');
disp('(NOTE: All pressure levels relative to 0 dB source level.');

disp('1)  Output starting field data;');
disp('2)  Compute data for single radial;');
disp('3)  Compute data for single range;');
disp('4)  Compute data for single depth;');
disp('5)  Compute data for single interface;');
disp('6)  Compute travel time data;');
if ~exist('opt1','var')
    opt1=input(' ');
else
    disp(opt1)
end
disp(' ');

if itype==1
    rhoc0=1.;
else
    rhoc0=1000.*c0;  % Note, all values normalized such that P0/(rho*c0)=1.
end
% =============================================================================
% OUTPUT STARTING FIELD DATA.
if opt1 == 1

% Skip past header only
  skip=0;
  fid=fopen(pefile,'r');
  fseek(fid,(header+skip)*4,0);

% Skip consists of size of data set for each freq except starting field, i.e.
% (number of radials)*(number of ranges)*(number of depths)
  skip=recl*nrad*nrout;
  for ifreq=1:nf
    data=fread(fid,2*nzout,'float32');
    if nf > 1 && ifreq < nf
      fseek(fid,(skip)*4,0);
    end
    psi0(:,ifreq)=(data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1))/rhoc0; %#ok<*SAGROW>
    clear data;
  end

  fclose(fid);

% Fill in any NaN gaps
  gaps=find(isnan(psi0));
  for i=1:size(gaps)
    psi0(gaps(i))=1.e-20;
  end  % end of 'for i=1:size(gaps)'

% Re-order freqs
  if nf > 1
    for iz=1:nzout
      psi0(iz,:)=fftshift(psi0(iz,:));
      if dep(iz)<surf0
          psi0(iz,:)=1.e-20;
      end      
    end
  end

  tlpsi0=-20*log10(max(abs(psi0),1.e-20));
  tlmin=0;tlmax=tlmin+100;

  if doPlots
      if nf == 1
        figure;plot(tlpsi0,dep);axis('ij');
        if itype==1
          xlabel('Source Level (dB re 1m) for Pressure');ylabel('Depth (m)');
        elseif itype==2
          xlabel('Source Level (dB re 1m) for Velocity (Radial)');ylabel('Depth (m)');
        else
          xlabel('Source Level (dB re 1m) for Vertical Velocity');ylabel('Depth (m)');
        end
      else
        figure;imagesc(freq,dep,tlpsi0);axis('ij');caxis([tlmin tlmax]);colormap(flipud(jet));
        colorbar;set(colorbar,'YDir','reverse');
        if itype==1
          title('Source Level (dB re 1m) for Pressure');xlabel('Freq (Hz)');ylabel('Depth (m)');
        elseif itype==2
          title('Source Level (dB re 1m) for Velocity (Radial)');xlabel('Freq (Hz)');ylabel('Depth (m)');
        else
          title('Source Level (dB re 1m) for Vertical Velocity');xlabel('Freq (Hz)');ylabel('Depth (m)');
        end
      end
  end

  if ~exist('savedat','var')
      savedat=input('Save data? (y or n) ','s');
  else
      disp(savedat)
  end
  disp(' ');

  if isempty(savedat)
  elseif savedat == 'y' || savedat == 'Y'
    if ~exist('fileout','var')
        fileout=input('Enter matlab output filename (no extension needed): ','s');
    else
        disp(fileout)
    end
    if itype==1
      eval(['save ' fileout ' psi0 tlpsi0 dep freq rng0 bath0 dbath0;']);
    elseif itype==2
      psivr0=psi0/rhoc0; tlpsivr0=tlpsi0;
      eval(['save ' fileout ' psivr0 tlpsivr0 dep freq rng0 bath0 dbath0 rhoc0;']);
    else
      psivz0=psi0/rhoc0; tlpsivz0=tlpsi0;
      eval(['save ' fileout ' psivz0 tlpsivz0 dep freq rng0 bath0 dbath0 rhoc0;']);
    end
  end  % end of 'if savedat == 'y' | savedat == 'Y''

end  % end of 'if opt1 == 1'

% =============================================================================
% COMPUTE DATA ALONG SINGLE RADIAL.
if opt1 == 2

  if nrad > 1
    disp(['Enter radial angle, ' num2str(rad(1)) ' deg to ' num2str(rad(nrad)) ' deg, to extract data:']);
    if ~exist('radout','var')
        radout=input(' ');
    else
        disp(radout)
    end
    nradout=find((rad-radout)>=0);
    if isempty(nradout)
      if radout<0
        nradout=1;
      else
        nradout=nrad;
      end
    else
      nradout=nradout(1);
    end
    radout=rad(nradout);
    disp(['Outputting data for radial ' num2str(rad(nradout)) 'deg.']);
    if nradout < nrad/2
      nradout=nradout+nrad/2+1;
    else
      nradout=nradout-nrad/2+1;
    end
    disp(' ');
  else
    radout=0;
    nradout=1;
  end  % end of 'if nrad > 1'

  if nf > 1 && nrout > 1
    disp(['Enter freq, ' num2str(freq(1)) ' Hz to ' num2str(freq(nf)) ' Hz, to extract data:']);
    if ~exist('freqout','var')
        freqout=input(' ');
    else
        disp(freqout)
    end
    nfout=find((freq-freqout)>=0);
    if isempty(nfout)
      if freqout<cfreq
        nfout=1;
      else
        nfout=nf;
      end
    else
      nfout=nfout(1);
    end
    freqout=freq(nfout);
    disp(['Outputting data for freq ' num2str(freq(nfout)) 'Hz.']);
    if nfout <= nf/2
      nfout=nfout+nf/2;
    else
      nfout=nfout-nf/2;
    end
    disp(' ');
  elseif nrout > 1
    nfout=1;
    freqout=cfreq;
  else
    nfout=nf;
    freqout=freq;
    rngout=rng(nrout);
  end  % end of 'if nf > 1 & nrout > 1'

% Skip consists of (out to proper freq)+(starting field)+(out to proper radial)
  if nrout > 1
    skip=(nfout-1)*(recl + recl*nrad*nrout) + recl + recl*(nradout-1);
  else
    skip=recl + recl*(nradout-1);
  end
  fid=fopen(pefile,'r');
  fseek(fid,(header+skip)*4,0);

  if nrout > 1
    skip=recl*(nrad-1);
    for irng=1:nrout
      data=fread(fid,2*nzout,'float32');
      if irng < nrout && nrad > 1
        fseek(fid,(skip)*4,0);
      end
      psi(:,irng)=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
      for iz=1:nzout
        if dep(iz)<surf(irng)
          psi(iz,irng)=1.e-20;
        end
      end      
      clear data;
    end
  else
    skip=recl*nrad*nrout;
    for ifreq=1:nf
      data=fread(fid,2*nzout,'float32');
      if ifreq < nf
        fseek(fid,(skip)*4,0);
      end
      psi(:,ifreq)=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
      for iz=1:nzout
        if dep(iz)<surf(1)
          psi(iz,ifreq)=1.e-20;
        end
      end      
      clear data;
    end
  end

  fclose(fid);

% Fill in any NaN gaps
  gaps=find(isnan(psi));
  for i=1:size(gaps)
    psi(gaps(i))=1.e-20;
  end  % end of 'for i=1:size(gaps)'

% Include cylindrical spreading and phase factors
  r0=1.;
  if nrout > 1
    k0=2*pi*freqout/c0;
    for irng=1:nrout
      fact=sqrt(r0/max(1000.*rng(irng),r0))*exp(1i*k0*1000.*rng(irng));
      press(:,irng)=fact*psi(:,irng);
    end
  else
% Re-order frequencies
    for iz=1:nzout
      psi(iz,:)=fftshift(psi(iz,:));
    end
    for ifreq=1:nf
      k0=2*pi*freq(ifreq)/c0;
      fact=sqrt(r0/max(1000.*rng(nrout),r0))*exp(1i*k0*1000.*rng(nrout));
      press(:,ifreq)=fact*psi(:,ifreq);
    end
  end

  tlpress=-20*log10(max(abs(press),1.e-20));
  tlmin=min(min(tlpress));tlmax=tlmin+100;
  bathout=bath(:,nradout);dbathout=dbath(:,nradout);surfout=surf(:,nradout);

  if doPlots
      if nrout > 1
        figure;imagesc(rng,dep,tlpress);axis('ij');caxis([tlmin tlmax]);colormap(flipud(jet));
        hold on;plot(rng,bathout,'w');plot(rng,dbathout,'w');plot(rng,surfout,'w');hold off
        colorbar;set(colorbar,'YDir','reverse');
        if itype==1
          title('Transmission Loss (dB re 1m) for Pressure');xlabel('Range (km)');ylabel('Depth (m)');
        elseif itype==2
          title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Range (km)');ylabel('Depth (m)');
        else
          title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Range (km)');ylabel('Depth (m)');
        end
        hold off;
      else
        figure;imagesc(freq,dep,tlpress);axis('ij');caxis([tlmin tlmax]);colormap(flipud(jet));
        hold on;plot(rng,bathout,'w');plot(rng,dbathout,'w');plot(rng,surfout,'w');hold off
        colorbar;set(colorbar,'YDir','reverse');
        if itype==1
          title('Transmission Loss (dB re 1m) for Pressure');xlabel('Freq (Hz)');ylabel('Depth (m)');
        elseif itype==2
          title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Freq (Hz)');ylabel('Depth (m)');
        else
          title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Freq (Hz)');ylabel('Depth (m)');
        end
      end
  end

  if ~exist('savedat','var')
      savedat=input('Save data? (y or n) ','s');
  else
      disp(savedat)
  end
  disp(' ');

  if isempty(savedat)
  elseif savedat == 'y' || savedat == 'Y'
    if ~exist('fileout','var')
        fileout=input('Enter matlab output filename (no extension needed): ','s');
    else
        disp(fileout)
    end
    if itype==1
      eval(['save ' fileout ' press tlpress dep rng radout freqout bathout dbathout;']);
    elseif itype==2
      apvr=press/rhoc0; tlapvr=tlpress;
      eval(['save ' fileout ' apvr tlapvr dep rng radout freqout bathout dbathout rhoc0;']);
    else
      apvz=press/rhoc0; tlapvz=tlpress;
      eval(['save ' fileout ' apvz tlapvz dep rng radout freqout bathout dbathout rhoc0;']);
    end
  end  % end of 'if savedat == 'y' || savedat == 'Y''

end  % end of 'if opt1 == 2'

% =============================================================================
% COMPUTE DATA AT SINGLE RANGE.
if opt1 == 3

  if nrout > 1
    disp(['Enter range, ' num2str(rng(1)) ' km to ' num2str(rng(nrout)) ' km, to extract data:']);
    if ~exist('rngout','var')
        rngout=input(' ');
    else
        disp(rngout)
    end
    nrngout=find((rng-rngout)>=0);
    if isempty(nrngout)
      if rngout<rng(1)
        nrngout=1;
      else
        nrngout=nrout;
      end
    else
      nrngout=nrngout(1);
    end
    rngout=rng(nrngout);
    disp(['Outputting data for range ' num2str(rng(nrngout)) 'km.']);
    disp(' ');
  else
    nrngout=1;
    rngout=rng(1);
  end  % end of 'if nrout > 1'
  irng=nrngout;
  
  if nf > 1 && nrad > 1

   disp(' ');
   disp('1)  Output data for single frequency;');
   disp('2)  Output data for single radial;');
   if ~exist('opt1t','var')
       opt1t=input(' ');
   else
       disp(opt1t)
   end
   disp(' ');
   
   if opt1t == 1

    disp(['Enter freq, ' num2str(freq(1)) ' Hz to ' num2str(freq(nf)) ' Hz, to extract data:']);
    if ~exist('freqout','var')
        freqout=input(' ');
    else
        disp(freqout)
    end
    nfout=find((freq-freqout)>=0);
    if isempty(nfout)
      if freqout<cfreq
        nfout=1;
      else
        nfout=nf;
      end
    else
      nfout=nfout(1);
    end
    freqout=freq(nfout);
    disp(['Outputting data for freq ' num2str(freq(nfout)) 'Hz.']);
    if nfout <= nf/2
      nfout=nfout+nf/2;
    else
      nfout=nfout-nf/2;
    end
    disp(' ');
    
   else
   
    disp(['Enter radial angle, ' num2str(rad(1)) ' deg to ' num2str(rad(nrad)) ' deg, to extract data:']);
    if ~exist('radout','var')
        radout=input(' ');
    else
        disp(radout)
    end
    nradout=find((rad-radout)>=0);
    if isempty(nradout)
      if radout<0
        nradout=1;
      else
        nradout=nrad;
      end
    else
      nradout=nradout(1);
    end
    radout=rad(nradout);
    disp(['Outputting data for radial ' num2str(rad(nradout)) 'deg.']);
    if nradout < nrad/2
      nradout=nradout+nrad/2+1;
    else
      nradout=nradout-nrad/2+1;
    end
    disp(' ');
   
   end  % end of 'if opt1t == 1'

  elseif nrad > 1
    nfout=1;
    freqout=cfreq;
    opt1t=1;
  else
    nradout=1;
    freqout=freq;
    radout=rad(nrad);
    opt1t=2;
  end  % end of 'if nf > 1 && nrad > 1'


  if opt1t == 1
  
% Skip consists of (out to proper freq)+(starting field)+(out to proper range)
    skip=(nfout-1)*(recl + recl*nrad*nrout) + recl + recl*nrad*(nrngout-1);
    fid=fopen(pefile,'r');
    fseek(fid,(header+skip)*4,0);

    skip=0; % don't need to use it since reading consecutive radials
    for irad=1:nrad
      data=fread(fid,2*nzout,'float32');
      psi(:,irad)=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
      for iz=1:nzout
        if dep(iz)<surf(irng)
          psi(iz,irad)=1.e-20;
        end
      end      
      clear data;
    end

  else % opt1t=2
  
% Skip consists of (starting field)+(out to proper range)+(out to proper radial)
    skip=recl+recl*nrad*(nrngout-1)+recl*(nradout-1);
    fid=fopen(pefile,'r');
    fseek(fid,(header+skip)*4,0);
    
    skip=recl*nrad*nrout;
    for ifreq=1:nf
      data=fread(fid,2*nzout,'float32');
      if ifreq < nf
        fseek(fid,(skip)*4,0);
      end
      psi(:,ifreq)=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
      for iz=1:nzout
        if dep(iz)<surf(irng)
          psi(iz,ifreq)=1.e-20;
        end
      end      
      clear data;
    end
    
  end

  fclose(fid);

% Fill in any NaN gaps
  gaps=find(isnan(psi));
  for i=1:size(gaps)
    psi(gaps(i))=1.e-20;
  end  % end of 'for i=1:size(gaps)'

% Include cylindrical spreading and phase factors
  r0=1.;
  if opt1t == 1
    k0=2*pi*freqout/c0;
    fact=sqrt(r0/max(1000.*rngout,r0))*exp(1i*k0*1000.*rngout);
    press=fact*psi;
    tmp=press(:,nrad/2+2:nrad);
    press(:,nrad/2:nrad)=press(:,1:nrad/2+1);press(:,1:nrad/2-1)=tmp;
  else
% Re-order frequencies
    for iz=1:nzout
      psi(iz,:)=fftshift(psi(iz,:));
    end
    for ifreq=1:nf
      k0=2*pi*freq(ifreq)/c0;
      fact=sqrt(r0/max(1000.*rngout,r0))*exp(1i*k0*1000.*rngout);
      press(:,ifreq)=fact*psi(:,ifreq);
    end
  end

  tlpress=-20*log10(max(abs(press),1.e-20));
  tlmin=min(min(tlpress));tlmax=tlmin+100;
  bathout=bath(nrngout,:);dbathout=dbath(nrngout,:);
  if nrad > 1
    tmp=bathout(nrad/2+2:nrad);
    bathout(nrad/2:nrad)=bathout(1:nrad/2+1);bathout(1:nrad/2-1)=tmp;
    tmp=dbathout(nrad/2+2:nrad);
    dbathout(nrad/2:nrad)=dbathout(1:nrad/2+1);dbathout(1:nrad/2-1)=tmp;
  end

  if doPlots
      if opt1t == 1
        figure;imagesc(rad,dep,tlpress);axis('ij');caxis([tlmin tlmax]);colormap(flipud(jet));
        hold on;plot(rad,bathout,'w');plot(rad,dbathout,'w');
        colorbar;set(colorbar,'YDir','reverse');
        if itype==1
          title('Transmission Loss (dB re 1m) for Pressure');xlabel('Radial (deg)');ylabel('Depth (m)');
        elseif itype==2
          title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Radial (deg)');ylabel('Depth (m)');
        else
          title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Radial (deg)');ylabel('Depth (m)');
        end
        hold off;
      else
        figure;imagesc(freq,dep,tlpress);axis('ij');caxis([tlmin tlmax]);colormap(flipud(jet));
        colorbar;set(colorbar,'YDir','reverse');
        if itype==1
          title('Transmission Loss (dB re 1m) for Pressure');xlabel('Freq (Hz)');ylabel('Depth (m)');
        elseif itype==2
          title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Freq (Hz)');ylabel('Depth (m)');
        else
          title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Freq (Hz)');ylabel('Depth (m)');
        end
      end
  end

  if ~exist('savedat','var')
      savedat=input('Save data? (y or n) ','s');
  else
      disp(savedat)
  end
  disp(' ');

  if isempty(savedat)
  elseif savedat == 'y' || savedat == 'Y'
    if ~exist('fileout','var')
        fileout=input('Enter matlab output filename (no extension needed): ','s');
    else
        disp(fileout)
    end
    if itype==1
      eval(['save ' fileout ' press tlpress dep rad rngout freqout bathout dbathout;']);
    elseif itype==2
      apvr=press/rhoc0; tlapvr=tlpress;
      eval(['save ' fileout ' apvr tlapvr dep rad rngout freqout bathout dbathout rhoc0;']);
    else
      apvz=press/rhoc0; tlapvz=tlpress;
      eval(['save ' fileout ' apvz tlapvz dep rad rngout freqout bathout dbathout rhoc0;']);
    end
  end  % end of 'if savedat == 'y' || savedat == 'Y''

end  % end of 'if opt1 == 3'

% =============================================================================
% COMPUTE DATA ALONG SINGLE DEPTH.
if opt1 == 4

  if ~exist('depout','var')
      depout=input('Enter depth (in meters positive downward) to extract data: ');
  else
      disp(depout)
  end
  disp(' ');

  if depout < dep(1) || depout > dep(nzout)
    disp('Requested depth not contained in input data file.');
    disp(' ');
  else

  if nf > 1 && nrad > 1

   disp(' ');
   disp('1)  Output data for single frequency;');
   disp('2)  Output data for single radial;');
   if ~exist('opt1','var')
       opt1t=input(' ');
   else
       disp(opt1)
   end
   disp(' ');
   
   if opt1t == 1

    disp(['Enter freq, ' num2str(freq(1)) ' Hz to ' num2str(freq(nf)) ' Hz, to extract data:']);
    if ~exist('freqout','var')
        freqout=input(' ');
    else
        disp(freqout)
    end
    nfout=find((freq-freqout)>=0);
    if isempty(nfout)
      if freqout<cfreq
        nfout=1;
      else
        nfout=nf;
      end
    else
      nfout=nfout(1);
    end
    freqout=freq(nfout);
    disp(['Outputting data for freq ' num2str(freq(nfout)) 'Hz.']);
    if nfout <= nf/2
      nfout=nfout+nf/2;
    else
      nfout=nfout-nf/2;
    end
    disp(' ');
    
   else
   
    disp(['Enter radial angle, ' num2str(rad(1)) ' deg to ' num2str(rad(nrad)) ' deg, to extract data:']);
    if ~exist('radout','var')
        radout=input(' ');
    else
        disp(radout)
    end
    nradout=find((rad-radout)>=0);
    if isempty(nradout)
      if radout<0
        nradout=1;
      else
        nradout=nrad;
      end
    else
      nradout=nradout(1);
    end
    radout=rad(nradout);
    disp(['Outputting data for radial ' num2str(rad(nradout)) 'deg.']);
    if nradout < nrad/2
      nradout=nradout+nrad/2+1;
    else
      nradout=nradout-nrad/2+1;
    end
    disp(' ');
   
   end  % end of 'if opt1t == 1'

  elseif nrad > 1
    nfout=1;
    freqout=cfreq;
    opt1t=1;
  else
    nradout=1;
    freqout=freq;
    radout=rad(nrad);
    opt1t=2;
  end  % end of 'if nf > 1 && nrad > 1'

  if opt1t == 1  % output for single freq
      
% First skip consists of (out to proper freq)

    skip=(nfout-1)*(recl + recl*nrad*nrout);
    fid=fopen(pefile,'r');
    fseek(fid,(header+skip)*4,0);

    skip=0; % don't need to use it since reading consecutive radials
    for irng=1:nrout
      for irad=1:nrad
        data=fread(fid,2*nzout,'float32');
        psi=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
        apsi=abs(psi); phipsi=angle(psi);
        atmp=interp1(dep',apsi,depout);phitmp=interp1(dep',phipsi,depout);
        psid(irng,irad)=atmp*cos(phitmp)+1i*atmp*sin(phitmp);
        clear data psi;
      end
    end

    fclose(fid);

% Fill in any NaN gaps
    gaps=find(isnan(psid));
    for i=1:size(gaps)
      psid(gaps(i))=1.e-20;
    end  % end of 'for i=1:size(gaps)'

% Include cylindrical spreading and phase factors
    r0=1.;
    k0=2*pi*freqout/c0;
    for irng=1:nrout
      fact=sqrt(r0/max(1000.*rng(irng),r0))*exp(1i*k0*1000.*rng(irng));
      pressd(irng,:)=fact*psid(irng,:);
      if nrad > 1
        tmp=pressd(irng,nrad/2+2:nrad);
        pressd(irng,nrad/2:nrad)=pressd(irng,1:nrad/2+1);
        pressd(irng,1:nrad/2-1)=tmp;
      end % end of 'if nrad > 1'
    end

    tlpressd=-20*log10(max(abs(pressd),1.e-20));
    tlmin=min(min(tlpressd));tlmax=tlmin+100;

    if doPlots
        if nrad == 1 % single radial
          figure;plot(rng,tlpressd);axis('ij');
          if itype==1
            ylabel('TL (dB re 1m) for Pressure');xlabel('Range (km)');
          elseif itype==2
            ylabel('TL (dB re 1m) for Radial Velocity');xlabel('Range (km)');
          else
            ylabel('TL (dB re 1m) for Vertical Velocity');xlabel('Range (km)');
          end
        elseif nrout == 1 % single range
          figure;plot(rad,tlpressd);axis('ij');
          if itype==1
            ylabel('TL (dB re 1m) for Pressure');xlabel('Radial (deg)');
          elseif itype==2
            ylabel('TL (dB re 1m) for Radial Velocity');xlabel('Radial (deg)');
          else
            ylabel('TL (dB re 1m) for Vertical Velocity');xlabel('Radial (deg)');
          end
        else
          figure;imagesc(rad,rng,tlpressd);caxis([tlmin tlmax]);colormap(flipud(jet));
          colorbar;set(colorbar,'YDir','reverse');
          if itype==1
            title('Transmission Loss (dB re 1m) for Pressure');xlabel('Radial (deg)');ylabel('Range (km)');
          elseif itype==2
            title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Radial (deg)');ylabel('Range (km)');
          else
            title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Radial (deg)');ylabel('Range (km)');
          end
        end % end of 'if nrad == 1'
    end

    if ~exist('savedat','var')
        savedat=input('Save data? (y or n) ','s');
    else
        disp(savedat)
    end
    disp(' ');

    if isempty(savedat)
    elseif savedat == 'y' || savedat == 'Y'
      if ~exist('fileout','var')
          fileout=input('Enter matlab output filename (no extension needed): ','s');
      else
          disp(fileout)
      end
      if itype==1
        eval(['save ' fileout ' pressd tlpressd depout rng rad freqout;']);
      elseif itype==2
        apvrd=pressd/rhoc0; tlapvrd=tlpressd;
        eval(['save ' fileout ' apvrd tlapvrd depout rng rad freqout rhoc0;']);
      else
        apvzd=pressd/rhoc0; tlapvzd=tlpressd;
        eval(['save ' fileout ' apvzd tlapvzd depout rng rad freqout rhoc0;']);
      end
    end  % end of 'if savedat == 'y' || savedat == 'Y''
    
  else  % opt1t == 2, output for single radial
      
% First skip consists of (starting field)+(out to proper radial)

    skip=recl + recl*(nradout-1);
    fid=fopen(pefile,'r');
    fseek(fid,(header+skip)*4,0);

    for ifreq=1:nf
      for irng=1:nrout
        data=fread(fid,2*nzout,'float32');
        psi=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
        apsi=abs(psi); phipsi=angle(psi);
        atmp=interp1(dep',apsi,depout);phitmp=interp1(dep',phipsi,depout);
        psid(irng,ifreq)=atmp*cos(phitmp)+1i*atmp*sin(phitmp);
        clear data psi;
        skip=recl*(nrad-1); % skip now consists of jumping to proper radial
        if irng < nrout && ifreq < nf
          fseek(fid,(skip)*4,0);
        end
      end
      skip=recl*(nrad-1)+recl; % skip now consists of jumping to next frequency (including starting field)
      if ifreq < nf
        fseek(fid,(skip)*4,0);
      end
    end

    fclose(fid);

% Fill in any NaN gaps
    gaps=find(isnan(psid));
    for i=1:size(gaps)
      psid(gaps(i))=1.e-20;
    end  % end of 'for i=1:size(gaps)'

% Include cylindrical spreading and phase factors
    r0=1.;
    for irng=1:nrout
      psid(irng,:)=fftshift(psid(irng,:));
      for ifreq=1:nf
        k0=2*pi*freq(ifreq)/c0;
        fact=sqrt(r0/max(1000.*rng(irng),r0))*exp(1i*k0*1000.*rng(irng));
        pressd(irng,:)=fact*psid(irng,:);
      end
    end

    tlpressd=-20*log10(max(abs(pressd),1.e-20));
    tlmin=min(min(tlpressd));tlmax=tlmin+100;

    if doPlots
        if nf == 1 % single frequency
          figure;plot(rng,tlpressd);axis('ij');
          if itype==1
            ylabel('TL (dB re 1m) for Pressure');xlabel('Range (km)');
          elseif itype==2
            ylabel('TL (dB re 1m) for Radial Velocity');xlabel('Range (km)');
          else
            ylabel('TL (dB re 1m) for Vertical Velocity');xlabel('Range (km)');
          end
        elseif nrout == 1 % single range
          figure;plot(freq,tlpressd);axis('ij');
          if itype==1
            ylabel('TL (dB re 1m) for Pressure');xlabel('Radial (deg)');
          elseif itype==2
            ylabel('TL (dB re 1m) for Radial Velocity');xlabel('Radial (deg)');
          else
            ylabel('TL (dB re 1m) for Vertical Velocity');xlabel('Radial (deg)');
          end
        else
          figure;imagesc(freq,rng,tlpressd);caxis([tlmin tlmax]);colormap(flipud(jet));
          colorbar;set(colorbar,'YDir','reverse');
          if itype==1
            title('Transmission Loss (dB re 1m) for Pressure');xlabel('Freq (Hz)');ylabel('Range (km)');
          elseif itype==2
            title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Freq (Hz)');ylabel('Range (km)');
          else
            title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Freq (Hz)');ylabel('Range (km)');
          end
        end % end of 'if nrad == 1'
    end

    if ~exist('savedat','var')
        savedat=input('Save data? (y or n) ','s');
    else
        disp(savedat)
    end
    disp(' ');

    if isempty(savedat)
    elseif savedat == 'y' || savedat == 'Y'
      if ~exist('fileout','var')
          fileout=input('Enter matlab output filename (no extension needed): ','s');
      else
          disp(fileout)
      end
      if itype==1
        eval(['save ' fileout ' pressd tlpressd depout rng freq radout;']);
      elseif itype==2
        apvrd=pressd/rhoc0; tlapvrd=tlpressd;
        eval(['save ' fileout ' apvrd tlapvrd depout rng freq radout rhoc0;']);
      else
        apvzd=pressd/rhoc0; tlapvzd=tlpressd;
        eval(['save ' fileout ' apvzd tlapvzd depout rng freq radout rhoc0;']);
      end
    end  % end of 'if savedat == 'y' || savedat == 'Y''
      
  end

  end  % end of 'if depout < dep(1) || depout > dep(nzout)'

end  % end of 'if opt1 == 4'

% =============================================================================
% COMPUTE DATA ALONG SINGLE INTERFACE.
if opt1 == 5

  disp('1)  Extract data at water/bottom interface.');
  disp('2)  Extract data at bottom/basement interface.');
  if ~exist('opt2','var')
      opt2=input(' ');
  else
      disp(opt2)
  end
  disp(' ');
  if opt2 == 1
    depout=bath;
  else
    depout=dbath;
  end

  if nf > 1 && nrad > 1

   disp(' ');
   disp('1)  Output data for single frequency;');
   disp('2)  Output data for single radial;');
   if ~exist('opt1t','var')
       opt1t=input(' ');
   else
       disp(opt1t)
   end
   disp(' ');
   
   if opt1t == 1

    disp(['Enter freq, ' num2str(freq(1)) ' Hz to ' num2str(freq(nf)) ' Hz, to extract data:']);
    if ~exist('freqout','var')
        freqout=input(' ');
    else
        disp(freqout)
    end
    nfout=find((freq-freqout)>=0);
    if isempty(nfout)
      if freqout<cfreq
        nfout=1;
      else
        nfout=nf;
      end
    else
      nfout=nfout(1);
    end
    freqout=freq(nfout);
    disp(['Outputting data for freq ' num2str(freq(nfout)) 'Hz.']);
    if nfout <= nf/2
      nfout=nfout+nf/2;
    else
      nfout=nfout-nf/2;
    end
    disp(' ');
    
   else
   
    disp(['Enter radial angle, ' num2str(rad(1)) ' deg to ' num2str(rad(nrad)) ' deg, to extract data:']);
    if ~exist('radout','var')
        radout=input(' ');
    end
    nradout=find((rad-radout)>=0);
    if isempty(nradout)
      if radout<0
        nradout=1;
      else
        nradout=nrad;
      end
    else
      nradout=nradout(1);
    end
    radout=rad(nradout);
    disp(['Outputting data for radial ' num2str(rad(nradout)) 'deg.']);
    if nradout < nrad/2
      nradout=nradout+nrad/2+1;
    else
      nradout=nradout-nrad/2+1;
    end
    disp(' ');
   
   end  % end of 'if opt1t == 1'

  elseif nrad > 1
    nfout=1;
    freqout=cfreq;
    opt1t=1;
  else
    nradout=1;
    freqout=freq;
    radout=rad(nrad);
    opt1t=2;
  end  % end of 'if nf > 1 && nrad > 1'

  if opt1t == 1  % output for single freq
      
% First skip consists of (out to proper freq)

    skip=(nfout-1)*(recl + recl*nrad*nrout);
    fid=fopen(pefile,'r');
    fseek(fid,(header+skip)*4,0);

    skip=0; % don't need to use it since reading consecutive radials
    for irng=1:nrout
      for irad=1:nrad
        data=fread(fid,2*nzout,'float32');
        psi=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
        apsi=abs(psi); phipsi=angle(psi);
        atmp=interp1(dep',apsi,depout(irng,irad));phitmp=interp1(dep',phipsi,depout(irng,irad));
        psid(irng,irad)=atmp*cos(phitmp)+1i*atmp*sin(phitmp);
        clear data psi;
      end
    end

    fclose(fid);

% Fill in any NaN gaps
    gaps=find(isnan(psid));
    for i=1:size(gaps)
      psid(gaps(i))=1.e-20;
    end  % end of 'for i=1:size(gaps)'

% Include cylindrical spreading and phase factors
    r0=1.;
    k0=2*pi*freqout/c0;
    for irng=1:nrout
      fact=sqrt(r0/max(1000.*rng(irng),r0))*exp(1i*k0*1000.*rng(irng));
      pressd(irng,:)=fact*psid(irng,:);
      if nrad > 1
        tmp=pressd(irng,nrad/2+2:nrad);
        pressd(irng,nrad/2:nrad)=pressd(irng,1:nrad/2+1);
        pressd(irng,1:nrad/2-1)=tmp;
      end % end of 'if nrad > 1'
    end

    tlpressd=-20*log10(max(abs(pressd),1.e-20));
    tlmin=min(min(tlpressd));tlmax=tlmin+100;

    if doPlots
        if nrad == 1 % single radial
          figure;plot(rng,tlpressd);axis('ij');
          if itype==1
            ylabel('TL (dB re 1m) for Pressure');xlabel('Range (km)');
          elseif itype==2
            ylabel('TL (dB re 1m) for Radial Velocity');xlabel('Range (km)');
          else
            ylabel('TL (dB re 1m) for Vertical Velocity');xlabel('Range (km)');
          end
        elseif nrout == 1 % single range
          figure;plot(rad,tlpressd);axis('ij');
          if itype==1
            ylabel('TL (dB re 1m) for Pressure');xlabel('Radial (deg)');
          elseif itype==2
            ylabel('TL (dB re 1m) for Radial Velocity');xlabel('Radial (deg)');
          else
            ylabel('TL (dB re 1m) for Vertical Velocity');xlabel('Radial (deg)');
          end
        else
          figure;imagesc(rad,rng,tlpressd);caxis([tlmin tlmax]);colormap(flipud(jet));
          colorbar;set(colorbar,'YDir','reverse');
          if itype==1
            title('Transmission Loss (dB re 1m) for Pressure');xlabel('Radial (deg)');ylabel('Range (km)');
          elseif itype==2
            title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Radial (deg)');ylabel('Range (km)');
          else
            title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Radial (deg)');ylabel('Range (km)');
          end
        end % end of 'if nrad == 1'
    end

    if ~exist('savedat','var')
        savedat=input('Save data? (y or n) ','s');
    else
        disp(savedat)
    end
    disp(' ');

    if isempty(savedat)
    elseif savedat == 'y' || savedat == 'Y'
      if ~exist('fileout','var')
          fileout=input('Enter matlab output filename (no extension needed): ','s');
      else
          disp(fileout)
      end
      if itype==1
        eval(['save ' fileout ' pressd tlpressd depout rng rad freqout;']);
      elseif itype==2
        apvrd=pressd/rhoc0; tlapvrd=tlpressd;
        eval(['save ' fileout ' apvrd tlapvrd depout rng rad freqout rhoc0;']);
      else
        apvzd=pressd/rhoc0; tlapvzd=tlpressd;
        eval(['save ' fileout ' apvzd tlapvzd depout rng rad freqout rhoc0;']);
      end
    end  % end of 'if savedat == 'y' || savedat == 'Y''
    
  else  % opt1t == 2, output for single radial
      
% First skip consists of (starting field)+(out to proper radial)

    skip=recl + recl*(nradout-1);
    fid=fopen(pefile,'r');
    fseek(fid,(header+skip)*4,0);

    for ifreq=1:nf
      for irng=1:nrout
        data=fread(fid,2*nzout,'float32');
        psi=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
        apsi=abs(psi); phipsi=angle(psi);
        atmp=interp1(dep',apsi,depout(irng,nradout));phitmp=interp1(dep',phipsi,depout(irng,nradout));
        psid(irng,ifreq)=atmp*cos(phitmp)+1i*atmp*sin(phitmp);
        clear data psi;
        skip=recl*(nrad-1); % skip now consists of jumping to proper radial
        if irng < nrout && ifreq < nf
          fseek(fid,(skip)*4,0);
        end
      end
      skip=recl*(nrad-1)+recl; % skip now consists of jumping to next frequency (including starting field)
      if ifreq < nf
        fseek(fid,(skip)*4,0);
      end
    end

    fclose(fid);

% Fill in any NaN gaps
    gaps=find(isnan(psid));
    for i=1:size(gaps)
      psid(gaps(i))=1.e-20;
    end  % end of 'for i=1:size(gaps)'

% Include cylindrical spreading and phase factors
    r0=1.;
    for irng=1:nrout
      psid(irng,:)=fftshift(psid(irng,:));
      for ifreq=1:nf
        k0=2*pi*freq(ifreq)/c0;
        fact=sqrt(r0/max(1000.*rng(irng),r0))*exp(1i*k0*1000.*rng(irng));
        pressd(irng,:)=fact*psid(irng,:);
      end
    end

    tlpressd=-20*log10(max(abs(pressd),1.e-20));
    tlmin=min(min(tlpressd));tlmax=tlmin+100;

    if doPlots
        if nf == 1 % single frequency
          figure;plot(rng,tlpressd);axis('ij');
          if itype==1
            ylabel('TL (dB re 1m) for Pressure');xlabel('Range (km)');
          elseif itype==2
            ylabel('TL (dB re 1m) for Radial Velocity');xlabel('Range (km)');
          else
            ylabel('TL (dB re 1m) for Vertical Velocity');xlabel('Range (km)');
          end
        elseif nrout == 1 % single range
          figure;plot(freq,tlpressd);axis('ij');
          if itype==1
            ylabel('TL (dB re 1m) for Pressure');xlabel('Radial (deg)');
          elseif itype==2
            ylabel('TL (dB re 1m) for Radial Velocity');xlabel('Radial (deg)');
          else
            ylabel('TL (dB re 1m) for Vertical Velocity');xlabel('Radial (deg)');
          end
        else
          figure;imagesc(freq,rng,tlpressd);caxis([tlmin tlmax]);colormap(flipud(jet));
          colorbar;set(colorbar,'YDir','reverse');
          if itype==1
            title('Transmission Loss (dB re 1m) for Pressure');xlabel('Freq (Hz)');ylabel('Range (km)');
          elseif itype==2
            title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Freq (Hz)');ylabel('Range (km)');
          else
            title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Freq (Hz)');ylabel('Range (km)');
          end
        end % end of 'if nrad == 1'
    end

    if ~exist('savedat','var')
        savedat=input('Save data? (y or n) ','s');
    else
        disp(savedat)
    end
    disp(' ');

    if isempty(savedat)
    elseif savedat == 'y' || savedat == 'Y'
      if ~exist('fileout','var')
          fileout=input('Enter matlab output filename (no extension needed): ','s');
      else
          disp(fileout)
      end
      depout=depout(:,nradout);
      if itype==1
        eval(['save ' fileout ' pressd tlpressd depout rng freq radout;']);
      elseif itype==2
        apvrd=pressd/rhoc0; tlapvrd=tlpressd;
        eval(['save ' fileout ' apvrd tlapvrd depout rng freq radout rhoc0;']);
      else
        apvzd=pressd/rhoc0; tlapvzd=tlpressd;
        eval(['save ' fileout ' apvzd tlapvzd depout rng freq radout rhoc0;']);
      end
    end  % end of 'if savedat == 'y' || savedat == 'Y''
      
  end


end  % end of 'if opt1 == 5'

% =============================================================================
% COMPUTE TRAVEL TIME DATA.
if opt1 == 6

if nf == 1
  disp('Input file contains only single frequency data.');
else

  disp('1)  Compute data for single radial/range;');
  disp('2)  Compute data for single depth;');
  disp('3)  Compute data for single interface;');
  if ~exist('optt','var')
      optt=input(' ');
  else
      disp(optt)
  end
  disp(' ');

% -----------------------------------------------------------------------------
% COMPUTE DATA ALONG SINGLE RADIAL/RANGE.
  if optt == 1

    if nrad > 1
      disp(['Enter radial angle, ' num2str(rad(1)) ' deg to ' num2str(rad(nrad)) ' deg, to extract data:']);
      if ~exist('radout','var')
          radout=input(' ');
      else
          disp(radout)
      end
      nradout=find((rad-radout)>=0);
      if isempty(nradout)
        if radout<0
          nradout=1;
        else
          nradout=nrad;
        end
      else
        nradout=nradout(1);
      end
      radout=rad(nradout);
      disp(['Outputting data for radial ' num2str(rad(nradout)) 'deg.']);
      if nradout < nrad/2
        nradout=nradout+nrad/2+1;
      else
        nradout=nradout-nrad/2+1;
      end
      disp(' ');
    else
      radout=0;
      nradout=1;
    end  % end of 'if nrad > 1'

    if nrout > 1
      disp(['Enter range, ' num2str(rng(1)) ' km to ' num2str(rng(nrout)) ' km, to extract data:']);
      if ~exist('rngout','var')
          rngout=input(' ');
      else
          disp(rngout)
      end
      nrngout=find((rng-rngout)>=0);
      if isempty(nrngout)
        if rngout<rng(1)
          nrngout=1;
        else
          nrngout=nrout;
        end
      else
        nrngout=nrngout(1);
      end
      rngout=rng(nrngout);
      disp(['Outputting data for range ' num2str(rng(nrngout)) 'km.']);
      disp(' ');
    else
      nrngout=1;
      rngout=rng(1);
    end  % end of 'if nrout > 1'

   disp('1)  Output data as time -vs- depth;');
   disp('2)  Output data as time -vs- arrival angle;');
   if ~exist('optt1','var')
       optt1=input(' ');
   else
       disp(optt1)
   end
   disp(' ');

% Skip consists of (starting field)+(out to proper radial)+(out to proper range)
    skip=recl + recl*(nradout-1) + recl*nrad*(nrngout-1);
    fid=fopen(pefile,'r');
    fseek(fid,(header+skip)*4,0);

    skip=recl+recl*(nrad*nrout-1);
    for ifreq=1:nf
      data=fread(fid,2*nzout,'float32');
      if ifreq < nf
        fseek(fid,(skip)*4,0);
      end
      psi(:,ifreq)=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
      for iz=1:nzout
        if dep(iz)<surf(nrngout)
          psi(iz,ifreq)=1.e-20;
        end
      end      
      clear data;
    end

    fclose(fid);

% Fill in any NaN gaps
    gaps=find(isnan(psi));
    for i=1:size(gaps)
      psi(gaps(i))=1.e-20;
    end  % end of 'for i=1:size(gaps)'

% Include cylindrical spreading
% NOTE:  Don't need to include phase factors since it results
%        only in time-shift relative to reduced time.  However,
%        phase factor does arise due to basebanding and is added
%        after FFT transform.
    r0=1.;
    press=psi*sqrt(r0/max(1000.*rngout,r0));

    hanwinf=hanning(nf+1);
    hanwinf=fftshift(hanwinf(1:nf)/sum(hanwinf(1:nf)));
    phs=2*pi*cfreq*time;

    if optt1 == 1

% Convert to time domain, re-order times and correct for basebanded phase
      for iz=1:nzout
        press(iz,:)=fftshift(fft(press(iz,:).*hanwinf')).*exp(-1i*phs);
      end
      timeout=time+rngout/(c0/1000.);

      tlpress=-20*log10(max(abs(press),1.e-20));
      tlmin=min(min(tlpress));tlmax=tlmin+100;
      bathout=bath(nrngout,nradout);dbathout=dbath(nrngout,nradout);

      if doPlots
          figure;imagesc(timeout,dep,tlpress);axis('ij');caxis([tlmin tlmax]);colormap(flipud(jet));
          colorbar;set(colorbar,'YDir','reverse');
          if itype==1
            title('Transmission Loss (dB re 1m) for Pressure');xlabel('Time (sec)');ylabel('Depth (m)');
          elseif itype==2
            title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Time (sec)');ylabel('Depth (m)');
          else
            title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Time (sec)');ylabel('Depth (m)');
          end
      end

      if ~exist('savedat','var')
          savedat=input('Save data? (y or n) ','s');
      else
          disp(savedat)
      end
      disp(' ');

      if isempty(savedat)
      elseif savedat == 'y' || savedat == 'Y'
        if ~exist('fileout','var')
            fileout=input('Enter matlab output filename (no extension needed): ','s');
        else
            disp(fileout)
        end
        if itype==1
          eval(['save ' fileout ' press tlpress dep timeout rngout radout;']);
        elseif itype==2
          apvr=press/rhoc0; tlapvr=tlpress;
          eval(['save ' fileout ' apvr tlapvr dep timeout rngout radout rhoc0;']);
        else
          apvz=press/rhoc0; tlapvz=tlpress;
          eval(['save ' fileout ' apvz tlapvz dep timeout rngout radout rhoc0;']);
        end
      end  % end of 'if savedat == 'y' || savedat == 'Y''

    elseif optt1 == 2
    
      disp(['Current depths stored from ' num2str(dep(1)) ' m to ' num2str(dep(nzout)) ' m.']);
      depmin=0;depmax=0;
      while (depmax<=depmin)
        if ~exist('depmin','var')
            depmin=input('Enter min depth [m] to use in beamformer: ');
        else
            disp(depmin)
        end
        if ~exist('depmax','var')
            depmax=input('Enter max depth [m] to use in beamformer: ');
        else
            disp(depmax)
        end
        % Error checking of scripted inputs to break out of while loop
        if depmin==depmax
            depmax = depmin + 1;
        end
        if depmin>depmax
            tmpdep = depmax;
            depmax = depmin;
            depmin = tmpdep;
            clear tmpdep;
        end
      end

% Perform FFT beamforming, re-order times and correct for basebanded phase

%  Set parameters in depth domain
      delz=dep(2)-dep(1);    % depth spacing
      idx_strt=find((dep-depmin)>=0);idx_strt=idx_strt(1);
      idx_end=find((dep-depmax)>=0);
      if isempty(idx_end)
        idx_end=nzout;
      else
        idx_end=idx_end(1);
      end
      nzbeam=idx_end-idx_strt+1;      
      pow_z=ceil(log10(nzbeam)/log10(2));
      mzt=2^(pow_z);       %  Total depth domain FFT size

      zvsfreq=zeros(mzt,nf);
      for iz=1:nzbeam
        zvsfreq(iz,:)=press(iz+idx_strt-1,:).*hanwinf';
      end
      
      kvsfreq=zeros(mzt,nf);

      hanwinz=hanning(mzt+1);
      hanwinz=hanwinz(1:mzt)/sum(hanwinz(1:mzt));
      theta=(-90:0.5:90)*pi/180;
      angvsfreq=zeros(length(theta),nf);
      
      freq=fftshift(freq);
      for n=1:nf
        kvsfreq(:,n)=fft(fftshift(zvsfreq(:,n).*hanwinz));
%  Transform K --> angle
        kz=mzt*delz*freq(n)*sin(theta)/c0;
        posk=find(kz>=0);
        for i=posk
          kz(i)=kz(i)+1;
        end
        negk=find(kz<0);
        for i=negk
          kz(i)=mzt+kz(i)+1;
        end
        for i=negk
          if kz(i) >= mzt/2+1
            kzi_1=ceil(kz(i));
            if kzi_1 > mzt
              kzi_1=1;
            end
            kzi=floor(kz(i));
            angvsfreq(i,n)=kvsfreq(kzi,n)+(kz(i)-kzi).*(kvsfreq(kzi_1,n)-kvsfreq(kzi,n));
          end
        end
        for i=posk
          if kz(i) < mzt/2
            kzi_1=ceil(kz(i));
            kzi=floor(kz(i));
            angvsfreq(i,n)=kvsfreq(kzi,n)+(kz(i)-kzi).*(kvsfreq(kzi_1,n)-kvsfreq(kzi,n));
          end
        end
      end
      freq=fftshift(freq);

%  Transform to time domain
      for m=1:length(theta)
        pressbeam(m,:)=fftshift(fft(angvsfreq(m,:).*hanwinf')).*exp(-1i*phs);
      end

      clear delz idx_strt idx_end nzbeam pow_z mzt hanwinz hanwinf
      clear zvsfreq kvsfreq angvsfreq posk negk phs kz kzi kzi_1

      timeout=time+rngout/(c0/1000.);
      theta=theta*180/pi;

      tlpressbeam=-20*log10(max(abs(pressbeam),1.e-20));
      tlmin=min(min(tlpressbeam));tlmax=tlmin+100;

      if doPlots
          figure;imagesc(timeout,theta,tlpressbeam);caxis([tlmin tlmax]);colormap(flipud(jet));
          colorbar;set(colorbar,'YDir','reverse');
          if itype==1
            title('Transmission Loss (dB re 1m) for Pressure');xlabel('Time (sec)');ylabel('Arrival Angle (deg)');
          elseif itype==2
            title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Time (sec)');ylabel('Arrival Angle (deg)');
          else
            title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Time (sec)');ylabel('Arrival Angle (deg)');
          end
      end

      if ~exist('savedat','var')
          savedat=input('Save data? (y or n) ','s');
      else
          disp(savedat)
      end
      disp(' ');

      if isempty(savedat)
      elseif savedat == 'y' || savedat == 'Y'
        if ~exist('fileout','var')
            fileout=input('Enter matlab output filename (no extension needed): ','s');
        else
            disp(fileout)
        end
        if itype==1
          eval(['save ' fileout ' pressbeam tlpressbeam theta timeout rngout radout;']);
        elseif itype==2
          apvrbeam=pressbeam/rhoc0; tlapvrbeam=tlpressbeam;
          eval(['save ' fileout ' apvrbeam tlapvrbeam theta timeout rngout radout rhoc0;']);
        else
          apvzbeam=pressbeam/rhoc0; tlapvzbeam=tlpressbeam;
          eval(['save ' fileout ' apvzbeam tlapvzbeam theta timeout rngout radout rhoc0;']);
        end
      end  % end of 'if savedat == 'y' || savedat == 'Y''

    end  % end of 'if optt1 == 1'

  end  % end of 'if optt == 1'

% -----------------------------------------------------------------------------
% COMPUTE DATA ALONG SINGLE DEPTH.
  if optt == 2

  if ~exist('depout','var')
      depout=input('Enter depth (in meters positive downward) to extract data: ');
  else
      disp(depout)
  end
  disp(' ');

  if depout < dep(1) || depout > dep(nzout)
    disp('Requested depth not contained in input data file.');
    disp(' ');
  else

    if nrad > 1 && nrout > 1
      optout=0;
      while optout ~= 1 && optout ~= 2
        disp('1)  Compute data for single radial.');
        disp('2)  Compute data for single range.');
        if ~exist('optout','var')
            optout=input(' ');
        else
            disp(optout);
        end
        disp(' ');
      end % end of 'while optout ~= 1 && optout ~= 2,'
      if optout == 1
        optt1 = 1;
        disp(['Enter radial angle, ' num2str(rad(1)) ' deg to ' num2str(rad(nrad)) ' deg, to extract data:']);
        if ~exist('radout','var')
            radout=input(' ');
        else
            disp(radout)
        end
        nradout=find((rad-radout)>=0);
        if isempty(nradout)
          if radout<0
            nradout=1;
          else
            nradout=nrad;
          end
        else
          nradout=nradout(1);
        end
        radout=rad(nradout);
        disp(['Outputting data for radial ' num2str(rad(nradout)) 'deg.']);
        if nradout < nrad/2
          nradout=nradout+nrad/2+1;
        else
          nradout=nradout-nrad/2+1;
        end
        rngout=rng;
        disp(' ');
      else % optout == 2
        disp(['Enter range, ' num2str(rng(1)) ' km to ' num2str(rng(nrout)) ' km, to extract data:']);
        if ~exist('rngout','var')
            rngout=input(' ');
        else
            disp(rngout)
        end
        nrngout=find((rng-rngout)>=0);
        if isempty(nrngout)
          if rngout<rng(1)
            nrngout=1;
          else
            nrngout=nrout;
          end
        else
          nrngout=nrngout(1);
        end
        rngout=rng(nrngout);
        disp(['Outputting data for range ' num2str(rng(nrngout)) 'km.']);
        radout=rad;
        disp(' ');
      end % end of 'if optout == 1'
    elseif nrad == 1
      optout=1;
      optt1=1;
      nradout=1;
      radout=rad;
      nrngout=nrout;
      rngout=rng;
    else
      optout=2;
      nrngout=1;
      rngout=rng;
      nradout=nrad;
      radout=rad;
    end  % end of 'if nrad > 1 && nrout > 1'

    for ifreq=1:nf

% First skip consists of (out to proper freq)+(starting field)+
%                          (out to proper radial/range)

      if optout == 1 % extract along single radial

        skip=(ifreq-1)*(recl + recl*nrad*nrout) + recl + recl*(nradout-1);
        fid=fopen(pefile,'r');
        fseek(fid,(header+skip)*4,0);

        skip=recl*(nrad-1);
        for irng=1:nrout
          data=fread(fid,2*nzout,'float32');
          if irng < nrout && nrad > 1
            fseek(fid,(skip)*4,0);
          end
          psi=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
          apsi=abs(psi); phipsi=angle(psi);
          atmp=interp1(dep',apsi,depout);phitmp=interp1(dep',phipsi,depout);
          psid(irng,ifreq)=atmp*cos(phitmp)+1i*atmp*sin(phitmp);
          clear data psi;
        end

      else % extract at single range

        skip=(ifreq-1)*(recl + recl*nrad*nrout) + recl + recl*nrad*(nrngout-1);
        fid=fopen(pefile,'r');
        fseek(fid,(header+skip)*4,0);

        skip=0; % don't need to use it since reading consecutive radials
        for irad=1:nrad
          data=fread(fid,2*nzout,'float32');
          psi=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
          apsi=abs(psi); phipsi=angle(psi);
          atmp=interp1(dep',apsi,depout);phitmp=interp1(dep',phipsi,depout);
          psid(irad,ifreq)=atmp*cos(phitmp)+1i*atmp*sin(phitmp);
          clear data psi;
        end

      end % end of 'if optout == 1'
      fclose(fid);

    end % end of 'for ifreq=1:nf,'

% Fill in any NaN gaps
    gaps=find(isnan(psid));
    for i=1:size(gaps)
      psid(gaps(i))=1.e-20;
    end  % end of 'for i=1:size(gaps)'

% Include cylindrical spreading
% NOTE:  Don't need to include phase factors since it results
%        only in time-shift relative to reduced time.  However,
%        phase factor does arise due to basebanding and is added
%        after FFT transform.
    r0=1.;
    if optout == 1
      for irng=1:nrout
        pressd(irng,:)=psid(irng,:)*sqrt(r0/max(1000.*rng(irng),r0));
      end
    else % computing for single range, optout == 1
      pressd=psid*sqrt(r0/max(1000.*rng(nrngout),r0))/rhoc0;
      disp('1)  Output data as time -vs- radial;');
      disp('2)  Output data as time -vs- azimuthal arrival angle;');
      if ~exist('outt1','var')
          optt1=input(' ');
      else
          disp(optt1)
      end
      disp(' ');
    end  % end of 'if optout == 1'
    
% Convert to time domain, re-order freqs and correct for basebanded phase
    hanwin=hanning(nf+1);
    hanwin=fftshift(hanwin(1:nf)/sum(hanwin(1:nf)));
    phs=2*pi*cfreq*time;
    if optout == 1
      for irng=1:nrout
        pressd(irng,:)=fftshift(fft(pressd(irng,:).*hanwin')).*exp(-1i*phs);
        timeout(irng,:)=time+rngout(1,irng)/(c0/1000.);
      end
      
      tlpressd=-20*log10(max(abs(pressd),1.e-20));
      tlmin=min(min(tlpressd));tlmax=tlmin+100;

    else
      timeout=time+rngout/(c0/1000.);
      tmp=pressd(nrad/2+2:nrad,:);pressd(nrad/2:nrad,:)=pressd(1:nrad/2+1,:);pressd(1:nrad/2-1,:)=tmp;
      if optt1 == 1
        for irad=1:nrad
          pressd(irad,:)=fftshift(fft(pressd(irad,:).*hanwin')).*exp(-1i*phs);
        end
        tlpressd=-20*log10(max(abs(pressd),1.e-20));
        tlmin=min(min(tlpressd));tlmax=tlmin+100;
      else % beamform in horizontal
         
% Perform horizontal FFT beamforming, re-order times and correct for basebanded phase
        el=radout*pi/180*rngout*1000;
        nel=size(radout,2);
        disp(['Current elements stored from ' num2str(el(1)) ' m to ' num2str(el(nel)) ' m.']);
        elmin=0;elmax=0;
        while (elmax<=elmin)
          if ~exist('elmin','var')
              elmin=input('Enter min element location [m] to use in beamformer: ');
          else
              disp(elmin)
          end
          if ~exist('elmax','var')
              elmax=input('Enter max element location [m] to use in beamformer: ');
          else
              disp(elmax)
          end
          % Error checking for scripted inputs
          if elmin == elmax
              elmax = elmin + 1;
          end
          if elmin > elmax
              tmpel = elmax;
              elmax = elmin;
              elmin = tmpel;
              clear tmpel;
          end
        end

%  Set parameters in horizontal domain
        del=el(2)-el(1);    % element spacing
        idx_strt=find((el-elmin)>=0);idx_strt=idx_strt(1);
        idx_end=find((el-elmax)>=0);
        if isempty(idx_end)
          idx_end=nel;
        else
          idx_end=idx_end(1);
        end
        nbeam=idx_end-idx_strt+1;      
        pow_n=floor(log10(nbeam)/log10(2));
        mt=2^(pow_n);       %  Total depth domain FFT size
        nbeam=mt;
        idx_strt=nel/2-nbeam/2+1;idx_end=nbeam+idx_strt-1;
      
        hanwinf=hanning(nf+1);
        hanwinf=fftshift(hanwinf(1:nf)/sum(hanwinf(1:nf)));
      
%         xvsfreq=zeros(mt,nf);
        xvsfreq=pressd(idx_strt:idx_end,:);
      
        kvsfreq=zeros(mt,nf);

        hanwinx=hanning(mt+1);
        hanwinx=hanwinx(1:mt)/sum(hanwinx(1:mt));
        theta=(-90:0.5:90)*pi/180;
        angvsfreq=zeros(length(theta),nf);
        freq=fftshift(freq);

        for n=1:nf
          kvsfreq(:,n)=fft(fftshift(xvsfreq(:,n).*hanwinx));
%          kvsfreq(:,n)=fft(xvsfreq(:,n).*fftshift(hanwinx));
%  Transform K --> angle
          kx=mt*del*freq(n)*sin(theta)/c0;
          posk=find(kx>=0);
          for i=posk
            kx(i)=kx(i)+1;
          end
          negk=find(kx<0);
          for i=negk
            kx(i)=mt+kx(i)+1;
          end
          for i=negk
            if kx(i) >= mt/2+1
              kxi_1=ceil(kx(i));
              if kxi_1 > mt
                kxi_1=1;
              end
              kxi=floor(kx(i));
              angvsfreq(i,n)=kvsfreq(kxi,n)+(kx(i)-kxi).*(kvsfreq(kxi_1,n)-kvsfreq(kxi,n));
            end
          end
          for i=posk
            if kx(i) < mt/2
              kxi_1=ceil(kx(i));
              kxi=floor(kx(i));
              angvsfreq(i,n)=kvsfreq(kxi,n)+(kx(i)-kxi).*(kvsfreq(kxi_1,n)-kvsfreq(kxi,n));
            end
          end
        end
        freq=fftshift(freq);

%  Transform to time domain
        for m=1:length(theta)
          pressbeam(m,:)=fftshift(fft(angvsfreq(m,:).*hanwinf'));
        end

        theta=theta*180/pi;

      end % end of 'if optt == 1'

    end % end of 'if optout == 1'

    if doPlots
        if optout == 1 % output single radial

          if nrout == 1 % single range
            figure;plot(timeout,tlpressd);axis('ij');
            if itype==1
              ylabel('TL (dB re 1m) for Pressure');xlabel('Time (sec)');
            elseif itype==2
              ylabel('TL (dB re 1m) for Radial Velocity');xlabel('Time (sec)');
            else
              ylabel('TL (dB re 1m) for Vertical Velocity');xlabel('Time (sec)');
            end
          else
            figure;imagesc(time,rng,tlpressd);caxis([tlmin tlmax]);colormap(flipud(jet));
            set(colorbar,'YDir','reverse');
            if itype==1
              title('Transmission Loss (dB re 1m) for Pressure');xlabel('Reduced Time (sec)');ylabel('Range (km)');
            elseif itype==2
              title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Reduced Time (sec)');ylabel('Range (km)');
            else
              title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Reduced Time (sec)');ylabel('Range (km)');
            end
          end % end of 'if nrout == 1'

        else % output single range

          if optt1 == 1
            if nrad == 1 % single radial
              figure;plot(timeout,tlpressd);axis('ij');
              if itype==1
                ylabel('TL (dB re 1m) for Pressure');xlabel('Time (sec)');
              elseif itype==2
                ylabel('TL (dB re 1m) for Radial Velocity');xlabel('Time (sec)');
              else
                ylabel('TL (dB re 1m) for Vertical Velocity');xlabel('Time (sec)');
              end
            else
              figure;imagesc(timeout,rad,tlpressd);caxis([tlmin tlmax]);colormap(flipud(jet));
              set(colorbar,'YDir','reverse');
              if itype==1
                title('Transmission Loss (dB re 1m) for Pressure');xlabel('Time (sec)');ylabel('Radial (deg)');
              elseif itype==2
                title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Time (sec)');ylabel('Radial (deg)');
              else
                title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Time (sec)');ylabel('Radial (deg)');
              end
            end % end of 'if nrad == 1'
          else % output beamformed data
            tlpressbeam=-20*log10(max(abs(pressbeam),1.e-20));
            tlmin=min(min(tlpressbeam));tlmax=tlmin+100;
            figure;imagesc(timeout,theta,tlpressbeam);caxis([tlmin tlmax]);colormap(flipud(jet));
            set(colorbar,'YDir','reverse');
            if itype==1
              title('Transmission Loss (dB re 1m) for Pressure');xlabel('Time (sec)');ylabel('Azimuthal Arrival Angle (deg)');
            elseif itype==2
              title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Time (sec)');ylabel('Azimuthal Arrival Angle (deg)');
            else
              title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Time (sec)');ylabel('Azimuthal Arrival Angle (deg)');
            end
          end

        end % end of 'if optout == 1'
    end

    if ~exist('savedat','var')
        savedat=input('Save data? (y or n) ','s');
    else
        disp(savedat)
    end
    disp(' ');

    if isempty(savedat)
    elseif savedat == 'y' || savedat == 'Y'
      if ~exist('fileout','var')
          fileout=input('Enter matlab output filename (no extension needed): ','s');
      else
          disp(fileout)
      end
      if optt1 == 1
        if itype==1
          eval(['save ' fileout ' pressd tlpressd depout time timeout rngout radout;']);
        elseif itype==2
          apvrd=pressd/rhoc0; tlapvrd=tlpressd;
          eval(['save ' fileout ' apvrd tlapvrd depout time timeout rngout radout rhoc0;']);
        else
          apvzd=pressd/rhoc0; tlapvzd=tlpressd;
          eval(['save ' fileout ' apvzd tlapvzd depout time timeout rngout radout rhoc0;']);
        end
      else
        if itype==1
          eval(['save ' fileout ' pressbeam tlpressbeam depout time timeout rngout theta;']);
        elseif itype==2
          apvrbeam=pressbeam/rhoc0; tlapvrbeam=tlpressbeam;
          eval(['save ' fileout ' apvrbeam tlapvrbeam depout time timeout rngout theta rhoc0;']);
        else
          apvzbeam=pressbeam/rhoc0; tlapvzbeam=tlpressbeam;
          eval(['save ' fileout ' apvzbeam tlapvzbeam depout time timeout rngout theta rhoc0;']);
        end
      end
    end  % end of 'if savedat == 'y' || savedat == 'Y''

  end  % end of 'if depout < dep(1) || depout > dep(nzout)'

  end  % end of 'if optt == 2'

% -----------------------------------------------------------------------------
% COMPUTE DATA ALONG SINGLE INTERFACE.
  if optt == 3

    disp('1)  Extract data at water/bottom interface.');
    disp('2)  Extract data at bottom/basement interface.');
    if ~exist('opt2','var')
        opt2=input(' ');
    else
        disp(opt2)
    end
    disp(' ');
    if opt2 == 1
      depout=bath;
    else
      depout=dbath;
    end

    if nrad > 1 && nrout > 1
      optout=0;
      while optout ~= 1 && optout ~= 2
        disp('1)  Compute data for single radial.');
        disp('2)  Compute data for single range.');
        if ~exist('optout','var')
            optout=input(' ');
        else
            disp(optout)
        end
        disp(' ');
      end % end of 'while optout ~= 0 && optout ~= 1,'
      if optout == 1
        disp(['Enter radial angle, ' num2str(rad(1)) ' deg to ' num2str(rad(nrad)) ' deg, to extract data:']);
        if ~exist('radout','var')
            radout=input(' ');
        else
            disp(radout)
        end
        nradout=find((rad-radout)>=0);
        if isempty(nradout)
          if radout<0
            nradout=1;
          else
            nradout=nrad;
          end
        else
          nradout=nradout(1);
        end
        radout=rad(nradout);
        disp(['Outputting data for radial ' num2str(rad(nradout)) 'deg.']);
        if nradout < nrad/2
          nradout=nradout+nrad/2+1;
        else
          nradout=nradout-nrad/2+1;
        end
        rngout=rng;
        disp(' ');
      else % optout == 2
        disp(['Enter range, ' num2str(rng(1)) ' km to ' num2str(rng(nrout)) ' km, to extract data:']);
        if ~exist('rngout','var')
            rngout=input(' ');
        else
            disp(rngout)
        end
        nrngout=find((rng-rngout)>=0);
        if isempty(nrngout)
          if rngout<rng(1)
            nrngout=1;
          else
            nrngout=nrout;
          end
        else
          nrngout=nrngout(1);
        end
        rngout=rng(nrngout);
        disp(['Outputting data for range ' num2str(rng(nrngout)) 'km.']);
        radout=rad;
        disp(' ');
      end % end of 'if optout == 1'
    elseif nrad == 1
      optout=1;
      nradout=1;
      radout=rad;
      nrngout=nrout;
      rngout=rng;
    else
      optout=2;
      nrngout=1;
      rngout=rng;
      nradout=nrad;
      radout=rad;
    end  % end of 'if nrad > 1 && nrout > 1'

    for ifreq=1:nf

% First skip consists of (out to proper freq)+(starting field)+
%                          (out to proper radial/range)

      if optout == 1 % extract along single radial

        skip=(ifreq-1)*(recl + recl*nrad*nrout) + recl + recl*(nradout-1);
        fid=fopen(pefile,'r');
        fseek(fid,(header+skip)*4,0);

        skip=recl*(nrad-1);
        for irng=1:nrout
          data=fread(fid,2*nzout,'float32');
          if irng < nrout && nrad > 1
            fseek(fid,(skip)*4,0);
          end
          psi=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
          apsi=abs(psi); phipsi=angle(psi);
          atmp=interp1(dep',apsi,depout(irng,nradout));phitmp=interp1(dep',phipsi,depout(irng,nradout));
          psid(irng,ifreq)=atmp*cos(phitmp)+1i*atmp*sin(phitmp);
          clear data psi;
        end

      else % extract at single range

        skip=(ifreq-1)*(recl + recl*nrad*nrout) + recl + recl*nrad*(nrngout-1);
        fid=fopen(pefile,'r');
        fseek(fid,(header+skip)*4,0);

        skip=0; % don't need to use it since reading consecutive radials
        for irad=1:nrad
          data=fread(fid,2*nzout,'float32');
          psi=data(1:2:2*nzout-1,1)+1i*data(2:2:2*nzout,1);
          apsi=abs(psi); phipsi=angle(psi);
          atmp=interp1(dep',apsi,depout(nrngout,irad));phitmp=interp1(dep',phipsi,depout(nrngout,irad));
          psid(irad,ifreq)=atmp*cos(phitmp)+1i*atmp*sin(phitmp);
          clear data psi;
        end

      end % end of 'if optout == 1'
      fclose(fid);

    end % end of 'for ifreq=1:nf,'

% Fill in any NaN gaps
    gaps=find(isnan(psid));
    for i=1:size(gaps)
      psid(gaps(i))=1.e-20;
    end  % end of 'for i=1:size(gaps)'

% Include cylindrical spreading
% NOTE:  Don't need to include phase factors since it results
%        only in time-shift relative to reduced time.  However,
%        phase factor does arise due to basebanding and is added
%        after FFT transform.
    r0=1.;
    if optout == 1
      for irng=1:nrout
        pressd(irng,:)=psid(irng,:)*sqrt(r0/max(1000.*rng(irng),r0));
      end
    else
      pressd=psid*sqrt(r0/max(1000.*rng(nrngout),r0));
    end

% Convert to time domain, re-order freqs and correct for basebanded phase
    hanwin=hanning(nf+1);
    hanwin=fftshift(hanwin(1:nf)/sum(hanwin(1:nf)));
    phs=2*pi*cfreq*time;
    if optout == 1
      for irng=1:nrout
        pressd(irng,:)=fftshift(fft(pressd(irng,:).*hanwin')).*exp(-1i*phs);
        timeout(irng,:)=time+rngout(1,irng)/(c0/1000.);
      end
    else
      for irad=1:nrad
        pressd(irad,:)=fftshift(fft(pressd(irad,:).*hanwin')).*exp(-1i*phs);
      end
      tmp=pressd(nrad/2+2:nrad,:);
      pressd(nrad/2:nrad,:)=pressd(1:nrad/2+1,:);pressd(1:nrad/2-1,:)=tmp;
      timeout=time+rngout/(c0/1000.);
    end

    tlpressd=-20*log10(max(abs(pressd),1.e-20));
    tlmin=min(min(tlpressd));tlmax=tlmin+100;

    if doPlots
        if optout == 1 % output single radial

          if nrout == 1 % single range
            figure;plot(timeout,tlpressd);axis('ij');
            if itype==1
              ylabel('TL (dB re 1m) for Pressure');xlabel('Time (sec)');
            elseif itype==2
              ylabel('TL (dB re 1m) for Radial Velocity');xlabel('Time (sec)');
            else
              ylabel('TL (dB re 1m) for Vertical Velocity');xlabel('Time (sec)');
            end
          else
            figure;imagesc(time,rng,tlpressd);caxis([tlmin tlmax]);colormap(flipud(jet));
            set(colorbar,'YDir','reverse');
            if itype==1
              title('Transmission Loss (dB re 1m) for Pressure');xlabel('Reduced Time (sec)');ylabel('Range (km)');
            elseif itype==2
              title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Reduced Time (sec)');ylabel('Range (km)');
            else
              title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Reduced Time (sec)');ylabel('Range (km)');
            end
          end % end of 'if nrout == 1'

        else % output single range

          if nrad == 1 % single radial
            figure;plot(timeout,tlpressd);axis('ij');
            if itype==1
              ylabel('TL (dB re 1m) for Pressure');xlabel('Time (sec)');
            elseif itype==2
              ylabel('TL (dB re 1m) for Radial Velocity');xlabel('Time (sec)');
            else
              ylabel('TL (dB re 1m) for Vertical Velocity');xlabel('Time (sec)');
            end
          else
            figure;pcolor(timeout,rad,tlpressd);caxis([tlmin tlmax]);colormap(flipud(jet));
            set(colorbar,'YDir','reverse');
            if itype==1
              title('Transmission Loss (dB re 1m) for Pressure');xlabel('Time (sec)');ylabel('Radial (deg)');
            elseif itype==2
              title('Transmission Loss (dB re 1m) for Radial Velocity');xlabel('Time (sec)');ylabel('Radial (deg)');
            else
              title('Transmission Loss (dB re 1m) for Vertical Velocity');xlabel('Time (sec)');ylabel('Radial (deg)');
            end
          end % end of 'if nrad == 1'

        end % end of 'if optout == 1'
    end
    
    if ~exist('savedat','var')
        savedat=input('Save data? (y or n) ','s');
    else
        disp(savedat)
    end
    disp(' ');

    if isempty(savedat)
    elseif savedat == 'y' || savedat == 'Y'
      if ~exist('fileout','var')
          fileout=input('Enter matlab output filename (no extension needed): ','s');
      else
          disp(fileout)
      end
      if itype==1
        eval(['save ' fileout ' pressd tlpressd depout time timeout rngout radout;']);
      elseif itype==2
        apvrd=pressd/rhoc0; tlapvrd=tlpressd;
        eval(['save ' fileout ' apvrd tlapvrd depout time timeout rngout radout rhoc0;']);
      else
        apvzd=pressd/rhoc0; tlapvzd=tlpressd;
        eval(['save ' fileout ' apvzd tlapvzd depout time timeout rngout radout rhoc0;']);
      end
    end  % end of 'if savedat == 'y' || savedat == 'Y''

  end  % end of 'if optt == 3'

% -----------------------------------------------------------------------------


end  % end of 'if nf == 1'

end  % end of 'if opt1 == 6'

% =============================================================================

end  % end of 'if ~exist('pefile')'
