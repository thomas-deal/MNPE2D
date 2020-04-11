      PROGRAM PEMP

C  This is the main program of the PE split-step Fourier algorithm.  The
C  FFT convention is such that the PE field function psi(z)=FFT(psi(k)) and
C  psi(k)=conj(FFT(conj(psi(z)))) where FFT(f(x))=sum(dx*f(x)*exp(+i*x*y)).
C  The PE/SSF algorithm is then implemented with a centered step scheme, i.e.
C  psi(z,r+dr)=envprop2(r+dr)*FFT(phsprop*conj(FFT(conj(envprop2(r)*psi(z,r)))))
C  where envprop2=envprop(dr/2) and phsprop are the wide-angle PE propagators.

      REAL*4 psir(32768),psii(32768),psir0(32768),psii0(32768)
      REAL*4 psirs(32768),psiis(32768),depth(32768)
      REAL*4 apvrr(32768),apvri(32768),apvzr(32768),apvzi(32768)
      REAL*4 topkr(32768),topki(32768),filtk(16385)
      REAL*4 phsr(32768),phsi(32768),fk(32768),enz(16384)
      REAL*4 envr(32768),envi(32768),envr2(32768),envi2(32768)
      REAL*4 rngout(10000)
      REAL*4 bdout(10000,512),dbdout(10000,512)
      REAL*4 sdout(10000,512),filt(16385)
      REAL*4 depmin,depmax,rng,rngmin,rngmax,cfreq,freqbw,freq,zero,c0
      REAL*4 bdint,dbdint,Hdep
      INTEGER*4 nz,nzout,nrout,nf,lrec,ifiltyp
      CHARACTER*40 srcdata,sspdata,surfdata,botbathdata,botprofdata
      CHARACTER*40 seddepdata,sedprofdata,outpdata,outpvrdata,outpvzdata
      CHARACTER*80 junk

      COMMON nz,dz,dk,dr,rng,freq,c0,sd,arl,thc,depcalc,bdint,dbdint,
     &irad,nrad,drad,rngmax,rsint,rsdotint,Hdep

      open(9,file='mmpe2d.log',status='unknown')
C11    format(1x,'START DATE : ',i2,'/',i2,'/',i4)
C12    format(1x,'START TIME : ',i2,':',i2,':',i2,'.',i2)
C13    format(1x,'STOP DATE : ',i2,'/',i2,'/',i4)
C14    format(1x,'STOP TIME : ',i2,':',i2,':',i2,'.',i2)
C      CALL GETDAT(iyr,imon,iday)
C      CALL GETTIM(ihr,imin,isec,isec100)
C      write(*,11)imon,iday,iyr
C      write(9,11)imon,iday,iyr
C      write(*,12)ihr,imin,isec,isec100
C      write(9,12)ihr,imin,isec,isec100


20    format(33x,a40)
21    format(33x,a80)
22    format(a80)
23    format(27x,a80)

C     Read file names for initialization and range/depth sections for data output
      open(10,file='pefiles.inp',status='old')
      read(10,20,err=999)srcdata
      read(10,20,err=999)sspdata
      read(10,20,err=999)surfdata
      read(10,20,err=999)botbathdata
      read(10,20,err=999)botprofdata
      read(10,20,err=999)seddepdata
      read(10,20,err=999)sedprofdata
      read(10,20,err=999)outpdata
      read(10,20,err=999)outpvrdata
      read(10,20,err=999)outpvzdata
      read(10,21,err=999)junk
      open(11,file='tmpjunk.inp',status='new')
      write(11,22)junk
      close(11)
      open(11,file='tmpjunk.inp',status='old')
      read(11,*,err=999)nzout,depmin,depmax
      close(11,status='delete')
      read(10,21,err=999)junk
      open(11,file='tmpjunk.inp',status='new')
      write(11,22)junk
      close(11)
      open(11,file='tmpjunk.inp',status='old')
      read(11,*,err=999)nrout,rngmin,rngmax
      rngmin=rngmin*1000.
      rngmax=rngmax*1000.
      nrout=min(nrout,10000)
      close(11,status='delete')
C     Allow other parameters to be defined in this file
      read(10,21,end=31)junk
      open(11,file='tmpjunk.inp',status='new')
      write(11,22)junk
      close(11)
      open(11,file='tmpjunk.inp',status='old')
      read(11,*,end=30)nzt,drt,depcalc,c0
30    close(11,status='delete')
31    close(10)
      drt=drt*1000.
      if(c0.eq.0)c0=1500.

C     Read source paramter data
      open(10,file=srcdata,status='old')
      read(10,23,err=999)junk
      open(11,file='tmpjunk.inp',status='new')
      write(11,22)junk
      close(11)
      open(11,file='tmpjunk.inp',status='old')
      read(11,*,err=999)sd
      close(11,status='delete')
      read(10,23,err=999)junk
      open(11,file='tmpjunk.inp',status='new')
      write(11,22)junk
      close(11)
      open(11,file='tmpjunk.inp',status='old')
      read(11,*,err=999)arl
      close(11,status='delete')
      read(10,23,err=999)junk
      open(11,file='tmpjunk.inp',status='new')
      write(11,22)junk
      close(11)
      open(11,file='tmpjunk.inp',status='old')
      read(11,*,err=999)thc
      close(11,status='delete')
      read(10,23,err=999)junk
      open(11,file='tmpjunk.inp',status='new')
      write(11,22)junk
      close(11)
      open(11,file='tmpjunk.inp',status='old')
      read(11,*,err=999)cfreq
      close(11,status='delete')
      read(10,23,err=999)junk
      open(11,file='tmpjunk.inp',status='new')
      write(11,22)junk
      close(11)
      open(11,file='tmpjunk.inp',status='old')
      read(11,*,err=999)freqbw
      close(11,status='delete')
      read(10,23,err=999)junk
      open(11,file='tmpjunk.inp',status='new')
      write(11,22)junk
      close(11)
      open(11,file='tmpjunk.inp',status='old')
      read(11,*,err=999)nf
      close(11,status='delete')
      close(10)

      pi=acos(-1.)
      tpi=2.*pi
      fk0=tpi*cfreq/c0
      wvln=c0/cfreq
      slopmax=0.

      if(nzt.eq.0.)then
        ieff=0
        print*,'You have chosen to use default sampling in depth.'
        print*,' '
        print*,'Enter 1 for accuracy, 2 for efficiency:'
        read*,ieff
        if(ieff.ne.1.and.ieff.ne.2)then
          print*,'Invalid input.  Sampling will be based on efficiency.'
          ieff=2
        end if
        if(ieff.eq.1)then
            dz=wvln/10.
        else
            dz=wvln/2.
        end if
        nz=0
        else
          nz=nzt
        end if

        if(drt.eq.0.)then
          if(ieff.eq.1)then
            print*,'Range sampling will be based on accuracy.'
        elseif(ieff.eq.2)then
            print*,'Range sampling will be based on efficiency.'
        elseif(ieff.eq.0)then
            print*,'You have chosen to use default sampling in range.'
            print*,' '
            print*,'Enter 1 for accuracy, 2 for efficiency:'
            read*,ieff
            if(ieff.ne.1.and.ieff.ne.2)then
          print*,'Invalid input.  Sampling will be based on efficiency.'
          ieff=2
            end if
        end if
        if(ieff.eq.1)then
            dr=wvln
        else
            dr=3.*wvln
        end if
      else
       dr=drt
      end if

      if(nf.gt.1)then
        lnf=int(alog(float(nf))/alog(2.)+0.99)
        nf=2**lnf
      else
        lnf=0
      end if
      if(lnf.gt.13)then
        print*,'Error in number of freqs; exceeding max limit.'
        goto 999
      end if

      istrt=0
C     Define frequency of calculation.
C     Note:  Ordering of freqs follows Matlab's FFT convention.
C	itstart reinitializes after crashed run
      itstart=1
      DO ifreq=itstart,nf
      if(nf.eq.1)then
        freq=cfreq
      else
        if(ifreq.le.nf/2)then
          freq=cfreq+float(ifreq-1)*(freqbw/float(nf-1))
        else
          freq=cfreq-float(nf-ifreq+1)*(freqbw/float(nf-1))
        end if
      end if
      print*,' '
      print*,' '
      print*,'FREQ=',freq,'Hz'
      print*,' '

C     Read environmental data and initialize parameters
      if(ifreq.eq.itstart)then
        open(10,file=sspdata,status='old')
        open(11,file=botbathdata,status='old')
        open(12,file=botprofdata,status='old')
        open(13,file=seddepdata,status='old')
        open(14,file=sedprofdata,status='old')
        open(15,file=surfdata,status='old')
      end if
      rng=0.0
      irad=1
      CALL ENVPROP(envr,envi,envr2,envi2,enz,filt,istrt,ird)
      slopmax=max(slopmax,rsdotint)
      if(ifreq.eq.itstart)then
      if(istrt.eq.2)goto 999
        close(10)
        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
      end if
CCCCCCCCCCCCCCCCCCCCCC
C	open(21,file='envprop.dat')
C	do ii=nz/2+1,nz
C	depth=-(float(nz-ii)+0.5)*dz
C	iii=nz-ii+1
C	write(21,*)envr(iii),envi(iii),envr2(iii),envi2(iii)
C	end do
C	do ii=1,nz/2
C	depth=(float(ii)-.5)*dz
C	write(21,*)envr(ii),envi(ii),envr2(ii),envi2(ii)
C	end do
C	close(21)
CCCCCCCCCCCCCCCCCCCCCC

      if(ifreq.eq.itstart)then
        lnz=int(alog(float(nz))/alog(2.)+0.99)
        nz=2**lnz
        nz2=nz/2
        if(nrad.gt.1)then
          lnrad=int(alog(float(nrad))/alog(2.)+0.99)
          nrad=2**lnrad
        else
          lnrad=0
        end if
        if(lnz.gt.15.or.lnrad.gt.9)then
          print*,'Error in array size; exceeding max limit.'
          goto 999
        end if
      end if

C	Generate depth mesh vector
      DO iz=1,nz
        if(iz.le.nz/2)then
          depth(iz)=(float(iz-1)+0.5)*dz
        else
          depth(iz)=-(float(nz-iz)+0.5)*dz
        end if
      END DO

C     Generate phase space propagator (conjugate form)
      CALL PHSPROP(phsr,phsi,topkr,topki,fk)
C     Generate wavenumber filter for vertical velocity
      dk=pi/depcalc
      nz2=nz/2
      nz2p1=nz2+1
      fk0=2.*pi*freq/c0
      nzmax=int(fk0/dk)
      nz10=nzmax/30
      do iz=1,nz2p1
        fj=float(iz-(nzmax-nz10))/float(nz10)
        arg=(1.-amin1(1.,amax1(0.,fj)))*pi
        cosa=cos(arg)
        filtk(iz)= 0.5*(1.0-cosa)
      end do
CCCCCCCCCCCCCCCCCCCCCC
C	open(21,file='phsprop.dat')
C	do ii=nz/2+2,nz
C	write(21,*)phsr(ii),phsi(ii),topkr(ii),topki(ii)
C	end do
C	do ii=1,nz/2+1
C	write(21,*)phsr(ii),phsi(ii),topkr(ii),topki(ii)
C	end do
C	close(21)
CCCCCCCCCCCCCCCCCCCCCC

C     Generate starting field (result is true k-domain field)
      CALL PESRC(psir0,psii0,fk)
CCCCCCCCCCCCCCCCCCCCCC
C	open(21,file='pesrck.dat')
C	do ii=nz/2+2,nz
C	write(21,*)fk(ii),psir0(ii),psii0(ii)
C	end do
C	do ii=1,nz/2+1
C	write(21,*)fk(ii),psir0(ii),psii0(ii)
C	end do
C	close(21)
CCCCCCCCCCCCCCCCCCCCCC

C     Transform to z-domain (result is true z-domain field)
      CALL FFT(lnz,psir0,psii0)
C
      if((rsint.ne.0.).or.(rsdotint.ne.0.))then
C	Apply phase correction for rough surface image field.
        DO iz=1,nz
          retasdot=2.0*fk0*rsdotint*(depth(iz)-rsint)
          cms=cos(retasdot)
          sms=sin(retasdot)
          psirs(iz)=cms*psir0(iz)-sms*psii0(iz)
          psiis(iz)=sms*psir0(iz)+cms*psii0(iz)
        END DO

        if(rsint.gt.0.)then
          ns=int(rsint/dz+0.5)
        else
          ns=int(rsint/dz-0.5)
        end if
C	Correct field for rough surface effect
        if(ns.gt.0)then
          do iz=1,ns
            psir0(iz)=psirs(iz)
            psii0(iz)=psiis(iz)
          end do
        end if
        ns=amin0(0,ns)
        do iz=nz/2+1,nz+ns
          psir0(iz)=psirs(iz)
          psii0(iz)=psiis(iz)
        end do
      end if
C
C  Output at this point is tilde field (equiv to normal field if rsint=rsdotint=0)
C
CCCCCCCCCCCCCCCCCCCCCC
C	open(21,file='pesrcz.dat')
C	do ii=nz/2+1,nz
C	depth=-(float(nz-ii)+0.5)*dz
C	write(21,*)depth,psir0(ii),psii0(ii)
C	end do
C	do ii=1,nz
C	depth=(float(ii)-.5)*dz
C	write(21,*)depth,psir0(ii),psii0(ii)
C	end do
C	close(21)
CCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCC
      fk0=2.*pi*freq/c0
C	Compute acoustic particle velocity field for output
C	Store conjugate
      DO iz=1,nz
        apvrr(iz)=psir0(iz)
        apvri(iz)=-psii0(iz)
      END DO
C     Transform to Fourier domain (result is conjugate k-domain field)
      CALL FFT(lnz,apvrr,apvri)
C     Multiply conjugate by k-space operators (result is true k-domain field)
      DO iz=1,nz
        apvzr(iz)=(fk(iz)*apvrr(iz))/float(nz)
        apvzi(iz)=(-fk(iz)*apvri(iz))/float(nz)
        tmp=(topkr(iz)*apvrr(iz)+topki(iz)*apvri(iz))/float(nz)
        apvri(iz)=(-topkr(iz)*apvri(iz)+topki(iz)*apvrr(iz))/float(nz)
        apvrr(iz)=tmp
      END DO
C     Transform to physical space domain (result is true z-domain field)
      CALL FFT(lnz,apvrr,apvri)
      CALL FFT(lnz,apvzr,apvzi)
C	Compute additional radial component factors
      DO iz=1,nz2
        apvrr(iz)=enz(iz)*psir0(iz)-1./(2.*fk0)*psii0(iz)-apvrr(iz)
        apvri(iz)=enz(iz)*psii0(iz)+1./(2.*fk0)*psir0(iz)-apvri(iz)
      END DO
CCCCCCCCCCCCCCCCCCCCCC

      if(ifreq.eq.itstart)then
C     Output header data at initial range
        zero=0.0
        rng=0.0
        izmin=max(int(depmin/dz+0.5),1)
        izmax=min(int(depmax/dz+0.5),nz2)
        nzout=min(nzout,izmax-izmin+1)
        izskip=int(float(izmax-izmin)/float(nzout))+1
        nzout=int(float(izmax-izmin)/float(izskip))+1
        izmax=izmin+(nzout-1)*izskip
        depmin=(float(izmin)-0.5)*dz
        depmax=(float(izmax)-0.5)*dz
        nrout=min(nrout,int((rngmax-rngmin)/dr+1))
        if(nrout.gt.1)then
          drout=(rngmax-rngmin)/float(nrout-1)
          if(drout.eq.dr.and.rngmin.eq.0)then
            nrout=nrout-1
            rngmin=dr
          end if
          do irout=1,nrout
            rngout(irout)=rngmin+float(irout-1)*drout
          end do
        else
          rngout(1)=rngmax
        end if
        do irad=1,nrad
          bdout(1,irad)=bdint
          dbdout(1,irad)=dbdint
          sdout(1,irad)=rsint
        end do

      open(unit=10,file=outpdata,status='unknown',access='direct',
     & err=999,form='unformatted',recl=4*2*nzout)
      open(unit=11,file=outpvrdata,status='unknown',access='direct',
     & err=999,form='unformatted',recl=4*2*nzout)
      open(unit=12,file=outpvzdata,status='unknown',access='direct',
     & err=999,form='unformatted',recl=4*2*nzout)

      if(ifreq.eq.1)then
        write(9,*)' '
        write(9,*)'FFT size, nz = ',nz,' and depth mesh, dz = ',dz,' m.'
        write(9,*)'Range step, dr = ',(dr/1000.),' km.'
        write(9,*)'Outputting ',nrout,' points from ',(rngmin/1000.),
     &' km to ',(rngmax/1000.),' km in range.'
        write(9,*)'Outputting ',nzout,' points from ',depmin,
     &' m to ',depmax,' m in depth.'
        write(9,*)'Outputting ',nrad,' radials separated by ',
     &drad,' deg.'
        write(9,*)' '

      ifiltyp=1
      write(10,rec=1)zero,c0,nf,cfreq,freqbw,nrout,rngmin,rngmax,dr,
     & nzout,depmin,depmax,nrad,drad,bdint,dbdint,sd,ifiltyp,rsint
      ifiltyp=2
      write(11,rec=1)zero,c0,nf,cfreq,freqbw,nrout,rngmin,rngmax,dr,
     & nzout,depmin,depmax,nrad,drad,bdint,dbdint,sd,ifiltyp,rsint
      ifiltyp=3
      write(12,rec=1)zero,c0,nf,cfreq,freqbw,nrout,rngmin,rngmax,dr,
     & nzout,depmin,depmax,nrad,drad,bdint,dbdint,sd,ifiltyp,rsint
      end if

      end if

C	If itstart <> 1, program is recovering, needs to know correct lrec
C	Last record included header + all previous freq data (starting field and
C	nrad*nrout records)
      lrec=1 + (ifreq-1)*(1 + nrad*nrout)

      write(10,rec=lrec+1)(psir0(iz),psii0(iz),iz=izmin,izmax,izskip)
      write(11,rec=lrec+1)(apvrr(iz),apvri(iz),iz=izmin,izmax,izskip)
      write(12,rec=lrec+1)(apvzr(iz),apvzi(iz),iz=izmin,izmax,izskip)

C     Multiply by 1/2 z-space propagator, exp(-i*0.5*dr*k0*Uop(z)), enforce
C     symmetry, and conjugate (result is conjugate z-domain field)
C     NOTE:  At zero range, z-space prop is same for all radials in 3-D
      DO iz=1,nz
        tmp=envr2(iz)*psir0(iz)-envi2(iz)*psii0(iz)
        psii0(iz)=-(envr2(iz)*psii0(iz)+envi2(iz)*psir0(iz))
        psir0(iz)=tmp
      END DO

C     Begin radial loop
      DO irad=1,nrad
      if(nrad.gt.1)then
        if(amod(float(irad),float(max(2**(lnrad-4),1))).eq.0.)
     &  write(6,100)irad
      end if
100   format(1x,'Radial # ',i4)

C     Copy starting field for 3-D calculations
      DO iz=1,nz
        psir(iz)=psir0(iz)
        psii(iz)=psii0(iz)
      END DO

C     For each range step, do...
      irout=1
      rng=0.
      DO WHILE (rng.lt.rngmax)
        if((amod(float(int(rng/dr)),float(max(1,int((rngmax/dr)/20.))))
     &.eq.0.).and.(nrad.eq.1))write(6,101)rng/1000.
101       format(1x,' Range = ',f12.4)

C     Update range
        rng=rng+dr

C	Convert back to non-tilde field for k-space operation
      if((rsint.ne.0.).or.(rsdotint.ne.0.))then
C	Apply phase correction for rough surface image field.
        if(rsint.gt.0.)then
          ns=int(rsint/dz+0.5)
        else
          ns=int(rsint/dz-0.5)
        end if
C	Correct field for rough surface effect
        if(ns.gt.0)then
          do iz=1,ns
            retasdot=-2.0*fk0*rsdotint*(depth(iz)-rsint)
            cms=cos(retasdot)
            sms=sin(retasdot)
            tmp=cms*psir(iz)-sms*psii(iz)
            psii(iz)=sms*psir(iz)+cms*psii(iz)
            psir(iz)=tmp
          end do
        end if
        ns=amin0(0,ns)
        do iz=nz/2+1,nz+ns
          retasdot=-2.0*fk0*rsdotint*(depth(iz)-rsint)
          cms=cos(retasdot)
          sms=sin(retasdot)
          tmp=cms*psir(iz)-sms*psii(iz)
          psii(iz)=sms*psir(iz)+cms*psii(iz)
          psir(iz)=tmp
        end do
      end if
C
C  Output at this point is normal field (equiv to tilde field if rsint=rsdotint=0)
C     Transform to Fourier domain (result is conjugate (k,s)-domain field)
      CALL FFT(lnz,psir,psii)

C     Multiply conjugate by k-space propagator, exp(-i*dr*k0*Top(k))
C     (result is true k-domain field)
      DO iz=1,nz
        tmp=phsr(iz)*psir(iz)+phsi(iz)*psii(iz)
        psii(iz)=-phsr(iz)*psii(iz)+phsi(iz)*psir(iz)
        psir(iz)=tmp
      END DO
C
C	Impose symmetry for exact rough surface scatter
      if((rsint.ne.0.).or.(rsdotint.ne.0.))then
        psirs(1)=psir(1)
        psiis(1)=psii(1)
        DO iz=2,nz
          retas=2.0*rsint*fk0*fk(iz)
          jz=nz-iz+2
          cks=cos(retas)
          sks=sin(retas)
          psirs(jz)=psir(iz)*cks-psii(iz)*sks
          psiis(jz)=psii(iz)*cks+psir(iz)*sks
        END DO
      end if
C

C     Transform to physical space domain (result is true z-domain field)
      CALL FFT(lnz,psir,psii)
C
      if((rsint.ne.0.).or.(rsdotint.ne.0.))then
        CALL FFT(lnz,psirs,psiis)
C	Apply phase correction for exact surface scatter symmetry
        DO iz=1,nz
          retasdot=2.0*fk0*rsdotint*(depth(iz)-rsint)
          cms=cos(retasdot)
          sms=sin(retasdot)
          tmp=-(cms*psirs(iz)-sms*psiis(iz))
          psiis(iz)=-(sms*psirs(iz)+cms*psiis(iz))
          psirs(iz)=tmp
        END DO

        if(rsint.gt.0.)then
          ns=int(rsint/dz+0.5)
        else
          ns=int(rsint/dz-0.5)
        end if
C	Correct field exactly for rough surface symmetry
        if(ns.gt.0)then
          do iz=1,ns
            psir(iz)=psirs(iz)
            psii(iz)=psiis(iz)
          end do
        end if
        ns=amin0(0,ns)
        do iz=nz/2+1,nz+ns
          psir(iz)=psirs(iz)
          psii(iz)=psiis(iz)
        end do
      end if
C

C     Read current environmental propagator
      if(ird.ne.1)then
        CALL ENVPROP(envr,envi,envr2,envi2,enz,filt,istrt,ird)
        slopmax=max(slopmax,rsdotint)
      end if

      IF(rng.ge.rngout(irout))THEN
C     Output psi data at current range if requested

C     Multiply by 1/2 z-space propagator, exp(-i*0.5*dr*k0*Uop(z)) (result is true z-domain field)
      DO iz=1,nz
        tmp=envr2(iz)*psir(iz)-envi2(iz)*psii(iz)
        psii(iz)=envr2(iz)*psii(iz)+envi2(iz)*psir(iz)
        psir(iz)=tmp
      END DO

      bdout(irout+1,irad)=bdint
      dbdout(irout+1,irad)=dbdint
      sdout(irout+1,irad)=rsint

CCCCCCCCCCCCCCCCCCCCCC
C	Compute acoustic particle velocity field for output
C	Store conjugate
      DO iz=1,nz
        apvrr(iz)=psir(iz)
        apvri(iz)=-psii(iz)
      END DO
C     Transform to Fourier domain (result is conjugate k-domain field)
      CALL FFT(lnz,apvrr,apvri)
C     Multiply conjugate by k-space operators (result is true k-domain field)
      DO iz=1,nz
        jz=min0(iz+1,nz+2-iz)
        apvzr(iz)=filtk(jz)*(fk(iz)*apvrr(iz))/float(nz)
        apvzi(iz)=filtk(jz)*(-fk(iz)*apvri(iz))/float(nz)
        tmp=(topkr(iz)*apvrr(iz)+topki(iz)*apvri(iz))/float(nz)
        apvri(iz)=(-topkr(iz)*apvri(iz)+topki(iz)*apvrr(iz))/float(nz)
        apvrr(iz)=tmp
      END DO
C     Transform to physical space domain (result is true z-domain field)
      CALL FFT(lnz,apvrr,apvri)
      CALL FFT(lnz,apvzr,apvzi)
C	Compute additional radial component factors
      DO iz=1,nz2
        apvrr(iz)=enz(iz)*psir(iz)-1./(2.*fk0*rng)*psii(iz)-apvrr(iz)
        apvri(iz)=enz(iz)*psii(iz)+1./(2.*fk0*rng)*psir(iz)-apvri(iz)
      END DO
CCCCCCCCCCCCCCCCCCCCCC

C	Convert tilde field to standard space
C	Note, at this time, no need to convert since only field below surface will be saved and they are equivalent there

C	Last record included header + all previous freq data (starting field and
C	nrad*nrout records) + current freq starting field + previous range's radials
C	+ previous radial
      lrec=1+(ifreq-1)*(1+nrad*nrout)+1+(irout-1)*nrad+(irad-1)
        write(10,rec=lrec+1)(psir(iz),psii(iz),
     &                         iz=izmin,izmax,izskip)
        write(11,rec=lrec+1)(apvrr(iz),apvri(iz),
     &                         iz=izmin,izmax,izskip)
        write(12,rec=lrec+1)(apvzr(iz),apvzi(iz),
     &                         iz=izmin,izmax,izskip)
        irout=irout+1

C     Multiply by 1/2 z-space propagator, exp(-i*0.5*dr*k0*Uop(z)) and
C     conjugate (result is conjugate z-domain field)
        DO iz=1,nz
          tmp=envr2(iz)*psir(iz)-envi2(iz)*psii(iz)
          psii(iz)=-(envr2(iz)*psii(iz)+
     &                     envi2(iz)*psir(iz))
          psir(iz)=tmp
        END DO

       ELSE

C     Multiply by z-space propagator, exp(-i*dr*k0*Uop(z)) and
C     conjugate (result is conjugate z-domain field)
        DO iz=1,nz
          tmp=envr(iz)*psir(iz)-envi(iz)*psii(iz)
          psii(iz)=-(envr(iz)*psii(iz)+envi(iz)*psir(iz))
          psir(iz)=tmp
        END DO

C  Output at this point is tilde field (equiv to normal field if rsint=rsdotint=0)

       END IF

C     End of main range loop
      END DO

C     End of radial loop
      END DO

C     End of frequency loop
      END DO

C     Output bathymetry to end of file and size of data file to first record
      close(10)
      close(11)
      close(12)
      open(unit=10,file=outpdata,status='unknown',access='direct',
     & err=999,form='unformatted',recl=4)
      open(unit=11,file=outpvrdata,status='unknown',access='direct',
     & err=999,form='unformatted',recl=4)
      open(unit=12,file=outpvzdata,status='unknown',access='direct',
     & err=999,form='unformatted',recl=4)

C	Total number of records output includes header + all freq * all ranges * all radials
C	Convert into bytes for final output
      nrec=1+nf*(1+nrad*nrout)
      lrec=nrec*2*nzout
      DO irad=1,nrad
      do ir=1,nrout+1
        write(10,rec=lrec+1)bdout(ir,irad)
        write(11,rec=lrec+1)bdout(ir,irad)
        write(12,rec=lrec+1)bdout(ir,irad)
        lrec=lrec+1
      end do
      do ir=1,nrout+1
        write(10,rec=lrec+1)dbdout(ir,irad)
        write(11,rec=lrec+1)dbdout(ir,irad)
        write(12,rec=lrec+1)dbdout(ir,irad)
        lrec=lrec+1
      end do
      do ir=1,nrout+1
        write(10,rec=lrec+1)sdout(ir,irad)
        write(11,rec=lrec+1)sdout(ir,irad)
        write(12,rec=lrec+1)sdout(ir,irad)
        lrec=lrec+1
      end do
      END DO
      write(10,rec=1)lrec
      write(11,rec=1)lrec
      write(12,rec=1)lrec
      close(10)
      close(11)
      close(12)

C	Evaluate maximum surface slope and print warning if exceeded 20deg
      if(slopmax.gt.0.35)then
        print*,'WARNING: SURFACE SLOPE EXCEEDED MAX RECOMMENDED 20DEG'
        open(9,file='surfslopwarn.log',status='unknown')
        write(9,*)'Max surface slope = ',slopmax*180/pi,' deg.'
        write(9,*)'Surface slopes greater than 20 deg'
        write(9,*)'  may introduce errors into solution.'
        close(9)
      end if

      goto 1000
999   print*,'ERROR IN RUN'

1000  continue
C1000  CALL GETDAT(iyr,imon,iday)
C      CALL GETTIM(ihr,imin,isec,isec100)
C      write(*,13)imon,iday,iyr
C      write(9,13)imon,iday,iyr
C      write(*,14)ihr,imin,isec,isec100
C      write(9,14)ihr,imin,isec,isec100
      close(9)

      STOP
      END
