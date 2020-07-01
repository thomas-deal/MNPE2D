      SUBROUTINE ZSGEN(etas,detas,d2etas,sigma,istrt)
C
C     This subroutine generates a rough surface realization.
C
      REAL*4 etas(1),detas(1),d2etas(1)
      REAL*4 fk(8192),filter(8192),omega(8192),domegadk(8192)
      REAL*4 etasi(8192),detasi(8192),d2etasi(8192)
      INTEGER*4 rstyp,nspec
      REAL*4 rmsrough,corrl,wspd,fetch,fspec(101),sspec(101)
      REAL*4 w1(8192)
      REAL*4 alpha,beta,omega0,g
      REAL*4 delta,sigma0,math1,math2,gamma0,tmpf,tmpw1,tmpw2

      integer n,seed,wvdir
      integer, dimension(:), allocatable :: iseed
      REAL*4 ran(4),mag,phs,ranr,rani,surftime,wt,cwt,swt,wp,wn,tide
C
      COMMON nz,dz,dk,dr,rng,freq,c0,sd,arl,thc,depcalc,bdint,dbdint,
     &irad,nrad,drad,rngmax,rsint,rsdotint,Hdep

C     Constants/Ranges
      nrng=rngmax/dr+1
      lnr=int(alog(float(nrng))/alog(2.)+0.99)
      if(lnr.gt.13)print*,'WARNING: lnr > 13'
	  lnr=min(lnr,13)
      nrng=2.**lnr
      rngmax2=(nrng-1)*dr

      pi=acos(-1.)
      dkr=2.*pi/rngmax2
      fkmax=nrng*dkr/2.
      Lmin=2.*pi/fkmax
      g=9.806

C     Filter parameters
      filtlow=2*pi/20000.
      filthigh=2*pi/2.
      iflow=0
      ifhigh=nrng/2

C     Define wavenumbers
      do ir=1,nrng/2+1
        fk(ir)=(ir-1)*dkr
        omega(ir)=sqrt(g*fk(ir)*tanh(fk(ir)*Hdep))
        if(ir.eq.1)then
          domegadk(ir)=0.
        else
        domegadk(ir)=g*tanh(fk(ir)*Hdep)
        domegadk(ir)=domegadk(ir)+g*fk(ir)*Hdep/cosh(fk(ir)*Hdep)**2.
        domegadk(ir)=domegadk(ir)/(2*sqrt(g*fk(ir)*tanh(fk(ir)*Hdep)))
        end if
        if(fk(ir).gt.filtlow.and.iflow.eq.0)then
          iflow=ir
        end if
        if(fk(ir).gt.filthigh.and.ifhigh.eq.0)then
          ifhigh=ir-1
        end if
        filter(ir)=1.
      end do

C     Define filter
      do ir=1,iflow
         filter(ir)=0.5*(sin(-pi/2+pi*(ir-1)/iflow)+1)
      end do
      do ir=ifhigh,nrng/2+1
         filter(ir)=0.5*(sin(pi/2-pi*(ir-ifhigh)/(nrng/2-ifhigh+1))+1)
      end do

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Read surface spectrum file data
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      write(9,*)'Rough Surface Parameters'
C Read spectrum type code
      read(15,*)rstyp
      if(rstyp.lt.-1.or.rstyp.gt.4)then
        print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: RSTYP'
        istrt=2
        goto 9999
      end if

C Read tide height relative to MLLW
      read(15,*)tide

C Case -1: Sinusoidal Surface
      if(rstyp.eq.-1)then
        write(9,*)'Sinusoidal Surface'
        read(15,*)rmsrough,corrl
        write(9,*)'Amp = ',rmsrough,' SurfWvln = ',corrl
        fk0s=2.*pi/corrl
        do ir=1,nrng
            etas(ir)=-tide+rmsrough*(sin(fk0s*ir*dr))
            detas(ir)=fk0s*rmsrough*cos(fk0s*ir*dr)
            d2etas(ir)=-fk0s**2.*rmsrough*sin(fk0s*ir*dr)
        end do
        sigma=1.
        goto 9999
      end if

C Case 0: Flat Surface
      if(rstyp.eq.0)then
        write(9,*)'Flat Surface'
        do ir=1,nrng
            etas(ir)=-tide
            detas(ir)=0.
            d2etas(ir)=0.
        end do
        if(tide.eq.0)then
          sigma=0.
        else
          sigma=1.
        end if
        goto 9999
      end if

C Read random number generator seed and time for all other surface types
      read(15,*)seed,surftime,wvdir
      write(9,*)'Random Seed = ',seed
      write(9,*)'Time Offset = ',surftime
      write(9,*)'Direction = ',wvdir
      n=12
      allocate(iseed(n))
      call random_seed(size = n)
      iseed = seed + 37 * [(i, i = 0,n-1)]
      call random_seed(put=iseed)

C Case 1: RMS Roughness and Correlation Length
      if(rstyp.eq.1)then
        write(9,*)'RMS Roughness Model'
        read(15,*)rmsrough,corrl
        write(9,*)'RMS = ',rmsrough,' Lcorr = ',corrl
        do ir=1,nrng/2+1
            w1(ir)=(1+((corrl**2)*(fk(ir)**2)))**(-(rmsrough/2)+0.5)
        end do
        sigma=1.

C Case 2: Pierson-Moskowitz Spectrum
      elseif(rstyp.eq.2)then
        write(9,*)'Pierson-Moskowitz Spectrum'
        read(15,*)wspd
        write(9,*)'Wind Speed = ',wspd
        alpha=0.0081
        beta=0.74
        omega0=g/wspd
        do ir=1,nrng/2+1
          if(ir.gt.1)then
            w1(ir)=abs(alpha*g**2./omega(ir)**5.)
            w1(ir)=w1(ir)*exp(-beta*omega0**4./omega(ir)**4.)
            w1(ir)=w1(ir)*domegadk(ir)
          else
            w1(ir)=0.
          end if
        end do
        sigma=1.

C Case 3: JONSWAP Spectrum
      elseif(rstyp.eq.3)then
        write(9,*)'JONSWAP Spectrum'
        read(15,*)wspd,fetch
        write(9,*)'Wind Speed = ',wspd,' Fetch = ',fetch
C     Compute JONSWAP coefficients
        alpha=0.076*(g*fetch/wspd)**(-0.22)
        omega0=7.*pi*(g/wspd)*(g*fetch/wspd**2.)**(-0.33)
        gamma0=3.3
        do ir=1,nrng/2+1
C     Determine value for sigma0
            if(omega(ir).gt.omega0) then
                sigma0=0.09
            else
                sigma0=0.07
            end if
C     delta computation
            math1=(omega(ir)-omega0)**2.
            math2=2.*(sigma0**2.)*(omega0**2.)
            delta=-math1/math2
            if(delta.lt.-50.)then
               delta=0.
            else
               delta=exp(-math1/math2)
            end if
C     w1 calculation
            if (ir.gt.1) then
                math1=abs(alpha*g**2./omega(ir)**5.)
                math2=-1.25*omega0**4./omega(ir)**4.
                if(math2.lt.-50.)then
                   math2=0.
                else
                   math2=exp(math2)
                end if
                w1(ir)=math1*math2*gamma0**delta
                w1(ir)=w1(ir)*domegadk(ir)
            else
                w1(ir)=0.
            end if
        end do
        sigma=1.

C Case 4: User-Defined Spectrum
      elseif(rstyp.eq.4)then
        write(9,*)'User-Defined Spectrum'
        read(15,*)nspec
        if(nspec.le.0.or.nspec.gt.100)then
            print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: NSPEC'
            istrt=2
            goto 9999
        end if
C Read in spectrum
        DO j=1,nspec
          read(15,*)fspec(j),sspec(j)
          if(j.ne.1)then
            if(fspec(j).le.fspec(j-1))then
             print*,'ERROR IN ENVIRONMENTAL INPUT: FSPEC(J)<=FSPEC(J-1)'
              istrt=2
              goto 9999
            end if
          end if
        END DO
C Interpolate spectrum onto fk vector - see ssi.f
        j=2
        do ir=1,nrng/2+1
            tmpf=1./(2.*PI)*omega(ir)
400         if(tmpf.le.fspec(j).or.j.ge.nspec) goto 401
            j=j+1
            goto 400
401         frac=(tmpf-fspec(j-1))/(fspec(j)-fspec(j-1))
            w1(ir)=sspec(j-1)+frac*(sspec(j)-sspec(j-1))
            w1(ir)=w1(ir)*domegadk(ir)/(2.*PI)
        end do
        tmpw1 = w1(1)
        do ir=2,nrng/2
            tmpw2=0.25*(w1(ir-1)+w1(ir)+w1(ir)+w1(ir+1))
            w1(ir-1)=tmpw1
            tmpw1=tmpw2
        end do
        w1(nrng/2)=tmpw1
        sigma=1.

      end if

C Write spectrum to log file
      do ir=1,nrng/2+1
        write(9,*)fk(ir),w1(ir)
      end do

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C One-sided spectrum has been saved in w1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Compute realization from spectrum
      do ir=2,nrng/2
C Gaussian random variables
      CALL RANDOM_NUMBER(ran)
      mag=(-2.0*log(ran(1)))**0.5
      phs = 2.0*PI*ran(2)
      ranr= mag*sin(phs)/sqrt(2.)
      mag=(-2.0*log(ran(3)))**0.5
      phs = 2.0*PI*ran(4)
      rani= mag*sin(phs)/sqrt(2.)
C Time dependence
      if(surftime.ne.0.)then
        wt = -omega(ir)*surftime
        cwt = cos(wt)
        swt = sin(wt)
      else
        cwt = 1.
        swt = 0.
      end if
C Positive wavenumbers
      if(wvdir.lt.0)then
        wp=0.
        wn=sqrt(w1(ir)*dkr)
      elseif(wvdir.eq.0)then
        wp=sqrt(w1(ir)*dkr/2.)
        wn=wp
      else
        wp=sqrt(w1(ir)*dkr)
        wn=0.
      end if
      etas(ir)=1./sqrt(2.)*(wp+wn)*(ranr*cwt-rani*swt)*filter(ir)
      etasi(ir)=wp*(ranr*swt+rani*cwt)-wn*(ranr*swt+rani*cwt)
      etasi(ir)=1./sqrt(2.)*etasi(ir)*filter(ir)
      detas(ir)=-fk(ir)*etasi(ir)
      detasi(ir)=fk(ir)*etas(ir)
      d2etas(ir)=-fk(ir)*fk(ir)*etas(ir)
      d2etasi(ir)=-fk(ir)*fk(ir)*etasi(ir)
C Negative wavenumbers
      etas(nrng-ir+2)=etas(ir)
      etasi(nrng-ir+2)=-etasi(ir)
      detas(nrng-ir+2)=detas(ir)
      detasi(nrng-ir+2)=-detasi(ir)
      d2etas(nrng-ir+2)=d2etas(ir)
      d2etasi(nrng-ir+2)=-d2etasi(ir)
      end do
C DC and Nyquist Values
      etas(1)=-tide
      etasi(1)=0.
      detas(1)=0.
      detasi(1)=0.
      d2etas(1)=0.
      d2etasi(1)=0.
      etas(nrng/2+1)=0.
      etasi(nrng/2+1)=0.
      detas(nrng/2+1)=0.
      detasi(nrng/2+1)=0.
      d2etas(nrng/2+1)=0.
      d2etasi(nrng/2+1)=0.

      CALL FFT(lnr,etas,etasi)
      CALL FFT(lnr,detas,detasi)
      CALL FFT(lnr,d2etas,d2etasi)

C
9999  RETURN
      END
