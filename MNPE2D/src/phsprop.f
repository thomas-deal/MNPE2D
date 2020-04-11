      SUBROUTINE PHSPROP(phsr,phsi,topkr,topki,fk)

C     This subroutine generates the complex k-space propagator function
C     defined by phsprop(k,dr)=exp(-i*dr*k0*T(k)) where T(k)=1-sqrt(1-(k/k0)^2)
C     and k=k0*sin(theta) is the vertical wavenumber.  It also defines the
C	k-space wavenumber filter.

      REAL*4 phsr(1),phsi(1),topkr(1),topki(1)
      REAL*4 fk(1),filtk(16385)

      COMMON nz,dz,dk,dr,rng,freq,c0,sd,arl,thc,depcalc,bdint,dbdint,
     &irad,nrad,drad

      pi=acos(-1.)
      pi2=pi/2.
      degrad=180./pi
      dk=pi/depcalc
      nz2=nz/2
      nz2p1=nz2+1
      fk0=2.*pi*freq/c0

C	If phase space limited to <90 deg, then filter from 1 down to 1/2 over outer 1/10 of array
C	Otherwise, filter from 65 deg to 90 deg
      if(nz2*dk.lt.fk0)then
        nzmax=nz2
        nz10=nzmax/10
        thrms=float(nzmax-nz10)*dk/fk0
        tprw=float(nzmax)*dk/fk0
        thrms=asin(thrms)*degrad
        tprw=asin(tprw)*degrad
110     format(//,' Phase space filtering of field from ',f5.1,' to ',
     &   f5.1,' deg due to small FFT size.')
        write(*,110)thrms,tprw
      else
        nzmax=int(fk0/dk)
        nz10=nzmax/30
C        nz10=nzmax/10
        thrms=float(nzmax-nz10)*dk/fk0
        tprw=float(nzmax)*dk/fk0
        thrms=asin(thrms)*degrad
        tprw=asin(tprw)*degrad
      end if
      do iz=1,nz2p1
        fj=float(iz-(nzmax-nz10))/float(nz10)
        arg=(1.-amin1(1.,amax1(0.,fj)))*pi
        cosa=cos(arg)
        filtk(iz)= 0.25*(1.0-cosa)+0.5
      end do

C     Generate propagator function based on WAPE approximation
      DO iz=1,nz

        if(iz.le.nz2+1)then
          fk(iz)=(iz-1)*dk/fk0
        else
          fk(iz)=-(nz-iz+1)*dk/fk0
        end if
        fksq=fk(iz)*fk(iz)
        fc=0.
C
        if(fksq.lt.1.)then
          fc=1
          ang=dr*fk0*(1.-sqrt(1.-fksq))
          topkr(iz)=1.-sqrt(1.-fksq)
          topki(iz)=0.
        else
          fc=exp(amax1(-dr*fk0*sqrt(fksq-1.),-30.))
          ang=dr*fk0
          topkr(iz)=1.
          topki(iz)=-sqrt(fksq-1.)
        endif
C
        jz=min0(iz+1,nz+2-iz)
        phsr(iz)=filtk(jz)*fc*cos(ang)/float(nz)
        phsi(iz)=-filtk(jz)*fc*sin(ang)/float(nz)

      END DO
C
      RETURN
      END
