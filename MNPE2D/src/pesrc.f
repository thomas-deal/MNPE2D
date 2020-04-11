      SUBROUTINE PESRC(psir,psii,fk)

C     This subroutine generates the complex source function in k-space.  The
C     form of this source function corresponds to a wide-angle point source with
C     the near-field correction of Thomson and Bohun.  If a non-zero array
C     length is specified, an ideal, continuous line array (sinc function) is
C     used to define the source function with a D/E steering angle as input.

      REAL*4 psir(1),psii(1),fk(1)
      REAL*4 dsdk,A0,A2,A4,Az,arg1,arg2

      COMMON nz,dz,dk,dr,rng,freq,c0,sd,arl,thc,depcalc,bdint,dbdint,
     &irad,nrad,drad

C     A wide angle source is used unless a >zero array length is specified.
      srcline(xx)=sin(xx/wdth)/(xx/wdth)
      srcwide(xx)=1./sqrt(sqrt(1.-xx*xx))

C     Define parameters
      r0=1.0
      pi=acos(-1.)
      tpi=2.*pi
      pi2=pi/2.
      dzdk2=dz*dk/2.
      nzp2=nz+2
      nz2=nz/2
      fk0=tpi*freq/c0
      deg2rad=pi/180.

C     Define default source and taper width (used for filter) and initial
C     D/E angle
      thrms=60.*deg2rad
      tprw=30.*deg2rad
      thc=thc*deg2rad

C     Define monopole source depth
      sdk=fk0*sd
C     Define horizontal dipole vertical spacing and strength.
      dsdk=0.1
      A0=1./(2.*dsdk)
      A2=-A0*A0
      A4=A2*A2
C     Define vertical dipole vertical spacing and strength.
C     Strength term includes +90 degree phase shift that is applied to
C     xamp, yamp below.
      dsdk=0.1
      Az=1./(2.*dsdk);

C     Define source gain
      sgain=.5*dk*sqrt(r0/(pi2*fk0))
      fkc=0.0

      if(arl.gt.0.0)then
C     Define width of array type source and correction to array gain
        wdth=2./(fk0*arl)
        sgain=sgain*sqrt(amax1(1.,fk0*arl/pi))
        fkc=sin(thc)
      end if

      DO iz=1,nz2+1
C Point Monopole
        if(arl.eq.0.0)then
          arg1=amin1(.9999,amax1(-.9999,fk(iz)))
          xamp=0.0
          yamp=-2.*sgain*(srcwide(arg1))*sin(sdk*fk(iz))
C Point Hoprizontal Dipole
        elseif(arl.eq.-1.0)then
          arg1=amin1(.9999,amax1(-.9999,fk(iz)))
          arg2=-1./8.*A4*sin((sdk+4.*dsdk)*fk(iz))
          arg2=arg2+(-1./2.*A2+1./2.*A4)*sin((sdk+2.*dsdk)*fk(iz))
          arg2=arg2+(1.+A2-3./4.*A4)*sin((sdk)*fk(iz))
          arg2=arg2+(-1./2.*A2+1./2.*A4)*sin((sdk-2.*dsdk)*fk(iz))
          arg2=arg2-1./8.*A4*sin((sdk-4.*dsdk)*fk(iz))
          xamp=0.0
          yamp=-2.*sgain*(srcwide(arg1))*arg2
C Point Vertical Dipole
        elseif(arl.eq.-2.0)then
          arg1=amin1(.9999,amax1(-.9999,fk(iz)))
          arg2=Az*sin((sdk+dsdk)*fk(iz))
          arg2=arg2-Az*sin((sdk-dsdk)*fk(iz))
          xamp=2.*sgain*(srcwide(arg1))*arg2
          yamp=0.0

        else

          arg1=fk(iz)-fkc
          arg1=amin1(.9999,amax1(-.9999,arg1))
          arg2=fk(iz)+fkc
          arg2=amin1(.9999,amax1(-.9999,arg2))
          if(arg1.ne.0.)then
           src1=srcline(arg1)
          else
           src1=1.
          end if
          if(arg2.ne.0.)then
           src2=srcline(arg2)
          else
           src2=1.
          end if
          xamp=+sgain*(src1-src2)*cos(sdk*fk(iz))
          yamp=-sgain*(src1+src2)*sin(sdk*fk(iz))

        end if

C     Correct k-space source for 1/2-mesh symmetry in z-space
        phs=pi*((float(iz)-1.)/float(nz))
        psir(iz)=xamp*cos(phs)-yamp*sin(phs)
        psii(iz)=xamp*sin(phs)+yamp*cos(phs)
        if(iz.gt.1.and.iz.lt.nz2+1)then
          jz=nzp2-iz
          psir(jz)=-xamp*cos(phs)-yamp*sin(phs)
          psii(jz)=xamp*sin(phs)-yamp*cos(phs)
        end if

      END DO

C     Perform k-space taper of source function
      DO iz=1,nz2+1
        tprmin=amin1(thc+thrms,pi2-tprw)
        tprmax=amin1(thc+thrms+tprw,pi2)
        if(fk(iz).gt.sin(tprmin))then
          if(fk(iz).lt.sin(tprmax))then
            taper=(fk(iz)-sin(tprmin))/(sin(tprmax)-sin(tprmin))
            taper=0.5*(cos(pi*taper)+1.)
          else
            taper=0.0
          end if
        else
          taper=1.
        end if

        psir(iz)=psir(iz)*taper
        psii(iz)=psii(iz)*taper

        if(iz.gt.1.and.iz.lt.nz2+1)then
          jz=nzp2-iz
          psir(jz)=psir(jz)*taper
          psii(jz)=psii(jz)*taper
        end if

      END DO
      if(arl.lt.0.)then
        arl=0.
      endif
      RETURN
      END
