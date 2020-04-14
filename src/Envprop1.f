      SUBROUTINE ENVPROP(envr,envi,envr2,envi2,enz,filt,istrt,ird)
C
C     This subroutine generates the complex environmental propagator function,
C     envprop(z,r)=exp(-i*dr*k0*(1-n(z,r))) where n(z,r) is the index of
C     refraction computed from the environmental inputs (and effective values
C     due to density contrasts).  The z-space filter is also defined.

C     Note envprop2=envprop(z,dr/2).

C     Environmental inputs are as follows:  rss(nss) are the ranges [km] of nss
C     sound speed profiles defined at nssd(nss) depths ssd(nssd(nss),nss) [m]
C     with sound speed values ss(nssd(nss),nss) [m/s].  If only one value is
C     given, the medium is homogeneous.  If only one value in depth is given,
C     that profile has a constant sound speed.  If more than one value in depth
C     is given, the profile is extended to the maximum computational depth (if
C     necessary) using the gradient between the last two values in depth.
C
C     Two bottom types may be input with at least one bottom interface required.
C     The ranges and depths of the first bottom (assumed to be the water/bottom
C     interface) are defined by rb(nb) [km] and bd(nb) [m], respectively.  The
C     optional deep bottom interface ranges and depths are rdb(ndb) [km] and
C     dbd(ndb) [m], respectively.  Each bottom type may have its own range-
C     independent or range- dependent parameters specified.  Parameters to be
C     defined are range, rbp(nbp) [km], sound speed, bss(nbp) [m/s], sound speed
C     gradient, bg(nbp) [1/s], compressional loss, blkm(nbp) [dB/km/Hz],
C     density, bden(nbp) [g/cc] (assumes water has density 1.0 g/cc), shear wave
C     speed, bsws(nbp) [m/s], and shear wave loss, bswlkm(nbp) [dB/km/Hz].
C     Similar parameters may be defined for the deep bottom layer at ranges
C     rdbp(ndbp) [km].

C     NOTE:  If more than one range profile of any parameter is input, even if
C     the profiles are identical, the environment will be treated as range-
C     dependent and this propagator function will be recomputed at every range
C     step.

C     NOTE:  Current restrictions on the number of environmental input data may
C     be easily changed by altering the array size specifications.

C     NOTE:  This version only reads a single radial from the environmental
C     data files.  To use more than one radial, use other version of MMPE.

C
C
      REAL*4 envr(1),envi(1),envr2(1),envi2(1),enz(1),filt(1)
      REAL*4 ss(512,251),ssd(512,251),rss(251)
      REAL*4 rb(2501),bd(2501),rbp(101),bss(101),bg(101),
     + blkm(101),bl(101),bden(101),bsws(101),bswlkm(101),bswl(101)
      REAL*4 rdb(101),dbd(101),rdbp(101),dbss(101),dbg(101),
     + dblkm(101),dbl(101),dbden(101),dbsws(101),dbswlkm(101),dbswl(101)
      REAL*4 ss1(16384),ss2(16384),ssint(16384),aloss(16384)
      REAL*4 etas(8192),detas(8192),d2etas(8192)
      COMPLEX ang(16384)
      COMPLEX ic,term1,term2,term3,term4,bdenint,dbdenint,eneff
      COMPLEX	brod,dbrod,ub0,udb0,env,env2,en,en2,f2,f3,avgden
      INTEGER*4 nssda(251)
C
      COMMON nz,dz,dk,dr,rng,freq,c0,sd,arl,thc,depcalc,bdint,dbdint,
     &irad,nrad,drad,rngmax,rsint,rsdotint,Hdep

      SAVE nradb,bdmax,nb,nbmax,rb,bd,nbp,iden,rbp,bss,bg,bden,bl,blkm
      SAVE bsws,bswl,bswlkm,ndb,rdb,dbd,ndbp,rdbp,dbss,dbg,dbden,dbl,dblkm
      SAVE dbsws,dbswl,dbswlkm,nradss,nvloss
      SAVE nss,rss,nssda,nssmax,nssd,ssd,ss,ss1,ss2
      SAVE nz2,nz3,ny2,ny4,nz2p1,ny2p1,alpha0
      SAVE iss,rss1,rss2,ib,rb1,rb2,ibp,rbp1,rbp2
      SAVE idb,rdb1,rdb2,idbp,rdbp1,rdbp2
      SAVE etas,detas,d2etas,sigsmax
C

      fctn1(xx)=1./(1.+exp(-xx))
      fctn2(xx)=exp(-xx)/(1.+exp(-xx))**2.
      fctn3(xx)=-exp(-xx)*(1.-exp(-xx))/(1.+exp(-xx))**3.
C
      pi=acos(-1.)
      tpi=2.*pi
      fk0=tpi*freq/c0
      ic=csqrt((-1.,0))

C     Formats for environmental outputs to confirm proper input structure

100   format(1x,///2X,'SOUND SPEED PROFILE AT RANGE =',
     &        F8.2,' km'/6X,'DEPTH(m)',3X,'SSPD(m/s)')
101   format(1x,2F12.2)
110   format(1x///2X,'BOTTOM DEPTH')
111   format(6X,'RANGE(km)',2X,'DEPTH(m)')
112   format(1x,2F12.2)
120   format(1x///2X,'BOTTOM PROPERTIES')
121   format(2X,'RNG(km)',1X,'VEL(m/s)',
     &1X,'GRAD(1/s)',1X,'DENS',1X,'LOSS(dB/km/Hz)',
     &1X,'SHR(m/s)',1X,'SHRLOSS(dB/km/Hz)')
122   format(F7.2,3X,F7.2,3X,F5.2,4X,F4.2,5X,F6.4,
     &2X,F7.2,4X,F6.4)
130   format(1x///2X,'DEEP BOTTOM DEPTH')
140   format(1x///2X,'DEEP BOTTOM PROPERTIES')
150   format(///)
C
C     Read environmental inputs if first pass, otherwise begin at 1000

      IF(istrt.ne.0) goto 1000
      istrt=1
      irtmpreset=0
C
C     Initialize arrays
      DO k=1,101
        rss(k)=0.
        nssda(k)=0
        DO i=1,512
          ss(i,k)=0.
        END DO
        rdb(k)=0.
        dbd(k)=0.
      END DO
      DO k=1,251
          rb(k)=0.
          bd(k)=0.
      END DO
      nssmax=0
      nbmax=0
      nbp=0
      ndb=0
      ndbp=0

C     Read bathymetry data (upper layer)
C     Use only single radial
      nradb=1
      drad=0.

      bdmax=0.
      read(11,*)nb
      if(nb.le.0.or.nb.gt.2500)then
       print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: NB'
       istrt=2
       goto 9999
      end if
      nbmax=max0(nbmax,nb)

      DO j=1,nb
        read(11,*)rb(j),bd(j)
        bdmax=amax1(bdmax,bd(j))
        if(j.eq.1)then
         if(rb(1).ne.0.)then
          print *,'ERROR IN ENVIRONMENTAL INPUT FORMAT: RB(1)<>0'
          istrt=2
          goto 9999
         end if
        else
         if(rb(j).le.rb(j-1))then
          print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: RB(J)<=RB(J-1)'
          istrt=2
          goto 9999
         end if
        end if
      END DO

      if(depcalc.eq.0.)depcalc=3.*bdmax
      if(nz.eq.0)then
        nz=int(2.*depcalc/dz)
        lnz=int(alog(float(nz))/alog(2.)+0.99)
      else
        lnz=int(alog(float(nz))/alog(2.)+0.99)
      end if
      if(lnz.gt.15)then
        print*,'Error in nz array size; exceeding max limit.'
        print*,'Setting nz = 32768 (max allowable).'
        print*,' '
        lnz=15
      end if
      nz=2**lnz
      dz=2.*depcalc/nz

C     Read bottom properties
      read(12,*)nbp
      if(nbp.le.0.or.nbp.gt.100)then
       print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: NBP'
       istrt=2
       goto 9999
      end if

      iden=0
      DO j=1,nbp
        read(12,*)rbp(j),bss(j),bg(j),bden(j),blkm(j),bsws(j),bswlkm(j)
        if(bden(j).ne.1.0)iden=iden+1
        if(j.eq.1)then
         if(rbp(1).ne.0.)then
          print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: RBP(1)<>0'
          istrt=2
          goto 9999
         end if
        else
         if(rbp(j).le.rbp(j-1))then
          print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: RBP(J)<=RBP(J-1)'
          istrt=2
          goto 9999
         end if
        end if
      END DO

C     Read bathymetry data (lower layer)
      read(13,*)ndb
      if(ndb.lt.0.or.ndb.gt.100)then
       print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: NDB'
       istrt=2
       goto 9999
      end if

      DO j=1,ndb
        read(13,*)rdb(j),dbd(j)
        if(j.eq.1)then
         if(rdb(1).ne.0.)then
          print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: RDB(1)<>0'
          istrt=2
          goto 9999
         end if
        else
         if(rdb(j).le.rdb(j-1))then
          print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: RDB(J)<=RDB(J-1)'
          istrt=2
          goto 9999
         end if
        end if
      END DO

C     Read lower layer bottom properties
      IF(ndb.gt.0)THEN
      read(14,*)ndbp
      if(ndbp.lt.0.or.ndbp.gt.100)then
       print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: NDBP'
       istrt=2
       goto 9999
      end if

      DO j=1,ndbp
        read(14,*)rdbp(j),dbss(j),dbg(j),dbden(j),dblkm(j),dbsws(j),
     &  dbswlkm(j)
        if(dbden(j).ne.1.0)iden=iden+1
        if(j.eq.1)then
         if(rdbp(1).ne.0.)then
          print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: RDBP(1)<>0'
          istrt=2
          goto 9999
         end if
        else
         if(rdbp(j).le.rdbp(j-1))then
        print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: RDBP(J)<=RDBP(J-1)'
          istrt=2
          goto 9999
         end if
        end if
      END DO
      END IF
C
200   if(ndb.eq.0)then
        ndb=1
        rdb(1)=0.
        dbd(1)=depcalc
      end if

      if(bd(1).lt.dbd(1))then
        Hdep = bd(1)
      else
        Hdep = dbd(1)
      endif

C     Read sound speed data
C     Use only single radial
      nradss=1

C	Use volume loss in water column? (1=yes,0=no)
      read(10,*)nvloss
      if(nvloss.ne.0)then
        if(nvloss.ne.1)then
          print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: NSS'
          istrt=2
          goto 9999
        end if
      end if


      read(10,*)nss
      if(nss.le.0.or.nss.gt.250)then
       print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: NSS'
       istrt=2
       goto 9999
      end if
C
      DO j=1,nss
        read(10,*)rss(j),nssda(j)
        if(j.eq.1)then
         if(rss(1).ne.0.)then
          print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: RSS(1)<>0'
          istrt=2
          goto 9999
         end if
        else
         if(rss(j).le.rss(j-1))then
          print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: RSS(J)<=RSS(J-1)'
          istrt=2
          goto 9999
         end if
        end if
        nssmax=max0(nssmax,nss)
        nssd=nssda(j)
        if(nssd.gt.512)then
         print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: NSSD'
         istrt=2
         goto 9999
        end if
        DO k=1,nssd
          read(10,*)ssd(k,j),ss(k,j)
          if(k.eq.1)then
           if(ssd(1,j).ne.0.)then
            print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: SSD(1,J)<>0'
            istrt=2
            goto 9999
           end if
          else
           if(ssd(k,j).le.ssd(k-1,j))then
      print*,'ERROR IN ENVIRONMENTAL INPUT FORMAT: SSD(K,J)<=SSD(K-1,J)'
            istrt=2
            goto 9999
           end if
          end if
        END DO
        if(ssd(nssd,j).lt.depcalc)then
          nssda(j)=nssda(j)+1
          ssd(nssd+1,j)=depcalc
          if(nssd.gt.1)then
            ss(nssd+1,j)=ss(nssd,j)+(ssd(nssd+1,j)-
     &ssd(nssd,j))*(ss(nssd,j)-ss(nssd-1,j))/
     &(ssd(nssd,j)-ssd(nssd-1,j))
          else
            ss(nssd+1,j)=ss(nssd,j)
          end if
        end if
      END DO

C     Check for any density contrasts.  Verify depth mesh is adequate.
      if((iden.ne.0).and.(dz.gt.(c0/freq)))then
        print*,'WARNING:  DEPTH MESH MAY NOT BE SMALL ENOUGH'
        print*,'TO ADEQUATELY MODEL DENSITY DISCONTINUITIES.'
      end if

C     Check if range-independent, then ird=1.
C      ird=max0(nssmax,nbmax,nbp,ndb,ndbp)
C	Force range-dependence for perturbations
      ird=2
C     Determine number of radials
      nrad=max0(nradss,nradb)

C
C     Output first and last environmental profiles to log file for checking
C     and extend all data beyond final ranges and depths
C
      DO k=1,nss
        if(k.eq.1.or.k.eq.nss)then
          nssd=nssda(k)
          write(9,100)rss(k)
          DO j=1,nssd
            if(amod(float(j),amax0(int(float(nssd)/20.),1)).eq.0.)
     &        write(9,101)ssd(j,k),ss(j,k)
          END DO
        end if
        rss(k)=1000.*rss(k)
      END DO

      DO j=1,nssda(nss)
        ss(j,nss+1)=ss(j,nss)
        ssd(j,nss+1)=ssd(j,nss)
      END DO
      nssda(nss+1)=nssda(nss)
      rss(nss+1)=1.e20
      nss=nss+1

C
      write(9,110)
      write(9,111)
       DO k=1,nb
        if((k.eq.1.or.k.eq.nb))then
          write(9,112)rb(k),bd(k)
        end if
        rb(k)=1000.*rb(k)
       END DO
       rb(nb+1)=1.e20
       bd(nb+1)=bd(nb)
       nb=nb+1

C
      write(9,120)
      write(9,121)
      DO k=1,nbp
        if(k.eq.1.or.k.eq.nbp)write(9,122)rbp(k),bss(k),bg(k),
     &   bden(k),blkm(k),bsws(k),bswlkm(k)
        rbp(k)=1000.*rbp(k)
      END DO
      rbp(nbp+1)=1.e20
      bss(nbp+1)=bss(nbp)
      bden(nbp+1)=bden(nbp)
      blkm(nbp+1)=blkm(nbp)
      bg(nbp+1)=bg(nbp)
      bsws(nbp+1)=bsws(nbp)
      bswlkm(nbp+1)=bswlkm(nbp)
      nbp=nbp+1

C
      if(ndb.gt.1.or.dbd(1).lt.depcalc)then
       write(9,130)
       write(9,111)
      end if
      DO k=1,ndb
        if(ndb.gt.1.or.dbd(1).lt.depcalc)then
         if(k.eq.1.or.k.eq.ndb)write(9,112)rdb(k),dbd(k)
        end if
        rdb(k)=1000.*rdb(k)
      END DO

C
      if(ndb.gt.1.or.dbd(1).lt.depcalc)then
       write(9,140)
       write(9,121)
      end if
      DO k=1,ndbp
        if(ndb.gt.1.or.dbd(1).lt.depcalc)then
         if(k.eq.1.or.k.eq.ndbp)write(9,122)rdbp(k),dbss(k),dbg(k),
     &    dbden(k),dblkm(k),dbsws(k),dbswlkm(k)
        end if
        rdbp(k)=1000.*rdbp(k)
      END DO

C
      rdb(ndb+1)=1.e20
      dbd(ndb+1)=dbd(ndb)
      ndb=ndb+1
      rdbp(ndbp+1)=1.e20
      dbss(ndbp+1)=dbss(ndbp)
      dbden(ndbp+1)=dbden(ndbp)
      dblkm(ndbp+1)=dblkm(ndbp)
      dbg(ndbp+1)=dbg(ndbp)
      dbsws(ndbp+1)=dbsws(ndbp)
      dbswlkm(ndbp+1)=dbswlkm(ndbp)
      ndbp=ndbp+1

C
      write(9,150)

C     Create filter over 1/3 end of vector - need to make maxdep>2*botdep
      nz2=nz/2
      nz2p1=nz2+1
      nz3=nz2/3
      do iz=1,nz2p1
        fj=float(nz2p1-iz)/float(nz3)
        arg=amin1(1.,fj)*pi
        cosa=cos(arg)
C        filt(iz)= 0.125*(1.0-cosa)+0.75
        filt(iz)= 0.25*(1.0-cosa)+0.5
      end do

C
C	Call subroutines to compute surface roughness.
      CALL ZSGEN(etas,detas,d2etas,sigsmax,istrt)
      if(istrt.eq.2) goto 9999


C
1000  CONTINUE
C      NOTE: After initialization, this subroutine begins here.

C      Define mixing lengths
      wvln=tpi/fk0
C      xld=amax1(xldmin,2.*wvln)
      xldmin=2./fk0
      xld=amax1(xldmin,5.*dz)
C      xldmin=dz*5.0
      xldsq=xld*xld
      xlmin=dz*1.0
C      xl=amax1(xlmin,0.1*wvln)
      xl=dz/100.

C     Compute water volume loss coefficients
      conv=1.094/(1000.0*8.686)
      fsq=(freq/1.e3)**2
      ALPHA0=0.003+0.1*fsq/(1.+fsq)+40.*fsq/(4100.+fsq)+2.75e-4*fsq
      ALPHA0=conv*ALPHA0
      if(nvloss.eq.0)then
        ALPHA0=0.
      end if
C
C	NOTE:  If rng=0, subroutine only called for irad=1
      if(rng.eq.0)then
        iss=1
        nssd=nssda(iss)
        rss2=rss(iss)
        CALL SSI(ss(1,iss),ssd(1,iss),nssd,ss2(1))
        ib=1
        rb2=rb(ib)
       ibp=1
       idb=1
       idbp=1
       rbp2=rbp(ibp)
       rdb2=rdb(idb)
       rdbp2=rdbp(idbp)
       do j=1,nbp
         bl(j)=blkm(j)/(1000.0*8.686)*freq
         bswl(j)=bswlkm(j)/(1000.0*8.686)*freq
       end do
       do j=1,ndbp
         dbl(j)=dblkm(j)/(1000.0*8.686)*freq
         dbswl(j)=dbswlkm(j)/(1000.0*8.686)*freq
       end do
      end if
C
C
1009  continue
      if(rng.lt.rss2) goto 1020
      iss=iss+1
C
C     Update profile values for water column
      rss1=rss2
      rss2=rss(iss)
      nssd=nssda(iss)
      DO 1010 i=1,nz2
1010  ss1(i)=ss2(i)
      CALL SSI(ss(1,iss),ssd(1,iss),nssd,ss2(1))
      goto 1009
C
1020  continue
      if(rng.lt.rb2) goto 1030
      ib=ib+1
C
C     Update values for bathymetry (upper)
      rb1=rb2
      rb2=rb(ib)
      goto 1020
C
1030  continue
      if(rng.lt.rbp2) goto 1040
      ibp=ibp+1
C
C     Update values for bottom (upper) properties
      rbp1=rbp2
      rbp2=rbp(ibp)
      goto 1030
C
1040  continue
      if(rng.lt.rdb2) goto 1050
      idb=idb+1
C
C     Update values for deep bottom layer
      rdb1=rdb2
      rdb2=rdb(idb)
      goto 1040
C
1050  continue
      if(rng.lt.rdbp2) goto 1060
      idbp=idbp+1
C
C     Update values for deep bottom properties
      rdbp1=rdbp2
      rdbp2=rdbp(idbp)
C
C     Compute interpolation values
      goto 1050
C
1060  continue
      if(sigsmax.ne.0.)then
C     Update surface displacement from random rough surface realization
        irtmp=max(int(rng/dr),1)
        if(irtmp.gt.8192)then
          irtmpreset=irtmpreset+1
          irtmp=irtmp-irtmpreset*8192
        end if
        rsint=etas(irtmp)
        rsdotint=detas(irtmp)
        rsdot2int=d2etas(irtmp)
      else
        rsint=0.
        rsdotint=0.
        rsdot2int=0.
      end if
C
      fss=(rng-rss1)/(rss2-rss1)
      fb=(rng-rb1)/(rb2-rb1)
      fbp=(rng-rbp1)/(rbp2-rbp1)
      fdb=(rng-rdb1)/(rdb2-rdb1)
      fdbp=(rng-rdbp1)/(rdbp2-rdbp1)
C
      bdint=bd(ib-1)+fb*(bd(ib)-bd(ib-1))
C
      bssint=bss(ibp-1)+fbp*(bss(ibp)-bss(ibp-1))
      bgint=bg(ibp-1)+fbp*(bg(ibp)-bg(ibp-1))
      blint=bl(ibp-1)+fbp*(bl(ibp)-bl(ibp-1))
      bdenint=bden(ibp-1)+fbp*(bden(ibp)-bden(ibp-1))
      bswsint=bsws(ibp-1)+fbp*(bsws(ibp)-bsws(ibp-1))
      bswlint=bswl(ibp-1)+fbp*(bswl(ibp)-bswl(ibp-1))
C
      dbdint=dbd(idb-1)+fdb*(dbd(idb)-dbd(idb-1))
C
      dbssint=dbss(idbp-1)+fdbp*(dbss(idbp)-dbss(idbp-1))
      dbgint=dbg(idbp-1)+fdbp*(dbg(idbp)-dbg(idbp-1))
      dblint=dbl(idbp-1)+fdbp*(dbl(idbp)-dbl(idbp-1))
      dbdenint=dbden(idbp-1)+fdbp*(dbden(idbp)-dbden(idbp-1))
      dbswsint=dbsws(idbp-1)+fdbp*(dbsws(idbp)-dbsws(idbp-1))
      dbswlint=dbswl(idbp-1)+fdbp*(dbswl(idbp)-dbswl(idbp-1))
C
C     Compute updated sound speed profile in water column
      DO 2000 i=1,nz2
2000  ssint(i)=ss1(i)+fss*(ss2(i)-ss1(i))

C     Convert shear wave effect into complex bottom density and define
C	density ratios at interfaces
C
C	if deep bottom rises above sediment bottom and interfaces with water
C	(such as a rock outcropping)
      if(dbdint.le.bdint)then
        ndbi=int(1.+dbdint/dz)
        ssdbw=ssint(ndbi)
C	if shear exists, require that bottom sound speed exceeds water
C	column; otherwise, just leave parameters alone as sound will
C	penetrate into slower bottom anyway (not likely to have high
C	shear in such cases)
       if((dbswsint.gt.0.).and.(dbssint.gt.ssdbw))then
         term1=ssdbw/dbswsint+ic*dbswlint*ssdbw/(tpi*freq)
         term1=1.-2./(term1*term1)
         term1=term1*term1
         term2=ssdbw/dbssint+ic*dblint*ssdbw/(tpi*freq)
         term2=csqrt(1.-term2*term2)
         term3=ssdbw/dbswsint+ic*dbswlint*ssdbw/(tpi*freq)
         term3=csqrt(term3*term3-1.)
         term4=ssdbw/dbswsint+ic*dbswlint*ssdbw/(tpi*freq)
         term4=term4*term4*term4*term4
         dbdenint=dbdenint*(term1+ic*4.*term2*term3/term4)
       end if
C        dbrod=(csqrt(dbdenint)-1.)/(csqrt(dbdenint)+1.)
C	else deep bottom is below water/sediment interface
      else
        nbi=int(1.+bdint/dz)
        ssbw=ssint(nbi)
       if((bswsint.gt.0.).and.(bssint.gt.ssbw))then
         term1=ssbw/bswsint+ic*bswlint*ssbw/(tpi*freq)
         term1=1.-2./(term1*term1)
         term1=term1*term1
         term2=ssbw/bssint+ic*blint*ssbw/(tpi*freq)
         term2=csqrt(1.-term2*term2)
         term3=ssbw/bswsint+ic*bswlint*ssbw/(tpi*freq)
         term3=csqrt(term3*term3-1.)
         term4=ssbw/bswsint+ic*bswlint*ssbw/(tpi*freq)
         term4=term4*term4*term4*term4
         bdenint=bdenint*(term1+ic*4.*term2*term3/term4)
       end if
C        brod=(csqrt(bdenint)-1.)/(csqrt(bdenint)+1.)
C
C	now check to see if deep bottom exists above computational depth
       if(dbdint.lt.depcalc)then
        bdbssint=bssint+bgint*(dbdint-bdint)
       if((dbswsint.gt.0.).and.(dbssint.gt.bdbssint))then
         term1=bdbssint/dbswsint+ic*dbswlint*bdbssint/(tpi*freq)
         term1=1.-2./(term1*term1)
         term1=term1*term1
         term2=bdbssint/dbssint+ic*dblint*bdbssint/(tpi*freq)
         term2=csqrt(1.-term2*term2)
         term3=bdbssint/dbswsint+ic*dbswlint*bdbssint/(tpi*freq)
         term3=csqrt(term3*term3-1.)
         term4=bdbssint/dbswsint+ic*dbswlint*bdbssint/(tpi*freq)
         term4=term4*term4*term4*term4
         dbdenint=dbdenint*(term1+ic*4.*term2*term3/term4)
       end if
C        dbrod=(csqrt(dbdenint)-1.)/(csqrt(dbdenint)+1.)
       end if
      end if
C
C      ub0=-brod/(fk0*fk0)
C      udb0=-dbrod/(fk0*fk0)
C
CCCCCCCCCCCCCCCCCCCCCC
C      if(rng.eq.0)then
C	open(21,file='en.dat')
C	open(22,file='en2.dat')
C	end if
CCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCC	Define rough surface displacement for test
C	swvln=20.
C	surfint=5.*sin(2*pi*rng/swvln)
C	rsint=0.
C	rsdotint=0.
C	rsdot2int=0.
CCCCCCCCCCC

      DO 5000 i=1,nz2

      depth=(float(i)-.5)*dz

      en2=0
C     Compute volume loss coefficient in water column
      vloss=ALPHA0*(1.-6.46E-5*depth)

      IF(bdint.lt.dbdint)THEN

C     Mix water and bottom sound speed with hyperbolic function at interface
      rdepth=depth-bdint
      if(depth.gt.bdint)then
       if(depth.gt.dbdint)then
        ssintb=bssint+bgint*(dbdint-bdint)
       else
        ssintb=bssint+bgint*rdepth
       end if
      else
       ssintb=bssint
      end if
      r=amin1(amax1(rdepth/xl,-32.),32.)
      f1=fctn1(r)
      ssint(i)=ssint(i)+f1*(ssintb-ssint(i))
      vloss=vloss+f1*(blint-vloss)
C
C      ibot=1

C     Compute mixing function for density discontinuity at interface
      if(xld.eq.xldmin)then
        rdepthd=depth-dz*(int(bdint/dz+1.)-0.5)
      else
        rdepthd=depth-bdint
      end if
C	Old method with hyperbolic mixing
      r=amin1(amax1(rdepthd/xld,-32.),32.)
      upden=1.
      avgden=(upden+bdenint)/2.
      f2=1.5*(bdenint-upden)*(bdenint-upden)*fctn2(r)*fctn2(r)
      f2=f2/(avgden*avgden)
      f3=(bdenint-upden)*fctn3(r)/avgden
      en2=(f3-f2)/(2.*(xld*fk0)**2.)
C	Tappert's method with hyperbolic mixing
C	  denfact=(1.-sqrt(upden/bdenint))/(1+sqrt(upden/bdenint))
C	  en2=-denfact*f3/(xld*fk0)**2.
C	Tappert's method with cubic spline mixing
C      goto 3100
C3000  en2=ub0*u

      END IF

C     Mix water and deep bottom sound speed with hyperbolic function at interface
      IF(dbdint.lt.depcalc)THEN
       rdepth=depth-dbdint
       if(depth.gt.dbdint)then
        ssintb=dbssint+dbgint*rdepth
       else
        ssintb=dbssint
       end if
       r=amin1(amax1(rdepth/xl,-32.),32.)
       f1=fctn1(r)
       ssint(i)=ssint(i)+f1*(ssintb-ssint(i))
       vloss=vloss+f1*(dblint-vloss)
C
C       ibot=2

C     Compute mixing function for density discontinuity at interface
       if(xld.eq.xldmin)then
         rdepthd=depth-dz*(int(dbdint/dz+1.)-0.5)
       else
         rdepthd=depth-dbdint
       end if
C	Old method with hyperbolic mixing
      r=amin1(amax1(rdepthd/xld,-32.),32.)
      if(bdint.lt.dbdint)then
        upden=bdenint
      else
        upden=1.
      end if
      avgden=(upden+dbdenint)/2.
      f2=1.5*(dbdenint-upden)*(dbdenint-upden)*fctn2(r)*fctn2(r)
      f2=f2/(avgden*avgden)
      f3=(dbdenint-upden)*fctn3(r)/avgden
      en2=(f3-f2)/(2.*(xld*fk0)**2.)
C	Tappert's method with hyperbolic mixing
C	  denfact=(1.-sqrt(upden/dbdenint))/(1+sqrt(upden/dbdenint))
C	  en2=-denfact*f3/(xld*fk0)**2.
C	Tappert's method with cubic spline mixing
C       goto 3100
C3010   en2=udb0*u+en2
      END IF

C     Mixing of sound speeds and densities complete
C      GOTO 4000
C

C3100  if (rdepthd.le.-xld) m=1
C      if (rdepthd.ge.-xld.and.rdepthd.le.-0.5*xld) m=2
C      if (rdepthd.ge.-0.5*xld.and.rdepthd.le.0.5*xld) m=3
C      if (rdepthd.ge.0.5*xld.and.rdepthd.le.xld) m=4
C      if (rdepthd.ge.xld) m=5

C      go to (3200,3210,3220,3230,3240),m

C3200  u=0.0
C      goto 3300

C3210  u=4.0*(1.+rdepthd/xld)/xldsq
C      goto 3300

C3220  u=-4.0*rdepthd/(xld*xldsq)
C      goto 3300

C3230  u=-4.0*(1.-rdepthd/xld)/xldsq
C      goto 3300

C3240  u=0.0

C3300  if(ibot.eq.1)then
C        goto 3000
C       else
C        goto 3010
C      end if


C     Compute index of refraction
4000  en=c0/ssint(i)

CCCCCCCCCCCCCCCCCCCCCC
C      if(rng.eq.0)then
C	write(21,*)depth,real(en)
C	write(22,*)real(en2),imag(en2),real(f2),real(f3)
C	end if
CCCCCCCCCCCCCCCCCCCCCC
C
C     Compute z-space operator for defined potential function
C
      eneff=csqrt(en*en+en2)
      ang(i)=dr*fk0*(1.-eneff)
C      ang(i)=dr*fk0*(1.-en+en2)
C
C     Compute volume and bottom loss contributions and spatial filter
      ip1=i+1
      aloss(i)=filt(ip1)*exp(-dr*vloss)
C
      env=aloss(i)*cexp(-ic*ang(i))
      env2=sqrt(aloss(i))*cexp(-ic*ang(i)/2)

      envr(i)=real(env)
      envi(i)=imag(env)
      envr2(i)=real(env2)
      envi2(i)=imag(env2)
      enz(i)=eneff
C     Populate image ocean propagator for flat surface
      envr(nz-i+1)=envr(i)
      envi(nz-i+1)=envi(i)
      envr2(nz-i+1)=envr2(i)
      envi2(nz-i+1)=envi2(i)
C	enz(i)=en-en2
C
5000  CONTINUE
C
C     Create image ocean potential function for surface scatter
      if(rsint.ne.0.)then
        if(rsint.gt.0.)then
          ns=int(rsint/dz+0.5)
        else
          ns=int(rsint/dz-0.5)
        end if
        if(ns.gt.0)then
          do iz=1,ns
            dep=(float(iz)-0.5)*dz
            angs=ang(2*ns+1-iz)-dr*fk0*2.0*(dep-rsint)*rsdot2int
            envr(iz)=aloss(2*ns+1-iz)*cos(angs)
            envi(iz)=-aloss(2*ns+1-iz)*sin(angs)
            envr2(iz)=sqrt(aloss(2*ns+1-iz))*cos(angs/2.)
            envi2(iz)=-sqrt(aloss(2*ns+1-iz))*sin(angs/2.)
          end do
          do iz=2*ns+1,nz/2
            jz=nz-iz+2*ns+1
            dep=(float(jz-nz)-0.5)*dz
            angs=ang(iz)-dr*fk0*2.0*(dep-rsint)*rsdot2int
            envr(jz)=aloss(iz)*cos(angs)
            envi(jz)=-aloss(iz)*sin(angs)
            envr2(jz)=sqrt(aloss(iz))*cos(angs/2.)
            envi2(jz)=-sqrt(aloss(iz))*sin(angs/2.)
          end do
          do iz=nz/2+1,nz/2+2*ns
            dep=(float(iz-nz)-0.5)*dz
            angs=ang(nz/2)-dr*fk0*2.0*(dep-rsint)*rsdot2int
            envr(iz)=aloss(nz/2)*cos(angs)
            envi(iz)=-aloss(nz/2)*sin(angs)
            envr2(iz)=sqrt(aloss(nz/2))*cos(angs/2.)
            envi2(iz)=-sqrt(aloss(nz/2))*sin(angs/2.)
          end do
        end if
        if(ns.le.0)then
          ns=-ns
          sspd1=ssint(1)
          sspd2=ssint(2)
          sspdint=sspd1-sspd2
          do iz=1,ns
            sspd=sspd1+iz*sspdint
            en=c0/sspd
            angs=dr*fk0*(1.-en)
            jz=nz-iz+1
            envr(jz)=aloss(1)*cos(angs)
            envi(jz)=-aloss(1)*sin(angs)
            envr2(jz)=sqrt(aloss(1))*cos(angs/2.)
            envi2(jz)=-sqrt(aloss(1))*sin(angs/2.)
            mz=nz-2*ns+iz
            dep=(float(mz-nz)-0.5)*dz
            angs=angs-dr*fk0*2.0*(dep-rsint)*rsdot2int
            envr(mz)=aloss(1)*cos(angs)
            envi(mz)=-aloss(1)*sin(angs)
            envr2(mz)=sqrt(aloss(1))*cos(angs/2.)
            envi2(mz)=-sqrt(aloss(1))*sin(angs/2.)
          end do
C	    do iz=1,ns
C	      jz=nz-iz+1
C	      mz=nz-2*ns+iz
C	      dep=(float(mz-nz)-0.5)*dz
C	      angs=ang(jz)-dr*fk0*2.0*(dep-rsint)*rsdot2int
C		  envr(mz)=aloss(jz)*cos(angs)
C		  envi(mz)=-aloss(jz)*sin(angs)
C		  envr2(mz)=sqrt(aloss(jz))*cos(angs/2.)
C		  envi2(mz)=-sqrt(aloss(jz))*sin(angs/2.)
C	    end do
          do iz=1,nz/2-2*ns
            jz=nz-2*ns+1-iz
            dep=(float(jz-nz)-0.5)*dz
            angs=ang(iz)-dr*fk0*2.0*(dep-rsint)*rsdot2int
            envr(jz)=aloss(iz)*cos(angs)
            envi(jz)=-aloss(iz)*sin(angs)
            envr2(jz)=sqrt(aloss(iz))*cos(angs/2.)
            envi2(jz)=-sqrt(aloss(iz))*sin(angs/2.)
          end do
        end if
      end if
C
CCCCCCCCCCCCCCCCCCCCCC
C      if(rng.eq.0)then
C	close(21)
C	close(22)
C	end if
CCCCCCCCCCCCCCCCCCCCCC
C
9999  RETURN
C
      END
