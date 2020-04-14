      SUBROUTINE SSI(ss,ssd,nssd,ss2)
C
C     This subroutine performs a linear interpolation in depth of given
C     sound speed profiles to the computational mesh.
C
      REAL*4 ss(1),ssd(1),ss2(1)
C
      COMMON nz,dz,dk,dr,rng,freq,c0,sd,arl,thc,depcalc,bdint,dbdint,
     &irad,nrad,drad
C
      nz2=nz/2
      j=2
      DO i=1,nz2
        depth=(float(i)-.5)*dz
40      if(depth.le.ssd(j).or.j.ge.nssd) goto 50
        j=j+1
        goto 40
50      frac=(depth-ssd(j-1))/(ssd(j)-ssd(j-1))
        ss2(i)=ss(j-1)+frac*(ss(j)-ss(j-1))
      END DO
C
C     Perform 1-2-1 filter over NIT iterations
      NIT=1
      DO 200 k=1,NIT
        tss1=ss2(1)
        DO 200 i=2,nz2-1
          tss2=.25*(ss2(i-1)+ss2(i)+ss2(i)+ss2(i+1))
          ss2(i-1)=tss1
          tss1=tss2
200   CONTINUE
      ss2(nz2-1)=tss1
C
      RETURN
      END
