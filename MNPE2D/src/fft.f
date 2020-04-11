c
c                       SUBROUTINE FFT
c
c This program replaces the vector Z=X+iY by its finite discrete complex
c Fourier transform. It performs as many base and iterations as possible
c and then finishes with a base 4 iteration or a base 2 iteration if needed.
c
c The subroutine is called as FFT(M,X,Y). The integer M (where N = 2**M),
c the N real location array X, and the N real location array Y must be 
c supplied to the subroutine.
c
c The execution time of this program on the G.E. 535 computer for N=1024 is
c approximately 0.6 seconds. The program extends the ideas of A.B. Langdons 
c 4+2program to reduce the computation time
c
c
c       
        SUBROUTINE FFT(N2POW,X,Y)
        DOUBLE PRECISION DT,TH
        DIMENSION X(1),Y(1),L(15),CS(32769),SS(32769)
        EQUIVALENCE (L15,L(1)), (L14,L(2)), (L13,L(3)), (L12,L(4)),
     1  (L11,L(5)), (L10,L(6)), (L9,L(7)), (L8,L(8)), (L7,L(9)),
     1  (L6,L(10)), (L5,L(11)), (L4,L(12)), (L3,L(13)), (L2,L(14)),
     1  (L1,L(15))
        NTHPO = 2**N2POW
        N8POW = N2POW/3
        IF (N8POW) 17,3,17
C       radix 8 passes, if any
17      LNG = NTHPO
        FLNG = LNG
        IF (LNG-LST) 171,173,171
171     LST = LNG
        DT = 6.28318530717958648D0/FLNG
        CS(1) = 1.
        SS(1) = 0.
        TH = 0.D0
        ND = LST/8
        DO 172  I=2,ND
           TH = TH+DT
           CS(I) = DCOS(TH)
172     SS(I) = DSIN(TH)
173     DO 1 IPASS=1,N8POW
           NXTLT = 2**(N2POW-3*IPASS)
           LENGT = 8*NXTLT
           SC = FLNG/FLOAT(LENGT)
           DO 1 J=1,NXTLT
              NRG = FLOAT(J-1)*SC+1.5
              C1 = CS(NRG)
              S1 = SS(NRG)
              C2 = C1**2-S1**2
              S2 = C1*S1+C1*S1
              C3 = C1*C2-S1*S2
              S3 = C2*S1+S2*C1
              C4 = C2**2-S2**2
              S4 = C2*S2+C2*S2
              C5 = C2*C3-S2*S3
              S5 = C3*S2+S3*C2
              C6 = C3**2-S3**2
              S6 = C3*S3+C3*S3
              C7 = C3*C4-S3*S4
              S7 = C4*S3+S4*C3
              DO 1 ISQLO = LENGT, NTHPO, LENGT
                 J0 = ISQLO-LENGT+J
                 J1 = J0+NXTLT
                 J2 = J1+NXTLT
                 J3 = J2+NXTLT
                 J4 = J3+NXTLT
                 J5 = J4+NXTLT
                 J6 = J5+NXTLT
                 J7 = J6+NXTLT
                 AR0 = X(J0)+X(J4)
                 AR1 = X(J1)+X(J5)
                 AR2 = X(J2)+X(J6)
                 AR3 = X(J3)+X(J7)
                 AR4 = X(J0)-X(J4)
                 AR5 = X(J1)-X(J5)
                 AR6 = X(J2)-X(J6)
                 AR7 = X(J3)-X(J7)
                 AI0 = Y(J0)+Y(J4)
                 AI1 = Y(J1)+Y(J5)
                 AI2 = Y(J2)+Y(J6)
                 AI3 = Y(J3)+Y(J7)
                 AI4 = Y(J0)-Y(J4)
                 AI5 = Y(J1)-Y(J5)
                 AI6 = Y(J2)-Y(J6)
                 AI7 = Y(J3)-Y(J7)
                 BR0 = AR0+AR2
                 BR1 = AR1+AR3
                 BR2 = AR0-AR2
                 BR3 = AR1-AR3
                 BR4 = AR4-AI6
                 BR5 = AR5-AI7
                 BR6 = AR4+AI6
                 BR7 = AR5+AI7
                 BI0 = AI0+AI2  
                 BI1 = AI1+AI3  
                 BI2 = AI0-AI2  
                 BI3 = AI1-AI3  
                 BI4 = AI4+AR6  
                 BI5 = AI5+AR7  
                 BI6 = AI4-AR6  
                 BI7 = AI5-AR7  
                 X(J0) = BR0+BR1
                 Y(J0) = BI0+BI1
                 IF(J-1) 2,2,18
18               X(J1) = C4*(BR0-BR1)-S4*(BI0-BI1)
                 Y(J1) = C4*(BI0-BI1)+S4*(BR0-BR1)
                 X(J2) = C2*(BR2-BI3)-S2*(BI2+BR3)
                 Y(J2) = C2*(BI2+BR3)+S2*(BR2-BI3)
                 X(J3) = C6*(BR2+BI3)-S6*(BI2-BR3)
                 Y(J3) = C6*(BI2-BR3)+S6*(BR2+BI3)
                 TR = 0.7071067812*(BR5-BI5)
                 TI = 0.7071067812*(BR5+BI5)
                 X(J4) = C1*(BR4+TR)-S1*(BI4+TI)
                 Y(J4) = C1*(BI4+TI)+S1*(BR4+TR)
                 X(J5) = C5*(BR4-TR)-S5*(BI4-TI)
                 Y(J5) = C5*(BI4-TI)+S5*(BR4-TR)
                 TR = -0.7071067812*(BR7+BI7)
                 TI = 0.7071067812*(BR7-BI7)
                 X(J6) = C3*(BR6+TR)-S3*(BI6+TI)
                 Y(J6) = C3*(BI6+TI)+S3*(BR6+TR)
                 X(J7) = C7*(BR6-TR)-S7*(BI6-TI)
                 Y(J7) = C7*(BI6-TI)+S7*(BR6-TR)
                 GO TO 1

2                X(J1) = BR0-BR1
                 Y(J1) = BI0-BI1
                 X(J2) = BR2-BI3
                 Y(J2) = BI2+BR3
                 X(J3) = BR2+BI3
                 Y(J3) = BI2-BR3
                 TR = 0.7071067812*(BR5-BI5)
                 TI = 0.7071067812*(BR5+BI5)
                 X(J4) = BR4+TR
                 Y(J4) = BI4+TI
                 X(J5) = BR4-TR
                 Y(J5) = BI4-TI
                 TR = -0.7071067812*(BR7+BI7)   
                 TI = 0.7071067812*(BR7-BI7)
                 X(J6) = BR6+TR
                 Y(J6) = BI6+TI
                 X(J7) = BR6-TR
                 Y(J7) = BI6-TI
1       CONTINUE
        
C       is there a four factor left

3       IF(N2POW-3*N8POW-1) 5,6,7
        
C       go through the base 2 iteration
6       DO 61 J=1,NTHPO,2
           R1 = X(J) +X(J+1)
           X(J+1) = X(J)-X(J+1)
           X(J) = R1
           FI1 = Y(J)+Y(J+1)
           Y(J+1) = Y(J)-Y(J+1)
61         Y(J) = FI1
        GO TO 5

C       go through the base 4 iteration

7       DO 71 J1=1,NTHPO,4
           J2 = J1+1
           J3 = J2+1
           J4 = J3+1
           R1 = X(J1)+X(J3)
           R2 = X(J1)-X(J3)
           R3 = X(J2)+X(J4)
           R4 = X(J2)-X(J4)
           FI1 = Y(J1)+Y(J3)
           FI2 = Y(J1)-Y(J3)
           FI3 = Y(J2)+Y(J4)
           FI4 = Y(J2)-Y(J4)
           X(J1) = R1+R3
           Y(J1) = FI1+FI3
           X(J2) = R1-R3
           Y(J2) = FI1-FI3
           X(J3) = R2-FI4
           Y(J3) = FI2+R4
           X(J4) = R2+FI4
71         Y(J4) = FI2-R4
5       DO 51 J=1,15
           L(J) = 1
           IF(J-N2POW) 19,19,51
19         L(J) = 2**(N2POW+1-J)
51      CONTINUE
        IJ = 1
        DO 8 J1=1,L1
         DO 8 J2=J1,L2,L1
          DO 8 J3=J2,L3,L2
           DO 8 J4=J3,L4,L3
            DO 8 J5=J4,L5,L4
             DO 8 J6=J5,L6,L5
              DO 8 J7=J6,L7,L6
               DO 8 J8=J7,L8,L7
                DO 8 J9=J8,L9,L8
                 DO 8 J10=J9,L10,L9
                  DO 8 J11=J10,L11,L10
                   DO 8 J12=J11,L12,L11
                    DO 8 J13=J12,L13,L12
                     DO 8 J14=J13,L14,L13
                      DO 8 JI=J14,L15,L14
                         IF (IJ-JI) 20,8,8
20                       R = X(IJ)
                         X(IJ)=X(JI)
                         X(JI)=R
                         FI=Y(IJ)
                         Y(IJ)=Y(JI)
                         Y(JI)=FI
8                        IJ=IJ+1

        RETURN
        END
