      SUBROUTINE NSCALE(N,EPSA,ALP,SCINIT,IFLAG,RANK,RESIDS,JRANK,
     .       SCORER,SCR1,SCR2,TAUHAT,SCOPAR,ISCOFN,MAXNO,ECDFAC,
     .         WORK5,HCPARM)
C
C       SUBROUTINE TO OBTAIN AN R-ESTIMATE OF SCALE PROPOSED
C       BY KOUL, SIEVERS, AND MCKEAN (1987,SCAND. J. STATIST.)
C       THE ESTIMATE IS OF THE FORM:
C                 GAMMA = HFUN(THAT/SQRT(N))/(2*(THAT/SQRT(N)))
C       WHERE HFUN(T) = (N**-2)*SUM(PHIPR(R(J))*I(ABS(RESID(J) - RESID(I))
C                                                        <T))
C       AND THE VALUE THAT SOLVES: HFUN(THAT) = ALP.
C       THE ACTUAL ESTIMATE RETURNED IS:
C                     TAUHAT = A/GAMMA.
C       NOTE THE R-SCORES IN RSCORER ARE ASSUMED TO SATISFY INTEG(PHI) = 0
C       AND INTEG(PHI**2) = 1.
C
C	ON ENTRY
C
C          N        INTEGER. SAMPLE SIZE.
C          JRANK    INTEGER(N). VECTOR OF ANTI-RANKS OF RESIDUALS.
C          RANK     INTEGER. RANK OF THE DESIGN MATRIX.
C          MAXNO    INTEGER. MAXIMUM NUMBER OF ITERATIONS IN LIN2.
C          EPSA     REAL. TOLERANCE.
C          ALP      REAL. ALPHA LEVEL OF ESTIMATE.
C          SCINIT   REAL. INITIAL ESTIMATE OF SCALE.
C          RESIDS   REAL(N). VECTOR OF RESIDUALS.
C          SCORER   REAL(N). VECTOR OF SCORES.
C          SCOPAR   REAL(5). VECTOR OF SCORE PARAMETERS.
C          SCR1     REAL(N). SCRATCH.
C          SCR2     REAL(N). SCRATCH.
C          HCPARM   REAL. DF CORRECTION FOR SCALE ESTIMATE.
C
C	ON RETURN
C
C          IFLAG    INTEGER. HAS THE VALUE ZERO IF LIN2 CONVERGED WITHIN MAXNO
C                   OF STEPS; OTHERWISE IT'S > 0.
C          TAUHAT   REAL. DESIRED ESTIMATE OF SCALE.
C
C        FUNCTIONS AND SUBROUTINES CALLED
C
C          HUBCOR
C          HFUNCD
C          LIN3
C
      INTEGER N,IFLAG,MAXNO,I,ISCOFN
c     double precision txsee
      double precision ZN,ZP,PHI1,PHI0,ECDFAC,EPSA,ALP,SCINIT,TAUHAT
      double precision CONST,GAM0,SOL,TEEH,HATT,TEMP,CORR,HUBCOR,XMAX
      double precision RESIDS(N),SCOPAR(5),SCORER(N),SCR1(N),TX
      double precision SCR2(N),WORK5(N),HCPARM
      INTEGER RANK,JRANK(N)
      ZN = DFLOAT(N)
      ZP = DFLOAT(RANK)
      IFLAG = 0
      PHI1 = SCORER(N)
      PHI0 = SCORER(1)
      ECDFAC = 0.0D0
C
C                        VECTOR OF DERIVATIVES OF SCORE FUNCTION,
C                            (STANDARDIZED).
C
      DO 10 I = 1, N
        SCR2(I) = WORK5(I)/(PHI1 - PHI0)
        ECDFAC = ECDFAC + SCR2(I)
10    CONTINUE
      ECDFAC = (ZN**2)/((ZN - 1.0D0)*ECDFAC)
C
C                         DERIVATIVES AT THE ANTI-RANKS.
C
      DO 20 I = 1, N
        SCR1(JRANK(I)) = SCR2(I)
20    CONTINUE
C
C                 SOLVING 1 - HFUN(T) = 1 - ALP
C
      CONST = ALP
      TX = SCINIT
      GAM0 = 1.0D0/(TX*(PHI1 - PHI0))
c      write(6,983) (resids(i),i=1,n)
c983   format('before lin3',5e16.6)
      CALL LIN3(N,EPSA,GAM0,CONST,MAXNO,RESIDS,SCR1,ECDFAC,IFLAG,SOL)
c      write(6,984) (resids(i),i=1,n)
c984   format('after lin3',5e16.6)
      IF(IFLAG .NE. 0) THEN
        IF(IFLAG .EQ. 1) THEN
          IFLAG = 21
        ELSE
          IFLAG = 22
        ENDIF
      ENDIF
      TEEH = SOL/DSQRT(DFLOAT(N))
      CALL HFUNCD(N,RESIDS,SCR1,ECDFAC,TEEH,HATT,0,XMAX)
C
C
      TX = HATT
      TAUHAT = ((2.0D0*TEEH)/TX)/(PHI1 - PHI0)
c      write(6,*) 'after sol: sol tauhat',sol,tauhat
C
C                 DEGREE OF FREEDOM CORRECTIONS
C
C                    STANDARD DEGREE OF FREEDOM CORRECTION
C
      TAUHAT = TAUHAT*DSQRT(DFLOAT(N)/DFLOAT(N - RANK))
c      write(6,*) 'after stan df: tauhat',tauhat
C
C                    A HUBER TYPE CORRECTION
C
C                        FOR POLICELLO SCORES
C
      IF((ISCOFN .EQ. 2) .OR. (ISCOFN .EQ. 3)) THEN
        TEMP = SCOPAR(2) - SCOPAR(1)
        CORR = 1.0D0 + ((ZP/ZN)*((1.0D0 - TEMP)/TEMP))
C
C                        FOR OTHER SCORES
C
      ELSE
c      write(6,985) (resids(i),i=1,n)
c985   format('before hub',5e16.6)
        CORR = HUBCOR(N,RESIDS,SCINIT,HCPARM)
c      write(6,986) (resids(i),i=1,n)
c986   format('after hub',5e16.6)
c       txsee = corr
        TX = CORR
        IF(TX .LT. EPSA) TX = EPSA
        CORR = 1.0D0 + ((ZP/ZN)*((1.0D0 - TX)/TX))
      ENDIF
      TAUHAT = TAUHAT*CORR
c      write(6,*) 'after huber df: tauhat tx scinit hcparm txsee',
c     .            tauhat,tx,scinit,hcparm,txsee
c      write(6,987) (resids(i),i=1,n)
c987   format(5e16.6)
      RETURN
      END
