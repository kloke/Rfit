      SUBROUTINE HFUNCD(N,Z,DERSC,CORR,X,Y,IWORK,XMAX)
C
C	EVALUATES THE FUNCTION H THAT IS USED FOR THE SCALE
C	ESTIAMTOR PROPOSED BY KOUL, SIEVERS, AND MCKEAN (1977,
C	SCAND. J. STAT.)
C
C	ON ENTRY
C
C          N         INTEGER. SAMPLE SIZE.
C          Z         REAL(N). VECTOR OF RESIDUALS.
C          DERSC     REAL(N). DERIVATIVE OF SCORE FUNCTION.
C                    EVALUATED AT THE ANTIRANKS OF RESIDS.
C          CORR      REAL. CONSTANT TO ENSURE H IS A CDF.
C          X         REAL. DOMAIN VALUE OF FUNCTION.
C          IWORK     INTEGER. 1 IF XMAX IS DESIRED.
C
C	ON RETURN
C
C	   Y	     double precision. VALUE OF THE FUNCTION.
C	   XMAX      double precision. MAXIMUM ABSOLUTE DIFFERENCE OF THE Z'S.
C
      INTEGER N,I,J,IP1,NM1,IWORK
      double precision ZN,CORR,X,XMAX,Y
      double precision Z(N),DERSC(N),S
      ZN = DBLE(N)
      Y = 0.0D0
      XMAX = 0.0D0
      NM1 = N - 1
      DO 10 I = 1, NM1
        IP1 = I + 1
        DO 5 J = IP1, N
          S = DABS(Z(I) - Z(J))
          IF(IWORK .EQ. 1) THEN
            IF(S .GT. XMAX) XMAX = S
          ENDIF
          IF(S .LE. X) Y = Y + DERSC(I) + DERSC(J)
5       CONTINUE
10    CONTINUE
      Y = (CORR/(ZN**2))*Y
      RETURN
      END
