      double precision FUNCTION HUBCOR(N,RESIDS,SCINIT,PARAM)
C
C	RETURNS A DEGREE OF FREEDOM CORRECTION FOR A SCALE ESTIMATE
C	WHICH WAS SUGGESTED BY HUBER (1973,ANN. STAT.).
C
C	ON ENTRY
C
C	   N        INTEGER. SAMPLE SIZE.
C          RESIDS   REAL(N). VECTOR OF RESIDUALS.
C          SCINIT   REAL. ESTIMATE OF SCALE.
C	   PARAM    REAL. PARAMETER THAT DETERMINES THE CORRECTION.
C
      INTEGER N,I
      double precision SCINIT
      double precision RESIDS(N),PARAM
      HUBCOR = 0.0D0
      DO 10 I = 1, N
        IF(DABS(RESIDS(I)/SCINIT) .LT. PARAM) HUBCOR = HUBCOR + 1.0D0
10    CONTINUE
      HUBCOR = HUBCOR/DBLE(N)
      RETURN
      END
