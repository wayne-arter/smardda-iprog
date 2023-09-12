C/ MODULE ISNUMB
C
         LOGICAL FUNCTION ISNUMB(KCSTR,KIN,KOFF)
C
C     DOES STRING BEGIN WITH AN INTEGER
       CHARACTER   *(*) KCSTR
       DIMENSION   IDIG(9)
         ISNUMB=.FALSE.
         KIN=0
C
C-----------------------------------------------------------------------
CL              1.         Check argument
C
         IOFF=KOFF-1
         ISLEN=LEN(KCSTR)-IOFF
         IF (ISLEN.LE.0) RETURN
C
C-----------------------------------------------------------------------
CL              2.         Find first non-blank
C
         I0=ICHAR('0')
         I9=ICHAR('9')
         DO 200 J=1,ISLEN
         IJ=J+IOFF
         IF (KCSTR(IJ:IJ).EQ.' '.OR.KCSTR(IJ:IJ).EQ.',') GO TO 200
         ICH=ICHAR(KCSTR(IJ:IJ))
         IF (ICH.GE.I0.AND.ICH.LE.I9) GO TO 201
         GO TO 205
  200    CONTINUE
         RETURN
C
  201    CONTINUE
         I=1
         IDIG(I)=ICH-I0
         DO 202 J=IJ+1,ISLEN
         ICH=ICHAR(KCSTR(J:J))
         IF (ICH.LT.I0.OR.ICH.GT.I9) GO TO 203
         I=I+1
         IDIG(I)=ICH-I0
  202    CONTINUE
         I=MIN(9,I)
  203    CONTINUE
         KIN=IDIG(1)
         DO 204 J=2,I
         KIN=KIN*10+IDIG(J)
  204    CONTINUE
         ISNUMB=.TRUE.
         RETURN
C
  205    CONTINUE
         RETURN
         END
