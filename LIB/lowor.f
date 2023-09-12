C
         SUBROUTINE LOWOR(KCSTR,KSTA,KEND)
C
C     Modify for testing
       CHARACTER   *(*) KCSTR
C     Check arguments
         ISLEN=LEN(KCSTR)
         IEND=MAX(KSTA,KEND)
         IF (KSTA.LT.1.OR.IEND.GT.ISLEN) THEN
            WRITE(*,*) KSTA,KEND,ISLEN,'*** RUN ABANDONED IN LOWOR'
            STOP
         END IF
C
         IALC=ICHAR('a')
         IAUC=ICHAR('A')
         IZUC=ICHAR('Z')
         IAOFF=IALC-IAUC
         DO 1 J=KSTA,KEND
         INCH=ICHAR(KCSTR(J:J))
         IF (INCH.GE.IAUC.AND.INCH.LE.IZUC)
     +   KCSTR(J:J)=CHAR(INCH+IAOFF)
    1    CONTINUE
         RETURN
         END
