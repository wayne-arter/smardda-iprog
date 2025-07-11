C
         FUNCTION RANDLX(PDUM)
C
       DIMENSION   ZR(1)
         CALL RANLUX(ZR,1)
         RANDLX=ZR(1)
         RETURN
         END
C
         SUBROUTINE RANLUX(RVEC,LENV)
C
C.        Subtract-and-borrow random number generator proposed by
C.        Marsaglia and Zaman, implemented by F. James with the name
C.        RCARRY in 1991, and later improved by Martin Luescher
C.        in 1993 to produce "Luxury Pseudorandom Numbers".
C.    Fortran 77 coded by F. James, 1993
C.  LUXURY LEVELS.
C.  ------ ------      The available luxury levels are:
C. level 0  (p=24): equivalent to the original RCARRY of Marsaglia
C.          and Zaman, very long period, but fails many tests.
C. level 1  (p=48): considerable improvement in quality over level 0,
C.          now passes the gap test, but still fails spectral test.
C. level 2  (p=97): passes all known tests, but theoretically still
C.          defective.
C. level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
C.          correlations have very small chance of being observed.
C. level 4  (p=389): highest possible luxury, all 24 bits chaotic.
C.!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C.!!  Calling sequences for RANLUX:                                  ++
C.!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
C.!!                   32-bit random floating point numbers between  ++
C.!!                   zero (not included) and one (also not incl.). ++
C.!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
C.!!               one 32-bit integer INT and sets Luxury Level LUX  ++
C.!!               which is integer between zero and MAXLEV, or if   ++
C.!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
C.!!               should be set to zero unless restarting at a break++
C.!!               point given by output of RLUXAT (see RLUXAT).     ++
C.!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
C.!!               which can be used to restart the RANLUX generator ++
C.!!               at the current point by calling RLUXGO.  K1 and K2++
C.!!               specify how many numbers were generated since the ++
C.!!               initialization with LUX and INT.  The restarting  ++
C.!!               skips over  K1+K2*E9   numbers, so it can be long.++
C.!!   A more efficient but less convenient way of restarting is by: ++
C.!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
C.!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
C.!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
C.!!                 32-bit integer seeds, to be used for restarting ++
C.!!      ISVEC must be dimensioned 25 in the calling program        ++
C.!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       DIMENSION   RVEC(LENV)
       DIMENSION   SEEDS(24),          ISEEDS(24),         ISDEXT(25)
       PARAMETER   (MAXLEV=4, LXDFLT=3)
       DIMENSION   NDSKIP(0:MAXLEV)
       DIMENSION   NEXT(24)
       PARAMETER   (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
       PARAMETER   (ITWO24=2**24, ICONS=2147483563)
       SAVE
     +   NOTYET,   I24,      J24,      CARRY,    SEEDS,    TWOM24,
     +   TWOM12,   LUXLEV
       SAVE
     +   NSKIP,    NDSKIP,   IN24,     NEXT,     KOUNT,    MKOUNT,
     +   INSEED
       INTEGER     LUXLEV
       LOGICAL     NOTYET
       DATA
     +   NOTYET,   LUXLEV,   IN24,     KOUNT,    MKOUNT /.TRUE.,
     +   LXDFLT,   0,        0,        0/
       DATA        I24,      J24,      CARRY/24, 10,       0./
C     default
C     Luxury Level   0     1     2   *3*    4
       DATA        NDSKIP/0, 24,       73,       199,      365 /
C.Corresponds to p=24    48    97   223   389
C     time factor 1     2     3     6    10   on slow workstation
C     1    1.5    2     3     5   on fast mainframe
C
C     NOTYET is .TRUE. if no initialization has been performed yet.
C     Default Initialization by Multiplicative Congruential
         IF (NOTYET) THEN
            NOTYET = .FALSE.
            JSEED = JSDFLT
            INSEED = JSEED
C.DONT            WRITE(*,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: '
C.DONT      +      ,JSEED
            LUXLEV = LXDFLT
            NSKIP = NDSKIP(LUXLEV)
            LP = NSKIP + 24
            IN24 = 0
            KOUNT = 0
            MKOUNT = 0
C.DONT            WRITE(*,'(A,I2,A,I4)')  
C.DONT      +      ' RANLUX DEFAULT LUXURY LEVEL =  ',
C.DONT      +      LUXLEV,'      p =',LP
            TWOM24 = 1.
            DO 1 I= 1, 24
            TWOM24 = TWOM24 * 0.5
            K = JSEED/53668
            JSEED = 40014*(JSEED-K*53668) -K*12211
            IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
            ISEEDS(I) = MOD(JSEED,ITWO24)
    1       CONTINUE
            TWOM12 = TWOM24 * 4096.
            DO 2 I= 1,24
            SEEDS(I) = REAL(ISEEDS(I))*TWOM24
            NEXT(I) = I-1
    2       CONTINUE
            NEXT(1) = 24
            I24 = 24
            J24 = 10
            CARRY = 0.
            IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
         ENDIF
C
C     The Generator proper: "Subtract-with-borrow",
C     as proposed by Marsaglia and Zaman,
C     Florida State University, March, 1989
C
         DO 4 IVEC= 1, LENV
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY
         IF (UNI .LT. 0.)  THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
         ELSE
            CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = NEXT(I24)
         J24 = NEXT(J24)
         RVEC(IVEC) = UNI
C     small numbers (with less than 12 "significant" bits) are "padded".
         IF (UNI .LT. TWOM12)  THEN
            RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
C     and zero is forbidden in case someone takes a logarithm
            IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
         ENDIF
C     Skipping to luxury.  As proposed by Martin Luscher.
         IN24 = IN24 + 1
         IF (IN24 .EQ. 24)  THEN
            IN24 = 0
            KOUNT = KOUNT + NSKIP
            DO 3 ISK= 1, NSKIP
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY
            IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
            ELSE
               CARRY = 0.
            ENDIF
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
    3       CONTINUE
         ENDIF
    4    CONTINUE
         KOUNT = KOUNT + LENV
         IF (KOUNT .GE. IGIGA)  THEN
            MKOUNT = MKOUNT + 1
            KOUNT = KOUNT - IGIGA
         ENDIF
         RETURN
C
C-----------------------------------------------------------------------
CL              1.         Entry input and float integer seeds from
C
         ENTRY RLUXIN(ISDEXT)
         NOTYET = .FALSE.
         TWOM24 = 1.
         DO 100 I= 1, 24
         NEXT(I) = I-1
  100    TWOM24 = TWOM24 * 0.5
         NEXT(1) = 24
         TWOM12 = TWOM24 * 4096.
C.DONT         WRITE(*,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 ',
C.DONT      +   'INTEGERS:'
C.DONT         WRITE(*,'(5X,5I12)') ISDEXT
         DO 101 I= 1, 24
         SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  101    CONTINUE
         CARRY = 0.
         IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
         ISD = IABS(ISDEXT(25))
         I24 = MOD(ISD,100)
         ISD = ISD/100
         J24 = MOD(ISD,100)
         ISD = ISD/100
         IN24 = MOD(ISD,100)
         ISD = ISD/100
         LUXLEV = ISD
         IF (LUXLEV .LE. MAXLEV) THEN
            NSKIP = NDSKIP(LUXLEV)
C.DONT            WRITE (42,'(A,I2)') 
C.DONT      +      ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',LUXLEV
         ELSE  IF (LUXLEV .GE. 24) THEN
            NSKIP = LUXLEV - 24
C.DONT            WRITE (42,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:'
C.DONT      +      ,LUXLEV
         ELSE
            NSKIP = NDSKIP(MAXLEV)
C.DONT            WRITE (42,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: '
C.DONT     +      ,LUXLEV
            LUXLEV = MAXLEV
         ENDIF
         INSEED = -1
         RETURN
C
C-----------------------------------------------------------------------
CL              2.         Entry to output seeds as integers
C
         ENTRY RLUXUT(ISDEXT)
         DO 200 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  200    CONTINUE
         ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
         IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
         RETURN
C
C-----------------------------------------------------------------------
CL              3.         Entry output the "convenient" restart points
C
         ENTRY RLUXAT(LOUT,INOUT,K1,K2)
         LOUT = LUXLEV
         INOUT = INSEED
         K1 = KOUNT
         K2 = MKOUNT
         RETURN
C
C-----------------------------------------------------------------------
CL              4.         Entry initialize from one or three integers
C
         ENTRY RLUXGO(LUX,INS,K1,K2)
         IF (LUX .LT. 0) THEN
            LUXLEV = LXDFLT
         ELSE IF (LUX .LE. MAXLEV) THEN
            LUXLEV = LUX
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
            LUXLEV = MAXLEV
C.DONT            WRITE (42,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
         ELSE
            LUXLEV = LUX
            DO 400 ILX= 0, MAXLEV
            IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
  400       CONTINUE
         ENDIF
         IF (LUXLEV .LE. MAXLEV)  THEN
            NSKIP = NDSKIP(LUXLEV)
C.DONT            WRITE(*,'(A,I2,A,I4)') 
C.DONT      +      ' RANLUX LUXURY LEVEL SET BY RLUXGO:',
C.DONT      +      LUXLEV,'     p=', NSKIP+24
         ELSE
            NSKIP = LUXLEV - 24
C.DONT            WRITE (42,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:'
C.DONT      +      ,LUXLEV
         ENDIF
         IN24 = 0
C.DONT         IF (INS .LT. 0)  WRITE (42,'(A)')
C.DONT      +   ' Illegal initialization by RLUXGO, negative input seed'
         IF (INS .GT. 0)  THEN
            JSEED = INS
C.DONT            WRITE(*,'(A,3I12)') 
C.DONT      +      ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
C.DONT      +      JSEED, K1,K2
         ELSE
            JSEED = JSDFLT
C.DONT            WRITE(*,'(A)')
C.DONT      +      ' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
         ENDIF
         INSEED = JSEED
         NOTYET = .FALSE.
         TWOM24 = 1.
         DO 401 I= 1, 24
         TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  401    CONTINUE
         TWOM12 = TWOM24 * 4096.
         DO 402 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
  402    CONTINUE
         NEXT(1) = 24
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
C     If restarting at a break point, skip K1 + IGIGA*K2
C     Note that this is the number of numbers delivered to
C     the user PLUS the number skipped (if luxury .GT. 0).
         KOUNT = K1
         MKOUNT = K2
         IF (K1+K2 .NE. 0)  THEN
            DO 404 IOUTER= 1, K2+1
            INNER = IGIGA
            IF (IOUTER .EQ. K2+1)  INNER = K1
            DO 403 ISK= 1, INNER
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY
            IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
            ELSE
               CARRY = 0.
            ENDIF
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
  403       CONTINUE
  404       CONTINUE
C     Get the right value of IN24 by direct calculation
            IN24 = MOD(KOUNT, NSKIP+24)
            IF (MKOUNT .GT. 0)  THEN
               IZIP = MOD(IGIGA, NSKIP+24)
               IZIP2 = MKOUNT*IZIP + IN24
               IN24 = MOD(IZIP2, NSKIP+24)
            ENDIF
C     Now IN24 had better be between zero and 23 inclusive
            IF (IN24 .GT. 23) THEN
C.DONT               WRITE (42,'(A/A,3I11,A,I5)')
C.DONT      +         '  Error in RESTARTING with RLUXGO:','  The values', INS,
C.DONT      +         K1, K2, ' cannot occur at luxury level', LUXLEV
               IN24 = 0
            ENDIF
         ENDIF
         RETURN
         END
