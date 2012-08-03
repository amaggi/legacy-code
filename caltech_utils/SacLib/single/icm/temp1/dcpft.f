      SUBROUTINE DCPFT(R, I, N, INCP, ISIGNP)
C
C DOUBLE PRECISION VERSION.
C  FORTRAN TRANSLITERATION OF SINGLETON'S 6600 ASSEMBLY-CODED FFT.
C  INTENDED TO BE OF ASSISTANCE IN UNDERSTANDING HIS CODE, AND IN
C  FUTURE WRITING OF AN FFT FOR ANOTHER MACHINE.
C  IT SHOULD BE TRANSLATED INTO MACHINE CODE RATHER THAN USED FOR
C  PRODUCTION AS IS BECAUSE- IT IS VERSATILE AND EFFICIENT ENOUGH TO
C  SEE LOTS OF USE, IT BENEFITS GREATLY FROM CAREFUL HAND CODING, AND IS
C  SHORT AND SIMPLE ENOUGH TO DO QUICKLY.
C  A. BRUCE LANGDON, M DIVISION, L.L.L., 1971.
C
C  COMMENTS BELOW ARE MOSTLY FROM 6600-7600 VERSION.
C  R      REAL PART OF DATA VECTOR.
C  I      IMAG PART OF DATA VECTOR.
C  N      NUMBER OF ELEMENTS (=1,2,4,8...32768).
C  INC    SPACING IN MEMORY OF DATA (USUALLY 1, BUT SEE BELOW).
C  SIGN   ITS SIGN WILL BE SIGN OF ARGUMENT IN TRANSFORM EXPONENTIAL.
C
C    ON ENTRY ARRAYS R AND I CONTAIN THE SEQUENCE TO BE TRANSFORMED.
C  ON EXIT THEY CONTAIN THE TRANSFORM. INPUT AND OUTPUT SEQUENCES ARE
C  BOTH IN NATURAL ORDER (I.E. NOT BIT-REVERSED SCRAMBLED).
C
C    A CALL TO CPFT WITH SIGN=+1, FOLLOWED BY ANOTHER CALL WITH THE
C  FIRST 4 PARAMETERS THE SAME AND SIGN=-1, WILL LEAVE R AND I WITH
C  THEIR ORIGINAL CONTENTS TIMES N. THE SAME IS TRUE IF FIRST SIGN=-1,
C  AND NEXT SIGN=+1.
C
C    THE USEFULNESS OF PARAMETER INC MAY BE ILLUSTRATED BY 2 EXAMPLES:
C    SUPPOSE THE COMPLEX SEQUENCE IS STORED AS A FORTRAN COMPLEX ARRAY
C  Z, I.E. REAL AND IMAGINARY PARTS IN ALTERNATE MEMORY CELLS. THE
C  SEPARATION BETWEEN CONSECUTIVE REAL (OR IMAGINARY) ELEMENTS IS 2
C  WORDS, SO INC=2. THE CALL MIGHT BE
C          CALL CPFT(REAL(Z), AIMAG(Z), N, 2, SIGN)
C  FOR MANY COMPILERS ONE WOULD INSTEAD HAVE TO DO SOMETHING LIKE
C          CALL CPFT(RI, RI(2), N, 2, SIGN)
C  WHERE RI IS A REAL ARRAY EQUIVALENCED TO Z.
C    SUPPOSE ONE HAD AN ARRAY C WITH DIMENSIONS N1, N2. ONE WANTS R TO
C  BE ROW I1 AND I TO BE ROW I2. THE SEPARATION OF CONSECUTIVE ELEMENTS
C  IS N1 AND STARTING ADDRESSES ARE C(I1,1) AND C(I2,1), SO USE
C          CALL CPFT(C(I1,1), C(I2,1), N2, N1, SIGN)
C
C  TIMING, ASSUMING MINIMAL MEMORY BANK CONFLICTS:
C    6400 TIME FOR N=1024 IS 220,000=21.5*N*LOG2N MICROSECONDS.
C    6600 TIME FOR N=1024 IS  44,500=4.35*N*LOG2N MICROSECONDS.
C    7600 TIME FOR N=1024 IS   8,300=0.81*N*LOG2N MICROSECONDS.
C
C    A RADIX 2 FFT PROVOKES MEMORY BANK CONFLICTS AT BEST, BUT TIMING
C  IS NOTICEABLY WORSENED WHEN LIKE ELEMENTS OF R AND I ARE IN THE SAME
C  BANK AND/OR INC IS A MULTIPLE OF A POWER OF 2. IN A WORST CASE ON
C  THE 7600 THE SPEED WAS DECREASED BY A FACTOR OF 3.
C    THUS IN THE EXAMPLE ABOVE, IF N1=MULTIPLE OF 32 ONE MIGHT
C  DECIDE TO WASTE A LITTLE MEMORY BY INCREASING N1 TO 33, THUS
C  DECREASING CONFLICTS FOR TRANSFORMS OVER ROWS OR OVER COLUMNS.
C
C  WRITTEN BY R. C. SINGLETON, STANFORD RESEARCH INSTITUTE, NOV. 1968.
C  COMMENTARY, LRL LINKAGE AND OTHER MINOR CHANGES BY A. BRUCE LANGDON
C  LAWRENCE RADIATION LABORATORY, LIVERMORE, APRIL 1971.
C
C  REFERENCES:
C    (1) R. C. SINGLETON, 'ON COMPUTING THE FAST FOURIER TRANSFORM',
C        COMM. ASSOC. COMP. MACH. VOL. 10, PP. 647-654 (1967).
C    (2) R. C. SINGLETON, ALGORITHM 345 'AN ALGOL CONVOLUTION PROCEDURE
C        BASED ON THE FAST FOURIER TRANSFORM', COMM. ACM VOL. 12,
C        PP. 179-184 (1969).
C    (3) W. M. GENTLEMAN AND G. SANDE, 'FAST FOURIER TRANSFORMS - FOR
C        FUN AND PROFIT', PROC. AFIPS 1966 FALL JOINT COMPUTER CONF.,
C        VOL. 29, PP. 563-578.
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 R(1), I(1)
      INTEGER SPAN, RC
      REAL*8 SINES(15), I0, I1
C
C  TABLE OF SINES.
C    THESE SHOULD BE GOOD TO THE VERY LAST BIT. THEY ARE GIVEN IN OCTAL
C  TO PREVENT AN ASSEMBLER FROM CONVERTING THEM POORLY. THEY MAY BE
C  OBTAINED BY EVALUATING THE INDICATED SINES IN DOUBLE PRECISION AND
C  PUNCHING THEM OUT IN OCTAL FORMAT (A SINGLE PRECISION SINE ROUTINE
C  IS NOT ACCURATE ENOUGH). USE THE MOST SIGNIFICANT WORD, ROUNDED
C  ACCORDING TO THE LEAST SIGNIFICANT WORD.
C  IN THIS VERSION I DO IT A LAZY WAY ON THE FIRST CALL.
      DATA SINES(1)/0.0D0/
      IF( SINES(1).EQ.1.0D0 ) GO TO 1
      SINES(1)=1.0D0
      T=DATAN(1.0D0)
      DO 2 IS=2,15
      SINES(IS)=DSIN(T)
    2 T=T/2.0D0
    1 CONTINUE
C
      IF( N.EQ.1 ) RETURN
C
C  SET UP VARIOUS INDICES.
C
      INC=INCP
      SGN=ISIGNP
      NINC=N*INC
      SPAN=NINC
      IT=N/2
      DO 3 IS=1,15
      IF( IT.EQ.1 ) GO TO 12
    3 IT=IT/2
C
C    THERE ARE 2 INNER LOOPS WHICH RUN OVER THE N/(2*SPAN) REPLICATIONS
C  OF TRANSFORMS OF LENGTH (2*SPAN). THESE LOOPS FIT INTO THE
C  INSTRUCTION STACK OF THE 6600 OR 7600. ONE LOOP IS FOR ARBITRARY
C  ROTATION FACTOR ANGLE. THE OTHER TAKES CARE OF THE SPECIAL CASE IN
C  WHICH THE ANGLE IS ZERO SO THAT NO COMPLEX MULTIPLICATION IS NEEDED.
C  THIS IS MORE EFFICIENT THAN TESTING AND BRANCHING INSIDE THE INNER
C  LOOP, AS IS OFTEN DONE. THE OTHER SPECIAL CASE IN WHICH NO COMPLEX
C  MULTIPLY IS NEEDED IS ANGLE=PI (I.E. FACTOR=I); THIS IS NOT HANDLED
C  SPECIALLY. THESE MEASURES ARE MOST HELPFUL FOR SMALL N.
C
C    THE ORGANIZATION OF THE RECURSION IS THAT OF SANDE (REF. (3),
C  PP. 566-568). THAT IS, THE DATA IS IN NORMAL ORDER TO START AND
C  SCRAMBLED AFTERWARD, AND THE EXPONENTIAL ROTATION ('TWIDDLE') FACTOR
C  ANGLES ARE USED IN ASCENDING ORDER DURING EACH RECURSION LEVEL.
C  ALL THE SINES AND COSINES NEEDED ARE GENERATED FROM A SHORT TABLE
C  USING A STABLE MULTIPLE-ANGLE RECURSION (REF. (1), P651 AND REF. (2),
C  PP. 179-180). THIS METHOD IS ECONOMICAL IN STORAGE AND TIME, AND
C  YIELDS ACCURACY COMPARABLE TO GOOD LIBRARY SIN-COS ROUTINES.
C  ANGLES BETWEEN 0 AND PI ARE NEEDED. THE RECURSION IS USED FOR
C  ANGLES UP TO PI/2; LARGER ANGLES ARE OBTAINED BY REFLECTION IN THE
C  IMAGINARY AXIS (ANGLE:=PI-ANGLE). THESE PAIRS OF ANGLES ARE USED
C  ONE RIGHT AFTER THE OTHER.
C
C    FOR SIMPLICITY, COMMENTARY BELOW APPLIES TO INC=1 CASE.
C
C  IF TRUNCATED RATHER THAN ROUNDED ARITHMETIC IS USED, SINGLETON'S
C  MAGNITUDE CORRECTION SHOULD BE APPLIED TO C AND S.
C
10    T=S+(S0*C-C0*S)
      C=C-(C0*C+S0*S)
      S=T
C  REPLICATION LOOP.
11    K1=K0+SPAN
      R0=R(K0+1)
      R1=R(K1+1)
      I0=I(K0+1)
      I1=I(K1+1)
      R(K0+1)=R0+R1
      I(K0+1)=I0+I1
      R0=R0-R1
      I0=I0-I1
      R(K1+1)=C*R0-S*I0
      I(K1+1)=S*R0+C*I0
      K0=K1+SPAN
      IF( K0.LT.NINC ) GO TO 11
      K1=K0-NINC
      C=-C
      K0=SPAN-K1
      IF( K1.LT.K0 ) GO TO 11
      K0=K0+INC
      IF( K0.LT.K1 ) GO TO 10
C  RECURSION TO NEXT LEVEL.
12    CONTINUE
      SPAN=SPAN/2
      K0=0
C  ANGLE=0 LOOP.
13    K1=K0+SPAN
      R0=R(K0+1)
      R1=R(K1+1)
      I0=I(K0+1)
      I1=I(K1+1)
      R(K0+1)=R0+R1
      I(K0+1)=I0+I1
      R(K1+1)=R0-R1
      I(K1+1)=I0-I1
      K0=K1+SPAN
      IF( K0.LT.NINC ) GO TO 13
C  ARE WE FINISHED...
      IF( SPAN.EQ.INC ) GO TO 20
C  NO. PREPARE NON-ZERO ANGLES.
      C0=2.0D0*SINES(IS)**2
      IS=IS-1
      S=DSIGN( SINES(IS),SGN )
      S0=S
      C=1.0D0-C0
      K0=INC
      GO TO 11
C
C    ARRAYS R AND I NOW CONTAIN TRANSFORM, BUT STORED IN 'REVERSE-
C    BINARY' ORDER. THE RE-ORDERING IS DONE BY PAIR EXCHANGES.
C    REFERENCE FOR SORTING PRINCIPLE IS P. 180 AND P. 182 OF REF. (2).
C
C    ONCE AGAIN, COMMENTARY APPLIES TO INC=1 CASE.
C  INDICES ARE:
C    IJ:=0,1,2...N/2-1 ( A SIMPLE COUNTER).
C    JI:=REVERSAL OF IJ.
C    RC:=REVERSAL OF 0,2,4...N/2 (INCREMENTED N/4 TIMES).
C  RC IS INCREMENTED THUSLY: STARTING WITH THE NEXT-TO-LEFTMOST BIT,
C  CHANGE EACH BIT UP TO AND INCLUDING FIRST 0. (THE ACTUAL CODING IS
C  DONE SO AS TO WORK FOR ANY INC>0 WITH EQUAL EFFICIENCY.)
C    FOR ALL EXCHANGES IJ FITS ONE OF THESE CASES:
C      (1) 1ST AND LAST BITS ARE 0 (IJ,JI EVEN AND <N/2), AND IJ<=JI.
C      (2) ONE'S COMPLEMENT OF CASE (1) (BOTH ODD AND >N/2).
C      (3) 1ST BIT 0, LAST BIT 1 (IJ ODD AND <N/2, JI>N/2).
C    THE CODE FROM LABEL EVEN DOWN TO ODD IS ENTERED WITH IJ EVEN AND
C  <=JI. FIRST TIME THRU THE COMPLEMENTS ARE DONE -CASE (2). SECOND
C  TIME THRU GETS CASE (1). THUS A PAIR OF ELEMENTS BOTH IN THE FIRST
C  HALF OF THE SEQUENCE, AND ANOTHER PAIR IN THE 2ND HALF, ARE
C  EXCHANGED. THE CONDITION IJ<JI PREVENTS A PAIR FROM BEING EXCHANGED
C  TWICE.
C    THE CODE FROM LABEL ODD DOWN TO INCREV DOES CASE (3).
C
20    N1=NINC-INC
      N2=NINC/2
      RC=0
      JI=RC
      IJ=JI
      IF( N2.EQ.INC ) RETURN
      GO TO 22
C
C  EVEN.
21    IJ=N1-IJ
      JI=N1-JI
      T=R(IJ+1)
      R(IJ+1)=R(JI+1)
      R(JI+1)=T
      T=I(IJ+1)
      I(IJ+1)=I(JI+1)
      I(JI+1)=T
      IF( IJ.GT.N2 ) GO TO 21
C  ODD.
22    IJ=IJ+INC
      JI=JI+N2
      T=R(IJ+1)
      R(IJ+1)=R(JI+1)
      R(JI+1)=T
      T=I(IJ+1)
      I(IJ+1)=I(JI+1)
      I(JI+1)=T
      IT=N2
C  INCREMENT REVERSED COUNTER.
23    IT=IT/2
      RC=RC-IT
      IF( RC.GE.0 ) GO TO 23
      RC=RC+2*IT
      JI=RC
      IJ=IJ+INC
      IF( IJ.LT.JI ) GO TO 21
      IF( IJ.LT.N2 ) GO TO 22
C
      RETURN
      END

