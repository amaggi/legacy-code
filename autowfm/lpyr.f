      INTEGER FUNCTION LPYR(YEAR)
C
C Function lpyr determines if year
C is a leap year.
C
C This function uses the intrinsic
C function mod. If your machine
C does not supply this function,
C make one -
C mod(i,j) = iabs(i - (i/j)*j)
C
C
C Calls:
C   mod - intrinsic funtion
C
C      Programmed by Madeleine Zirbes
C         September 15,1980
C
C YEAR - INPUT
      INTEGER YEAR
      IF (.NOT.(MOD(YEAR, 400) .EQ. 0)) GOTO 2000
        LPYR = (1)
        RETURN
2000  CONTINUE
      IF (.NOT.(MOD(YEAR, 4) .NE. 0)) GOTO 2020
        LPYR = (0)
        RETURN
2020  CONTINUE
      IF (.NOT.(MOD(YEAR, 100) .EQ. 0)) GOTO 2040
        LPYR = (0)
        RETURN
2040  CONTINUE
      LPYR = (1)
      RETURN
      END

