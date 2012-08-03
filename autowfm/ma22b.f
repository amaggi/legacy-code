C             @(#)ma22b.f	1.5      9/4/86
      SUBROUTINE MA22B(A,IA,N,W,E)
*
*     Sous-Programme d'Inversion de Matrice Symetrique
*                                           ----------
      DIMENSION A(1),W(1)
      CALL A22M(A,A,IA,N,W,E)
      RETURN
      END
*---------------------------------------------
      SUBROUTINE A22M(A,B,IA,N,W,E)
      DIMENSION A(IA,1),W(1),B(1)
*          le tableau W semble ne servir a rien dans ce ss-programme !
      L=0
      EPS=1.0E-15
C     PRINT 401
401   FORMAT('A22M')
C     CALL TIME
      DO 1 J=1,N
      DO 1 I=1,J
      L=L+1
      B(L)=A(I,J)
1     CONTINUE
C     PRINT 402
402   FORMAT('SINV')
C     CALL TIME
      CALL SINV(A,N,EPS,IERR)
C     PRINT 402
      E=IERR
      II=N
      JJ=N
20    II=JJ
24    A(II,JJ)=B(L)
      L=L-1
      II=II-1
      IF(II)22,22,24
22    JJ=JJ-1
      IF(JJ)30,30,20
30    DO 40  J=1,N
      DO 40 I=1,J
40    A(J,I)=A(I,J)
C     PRINT 401
C     CALL TIME
      RETURN
      END
*---------------------------------------------------
      SUBROUTINE SINV(A,N,EPS,IER)
      REAL*8 DIN,WORK
      DIMENSION A(1)
      CALL MFSD(A,N,EPS,IER)
C     CALL TIME
      IF (IER) 9,1,1
1     IPIV=N*(N+1)/2
      IND=IPIV
      DO 6 I=1,N
      DIN=1.D0/DBLE(A(IPIV))
      A(IPIV)=DIN
      MIN=N
      KEND=I-1
      LANF=N-KEND
      IF(KEND)5,5,2
2     J=IND
      DO 4 K=1,KEND
      WORK=0.D0
      MIN=MIN-1
      LHOR=IPIV
      LVER=J
      DO 3 L=LANF,MIN
      LVER=LVER+1
      LHOR=LHOR+L
3     WORK=WORK+DBLE(A(LVER)*A(LHOR))
      A(J)=-WORK*DIN
4     J=J-MIN
5     IPIV=IPIV-MIN
6     IND=IND-1
C     CALL TIME
      DO 8 I=1,N
      IPIV=IPIV+I
      J=IPIV
      DO 8 K=I,N
      WORK=0.D0
      LHOR=J
      DO 7 L=K,N
      LVER=LHOR+K-I
      WORK=WORK+DBLE(A(LHOR)*A(LVER))
7     LHOR=LHOR+L
      A(J)=WORK
8     J=J+K
9     RETURN
      END
*---------------------------------------------------
      SUBROUTINE MFSD(A,N,EPS,IER)
      REAL*8 DPIV,DSUM
      DIMENSION A(1)
      IF(N-1) 12,1,1
1     IER=0
      KPIV=0
      DO 11 K=1,N
      KPIV=KPIV+K
      IND=KPIV
      LEND=K-1
      TOL=ABS(EPS*A(KPIV))
      DO 11 I=K,N
      DSUM=0.D0
      IF(LEND) 2,4,2
2     DO 3 L=1,LEND
      LANF=KPIV-L
      LIND=IND-L
3     DSUM=DSUM+DBLE(A(LANF)*A(LIND))
4     DSUM=DBLE(A(IND))-DSUM
      IF(I-K) 10,5,10
5     IF(SNGL(DSUM)-TOL) 6,6,9
6     IF(DSUM) 12,12,7
7     IF(IER) 8,8,9
8     IER=K-1
9     DPIV=DSQRT(DSUM)
      A(KPIV)=DPIV
      DPIV=1.D0/DPIV
      GOTO 11
10    A(IND)=DSUM*DPIV
11    IND=IND+I
      IF(IER.NE.0) PRINT 21,IER
      RETURN
12    IER=-1
      PRINT 22
21    FORMAT (' RET',I4)
22    FORMAT (' BREAK')
      RETURN
      END
*----------------------------------------------------------
