
c @@@@@@@@@@@@@@@@@@ Liste des Subroutines @@@@@@@@@@@@@@@@@@@@
      SUBROUTINE MYGRT(ALAT1,ALON1,ALAT2,ALON2,DISD,AZ12,AZ21,GC)

c Great circle program. Given epicenter and station, finds distance,
c take-off and back- azimuths, and length of great circle

      PI = 4.*atan(1.0)
      ATH=6378.388
      BTH=6356.912
      RAD = PI/180.
      H = 1. - BTH*BTH/(ATH*ATH)
      P = H/(1. - H)
      GR = ALON1*RAD
      TR = ALAT1*RAD
      SINTR = SIN(TR)
      COSTR = COS(TR)
      IF (SINTR .EQ. 0.) SINTR = .000001
      IF (COSTR .EQ. 0.) COSTR = .000001
      R1 = ATH/SQRT(1. - H*SINTR*SINTR)
      Z1 = R1*(1. - H)*SINTR
      G = ALON2*RAD
      T = ALAT2*RAD
      IF (T .EQ. 0.) T = .00001
      SINT = SIN(T)
      COST = COS(T)
      R2 = ATH/SQRT(1. - H*SINT*SINT)
      DG = G - GR
      COSDG = COS(DG)
      COSDG = COS(DG)
      SINDG = SIN(DG)
      DGR = GR - G
      DT = T - TR
      Q = SINT*COSTR/((1. + P)*COST*SINTR) + H*R1*COSTR/(R2*COST)
      X = R2*COST*COSDG
      Y = R2*COST*SINDG
      Z = R2*(1. - H)*SINT
      AZ12 = ATAN2(SINDG,(Q - COSDG)*SINTR)
      Q = SINTR*COST/(COSTR*SINT*(1. + P)) + H*R2*COST/(R1*COSTR)
      AZ21 = ATAN2(SIN(DGR),SINT*(Q - COS(DGR)))
      COS12 = COS(AZ12)
      CTA2 = COSTR*COSTR*COS12*COS12
      P0 = P*(CTA2 + SINTR*SINTR)
      B0 = (R1/(1. + P0))*SQRT(1. + P*CTA2)
      E0 = P0/(1. + P0)
      GC = 2.*PI*B0*SQRT(1. + P0)*(1. - E0*(.25 + E0*(3./64.
     *                                          + 5.*E0/256.)))
      C0 = 1. + P0*(.25 - P0*(3./64 - 5.*P0/256.))
      C2 = P0*(-.125 + P0*(1./32. - 15.*P0/1024.))
      C4 = (-1./256. + 3.*P0/1024.)*P0*P0
      U0 = ATAN2(SINTR,COSTR*COS12*SQRT(1. + P0))
      U = ATAN2(R1*SINTR + (1. + P0)*(Z - Z1),(X*COS12 - Y*SINTR*
     *                                         SIN(AZ12))*SQRT(1. + P0))
      DISD = U - U0
      IF (U .LT. U0) DISD = PI + PI + DISD
      DIST = B0*(C0*( DISD ) + C2*(SIN(U + U) - SIN(U0 + U0))
     *                       + C4*(SIN(4.*U) - SIN(4.*U0)))
      DISD = DISD/RAD
      AZ12 = AZ12/RAD
      AZ21 = AZ21/RAD
      IF (AZ12 .LT. 0.) AZ12 = 360. + AZ12
      IF (AZ21 .LT. 0.) AZ21 = 360. + AZ21
        if(disd.gt.355.) disd=360.-disd
      RETURN
      END

c-------------------------------------------------------------

      SUBROUTINE MYGDS(ALAT,ALON,AZ,DIST,BLAT,BLON)

c Geodesic program. Given epicenter (alat,alon), azimuth of travel and distance
c traveled, finds arrival point (blat,blon);  dist in km

      ALT=ALAT
      ALN=ALON
      IF(ALAT.EQ.90.)ALAT=ALAT-0.000001
      IF(AZ.EQ.90.0.OR.AZ.EQ.270.0)AZ=AZ-0.000001
      A0=6378.388
      B0=6356.912
      F=(A0-B0)/A0
      PI=4.*ATAN(1.0)
      RAD=PI/180.0
      ALAT=ALAT*RAD
      ALON=ALON*RAD
      AZ=AZ*RAD
      E2=1.-(B0**2)/(A0**2)
      EPS=E2/(1.-E2)
      SIN1=SIN(ALAT)
      V=A0/(SQRT(1.-E2*(SIN1**2)))
      SINA=SIN(AZ)
      COSA=COS(AZ)
      COS1=COS(ALAT)
      TC2=(COSA*COSA)*(COS1*COS1)
      C2=TC2+(SIN1*SIN1)
      EPS0=C2*EPS
      TB1=SQRT(1.+EPS*TC2)
      B1=(V*TB1)/(1.+EPS0)
      TAN1=TAN(ALAT)
      IF(COSA.EQ.0.0)COSA=0.00001
      S1=SQRT(1.+EPS0)
      G0=1.-(EPS0/4.)+((7.*EPS0*EPS0)/64.)-((15.*(EPS0**3))/256.)
      G2=(EPS0/8.)-(.0625*EPS0*EPS0)+((145.*(EPS0**3))/2048.)
      G4=((5.*EPS0*EPS0)/256.)-((5.*(EPS0**3))/256.)
      G6=(29.*(EPS0**3))/6144.
      SIGM=(DIST*G0)/B1
      U1P=ATAN2(TAN1,(COSA*S1))
      SIN2P=SIN(2.*U1P)
      SIN4P=SIN(4.*U1P)
      TSS=((EPS0/4.)-((EPS0*EPS0)/8.))
      S12=2.*U1P-TSS*SIN2P-((EPS0*EPS0)/128.)*SIN4P
      SIGP=S12+SIGM
      T1=SIGM+(2.*G2*SIN(SIGM))*COS(SIGP)
      T2=(2.*G4*SIN(2.*SIGM))*COS(2.*SIGP)
      T3=(2.*G6*SIN(3.*SIGM))*COS(3.*SIGP)
      U2P=U1P+T1+T2+T3
      SINU1=TAN1/(SQRT(1.+EPS+(TAN1*TAN1)))
      C=SQRT(C2) 
      SINU2=(((B1*C)/B0)*SIN(U2P))-((EPS-EPS0)/(1.+EPS0))*SINU1
      U2=ASIN(SINU2)
      SINP1=SINU2/(SQRT(1.-E2*(COS(U2)*COS(U2))))
      BLAT=ASIN(SINP1)
      A1=B1*(SQRT(1.+EPS0))
      IF(COS(U2).EQ.0.0)U2=U2-0.00001
      Q1=(A1*COS(U2P))/(A0*COS(U2))
      if(q1.gt.-1..and.q1.lt.1.) Q2=ACOS(Q1)
      if(q1.eq.1) q2=0.
      if(q1.eq.-1.) q2=pi
      if(q1.gt.1.) go to 100
      if(q1.lt.-1.) go to 200
300   X1=SIN1*SINA
      AMU=ATAN2(X1,COSA)
      AZ=AZ/RAD
      U2P=U2P/RAD
      IF(AZ.GT.180.0)Q2=-Q2
      IF(U2P.GT.180.0.OR.U2P.LT.0.0)Q2=-Q2
      DLAMB=Q2-AMU
      BLON=DLAMB+ALON
      BLAT=BLAT/RAD
      BLON=BLON/RAD
      IF(ABS(BLON).GT.180.0)BLON=BLON-SIGN(360.,BLON)
      ALAT=ALT
      ALON=ALN
      RETURN
c100   write(6,*)'Flag in "Q2=ACOS(Q1)",Q1 = ',q1,';Q2 taken as 0'
100   q2=0.
      go to 300
c200   write(6,*)'     Flag in "Q2=ACOS(Q1)", Q1= ',q1,';Q2 taken as PI'
200   q2=pi
      go to 300
      END
c-------------------------------------------------------------
c----------------------------------------------------------
      subroutine trapez(npoint,f,ds,sint)
      dimension f(*)
      sint=0.
      nmax=npoint-1
      do 1 ns=2,nmax
  1   sint=sint+2.*f(ns)
      sint=sint+f(1)+f(npoint)
      sint=sint*ds*0.5
      return
      end
c----------------------------------------------------------
      SUBROUTINE  LISSE(N,X,Y,DY,S,A,B,C,D,R)
C-----------------------------------------------------------------------
C     CETTE SUBROUTINE PERMET LE CALCUL DE LA FONCTION SPLINE D'AJUSTEME
C        D'ORDRE 2 SUR LES POINTS X(I),Y(I),I=1,N .
C        LES X(I) SONT DONNES DANS L'ORDRE CROISSANT  . F EST CHOISIE TELLE
C        QUE SIGMA(((F(X(I))-Y(I))/DY(I))**2) )S .
C     VALEURS DE SORTIES : A,B,C,D (VECTEURS DE DIMENSION N+1)
C        F(T)=A(I)+B(I)*H+C(I)*H**2+D(I)*H**3  AVEC H=T-X(I-1) , X(I-1))
C        F(T)=A(1)+B(1)*(T-X(1)) SI T)X(1)
C        F(T)=A(N+1)+B(N+1)*(T-X(N)) SI T>X(N)
C     R : VECTEUR DE TRAVAIL EN DOUBLE PRECISION DE DIMENSION 8*(N+1)
C
C-----------------------------------------------------------------------
      DIMENSION X(N),Y(N),DY(N),A(N),B(N),C(N),D(N),R(N) 
      DOUBLE PRECISION BI,CI,DI,CI1,DI1,DI2,R,DH,DG,KA,KB,
     1DYJ,DP,AI,DV
      IF (S.LT.1.E-6) GO TO 999
C
C     INITIALISATIONS
C
      SS=S+1.E-06*S
      KA=2.D0/3.D0
      KB=1.D0/3.D0
      P=0.
      NIT=0
      NA=N + 1
      NB=NA + NA
      NC=NB+NA
      ND=NC+NA
      NE=ND+NA
      NF=NE+NA
      NG=NF+NA
      NH=NG+NA
      DO 10 J=1,NH
   10 R(J)=0.D0
      DO 15 I=1,NA
      A(I)=0.
      B(I)=0.
      C(I)=0.
   15 D(I)=0.
C
C     CALCUL DE C=Q*DDQ ET DE Q*D=Y
C
      H=X(2) -X(1)
      DH=DBLE(H)
      F=(Y(2)-Y(1))/H
      DO 20 I=3,N
      J=I-1
      G=H
      H=X(I)-X(J)
      DG=DH
      DH=DBLE(H)
      E=F
      F=(Y(I)-Y(J))/H
      A(I)=F-E
      R(NC+I)=KA*(DG+DH)
      R(ND+I)=KB*DH
      DYJ=DBLE(DY(J))
      R(I)=-DYJ*(1.D0/DG+1.D0/DH)
      R(NA+I)=DBLE(DY(J-1))/DG
   20 R(NB+I)=DBLE(DY(I))/DH
      IA=NA +2
      IB=NB+2
      DO 30 I=3,N
      IA=IA +1
      IB=IB +1
      R(NE+I)=R(IA)*R(IA) +R(I)*R(I) +R(IB)*R(IB)
      R(NF+I)=R(I)*R(IA+1)+R(I+1)*R(IB)
   30 R(IB)=R(IB)*R(IA+2)
   35 IF(NIT.GT.200) GO TO 999
      NIT=NIT +1
      DP=DBLE(P)
C
C     DECOMPOSITION CHOLESKI RR*=C
C
      DO 40 I=3,N
      I1=I-1
      I2=I-2
      AI=DBLE(A(I))
      BI=R(NE+I)+R(NC+I)*DP
      CI=R(NF+I)+R(ND+I)*DP
      DI=R(NB+I)
      TOL=1.E-16*ABS(SNGL(BI))
      DI1=DBLE(D(I1))
      CI1=DBLE(C(I1))
      DI2=DBLE(D(I2))
      BI=BI-DI2*DI2-CI1*CI1
      IF(SNGL(BI).LT.TOL) GO TO 999
      BI=1.D0/DSQRT(BI)
      R(I)=BI
      C(I)=SNGL(BI*(CI-CI1*DI1))
      D(I)=SNGL(BI*DI)
   40 R(NG+I)=(AI-CI1*R(NG+I1)-DI2*R(NG+I2))*BI
C
C     RESOLUTION CU=Y
C
      R(NH)=0.D0
      II=NH-1
      IJ=N
      R(II)=R(II)*R(IJ)
      DO 50 I=4,N
      II=II-1
      IJ=IJ-1
      CI=DBLE(C(IJ))
      DI=DBLE(D(IJ))
   50 R(II)=(R(II)-CI*R(II+1)-DI*R(II+2))*R(IJ)
C
C     CALCUL DE V=DQU ET E=V*V
C
      RES=0.
      H=0.
      F=0.
      IG=NG+1
      DO 60 I=2,N
      IG=IG+1
      I1=I-1
      G=H
      H=X(I)-X(I1)
      E=F
      F=(SNGL(R(IG+1)-R(IG)))/H
      B(I)=(F-E)*DY(I1)*DY(I1)
   60 RES=RES+B(I)*(F-E)
      B(NA)=-F*DY(N)*DY(N)
      RES=RES-B(NA)*F
C
C     TEST RES>S
C
      IF(RES.LT.SS) GO TO 80
C
C     CALCUL DE G=W*W ET F=U*TU
C
      G=0.
      F=0.
      IA =NA+2
      IC=NC+2
      ID=ND+2
      IG=NG+2
      DO 70 I=3,N
      IA=IA+1
      IC=IC+1
      ID=ID+1
      IG=IG+1
      DV=R(ID-1)*R(IG-1)+R(IC)*R(IG)+R(ID)*R(IG+1)
      CI1=DBLE(C(I-1))
      DI2=DBLE(D(I-2))
      R(IA)=(DV-CI1*R(IA-1)-DI2*R(IA-2))*R(I)
      G=G+SNGL(R(IA)*R(IA))
   70 F=F+SNGL(R(IG)*DV)
C
C     NOUVELLE VALEUR DE P
C
      P=P+(RES-SQRT(S*RES))/(F-P*G)
      GO TO 35
C
C     CALCUL DE A,B,C,D
C
   80 DO 90 I=2,NA
      C(I)=P*SNGL(R(NG+I))
   90 A(I)=Y(I-1)-B(I)
      DO 100 I=2,N
      H=X(I)-X(I-1)
      D(I)=(C(I+1)-C(I))/(3.*H)
  100 B(I)=(A(I+1)-A(I))/H-H*(C(I)+H*D(I))
      B(1)=B(2)
      A(1)=A(2)
      B(NA)=B(N)+(2.*C(N)+3.*D(N)*H)*H
      RETURN
  999 DO 2000 I=1,NA
      A(I)=0.
      B(I)=0.
      C(I)=0.
 2000 D(I)=0.
      RETURN
      END
c---------------------------------------------------------------
        FUNCTION SMOO(N,X,C0,C1,C2,C3,XX)
      DIMENSION X(1),C0(1),C1(1),C2(1),C3(1)
      SS=XX-X(N)
      S=XX-X(1)
      IF(S.LE.0.) GOTO 5
      IF(SS.GE.0.) GO TO 6
      DO 1 I=1,N
      T=XX-X(I)
      IF (T.LE.0.)GO TO 2
    1 CONTINUE
    2 D=XX-X(I-1)
      SMOO=C0(I)+D*C1(I)+D*D*C2(I)+D*D*D*C3(I)
      RETURN
    5 SMOO=C0(1)+C1(1)*S
      RETURN
    6 SMOO=C0(N+1)+C1(N+1)*SS
      RETURN
      END
c---------------------------------------
