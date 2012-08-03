c     10/31/86             @(#)tvinvma22b.f	1.4 
*-----------------------------------------------------------------------
      Subroutine INTV(P0,D0,CP0,CD0,G,Dinit,Pinit,n,m,nd,mp,P,CP,Resol,
     *                                                        Residu)
      Parameter (mmp=200, nnd=200)

*     Tableaux en INPUT:
      Dimension P0(mp), D0(nd), CP0(mp,mp), CD0(nd,nd), G(nd,mp),
     *          Pinit(mp), Dinit(nd)

*     Tableaux en OUTPUT:
      Dimension P(mp), CP(mp,mp), Resol(mp,mp), Residu(nd)

*     Tableaux de travail:
      Dimension CP0G(mmp,nnd), B(nnd,nnd), BB(mmp,nnd)
     *          , B1(nnd,nnd), B2(nnd,nnd)
 
*    Ce sous-programme calcul le modele estime par inversion, algorithme
*    Tarantola & Valette (1982) moindres carres, discret, explicite,
*                               sous-determine (formule 23).
*
*  INPUT:
*        P0, CP0 : modele a priori et covariance a priori du modele
*        D0, CD0 : donnees (mesures) et covariance des donnees
*        G       : matrice des derivees partielles (adaptee a la discretisation)
*        Pinit   : modele provenant de l'iteration precedente
*        Dinit  : donnees synthetiques pour Pinit
*        m       : nombre de parametres reellement inverses
*        n       : nombre de donnees reellement inversees
*
* OUTPUT:
*        P       : modele estime (final)
*        CP      : covariance a posteriori du modele estime (basee sur
*                  la formule lineaire, cf. remarque 2 du par. 2f de T&V)
*        Resol   : Resolution ( = Inverse.Direct)
*        Residu    : residus a posteriori (misfits) estimes lineairement
*
*  LOCAL:
*        CP0G    : CP0.G~
*        B       : (CD0+G.CP0.G~)-1
*        BB      : CP0.G~.(CD0+G.CP0.G~)-1 = CP0.G~.B
 
 
      if (mp.ne.mmp.or.m.gt.mp) stop 'INVTV:erreur de dimension m ou mp'
      write(*,*) nd
      if (nd.ne.nnd.or.n.gt.nd) stop 'INVTV:erreur de dimension n ou nd'
 

      do 1 i=1,m
        do 1 j=1,n
        CP0G(i,j)=0.
          do 1 k=1,m
          CP0G(i,j)=CP0G(i,j)+CP0(i,k)*G(j,k)
    1 continue
 
      do 2 i=1,n
        do 2 j=1,n
        B(i,j)=0.
          do 3 k=1,m
          B(i,j)=B(i,j) + G(i,k)*CP0G(k,j)
    3     continue
        B(i,j)= B(i,j) + CD0(i,j)
        B1(i,j)= B(i,j)
    2 continue
 
      e=0.
      call MA22B(B,nd,n,W,e)
      if(e.ne.0.) stop 'erreur dans MA22B (inversion de matrice)'
c     verification de la qualite de ce qui sort de MA22B
      do 22 i=1,n
      do 22 j=1,n
      B2(i,j)=0.
      do 22 k=1,n
      B2(i,j)=B2(i,j)+B1(i,k)*B(k,j)
   22 continue
      amax=0.
      do 23 i=1,n
      if(abs(B2(i,i)-1.).gt.amax) amax=abs(B2(i,i)-1.)
      B2(i,i)=0.
      do 23 j=1,n
      if(abs(B2(i,j)).gt.amax) amax=abs(B2(i,j))
   23 continue
      write(*,*) ' Ecart maxi % identite : ', amax
      if(amax.gt.1.e-2) write(0,*) ' B*B-1 differe de l''Identite'

      do 4 i=1,m
        do 4 j=1,n
        BB(i,j)=0.
          do 4 k=1,n
          BB(i,j)=BB(i,j) + CP0G(i,k)*B(k,j)
    4 continue
 
      do 5 i=1,m
        do 5 j=1,m
        Resol(i,j)=0.
          do 5 k=1,n
          Resol(i,j)=Resol(i,j) + BB(i,k)*G(k,j)
    5 continue
 
        do 77 j=1,n
c  le tableau Residu est employe ici comme tableau de travail
      Residu(j)=0.
        do 7 k=1,m
    7   Residu(j)=Residu(j) + G(j,k)*(Pinit(k)-P0(k))
      Residu(j)=Residu(j)+D0(j)-Dinit(j)
   77 continue
      do 6 i=1,m
      P(i)=0.
      do 66 j=1,n
   66 P(i)=P(i)+ BB(i,j)*Residu(j)
      P(i)=P(i)+P0(i)
    6 continue
 
c     CP0G utilise ici comme tableau de travail
C     do 8 i=1,n
C       do 8 j=1,m
C       CP0G(i,j)=0.
C         do 8 k=1,m
C   8     CP0G(i,j)=CP0G(i,j) + G(i,k)*CP0(k,j)

      do 88 i=1,m
      do 88 j=1,m
      CP(i,j)=0.
        do 89 k=1,n
   89   CP(i,j)=CP(i,j)+ BB(i,k)*CP0G(j,k)
   88   CP(i,j)=CP0(i,j)-CP(i,j)
 
      do 10 i=1,n
      Residu(i)=0.
        do 11 k=1,m
        Residu(i)=Residu(i) + G(i,k)*(P(k)-Pinit(k))
   11   continue
      Residu(i)=D0(i)-Dinit(i)-Residu(i)
   10 continue
 
      return
      end
*-------------------------------------------------------------------------
