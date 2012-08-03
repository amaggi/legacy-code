c     Pgm d'estimation de l'amplitude relative de chaque mode 
c     d'un sismograme synthetique filtres a differentes periodes

C Dimensions des tableaux:
C xr,xi,xrs .......................... 2**ip
C
c     parameter(ndim=4096, ndim2=ndim/2+1)
      parameter(ndim=8192, ndim2=ndim/2+1)
      dimension xr(ndim),xi(ndim),xrs(ndim)
      dimension per(10),som1(10,10),som2(10),ratio(10,10)
      character fi*80,ficpi*80
      real*8 pi, pi2

      pi=3.1415926536
      pi2=2*pi
      r0=6371.
      io=10
      ipi=11
      iout=12


9     write(*,*)'Nom du fichier contenant la TF du sismo:'
      read(*,'(a)')fi
      open(io,status='old',form='unformatted',err=9,file=fi)
      read(io)npp,dnu,t1,dt,dlt,km,ip,dist
      write(*,*)npp,dnu,t1,dt,dlt,km,ip,dist
      
      write(*,*)'Nom du fichier pick'
      read(*,'(a)') ficpi                 
      open(ipi,status='old',file=ficpi)
      read(ipi,*)nmode,nper,(per(i),i=1,nper)
      close(ipi)

      np=2*npp-1
      do 10 i=1,np
10    xrs(i)=0.
c  Boucle sur le rang k des modes
      do 20 k=1,km
      read(io) (xr(i),xi(i),i=1,npp)
c     write(*,*)(xr(i),xi(i),i=1,npp)
c  Boucle sur les periodes de filtrage des intercorrelogrames           
      do 30 iper=1,nper
      som1(iper,k)=0.
      if(k.eq.1) som2(iper)=0.
      p4=per(iper)-5.
      p3=per(iper)
      p2=per(iper)+10.
      p1=2*p2
      write(*,*)'PERIODE',p1,p2,p3,p4
      call fildes(xr,xi,npp,dnu,p1,p2,p3,p4)
c     write(*,*)(xr(i),xi(i),i=1,npp)
c  Passage en temporel
      call nlogn(ip,xr,xi,+1.)
      write(*,*)'IP',ip
c     write(*,*)(xr(i),xi(i),i=1,np)
      do 40 i=1,np
40    som1(iper,k)=som1(iper,k)+(2*xr(i))**2
      som2(iper)=som2(iper)+som1(iper,k)
30    continue
20    continue
      read(io) amoment
      write(*,*) ' amoment =',amoment
      close(io)

c     calcul de l'energie de chaque mode filtre/sismo somme filtre
      open(iout,file='energiemode')
      ratio(iper,k)=0.
      do 60 k=1,km
c     write(*,*)('som2 ',som2(iper),' som1 ',som1(iper,k),iper=1,nper)
      do 70 iper=1,nper
      ratio(iper,k)=som1(iper,k)/som2(iper)
70    continue
      write(iout,1000),k,(per(iper),ratio(iper,k)*100,iper=1,nper)
60    continue
      close(iout)
1000  format(i2(10f8.2))
      end

c---------------------------------------------------------------
      SUBROUTINE NLOGN(N,xr,xi,SIGN)
C ALGORITHME COOLEY TUKEY.
C      SIGN=-1.0  EXP (-2*I*PI*...)
C      SIGN =1.0 (1/Q) * EXP(2*I*PI*...)
C      NMAX=PLUS GRANDE VALEUR DE N AVEC
C      DIMENSION X(2**N) ,M(NMAX)
      DIMENSION XR(2),XI(2),M(20)
      LX=2**N
      DO 1 I=1,N
    1 M(I)=2**(N-I)
      DO 4 L=1,N
      NBLOC=2**(L-1)
      LBLOC=LX/NBLOC
      LBHAF=LBLOC/2
      K=0
      DO 4 IBLOC=1,NBLOC
      FK=K
      FLX=LX
      V=SIGN*6.2831853*FK/FLX
      WKR=COS(V)
      WKI=SIN(V)
      ISTAT=LBLOC*(IBLOC-1)
      DO 2 I=1,LBHAF
      J=ISTAT+I
      JH=J+LBHAF
      QR=XR(JH)*WKR-XI(JH)*WKI
      QI=XI(JH)*WKR+XR(JH)*WKI
      XR(JH)=XR(J)-QR
      XI(JH)=XI(J)-QI
      XR(J)=XR(J)+QR
      XI(J)=XI(J)+QI
    2 CONTINUE
      DO 3 I=2,N
      II=I
      IF (K-M(I)) 4,3,3
    3 K=K-M(I)
    4 K=K+M(II)
      K=0
      DO 8 J=1,LX
      IF (K-J) 5,6,6
    6 HOLDR=XR(J)
      HOLDI=XI(J)
      XR(J)=XR(K+1)
      XI(J)=XI(K+1)
      XR(K+1)=HOLDR
      XI(K+1)=HOLDI
    5 DO 7 I=1,N
      II=I
      IF (K-M(I)) 8,7,7
    7 K=K-M(I)
    8 K=K+M(II)
      IF (SIGN) 11,9,9
    9 DO 10 I=1,LX
      XR(I)=XR(I)/FLX
   10 XI(I)=XI(I)/FLX
   11 RETURN
      END
c---------------------------------------------------------------
      subroutine fildes(xr,xi,npp,dnu,p1,p2,p3,p4)
      dimension xr(1),xi(1)
c
c          Filtrage passe-bande d'un sismogramme sous forme TF
c
c                 npp  indice de la frequence de Nyquist du signal
c                 xr   partie reelle de la TF
c                 xi   partie imaginaire de la TF
c                 dnu  pas en frequence de la TF
c
c
c Filtrage passe bande
c
      per=1./dnu/(npp-1)
15    write(*,'("la periode de Nyquist doit etre superieure a",
     *      f10.4)')per
      if1=1./p1/dnu+1.5
      if2=1./p2/dnu+1.5
      if3=1./p3/dnu+1.5
      if4=1./p4/dnu+1.5
        if(if1.le.1) if1=2
        if(if2.le.1) if2=2
        if(if3.gt.npp) if3=npp
        if(if4.gt.npp) if4=npp
      f1=dnu*(if1-1)
      f2=dnu*(if2-1)
      f3=dnu*(if3-1)
      f4=dnu*(if4-1)
      write(*,*)'frequences reelles du filtre:',f1,f2,f3,f4
      write(*,*)'periodes reelles du filtre:', 1./f4,1./f3,1./f2,1./f1
      id2=if2-if1
      id3=if4-if3
      if(id2.lt.1) goto 88
      do 7 i=if1,if2
      u=float(i-if2)/float(id2)
      u=u*u
      u=(1.-u)*(1.-u)
      xr(i)=xr(i)*u
7     xi(i)=xi(i)*u
88    if(id3.lt.1) goto 19
      do 71 i=if3,if4
      u=float(i-if3)/float(id3)
      u=u*u
      u=(1.-u)*(1.-u)
      xr(i)=xr(i)*u
71    xi(i)=xi(i)*u
19    continue
      if(if1.le.1) goto 18
      do 17 i=1,if1-1
      xr(i)=0.
      xi(i)=0.
17    continue
18    continue
      if(if4.ge.npp) goto 172
      do 171 i=if4+1,npp
      xr(i)=0.
      xi(i)=0.
171   continue
172   continue
c     write(*,*)(xr(i),xi(i),i=1,npp)

      return
      end

