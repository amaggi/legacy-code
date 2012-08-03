      PROGRAM sfroid
      INTEGER NE,M,NB,NCI,NCJ,NCK,NSI,NSJ,NYJ,NYK
      COMMON /sfrcom/ x,h,mm,n,c2,anorm
      PARAMETER (NE=3,M=41,NB=1,NCI=NE,NCJ=NE-NB+1,NCK=M+1,NSI=NE,NSJ=2*
     *NE+1,NYJ=NE,NYK=M)
CU    USES plgndr,solvde
      INTEGER i,itmax,k,mm,n,indexv(NE)
      REAL anorm,c2,conv,deriv,fac1,fac2,h,q1,slowc,c(NCI,NCJ,NCK),
     *s(NSI,NSJ),scalv(NE),x(M),y(NE,M),plgndr
      itmax=100
      conv=5.e-6
      slowc=1.
      h=1./(M-1)
      c2=0.
      write(*,*)'ENTER M,N'
      read(*,*)mm,n
      if(mod(n+mm,2).eq.1)then
        indexv(1)=1
        indexv(2)=2
        indexv(3)=3
      else
        indexv(1)=2
        indexv(2)=1
        indexv(3)=3
      endif
      anorm=1.
      if(mm.NE.0)then
        q1=n
        do 11 i=1,mm
          anorm=-.5*anorm*(n+i)*(q1/i)
          q1=q1-1.
11      continue
      endif
      do 12 k=1,M-1
        x(k)=(k-1)*h
        fac1=1.-x(k)**2
        fac2=fac1**(-mm/2.)
        y(1,k)=plgndr(n,mm,x(k))*fac2
        deriv=-((n-mm+1)*plgndr(n+1,mm,x(k))-(n+1)*x(k)*plgndr(n,mm,
     *x(k)))/fac1
        y(2,k)=mm*x(k)*y(1,k)/fac1+deriv*fac2
        y(3,k)=n*(n+1)-mm*(mm+1)
12    continue
      x(M)=1.
      y(1,M)=anorm
      y(3,M)=n*(n+1)-mm*(mm+1)
      y(2,M)=(y(3,M)-c2)*y(1,M)/(2.*(mm+1.))
      scalv(1)=abs(anorm)
      scalv(2)=max(abs(anorm),y(2,M))
      scalv(3)=max(1.,y(3,M))
1     continue
      write (*,*) 'ENTER C**2 OR 999 TO END'
      read (*,*) c2
      if (c2.eq.999.) stop
      call solvde(itmax,conv,slowc,scalv,indexv,NE,NB,M,y,NYJ,NYK,c,NCI,
     *NCJ,NCK,s,NSI,NSJ)
      write (*,*) ' M = ',mm,'  N = ',n,'  C**2 = ',c2,'  LAMBDA = ',
     *y(3,1)+mm*(mm+1)
      goto 1
      END
