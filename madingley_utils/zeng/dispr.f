       PROGRAM DISPR
c.....................................................................
c
c  Program : dispr.f
c  Purpose : for calculating the Rayleigh wave dispersion curve in a 
c            layered medium over a half space.
c  Input :
c            nlay - number of layer
c            h - layer thickness
c            vp, vs - P and S wave velocity
c            Qp, Qs - P and S wave attenuation factor
c            dens - medium density
c            T - period
c            cmn,cmx,dc - min, max guess and step of phase velocity
c            eps - relative tolerance of phase velocity
c            ite - max number of iterations
c  Output :
c            T - period
c            c - phase velocity
c            u - group velocity
c            cdvs,cdvp,cden - partial derivatives of c with respect
c                             to vs,vp and dens, respectively
c
c                                           -- Yuehua Zeng --
c                                           Apr., 1993, at UNR
c  All copy rights reserved.
c.....................................................................
c
	real*4 h(100),vp(100),Qp(100),vs(100),Qs(100),dens(100),
     +         cden(100), cdvs(100),cdvp(100)
	common/input/h,dens,vp,Qp,vs,Qs,T,cmin,cmax,dc,eps,ite,nlay
        open(10,file='dispr.in',form='formatted')
	read(10,*)nlay
        if(nlay .gt. 100) then
		print *, 'Increasce nlay'
		stop
	endif
	read(10,*)(h(i),dens(i),vp(i),vs(i),i=1,nlay)
	h(nlay)=1.0e+5
	read(10,*)nc
	c=0.0
	do 20 n=1,nc
	  read(10,*)T,cmin,cmax,dc,eps,ite
c
	  call crootr(c,u,cden,cdvp,cdvs,ier)
	  if(ier.eq.1)goto 20
          write(6,"(/12x,'T',11x,'C',11x,'U')")
	  write(6,'(4x,3f12.5)')  T,c,u
	  write(3,'(4x,3f12.5)')  T,c,u
          write(6,"(/9x,'DC/DVS',6x,'DC/DVP',6x,'DC/DDEN')")
          do 10 i=1,nlay
	    write(6,'(i2,3x,3f12.8)')i,cdvs(i),cdvp(i),cden(i)
 10	  continue
 20	continue
	close(10)
	stop
	end
c
c.........................................................
c
	subroutine crootr(c,u,cden,cdvp,cdvs,ier)
	real*4 h(100),vp(100),Qp(100),vs(100),Qs(100),mu(100),lam(100),
     +	       dens(100),cden(100),cdvs(100),cdvp(100),w,f,c,u,c1,c2
	complex*16 ka(100),kb(100),nu(100),ga(100),kbeta(100),
     +         exu(0:100,2),
     +	       I11(100,2,2),I12(100,2,2),I21(100,2,2),I22(100,2,2),
     +	       Td(100,2,2),Tu(100,2,2),Rd(100,2,2),Ru(0:100,2,2),
     +	       HTd(100,2,2),HRd(100,2,2),exd(100,2),aux(2,2),
     +	       temp,tmpr,ci,one,zero,Ed(100,2),Eu(100,2),Es(2),dzd(4),
     +	       ex(4,4),Er(4,4),a1,a2,a3,a4,a5,a6,temp1,temp2,temp3,
     +	       temp4,temp5,temp6
	complex ar(2)
	real*8 ak,ak2
	common/input/h,dens,vp,Qp,vs,Qs,T,cmin,cmax,dc,eps,ite,nlay
	common/para/ci,one,zero,ga,nu,kbeta,kb,mu,w,nwater,nly
	common/matrx1/I11,I12,I21,I22
	common/matrx2/exd,exu,Td,Tu,Rd,Ru,HTd,HRd
c
	ci=cmplx(0.0,1.0)
	pi=3.1415926535897936
	pi2=pi+pi
	one=cmplx(1.0,0.0)
	zero=cmplx(0.0,0.0)
	nly=nlay
	nwater=0
	if(vs(1).le.0.01)nwater=1
c
	f=1.0/T
	w=f*pi2
	kc=1
	c1=cmin
	if(c.ne.0.0.and.t.gt.10.0.and.ier.ne.1)c1=c+dc*0.2
 	c2=c1+dc
 	ak=w/c1
	do 5 j=1,nlay	
	   mu(j)=vs(j)*vs(j)*dens(j)
	   lam(j)=vp(j)*vp(j)*dens(j)-2*mu(j)
	   ka(j)=w/vp(j)
	   if(nwater.eq.0.or.j.ne.1)then
	     kb(j)=w/vs(j)
	   else
	     kb(1)=zero
	   endif
  5	continue
c
 10	do 40 k=kc,2
	  ak2=ak*ak
	  do 20 j=1,nlay
	     nu(j)=sqrt(cmplx(ak2,0.0)-ka(j)*ka(j))
	     if(dreal(nu(j)).lt.0.)nu(j)=-nu(j)
	     ga(j)=sqrt(cmplx(ak2,0.0)-kb(j)*kb(j))
	     if(dreal(ga(j)).lt.0.)ga(j)=-ga(j)
	     kbeta(j)=2*ak2-kb(j)*kb(j)
	     exd(j,1)=exp(-nu(j)*h(j))
	     exd(j,2)=exp(-ga(j)*h(j))
	     exu(j-1,1)=exd(j,1)
	     exu(j-1,2)=exd(j,2)
 20	  continue
	  call gcoefra(ak,dens(1),dens(2))
	  if(nwater.eq.0)then
	    do 30 i=1,2
	      do 30 j=1,2
	        aux(i,j)=0.
	        do 30 l=1,2
	           aux(i,j)=aux(i,j)-Ru(0,i,l)*HRd(1,l,j)
 30	    continue
 	    ar(k)=(one+aux(1,1))*(one+aux(2,2))-aux(1,2)*aux(2,1)
	  else
	    ar(k)=one-Ru(0,1,1)*HRd(1,1,1)
	  endif
	  ak=w/c2
 40	continue
	kc=2
	if(real(ar(1))*real(ar(2)).lt.0.0)then
	   nroot=1
	   tmp1=real(ar(1))
	   tmp2=real(ar(2))
	elseif(aimag(ar(1))*aimag(ar(2)).lt.0.0)then
	   nroot=2
	   tmp1=aimag(ar(1))
	   tmp2=aimag(ar(2))
	else
	   c1=c2
	   ar(1)=ar(2)
	   c2=c2+dc
	   if(c2.gt.cmax)then
	      print*,' root not found!'
	      ier=1
	      return
	   endif
	   ak=w/c2
	   goto 10
	endif
c
	csave=c2
	tmpr=ar(2)
	do 70 k=1,ite
	  c=0.5*(c2+c1)
	  ak=w/c
	  ak2=ak*ak
	  do 50 j=1,nlay
	     nu(j)=sqrt(cmplx(ak2,0.0)-ka(j)*ka(j))
	     ga(j)=sqrt(cmplx(ak2,0.0)-kb(j)*kb(j))
	     if(dreal(nu(j)).lt.0.)nu(j)=-nu(j)
	     if(dreal(ga(j)).lt.0.)ga(j)=-ga(j)
	     kbeta(j)=2*ak2-kb(j)*kb(j)
	     exd(j,1)=exp(-nu(j)*h(j))
	     exd(j,2)=exp(-ga(j)*h(j))
	     exu(j-1,1)=exd(j,1)
	     exu(j-1,2)=exd(j,2)
 50	  continue
	  call gcoefra(ak,dens(1),dens(2))
	  if(nwater.eq.0)then
	    do 60 i=1,2
	      do 60 j=1,2
	        aux(i,j)=0.
	        do 60 l=1,2
	           aux(i,j)=aux(i,j)-Ru(0,i,l)*HRd(1,l,j)
 60	    continue
 	    temp=(one+aux(1,1))*(one+aux(2,2))-aux(1,2)*aux(2,1)
	  else
	    temp=one-Ru(0,1,1)*HRd(1,1,1)
	  endif
	  tmp2=dreal(temp)
	  if(nroot.eq.2)tmp2=dimag(temp)
	  if(abs(temp).lt.eps)then
	    goto 80
	  elseif(tmp1*tmp2.le.0.0)then
	    c2=c
	    ar(2)=temp
	  else
	    c1=c
	    ar(1)=temp
	    tmp1=tmp2
	  endif
 70	continue
	c1=csave
	ar(1)=tmpr
	c2=c1+dc
	ak=w/c2
	goto 10
c
 80	call coefra(ak,dens(1))
	if(nwater.eq.0)then
	  Ed(1,1)=-aux(1,2)
	  Ed(1,2)=(one+aux(1,1))
	else
	  Ed(1,1)=one
	  Ed(1,2)=zero
	endif
	do 110 j=1,2
	  Eu(1,j)=zero
	  do 110 l=1,2
	     Eu(1,j)=Eu(1,j)+HRd(1,j,l)*Ed(1,l)
110	continue
	temp=I11(1,2,1)*Ed(1,1)+I11(1,2,2)*Ed(1,2)+I12(1,2,1)
     +	     *Eu(1,1)*exd(1,1)+I12(1,2,2)*Eu(1,2)*exd(1,2)
	do 115 i=1,2
	  Ed(1,i)=Ed(1,i)/temp
	  Eu(1,i)=Eu(1,i)/temp
	  Es(i)=Ed(1,i)
115	continue
c
	do 140 i=2,nlay
	  do 120 j=1,2
	    Ed(i,j)=zero
	    do 120 l=1,2
	       Ed(i,j)=Ed(i,j)+HTd(i-1,j,l)*Es(l)
120	  continue
	  Es(1)=Ed(i,1)
	  Es(2)=Ed(i,2)
	  do 130 j=1,2
	    Eu(i,j)=zero
	    do 130 l=1,2
	       Eu(i,j)=Eu(i,j)+HRd(i,j,l)*Ed(i,l)
130	  continue
140	continue
c
c  compute the displacement and stress field
	en1=0.0
	en2=0.0
	en3=0.0
	en4=0.0
	ak2=ak*ak
	w2=w*w
	do 150 i=1,nlay
	  a1=0.5*(one-exd(i,1)*exd(i,1))/nu(i)
	  a2=0.5*(one-exd(i,2)*exd(i,2))/ga(i)
	  a3=h(i)*exd(i,1)
	  a4=h(i)*exd(i,2)
	  a5=(one-exd(i,1)*exd(i,2))/(nu(i)+ga(i))
	  a6=(exd(i,1)-exd(i,2))/(ga(i)-nu(i))
c
	  dzd(1)=-nu(i)
	  dzd(2)=-ga(i)
	  dzd(3)=nu(i)
	  dzd(4)=ga(i)
c
	  ex(1,1)=a1
	  ex(1,2)=a5
	  ex(1,3)=a3
	  ex(1,4)=a6
	  ex(2,1)=a5
	  ex(2,2)=a2
	  ex(2,3)=a6
	  ex(2,4)=a4
	  ex(3,1)=a3
	  ex(3,2)=a6
	  ex(3,3)=a1
	  ex(3,4)=a5
	  ex(4,1)=a6
	  ex(4,2)=a4
	  ex(4,3)=a5
	  ex(4,4)=a2
c
	  Er(1,1)=I11(i,1,1)*Ed(i,1)
	  Er(1,2)=I11(i,1,2)*Ed(i,2)
	  Er(1,3)=I12(i,1,1)*Eu(i,1)
	  Er(1,4)=I12(i,1,2)*Eu(i,2)
	  Er(2,1)=I11(i,2,1)*Ed(i,1)
	  Er(2,2)=I11(i,2,2)*Ed(i,2)
	  Er(2,3)=I12(i,2,1)*Eu(i,1)
	  Er(2,4)=I12(i,2,2)*Eu(i,2)
	  temp1=0.0
	  temp2=0.0
	  temp3=0.0
	  temp4=0.0
	  temp5=0.0
	  temp6=0.0
	  do 146 k=1,4
	    do 145 l=1,4
	       temp1=temp1+(Er(1,k)*ex(k,l)*Er(1,l))
	       temp2=temp2+(Er(2,k)*ex(k,l)*Er(2,l))
	       temp3=temp3+(dzd(k)*Er(1,k)*ex(k,l)*Er(2,l))
	       temp4=temp4+(Er(1,k)*ex(k,l)*dzd(l)*Er(2,l))
	       temp5=temp5+(dzd(k)*Er(1,k)*ex(k,l)*dzd(l)*Er(1,l))
	       temp6=temp6+(dzd(k)*Er(2,k)*ex(k,l)*dzd(l)*Er(2,l))
145	    continue
146	  continue
	  cdvs(i)=2.0*(ak2*temp1+temp6-ak*temp3)+ak2*temp2+temp5
	  cdvp(i)=ak2*temp1+2.0*ak*temp4+temp6
	  cden(i)=-(temp1+temp2)*w2
	  en1=en1+0.5*dens(i)*(temp1+temp2)
	  en2=en2+0.5*((lam(i)+2.0*mu(i))*temp1+mu(i)*temp2)
	  en3=en3+lam(i)*temp4-mu(i)*temp3
	  en4=en4+0.5*((lam(i)+2.0*mu(i))*temp6+mu(i)*temp5)
150	continue
	u=(en2+0.5*en3/ak)/(c*en1)
	c2=0.25/(ak2*u*en1)
	do 155 i=1,nlay
	  cden(i)=(cden(i)+(cdvs(i)*mu(i)+cdvp(i)*lam(i))/dens(i))*c2
	  cdvp(i)=cdvp(i)*dens(i)*2.0*vp(i)*c2
	  cdvs(i)=(cdvs(i)*dens(i)*c2-cdvp(i)/vp(i))*2.0*vs(i)
155	continue
	ier=0
c
	return
	end
c
c.........................................................
c
	subroutine coefra(ak,dens1)
	real*8 ak
	real*4 mu(100)
	complex*16 kb(100),nu(100),ga(100),kbeta(100),ci,one,zero,a1,
     +         a2,a3,
     +	       I11(100,2,2),I12(100,2,2),I21(100,2,2),I22(100,2,2),a4,
     +         a5,a6
	common/para/ci,one,zero,ga,nu,kbeta,kb,mu,w,nwater,nlay
	common/matrx1/I11,I12,I21,I22
c	
c  generate displacement-stress matrices
	do 10 n=1,nlay
	  if(nwater.eq.0.or.n.ne.1)then
	    a1=ak
	    a2=ga(n)
	    a3=nu(n)
	    a4=2*nu(n)*ak*mu(n)
	    a5=kbeta(n)*mu(n)
	    a6=2*ga(n)*ak*mu(n)
c
	    I11(n,1,1)=a1
	    I11(n,1,2)=a2
	    I11(n,2,1)=a3
	    I11(n,2,2)=a1
	    I12(n,1,1)=a1
	    I12(n,1,2)=a2
	    I12(n,2,1)=-a3
	    I12(n,2,2)=-a1
	    I21(n,1,1)=-a4
	    I21(n,1,2)=-a5
	    I21(n,2,1)=-a5
	    I21(n,2,2)=-a6
	    I22(n,1,1)=a4
	    I22(n,1,2)=a5
	    I22(n,2,1)=-a5
	    I22(n,2,2)=-a6
	  else
	    a1=ak
	    a2=nu(1)
	    a3=dens1*w*w
	    I11(1,1,1)=a1
	    I11(1,1,2)=zero
	    I11(1,2,1)=a2
	    I11(1,2,2)=zero
	    I12(1,1,1)=a1
	    I12(1,1,2)=zero
	    I12(1,2,1)=-a2
	    I12(1,2,2)=zero
	    I21(1,1,1)=zero
	    I21(1,1,2)=zero
	    I21(1,2,1)=a3
	    I21(1,2,2)=zero
	    I22(1,1,1)=zero
	    I22(1,1,2)=zero
	    I22(1,2,1)=a3
	    I22(1,2,2)=zero
	  endif
 10	continue
	return
	end
c
c.........................................................
c
	subroutine gcoefra(ak,dens1,dens2)
	real*8 ak,ak2
	real*4 mu(100)
	complex*16 kb(100),nu(100),ga(100),kbeta(100),aux(2,2),
     +         ainv(2,2),
     +	       Td(100,2,2),Tu(100,2,2),Rd(100,2,2),Ru(0:100,2,2),one,
     +         zero,
     +	       HTd(100,2,2),HRd(100,2,2),exd(100,2),exu(0:100,2),ci,det,
     +	       a1,a2,a3,a4,a5,a6,a,b,c,d,E1,E2,F1,F2,G1,G2,H1,H2
	common/para/ci,one,zero,ga,nu,kbeta,kb,mu,w,nwater,nlay
	common/matrx2/exd,exu,Td,Tu,Rd,Ru,HTd,HRd
c
c  modified reflection and transmission coefficients
	ak2=ak*ak
c
c  at the free surface
	if(nwater.eq.0)then
	  a=kbeta(1)*kbeta(1)
	  b=4.0*ak2*nu(1)*ga(1)
	  d=a-b
          Ru(0,1,1)=-(a+b)*(exu(0,1)/d)
          Ru(0,1,2)=-4.0*ak*ga(1)*kbeta(1)*(exu(0,2)/d)
          Ru(0,2,1)= 4.0*ak*nu(1)*kbeta(1)*(exu(0,1)/d)
          Ru(0,2,2)= (a+b)*(exu(0,2)/d)
	else
          Ru(0,1,1)=-exu(0,1)
          Ru(0,1,2)=zero
          Ru(0,2,1)=zero
          Ru(0,2,2)=zero
	endif
c
	do 20 n=1,nlay-1
	  if(nwater.eq.0.or.n.ne.1)then
	    a=mu(n+1)*kbeta(n+1)-mu(n)*kbeta(n)
	    b=mu(n+1)*kbeta(n+1)-2.*mu(n)*ak2
	    c=mu(n)*kbeta(n)-2.*mu(n+1)*ak2
	    d=2.*(mu(n+1)-mu(n))
	    E1=b*nu(n)+c*nu(n+1)
	    E2=b*nu(n)-c*nu(n+1)
	    F1=b*ga(n)+c*ga(n+1)
	    F2=b*ga(n)-c*ga(n+1)
	    G1=a+d*nu(n)*ga(n+1)
	    G2=a-d*nu(n)*ga(n+1)
	    H1=a+d*nu(n+1)*ga(n)
	    H2=a-d*nu(n+1)*ga(n)
	    det=-E1*F1+G2*H2*ak2
	    a1=mu(n)*kb(n)*kb(n)
	    a2=mu(n+1)*kb(n+1)*kb(n+1)
	    a3=-a*c-b*d*nu(n)*ga(n)
	    a4=a*b+c*d*nu(n+1)*ga(n+1)
c
	    Td(n,1,1)=2*a1*nu(n)*F1/det*exd(n,1)
	    Td(n,1,2)=2*ak*a1*ga(n)*G2/det*exd(n,2)
	    Td(n,2,1)=2*ak*a1*nu(n)*H2/det*exd(n,1)
	    Td(n,2,2)=2*a1*ga(n)*E1/det*exd(n,2)
	    Ru(n,1,1)=(E2*F1-G2*H1*ak2)/det*exu(n,1)
 	    Ru(n,1,2)=-2*ak*ga(n+1)*a3/det*exu(n,2)
 	    Ru(n,2,1)=2*ak*nu(n+1)*a3/det*exu(n,1)
 	    Ru(n,2,2)=(-F2*E1+H2*G1*ak2)/det*exu(n,2)
	    Rd(n,1,1)=(-E2*F1-G1*H2*ak2)/det*exd(n,1)
	    Rd(n,1,2)=-2*ak*ga(n)*a4/det*exd(n,2)
	    Rd(n,2,1)=2*ak*nu(n)*a4/det*exd(n,1)
	    Rd(n,2,2)=(F2*E1+H1*G2*ak2)/det*exd(n,2)
	    Tu(n,1,1)=2*a2*nu(n+1)*F1/det*exu(n,1)
	    Tu(n,1,2)=-2*ak*a2*ga(n+1)*H2/det*exu(n,2)
	    Tu(n,2,1)=-2*ak*a2*nu(n+1)*G2/det*exu(n,1)
	    Tu(n,2,2)=2*a2*ga(n+1)*E1/det*exu(n,2)
	  else
	    a=dens1*w*w
	    b=dens2*w*w
	    c=mu(2)*kbeta(2)
	    d=2.0*ak*mu(2)
	    det=-nu(2)*a*b+nu(1)*(2.0*ak*mu(2)*ga(2)*nu(2)*d-c*c)
	    a1=nu(1)*c
	    a2=-nu(1)*nu(2)*d
	    a3=nu(2)*b
	    a4=2.0*a/det
	    a5=2.0*c/det
	    a6=2.0*ga(2)*d/det
	    Td(n,1,1)=a1*a4*exd(n,1)
	    Td(n,1,2)=zero
	    Td(n,2,1)=a2*a4*exd(n,1)
	    Td(n,2,2)=zero
	    Ru(n,1,1)=(a1*a5+one)*exu(n,1)
 	    Ru(n,1,2)=a1*a6*exu(n,2)
 	    Ru(n,2,1)=a2*a5*exu(n,1)
 	    Ru(n,2,2)=(a2*a6+one)*exu(n,2)
	    Rd(n,1,1)=(a3*a4+one)*exd(n,1)
	    Rd(n,1,2)=zero
	    Rd(n,2,1)=zero
	    Rd(n,2,2)=zero
	    Tu(n,1,1)=a3*a5*exu(n,1)
	    Tu(n,1,2)=a3*a6*exu(n,2)
	    Tu(n,2,1)=zero
	    Tu(n,2,2)=zero
	  endif
 20	continue
c
c  generalized reflection and transmission coefficients for
c  layers below the source.
	HRd(nlay,1,1)=zero
	HRd(nlay,1,2)=zero
	HRd(nlay,2,1)=zero
	HRd(nlay,2,2)=zero
c
	do 60 n=nlay-1,1,-1
	do 30 i=1,2
	do 30 j=1,2
	  aux(i,j)=zero
	  do 30 k=1,2
	    aux(i,j)=aux(i,j)-Ru(n,i,k)*HRd(n+1,k,j)
 30	continue
	call cinv2(aux,ainv)
	do 40 i=1,2
	do 40 j=1,2
	  HTd(n,i,j)=zero
	  do 40 k=1,2
	    HTd(n,i,j)=HTd(n,i,j)+ainv(i,k)*Td(n,k,j)
 40	continue
	do 50 i=1,2
	do 50 j=1,2
	  HRd(n,i,j)=Rd(n,i,j) 
	  do 50 k=1,2
	  do 50 l=1,2
	    HRd(n,i,j)=HRd(n,i,j)+Tu(n,i,k)*HRd(n+1,k,l)*HTd(n,l,j)
 50	continue
 60	continue
c	
	return
	end
c
c..........................................................
c
        subroutine cinv2(a,ainv)
        complex*16 a(2,2),ainv(2,2),det,one
	one=cmplx(1.0,0.0)
	a(1,1)=one+a(1,1)
	a(2,2)=one+a(2,2)
        det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
        ainv(1,1)=a(2,2)/det
        ainv(1,2)=-a(1,2)/det
        ainv(2,1)=-a(2,1)/det
        ainv(2,2)=a(1,1)/det
        return
        end
