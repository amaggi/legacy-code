       PROGRAM DISPL
c.....................................................................
c
c  Program : displ.f
c  Purpose : for calculating the Love wave dispersion curve in a 
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
c            cdvs,cden - partial derivatives of c with respect to vs
c                        and dens, respectively 
c
c                                           -- Yuehua Zeng --
c                                           Apr., 1993, at UNR
c  All copy rights reserved.
c.....................................................................
c
	real*4 h(50),vp,Qp,vs(50),Qs(50),dens(50),cden(50),cdvs(50)
	common/input/h,dens,vs,Qs,nlay,T,cmin,cmax,dc,eps,ite
c
        open(10,file='displ.in',form='formatted')
	read(10,*)nlay
	read(10,*)(h(i),dens(i),vp,vs(i),i=1,nlay)
	if(vs(1).le.0.01)then
	  nlay=nlay-1
	  do 10 i=1,nlay
	    h(i)=h(i+1)
	    dens(i)=dens(i+1)
	    vs(i)=vs(i+1)
 10	  continue
	endif
	h(nlay)=1.0e+5
	read(10,*)nc
	c=0.0
c
	do 30 n=1,nc
	  read(10,*)T,cmin,cmax,dc,eps,ite
	  call crootl(c,u,cden,cdvs,ier)
	  if(ier.eq.1)goto 30
          write(6,"(/12x,'T',11x,'C',11x,'U')")
	  write(6,'(4x,3f12.5)')  T,c,u
          write(6,"(/9x,'DC/DVS',6x,'DC/DDEN')")
	  do 20 i=1,nlay
	    write(6,'(i2,3x,2f12.8)')i,cdvs(i),cden(i)
 20	  continue
 30	continue
c
	close(10)
	stop
	end
c
c....................................................................
c
	subroutine crootl(c,u,cden,cdvs,ier)
	real*4 h(50),vs(50),Qs(50),mu(50),dens(50),cden(50),cdvs(50),
     +         I1,I2,I3,dI1,dI2,dI3,f,w,c,c1,c2,u
	complex*16 kb(50),ga(50),Edsh(50),Eush(50),e1,e2,ci,one,zero,
     +	       Ish11(50),Ish12(50),Ish21(50),Ish22(50),exu(0:50,2),
     +	       Tdsh(50),Tush(50),Rdsh(50),Rush(0:50),exd(50,2),a1,a2,
     +	       HTdsh(50),HRdsh(50),temp,tmpl
	complex al(2)
	real*8 ak,ak2
	real tmp1,tmp2
	common/input/h,dens,vs,Qs,nlay,T,cmin,cmax,dc,eps,ite
	common/para/ci,one,zero,ga,mu,nly
	common/mtrx1/Ish11,Ish12,Ish21,Ish22
	common/mtrx2/exd,exu,Tdsh,Tush,Rdsh,Rush,HTdsh,HRdsh
c
	ci=cmplx(0.0,1.0)
	pi=3.1415926535897936
	pi2=pi+pi
	one=cmplx(1.0,0.0)
	zero=cmplx(0.0,0.0)
	nly=nlay
c
	f=1.0/T
	w=f*pi2
	kc=1
	c1=cmin
	if(c.ne.0.0.and.t.gt.10.0.and.ier.ne.1)c1=c
 	c2=c1+dc
 	ak=w/c1
	do 5 j=1,nlay	
	   mu(j)=vs(j)*vs(j)*dens(j)
	   kb(j)=w/vs(j)
  5	continue
c
 10	do 40 k=kc,2
	  ak2=ak*ak
	  do 20 j=1,nlay
	     ga(j)=sqrt(cmplx(ak2,0.0)-kb(j)*kb(j))
	     if(dreal(ga(j)).lt.0.)ga(j)=-ga(j)
	     exd(j,2)=exp(-ga(j)*h(j))
	     exu(j-1,2)=exd(j,2)
 20	  continue
	  call gcoeflo
	  al(k)=one-Rush(0)*HRdsh(1)
	  ak=w/c2
 40	continue
	kc=2
	if(real(al(1))*real(al(2)).lt.0.0)then
	   nroot=1
	   tmp1=real(al(1))
	   tmp2=real(al(2))
	elseif(aimag(al(1))*aimag(al(2)).lt.0.0)then
	   nroot=2
	   tmp1=aimag(al(1))
	   tmp2=aimag(al(2))
	else
	   c1=c2
	   al(1)=al(2)
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
	tmpl=al(2)
	do 70 k=1,ite
	  c=0.5*(c2+c1)
	  ak=w/c
	  ak2=ak*ak
	  do 50 j=1,nlay
	     ga(j)=sqrt(cmplx(ak2,0.0)-kb(j)*kb(j))
	     if(dreal(ga(j)).lt.0.)ga(j)=-ga(j)
	     exd(j,2)=exp(-ga(j)*h(j))
	     exu(j-1,2)=exd(j,2)
 50	  continue
	  call gcoeflo
	  temp=one-Rush(0)*HRdsh(1)
	  tmp2=dreal(temp)
	  if(nroot.eq.2)tmp2=dimag(temp)
	  if(abs(temp).lt.eps)then
	    goto 75
	  elseif(tmp1*tmp2.le.0.0)then
	    c2=c
	    al(2)=temp
	  else
	    c1=c
	    al(1)=temp
	    tmp1=tmp2
	  endif
 70	continue
	c1=csave
	al(1)=tmpl
	c2=c1+dc
	ak=w/c2
	goto 10
c
c  compute the stress and displacement field
 75	call coeflo
	Edsh(1)=one
	temp=one
	Eush(1)=HRdsh(1)
	do 80 i=2,nlay
	  Edsh(i)=HTdsh(i-1)*temp
	  temp=Edsh(i)
	  Eush(i)=HRdsh(i)*Edsh(i)
 80	continue
c
	I1=0.0
	I2=0.0
	I3=0.0
	w2=w*w
	ak2=ak*ak
	do 90 i=1,nlay
	  a1=0.5*(one-exd(i,2)*exd(i,2))/ga(i)
	  a2=h(i)*exd(i,2)
	  e1=Ish11(i)*Edsh(i)
	  e2=Ish12(i)*Eush(i)
	  dI1=dens(i)*(e1*(e1*a1+e2*a2)+e2*(e1*a2+e2*a1))
	  dI2=mu(i)*(e1*(e1*a1+e2*a2)+e2*(e1*a2+e2*a1))
	  dI3=mu(i)*ga(i)*ga(i)*(e1*(e1*a1-e2*a2)-e2*(e1*a2-e2*a1))
	  I1=I1+dI1
	  I2=I2+dI2
	  I3=I3+dI3
	  cdvs(i)=(ak2*dI2+dI3)/vs(i)
	  cden(i)=(cdvs(i)*vs(i)-w2*dI1)/dens(i)
 90	continue
	c1=c/(ak2*I2)
	c2=0.5*c1
	do 95 i=1,nlay
	  cdvs(i)=cdvs(i)*c1
	  cden(i)=cden(i)*c2
 95	continue
	u=I2/(I1*c)
	ier=0
c
	return
	end
c
c.........................................................
c
	subroutine coeflo
	real*4 mu(50)
	complex*16 ga(50),Ish11(50),Ish12(50),Ish21(50),Ish22(50),
     +	       ci,one,zero
	common/para/ci,one,zero,ga,mu,nlay
	common/mtrx1/Ish11,Ish12,Ish21,Ish22
c	
c  generate displacement-stress matrices
	do 10 n=1,nlay
	  Ish11(n)=one
	  Ish12(n)=one
	  Ish21(n)=-ga(n)*mu(n)
	  Ish22(n)=-Ish21(n)
 10	continue
	return
	end
c
c.........................................................
c
	subroutine gcoeflo
	real*4 mu(50)
	complex*16 ga(50),HTdsh(50),HRdsh(50),exd(50,2),exu(0:50,2),ci,
     +	       Tdsh(50),Tush(50),Rdsh(50),Rush(0:50),det,a1,a2,one,zero
	common/para/ci,one,zero,ga,mu,nlay
	common/mtrx2/exd,exu,Tdsh,Tush,Rdsh,Rush,HTdsh,HRdsh
c
c  modified reflection and transmission coefficients
        Rush(0)=exu(0,2)
	do 10 n=1,nlay-1
	  a1=mu(n)*ga(n)
	  a2=mu(n+1)*ga(n+1)
	  det=a1+a2
	  Tdsh(n)=2*a1*exd(n,2)/det
	  Rush(n)=(a2-a1)*exu(n,2)/det
	  Rdsh(n)=(a1-a2)*exd(n,2)/det
	  Tush(n)=2*a2*exu(n,2)/det
 10	continue
c
c  generalized reflection and transmission coefficients for
c  layers below the source.
	HRdsh(nlay)=zero
	do 20 n=nlay-1,1,-1
 	  HTdsh(n)=Tdsh(n)/(one-Rush(n)*HRdsh(n+1))
	  HRdsh(n)=Rdsh(n)+Tush(n)*HRdsh(n+1)*HTdsh(n)
 20	continue
c	
	return
	end
