      SUBROUTINE period(x,y,n,ofac,hifac,px,py,np,nout,jmax,prob)
      INTEGER jmax,n,nout,np,NMAX
      REAL hifac,ofac,prob,px(np),py(np),x(n),y(n)
      PARAMETER (NMAX=2000)
CU    USES avevar
      INTEGER i,j
      REAL ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,
     *sumsh,sumsy,swtau,var,wtau,xave,xdif,xmax,xmin,yy
      DOUBLE PRECISION arg,wtemp,wi(NMAX),wpi(NMAX),wpr(NMAX),wr(NMAX),
     *TWOPID
      PARAMETER (TWOPID=6.2831853071795865D0)
      nout=0.5*ofac*hifac*n
      if(nout.gt.np) pause 'output arrays too short in period'
      call avevar(y,n,ave,var)
      xmax=x(1)
      xmin=x(1)
      do 11 j=1,n
        if(x(j).gt.xmax)xmax=x(j)
        if(x(j).lt.xmin)xmin=x(j)
11    continue
      xdif=xmax-xmin
      xave=0.5*(xmax+xmin)
      pymax=0.
      pnow=1./(xdif*ofac)
      do 12 j=1,n
        arg=TWOPID*((x(j)-xave)*pnow)
        wpr(j)=-2.d0*sin(0.5d0*arg)**2
        wpi(j)=sin(arg)
        wr(j)=cos(arg)
        wi(j)=wpi(j)
12    continue
      do 15 i=1,nout
        px(i)=pnow
        sumsh=0.
        sumc=0.
        do 13 j=1,n
          c=wr(j)
          s=wi(j)
          sumsh=sumsh+s*c
          sumc=sumc+(c-s)*(c+s)
13      continue
        wtau=0.5*atan2(2.*sumsh,sumc)
        swtau=sin(wtau)
        cwtau=cos(wtau)
        sums=0.
        sumc=0.
        sumsy=0.
        sumcy=0.
        do 14 j=1,n
          s=wi(j)
          c=wr(j)
          ss=s*cwtau-c*swtau
          cc=c*cwtau+s*swtau
          sums=sums+ss**2
          sumc=sumc+cc**2
          yy=y(j)-ave
          sumsy=sumsy+yy*ss
          sumcy=sumcy+yy*cc
          wtemp=wr(j)
          wr(j)=(wr(j)*wpr(j)-wi(j)*wpi(j))+wr(j)
          wi(j)=(wi(j)*wpr(j)+wtemp*wpi(j))+wi(j)
14      continue
        py(i)=0.5*(sumcy**2/sumc+sumsy**2/sums)/var
        if (py(i).ge.pymax) then
          pymax=py(i)
          jmax=i
        endif
        pnow=pnow+1./(ofac*xdif)
15    continue
      expy=exp(-pymax)
      effm=2.*nout/ofac
      prob=effm*expy
      if(prob.gt.0.01)prob=1.-(1.-expy)**effm
      return
      END
