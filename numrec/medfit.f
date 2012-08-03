      SUBROUTINE medfit(x,y,ndata,a,b,abdev)
      INTEGER ndata,NMAX,ndatat
      PARAMETER (NMAX=1000)
      REAL a,abdev,b,x(ndata),y(ndata),arr(NMAX),xt(NMAX),yt(NMAX),aa,
     *abdevt
      COMMON /arrays/ xt,yt,arr,aa,abdevt,ndatat
CU    USES rofunc
      INTEGER j
      REAL b1,b2,bb,chisq,del,f,f1,f2,sigb,sx,sxx,sxy,sy,rofunc
      sx=0.
      sy=0.
      sxy=0.
      sxx=0.
      do 11 j=1,ndata
        xt(j)=x(j)
        yt(j)=y(j)
        sx=sx+x(j)
        sy=sy+y(j)
        sxy=sxy+x(j)*y(j)
        sxx=sxx+x(j)**2
11    continue
      ndatat=ndata
      del=ndata*sxx-sx**2
      aa=(sxx*sy-sx*sxy)/del
      bb=(ndata*sxy-sx*sy)/del
      chisq=0.
      do 12 j=1,ndata
        chisq=chisq+(y(j)-(aa+bb*x(j)))**2
12    continue
      sigb=sqrt(chisq/del)
      b1=bb
      f1=rofunc(b1)
      b2=bb+sign(3.*sigb,f1)
      f2=rofunc(b2)
1     if(f1*f2.gt.0.)then
        bb=2.*b2-b1
        b1=b2
        f1=f2
        b2=bb
        f2=rofunc(b2)
        goto 1
      endif
      sigb=0.01*sigb
2     if(abs(b2-b1).gt.sigb)then
        bb=0.5*(b1+b2)
        if(bb.eq.b1.or.bb.eq.b2)goto 3
        f=rofunc(bb)
        if(f*f1.ge.0.)then
          f1=f
          b1=bb
        else
          f2=f
          b2=bb
        endif
        goto 2
      endif
3     a=aa
      b=bb
      abdev=abdevt/ndata
      return
      END
