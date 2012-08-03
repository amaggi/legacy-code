      SUBROUTINE fasper(x,y,n,ofac,hifac,wk1,wk2,nwk,nout,jmax,prob)
      INTEGER jmax,n,nout,nwk,MACC
      REAL hifac,ofac,prob,wk1(nwk),wk2(nwk),x(n),y(n)
      PARAMETER (MACC=4)
CU    USES avevar,realft,spread
      INTEGER j,k,ndim,nfreq,nfreqt
      REAL ave,ck,ckk,cterm,cwt,den,df,effm,expy,fac,fndim,hc2wt,hs2wt,
     *hypo,pmax,sterm,swt,var,xdif,xmax,xmin
      nout=0.5*ofac*hifac*n
      nfreqt=ofac*hifac*n*MACC
      nfreq=64
1     if (nfreq.lt.nfreqt) then
        nfreq=nfreq*2
      goto 1
      endif
      ndim=2*nfreq
      if(ndim.gt.nwk) pause 'workspaces too small in fasper'
      call avevar(y,n,ave,var)
      xmin=x(1)
      xmax=xmin
      do 11 j=2,n
        if(x(j).lt.xmin)xmin=x(j)
        if(x(j).gt.xmax)xmax=x(j)
11    continue
      xdif=xmax-xmin
      do 12 j=1,ndim
        wk1(j)=0.
        wk2(j)=0.
12    continue
      fac=ndim/(xdif*ofac)
      fndim=ndim
      do 13 j=1,n
        ck=1.+mod((x(j)-xmin)*fac,fndim)
        ckk=1.+mod(2.*(ck-1.),fndim)
        call spread(y(j)-ave,wk1,ndim,ck,MACC)
        call spread(1.,wk2,ndim,ckk,MACC)
13    continue
      call realft(wk1,ndim,1)
      call realft(wk2,ndim,1)
      df=1./(xdif*ofac)
      k=3
      pmax=-1.
      do 14 j=1,nout
        hypo=sqrt(wk2(k)**2+wk2(k+1)**2)
        hc2wt=0.5*wk2(k)/hypo
        hs2wt=0.5*wk2(k+1)/hypo
        cwt=sqrt(0.5+hc2wt)
        swt=sign(sqrt(0.5-hc2wt),hs2wt)
        den=0.5*n+hc2wt*wk2(k)+hs2wt*wk2(k+1)
        cterm=(cwt*wk1(k)+swt*wk1(k+1))**2/den
        sterm=(cwt*wk1(k+1)-swt*wk1(k))**2/(n-den)
        wk1(j)=j*df
        wk2(j)=(cterm+sterm)/(2.*var)
        if (wk2(j).gt.pmax) then
          pmax=wk2(j)
          jmax=j
        endif
        k=k+2
14    continue
      expy=exp(-pmax)
      effm=2.*nout/ofac
      prob=effm*expy
      if(prob.gt.0.01)prob=1.-(1.-expy)**effm
      return
      END
