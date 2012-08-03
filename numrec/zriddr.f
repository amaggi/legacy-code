      FUNCTION zriddr(func,x1,x2,xacc)
      INTEGER MAXIT
      REAL zriddr,x1,x2,xacc,func,UNUSED
      PARAMETER (MAXIT=60,UNUSED=-1.11E30)
      EXTERNAL func
CU    USES func
      INTEGER j
      REAL fh,fl,fm,fnew,s,xh,xl,xm,xnew
      fl=func(x1)
      fh=func(x2)
      if((fl.gt.0..and.fh.lt.0.).or.(fl.lt.0..and.fh.gt.0.))then
        xl=x1
        xh=x2
        zriddr=UNUSED
        do 11 j=1,MAXIT
          xm=0.5*(xl+xh)
          fm=func(xm)
          s=sqrt(fm**2-fl*fh)
          if(s.eq.0.)return
          xnew=xm+(xm-xl)*(sign(1.,fl-fh)*fm/s)
          if (abs(xnew-zriddr).le.xacc) return
          zriddr=xnew
          fnew=func(zriddr)
          if (fnew.eq.0.) return
          if(sign(fm,fnew).ne.fm) then
            xl=xm
            fl=fm
            xh=zriddr
            fh=fnew
          else if(sign(fl,fnew).ne.fl) then
            xh=zriddr
            fh=fnew
          else if(sign(fh,fnew).ne.fh) then
            xl=zriddr
            fl=fnew
          else
            pause 'never get here in zriddr'
          endif
          if(abs(xh-xl).le.xacc) return
11      continue
        pause 'zriddr exceed maximum iterations'
      else if (fl.eq.0.) then
        zriddr=x1
      else if (fh.eq.0.) then
        zriddr=x2
      else
        pause 'root must be bracketed in zriddr'
      endif
      return
      END
