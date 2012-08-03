      SUBROUTINE mgfas(u,n,maxcyc)
      INTEGER maxcyc,n,NPRE,NPOST,NG,MEMLEN
      DOUBLE PRECISION u(n,n),ALPHA
      PARAMETER (NG=5,MEMLEN=17*2**(2*NG)/3+18*2**NG+10*NG-86/3)
      PARAMETER (NPRE=1,NPOST=1,ALPHA=.33d0)
CU    USES anorm2,copy,interp,lop,maloc,matadd,matsub,relax2,rstrct,slvsm2
      INTEGER j,jcycle,jj,jm1,jpost,jpre,mem,nf,ngrid,nn,irho(NG),
     *irhs(NG),itau(NG),itemp(NG),iu(NG),maloc
      DOUBLE PRECISION res,trerr,z,anorm2
      COMMON /memory/ z(MEMLEN),mem
      mem=0
      nn=n/2+1
      ngrid=NG-1
      irho(ngrid)=maloc(nn**2)
      call rstrct(z(irho(ngrid)),u,nn)
1     if (nn.gt.3) then
        nn=nn/2+1
        ngrid=ngrid-1
        irho(ngrid)=maloc(nn**2)
        call rstrct(z(irho(ngrid)),z(irho(ngrid+1)),nn)
      goto 1
      endif
      nn=3
      iu(1)=maloc(nn**2)
      irhs(1)=maloc(nn**2)
      itau(1)=maloc(nn**2)
      itemp(1)=maloc(nn**2)
      call slvsm2(z(iu(1)),z(irho(1)))
      ngrid=NG
      do 16 j=2,ngrid
        nn=2*nn-1
        iu(j)=maloc(nn**2)
        irhs(j)=maloc(nn**2)
        itau(j)=maloc(nn**2)
        itemp(j)=maloc(nn**2)
        call interp(z(iu(j)),z(iu(j-1)),nn)
        if (j.ne.ngrid) then
          call copy(z(irhs(j)),z(irho(j)),nn)
        else
          call copy(z(irhs(j)),u,nn)
        endif
        do 15 jcycle=1,maxcyc
            nf=nn
            do 12 jj=j,2,-1
          do 11 jpre=1,NPRE
              call relax2(z(iu(jj)),z(irhs(jj)),nf)
11        continue
          call lop(z(itemp(jj)),z(iu(jj)),nf)
          nf=nf/2+1
          jm1=jj-1
          call rstrct(z(itemp(jm1)),z(itemp(jj)),nf)
          call rstrct(z(iu(jm1)),z(iu(jj)),nf)
          call lop(z(itau(jm1)),z(iu(jm1)),nf)
          call matsub(z(itau(jm1)),z(itemp(jm1)),z(itau(jm1)),nf)
          if(jj.eq.j)trerr=ALPHA*anorm2(z(itau(jm1)),nf)
          call rstrct(z(irhs(jm1)),z(irhs(jj)),nf)
          call matadd(z(irhs(jm1)),z(itau(jm1)),z(irhs(jm1)),nf)
12          continue
            call slvsm2(z(iu(1)),z(irhs(1)))
            nf=3
            do 14 jj=2,j
          jm1=jj-1
          call rstrct(z(itemp(jm1)),z(iu(jj)),nf)
          call matsub(z(iu(jm1)),z(itemp(jm1)),z(itemp(jm1)),nf)
          nf=2*nf-1
          call interp(z(itau(jj)),z(itemp(jm1)),nf)
          call matadd(z(iu(jj)),z(itau(jj)),z(iu(jj)),nf)
          do 13 jpost=1,NPOST
              call relax2(z(iu(jj)),z(irhs(jj)),nf)
13        continue
14          continue
            call lop(z(itemp(j)),z(iu(j)),nf)
            call matsub(z(itemp(j)),z(irhs(j)),z(itemp(j)),nf)
            res=anorm2(z(itemp(j)),nf)
            if(res.lt.trerr)goto 2
15      continue
2       continue
16    continue
      call copy(u,z(iu(ngrid)),n)
      return
      END
