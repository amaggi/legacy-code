      SUBROUTINE mglin(u,n,ncycle)
      INTEGER n,ncycle,NPRE,NPOST,NG,MEMLEN
      DOUBLE PRECISION u(n,n)
      PARAMETER (NG=5,MEMLEN=13*2**(2*NG)/3+14*2**NG+8*NG-100/3)
      PARAMETER (NPRE=1,NPOST=1)
CU    USES addint,copy,fill0,interp,maloc,relax,resid,rstrct,slvsml
      INTEGER j,jcycle,jj,jpost,jpre,mem,nf,ngrid,nn,ires(NG),irho(NG),
     *irhs(NG),iu(NG),maloc
      DOUBLE PRECISION z
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
      call slvsml(z(iu(1)),z(irho(1)))
      ngrid=NG
      do 16 j=2,ngrid
        nn=2*nn-1
        iu(j)=maloc(nn**2)
        irhs(j)=maloc(nn**2)
        ires(j)=maloc(nn**2)
        call interp(z(iu(j)),z(iu(j-1)),nn)
        if (j.ne.ngrid) then
          call copy(z(irhs(j)),z(irho(j)),nn)
        else
          call copy(z(irhs(j)),u,nn)
        endif
        do 15 jcycle=1,ncycle
            nf=nn
            do 12 jj=j,2,-1
          do 11 jpre=1,NPRE
              call relax(z(iu(jj)),z(irhs(jj)),nf)
11        continue
          call resid(z(ires(jj)),z(iu(jj)),z(irhs(jj)),nf)
          nf=nf/2+1
          call rstrct(z(irhs(jj-1)),z(ires(jj)),nf)
          call fill0(z(iu(jj-1)),nf)
12          continue
            call slvsml(z(iu(1)),z(irhs(1)))
            nf=3
            do 14 jj=2,j
          nf=2*nf-1
          call addint(z(iu(jj)),z(iu(jj-1)),z(ires(jj)),nf)
          do 13 jpost=1,NPOST
              call relax(z(iu(jj)),z(irhs(jj)),nf)
13        continue
14          continue
15      continue
16    continue
      call copy(u,z(iu(ngrid)),n)
      return
      END
