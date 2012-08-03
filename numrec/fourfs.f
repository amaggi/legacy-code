      SUBROUTINE fourfs(iunit,nn,ndim,isign)
      INTEGER ndim,nn(ndim),isign,iunit(4),KBF
      PARAMETER (KBF=128)
CU    USES fourew
      INTEGER j,j12,jk,k,kk,n,mm,kc,kd,ks,kr,nr,ns,nv,jx,mate(4),na,nb,
     *nc,nd
      REAL tempr,tempi,afa(KBF),afb(KBF),afc(KBF)
      DOUBLE PRECISION wr,wi,wpr,wpi,wtemp,theta
      SAVE mate
      DATA mate /2,1,4,3/
      n=1
      do 11 j=1,ndim
        n=n*nn(j)
        if (nn(j).le.1)pause 'invalid dimension or wrong ndim in fourfs'
11    continue
      nv=ndim
      jk=nn(nv)
      mm=n
      ns=n/KBF
      nr=ns/2
      kc=0
      kd=KBF/2
      ks=n
      call fourew(iunit,na,nb,nc,nd)
1     continue
        theta=3.141592653589793d0/(isign*n/mm)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        mm=mm/2
        do 13 j12=1,2
          kr=0
2         continue
            read (iunit(na)) (afa(jx),jx=1,KBF)
            read (iunit(nb)) (afb(jx),jx=1,KBF)
            do 12 j=1,KBF,2
              tempr=sngl(wr)*afb(j)-sngl(wi)*afb(j+1)
              tempi=sngl(wi)*afb(j)+sngl(wr)*afb(j+1)
              afb(j)=afa(j)-tempr
              afa(j)=afa(j)+tempr
              afb(j+1)=afa(j+1)-tempi
              afa(j+1)=afa(j+1)+tempi
12          continue
            kc=kc+kd
            if (kc.eq.mm) then
              kc=0
              wtemp=wr
              wr=wr*wpr-wi*wpi+wr
              wi=wi*wpr+wtemp*wpi+wi
            endif
            write (iunit(nc)) (afa(jx),jx=1,KBF)
            write (iunit(nd)) (afb(jx),jx=1,KBF)
          kr=kr+1
          if (kr.lt.nr) goto 2
          if(j12.eq.1.and.ks.ne.n.and.ks.eq.KBF) then
            na=mate(na)
            nb=na
          endif
          if (nr.eq.0) goto 3
13      continue
3       call fourew(iunit,na,nb,nc,nd)
        jk=jk/2
4       if (jk.eq.1) then
          mm=n
          nv=nv-1
          jk=nn(nv)
        goto 4
        endif
        ks=ks/2
        if (ks.gt.KBF) then
          do 16 j12=1,2
            do 15 kr=1,ns,ks/KBF
              do 14 k=1,ks,KBF
                read (iunit(na)) (afa(jx),jx=1,KBF)
                write (iunit(nc)) (afa(jx),jx=1,KBF)
14            continue
              nc=mate(nc)
15          continue
            na=mate(na)
16        continue
          call fourew(iunit,na,nb,nc,nd)
          goto 1
        else if (ks.eq.KBF) then
          nb=na
          goto 1
        endif
      continue
      j=1
5     continue
        theta=3.141592653589793d0/(isign*n/mm)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        mm=mm/2
        ks=kd
        kd=kd/2
        do 18 j12=1,2
          do 17 kr=1,ns
            read (iunit(na)) (afc(jx),jx=1,KBF)
            kk=1
            k=ks+1
6           continue
              tempr=sngl(wr)*afc(kk+ks)-sngl(wi)*afc(kk+ks+1)
              tempi=sngl(wi)*afc(kk+ks)+sngl(wr)*afc(kk+ks+1)
              afa(j)=afc(kk)+tempr
              afb(j)=afc(kk)-tempr
              afa(j+1)=afc(kk+1)+tempi
              afb(j+1)=afc(kk+1)-tempi
              j=j+2
              kk=kk+2
            if (kk.lt.k) goto 6
            kc=kc+kd
            if (kc.eq.mm) then
              kc=0
              wtemp=wr
              wr=wr*wpr-wi*wpi+wr
              wi=wi*wpr+wtemp*wpi+wi
            endif
            kk=kk+ks
            if (kk.le.KBF) then
              k=kk+ks
              goto 6
            endif
            if (j.gt.KBF) then
              write (iunit(nc)) (afa(jx),jx=1,KBF)
              write (iunit(nd)) (afb(jx),jx=1,KBF)
              j=1
            endif
17        continue
          na=mate(na)
18      continue
        call fourew(iunit,na,nb,nc,nd)
        jk=jk/2
      if (jk.gt.1) goto 5
      mm=n
7     if (nv.gt.1) then
        nv=nv-1
        jk=nn(nv)
        if (jk.eq.1) goto 7
        goto 5
      endif
      return
      END
