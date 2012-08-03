      SUBROUTINE wtn(a,nn,ndim,isign,wtstep)
      INTEGER isign,ndim,nn(ndim),NMAX
      REAL a(*)
      EXTERNAL wtstep
      PARAMETER (NMAX=1024)
CU    USES wtstep
      INTEGER i1,i2,i3,idim,k,n,nnew,nprev,nt,ntot
      REAL wksp(NMAX)
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 16 idim=1,ndim
        n=nn(idim)
        nnew=n*nprev
        if (n.gt.4) then
          do 15 i2=0,ntot-1,nnew
            do 14 i1=1,nprev
              i3=i1+i2
              do 12 k=1,n
                wksp(k)=a(i3)
                i3=i3+nprev
12            continue
              if (isign.ge.0) then
                nt=n
1               if (nt.ge.4) then
                call wtstep(wksp,nt,isign)
                nt=nt/2
                goto 1
                endif
              else
                nt=4
2               if (nt.le.n) then
                call wtstep(wksp,nt,isign)
                nt=nt*2
                goto 2
                endif
              endif
              i3=i1+i2
              do 13 k=1,n
                a(i3)=wksp(k)
                i3=i3+nprev
13            continue
14          continue
15        continue
        endif
        nprev=nnew
16    continue
      return
      END
