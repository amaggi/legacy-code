      SUBROUTINE anneal(x,y,iorder,ncity)
      INTEGER ncity,iorder(ncity)
      REAL x(ncity),y(ncity)
CU    USES irbit1,metrop,ran3,revcst,revers,trncst,trnspt
      INTEGER i,i1,i2,idec,idum,iseed,j,k,nlimit,nn,nover,nsucc,n(6),
     *irbit1
      REAL de,path,t,tfactr,ran3,alen,x1,x2,y1,y2
      LOGICAL ans
      alen(x1,x2,y1,y2)=sqrt((x2-x1)**2+(y2-y1)**2)
      nover=100*ncity
      nlimit=10*ncity
      tfactr=0.9
      path=0.0
      t=0.5
      do 11 i=1,ncity-1
        i1=iorder(i)
        i2=iorder(i+1)
        path=path+alen(x(i1),x(i2),y(i1),y(i2))
11    continue
      i1=iorder(ncity)
      i2=iorder(1)
      path=path+alen(x(i1),x(i2),y(i1),y(i2))
      idum=-1
      iseed=111
      do 13 j=1,100
        nsucc=0
        do 12 k=1,nover
1         n(1)=1+int(ncity*ran3(idum))
          n(2)=1+int((ncity-1)*ran3(idum))
          if (n(2).ge.n(1)) n(2)=n(2)+1
          nn=1+mod((n(1)-n(2)+ncity-1),ncity)
          if (nn.lt.3) goto 1
          idec=irbit1(iseed)
          if (idec.eq.0) then
            n(3)=n(2)+int(abs(nn-2)*ran3(idum))+1
            n(3)=1+mod(n(3)-1,ncity)
            call trncst(x,y,iorder,ncity,n,de)
            call metrop(de,t,ans)
            if (ans) then
              nsucc=nsucc+1
              path=path+de
              call trnspt(iorder,ncity,n)
            endif
          else
            call revcst(x,y,iorder,ncity,n,de)
            call metrop(de,t,ans)
            if (ans) then
              nsucc=nsucc+1
              path=path+de
              call revers(iorder,ncity,n)
            endif
          endif
          if (nsucc.ge.nlimit) goto 2
12      continue
2       write(*,*)
        write(*,*) 'T =',t,' Path Length =',path
        write(*,*) 'Successful Moves: ',nsucc
        t=t*tfactr
        if (nsucc.eq.0) return
13    continue
      return
      END
