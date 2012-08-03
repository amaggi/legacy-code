      SUBROUTINE spctrm(p,m,k,ovrlap,w1,w2)
      INTEGER k,m
      REAL p(m),w1(4*m),w2(m)
      LOGICAL ovrlap
CU    USES four1
      INTEGER j,j2,joff,joffn,kk,m4,m43,m44,mm
      REAL den,facm,facp,sumw,w,window
      window(j)=(1.-abs(((j-1)-facm)*facp))
C     window(j)=1.
C     window(j)=(1.-(((j-1)-facm)*facp)**2)
      mm=m+m
      m4=mm+mm
      m44=m4+4
      m43=m4+3
      den=0.
      facm=m
      facp=1./m
      sumw=0.
      do 11 j=1,mm
        sumw=sumw+window(j)**2
11    continue
      do 12 j=1,m
        p(j)=0.
12    continue
      if(ovrlap)then
        read (9,*) (w2(j),j=1,m)
      endif
      do 18 kk=1,k
        do 15 joff=-1,0,1
          if (ovrlap) then
            do 13 j=1,m
              w1(joff+j+j)=w2(j)
13          continue
            read (9,*) (w2(j),j=1,m)
            joffn=joff+mm
            do 14 j=1,m
              w1(joffn+j+j)=w2(j)
14          continue
          else
            read (9,*) (w1(j),j=joff+2,m4,2)
          endif
15      continue
        do 16 j=1,mm
          j2=j+j
          w=window(j)
          w1(j2)=w1(j2)*w
          w1(j2-1)=w1(j2-1)*w
16      continue
        call four1(w1,mm,1)
        p(1)=p(1)+w1(1)**2+w1(2)**2
        do 17 j=2,m
          j2=j+j
          p(j)=p(j)+w1(j2)**2+w1(j2-1)**2+w1(m44-j2)**2+w1(m43-j2)**2
17      continue
        den=den+sumw
18    continue
      den=m4*den
      do 19 j=1,m
        p(j)=p(j)/den
19    continue
      return
      END
