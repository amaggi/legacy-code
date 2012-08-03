      SUBROUTINE miser(func,region,ndim,npts,dith,ave,var)
      INTEGER ndim,npts,MNPT,MNBS,MAXD,NSTACK
      REAL ave,dith,var,region(2*ndim),func,TINY,BIG,PFAC
      PARAMETER (MNPT=15,MNBS=4*MNPT,MAXD=10,TINY=1.e-30,BIG=1.e30,
     *NSTACK=1000,PFAC=0.1)
      EXTERNAL func
CU    USES func,ranpt
      INTEGER iran,j,jb,jstack,n,naddr,np,npre,nptl,nptr,nptt
      REAL avel,fracl,fval,rgl,rgm,rgr,s,sigl,siglb,sigr,sigrb,sum,sumb,
     *summ,summ2,varl,fmaxl(10),fmaxr(10),fminl(10),fminr(10),pt(10),
     *rmid(10),stack(NSTACK),stf(9)
      EQUIVALENCE (stf(1),avel),(stf(2),varl),(stf(3),jb),(stf(4),nptr),
     *(stf(5),naddr),(stf(6),rgl),(stf(7),rgm),(stf(8),rgr),(stf(9),
     *fracl)
      SAVE iran
      DATA iran /0/
      jstack=0
      nptt=npts
1     continue
      if (nptt.lt.MNBS) then
        np=abs(nptt)
        summ=0.
        summ2=0.
        do 11 n=1,np
          call ranpt(pt,region,ndim)
          fval=func(pt)
          summ=summ+fval
          summ2=summ2+fval**2
11      continue
        ave=summ/np
        var=max(TINY,(summ2-summ**2/np)/np**2)
      else
        npre=max(int(nptt*PFAC),MNPT)
        do 12 j=1,ndim
          iran=mod(iran*2661+36979,175000)
          s=sign(dith,float(iran-87500))
          rmid(j)=(0.5+s)*region(j)+(0.5-s)*region(j+ndim)
          fminl(j)=BIG
          fminr(j)=BIG
          fmaxl(j)=-BIG
          fmaxr(j)=-BIG
12      continue
        do 14 n=1,npre
          call ranpt(pt,region,ndim)
          fval=func(pt)
          do 13 j=1,ndim
            if(pt(j).le.rmid(j))then
              fminl(j)=min(fminl(j),fval)
              fmaxl(j)=max(fmaxl(j),fval)
            else
              fminr(j)=min(fminr(j),fval)
              fmaxr(j)=max(fmaxr(j),fval)
            endif
13        continue
14      continue
        sumb=BIG
        jb=0
        siglb=1.
        sigrb=1.
        do 15 j=1,ndim
          if(fmaxl(j).gt.fminl(j).and.fmaxr(j).gt.fminr(j))then
            sigl=max(TINY,(fmaxl(j)-fminl(j))**(2./3.))
            sigr=max(TINY,(fmaxr(j)-fminr(j))**(2./3.))
            sum=sigl+sigr
            if (sum.le.sumb) then
              sumb=sum
              jb=j
              siglb=sigl
              sigrb=sigr
            endif
          endif
15      continue
        if (jb.eq.0) jb=1+(ndim*iran)/175000
        rgl=region(jb)
        rgm=rmid(jb)
        rgr=region(jb+ndim)
        fracl=abs((rgm-rgl)/(rgr-rgl))
        nptl=MNPT+(nptt-npre-2*MNPT)*fracl*siglb/(fracl*siglb+
     *(1.-fracl)*sigrb)
        nptr=nptt-npre-nptl
        region(jb+ndim)=rgm
        naddr=1
        do 16 j=1,9
          stack(jstack+j)=stf(j)
16      continue
        jstack=jstack+9
        nptt=nptl
        goto 1
10      continue
        avel=ave
        varl=var
        region(jb)=rgm
        region(jb+ndim)=rgr
        naddr=2
        do 17 j=1,9
          stack(jstack+j)=stf(j)
17      continue
        jstack=jstack+9
        nptt=nptr
        goto 1
20      continue
        region(jb)=rgl
        ave=fracl*avel+(1.-fracl)*ave
        var=fracl**2*varl+(1.-fracl)**2*var
      endif
      if (jstack.ne.0) then
        jstack=jstack-9
        do 18 j=1,9
          stf(j)=stack(jstack+j)
18      continue
        goto (10,20),naddr
        pause 'miser: never get here'
      endif
      return
      END
