      PROGRAM badluk
      INTEGER ic,icon,idwk,ifrac,im,iybeg,iyend,iyyy,jd,jday,n,julday
      REAL TIMZON,frac
      PARAMETER (TIMZON=-5./24.)
      DATA iybeg,iyend /1900,2000/
CU    USES flmoon,julday
      write (*,'(1x,a,i5,a,i5)') 'Full moons on Friday the 13th from',
     *iybeg,' to',iyend
      do 12 iyyy=iybeg,iyend
        do 11 im=1,12
          jday=julday(im,13,iyyy)
          idwk=mod(jday+1,7)
          if(idwk.eq.5) then
            n=12.37*(iyyy-1900+(im-0.5)/12.)
            icon=0
1           call flmoon(n,2,jd,frac)
            ifrac=nint(24.*(frac+TIMZON))
            if(ifrac.lt.0)then
              jd=jd-1
              ifrac=ifrac+24
            endif
            if(ifrac.gt.12)then
              jd=jd+1
              ifrac=ifrac-12
            else
              ifrac=ifrac+12
            endif
            if(jd.eq.jday)then
              write (*,'(/1x,i2,a,i2,a,i4)') im,'/',13,'/',iyyy
              write (*,'(1x,a,i2,a)') 'Full moon ',ifrac,
     *' hrs after midnight (EST).'
              goto 2
            else
              ic=isign(1,jday-jd)
              if(ic.eq.-icon) goto 2
              icon=ic
              n=n+ic
            endif
            goto 1
2           continue
          endif
11      continue
12    continue
      END
