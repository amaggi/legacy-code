      SUBROUTINE pwtset(n)
      INTEGER n,NCMAX,ncof,ioff,joff
      PARAMETER (NCMAX=50)
      REAL cc(NCMAX),cr(NCMAX)
      COMMON /pwtcom/ cc,cr,ncof,ioff,joff
      INTEGER k
      REAL sig,c4(4),c12(12),c20(20)
      SAVE c4,c12,c20,/pwtcom/
      DATA c4/0.4829629131445341, 0.8365163037378079,0.2241438680420134,
     *-0.1294095225512604/
      DATA c12 /.111540743350, .494623890398, .751133908021,
     *.315250351709,-.226264693965,-.129766867567,.097501605587, 
     *.027522865530,-.031582039318,.000553842201, .004777257511,
     *-.001077301085/
      DATA c20 /.026670057901, .188176800078, .527201188932,
     *.688459039454, .281172343661,-.249846424327,-.195946274377, 
     *.127369340336, .093057364604,-.071394147166,-.029457536822, 
     *.033212674059,.003606553567,-.010733175483, .001395351747,
     *.001992405295,-.000685856695,-.000116466855,.000093588670,
     *-.000013264203 /
      ncof=n
      sig=-1.
      do 11 k=1,n
        if(n.eq.4)then
          cc(k)=c4(k)
        else if(n.eq.12)then
          cc(k)=c12(k)
        else if(n.eq.20)then
          cc(k)=c20(k)
        else
          pause 'unimplemented value n in pwtset'
        endif
        cr(ncof+1-k)=sig*cc(k)
        sig=-sig
11    continue
      ioff=-n/2
      joff=-n/2
      return
      END
