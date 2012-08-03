C         @(#)tv.cross.f	1.3       10/31/86  
*
*     Programme principal
 
      parameter (mp=200,nd=200,nt=10)
* Attention: instructions Parameter aussi dans INVTV et LECDGK 
 
*     Tableaux en INPUT:
      Dimension P0(mp), D0(nd), DP0(mp), CP0(mp,mp), CD0(nd,nd),
     *          G(nd,mp), Pinit(mp), Dinit(nd)

*     Tableaux en OUTPUT:
      Dimension P(mp), CP(mp,mp), Resol(mp,mp), Residu(nd), Prof(mp)
      Dimension Sig2r(nd), Sig2i(nd)

*     Tableaux divers:
      Dimension Coupl(7,7), Ndata(nt,nd,nd),nfiltr(nt)
      character*80 titre

*-------------------------------------------------

*     lecture du modele a priori et de sa covariance
*     lecture des donnees et de leur covariance
*     lecture de la matrice des derivees partielles
      call LECDGK(P0,CP0,DP0,D0,CD0,Pinit,Dinit,G,n,m,nd,mp,nvar,Prof,
     *            titre,distcor,Coupl,ncouch,nmode,nfiltr,ienv,Ndata,
     *            Sig2r,Sig2i,nbs,nbt,nt)
      write(*,*) n,m,nvar,ncouch,nbs

*     inversion Tarantola & Valette (1982)
      Call INTV(P0,D0,CP0,CD0,G,Dinit,Pinit,n,m,nd,mp,P,CP,Resol,
     *Residu)

*     impression des resultats
      Call SORDGK(P,CP,Resol,Residu,P0,CP0,DP0,D0,CD0,Dinit,n,m,nd,mp,
     *            titre,nvar, Prof,distcor,Coupl,ncouch,nmode,
     *            nfiltr,ienv,Ndata,Sig2r,Sig2i,nbs,nbt,nt)

      end
