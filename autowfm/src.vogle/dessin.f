c       program()
c Version faite pour tourner avec la bibliotheque graphique Vogle
C     compil: f77 dessin.f -o dessin -lcore2vogle -lvogle -lxview \
C                          -lolgx -lsunwindow -lpixrect -lX11 -lm

c Pour les besoins de la compilation Linux 07/02/2002
c tout les
c       call InqTextExtent2(carte(1:ncar),ddx,ddy)
c ont ete mis en commentaires (routine pas retrouvee
c dans les librairies linux

	common/fenetre/xmin,xmax,ymin,ymax
	parameter (NDI=500,NK=30)
      dimension x(NK,NDI), y(NK,NDI), xx(NDI), yy(NDI), nper(NK)
      character texte*70
      inicor=0
      ibid=-1
      texte='-----------------------------------------------'
  100 write(*,'(a,$)')' entrer le no de la colonne a dessiner (0=Fin, re
     &nversement d''axe si negatif): '
      read(*,*)  iy
	igroup=abs(iy)
      if(igroup.eq.0) goto 999
   99 write(*,*) texte
   10 write(*,'(a)') ' entrer le nom du fichier a dessiner:'
	write(*,'(a,$)') '     ( 0 pour modifier la colonne):'
      read(*,'(a)') texte
	if(texte(1:lnblnk(texte)).eq.'0') goto 100
*  LECTURE DES DONNEES
      open(10,file=texte,status='old',err=99)
      read(10,'(a)') texte
      write(*,*) texte
      read(10,*) km
	if(km.gt.nk) then
		write(*,'(a,i3)') ' Trop de courbes. Maxi=',NK
		stop
		endif
      do 1 i=1,km
      read(10,*) nper(i)
	if(nper(i).gt.NDI) then
		write(*,'(a,i3)') ' Trop de points dans la courbe', i
		write(*,'(a,i4)') ' Maxi=', NDI
		stop
		endif
      do 2 j=1,nper(i)
      read(10,*) x(i,j),(y(i,j),k=1,igroup-1),y(i,j)
	if(iy.lt.0) then
		xx(j)=x(i,j)
		x(i,j)=y(i,j)
		y(i,j)=-xx(j)
		endif
    2 continue
    1 continue
      close(10)
      if(ibid.lt.2) then
           tmin=x(1,1)
           tmax=tmin
           cmin=y(1,1)
           cmax=cmin
           do 3 i=1,km
           do 3 j=1,nper(i)
           if(x(i,j).lt.tmin) tmin=x(i,j)
           if(x(i,j).gt.tmax) tmax=x(i,j)
           if(y(i,j).lt.cmin) cmin=y(i,j)
           if(y(i,j).gt.cmax) cmax=y(i,j)
    3      continue
           xl=tmax-tmin
           yl=cmax-cmin
           xmin=tmin-xl*.2
           xmax=tmax+xl*.1
           ymin=cmin-yl*.1
           ymax=cmax+yl*.2
      endif
*  DESSIN
      if(inicor.eq.0) then
        write(*,*) ' initialisation du dessin: patience !'
CJJ        call desismo(xmin,xmax,ymin,ymax,0,1,0)
        call initdes(xmin,xmax,ymin,ymax,7,1,0,8)
        call SetCharPrecision(1)
        call SetFont(4)
        call SetCharSize((tmax-tmin)*.02,(cmax-cmin)*.02)
        inicor=1
      endif
      if(ibid.eq.0.or.ibid.eq.1) then
        call NewFrame
      call DelAllRetainSegs
        iseg=0
        call SetWindow(xmin,xmax,ymin,ymax)
	  call setfont(4)
        call SetCharSize((tmax-tmin)*.02,(cmax-cmin)*.02)
      endif
* DESSIN D'UN CADRE ET DES GRADUATIONS
      if(ibid.eq.2) call DelRetainSegment(iseg-1)
      iseg=iseg+1
	call SetLineWidth(0.03)
      call CreateRetainSeg(iseg)
	if(iy.gt.0) then
        call MoveAbs2(tmin,(2*ymax+cmax)/3)
	  else
        call MoveAbs2(tmin,ymin)
	  endif
      call Text(texte)
      call gradue(tmin,tmax,cmin,cmax,iy/igroup)
      call CloseRetainSeg(iseg)
*
*  DESSIN DES COURBES
	call SetLineWidth(0.05)
      do 4 i=1,km
	iseg=iseg+1
	call CreateRetainSeg(iseg)
      do 5 j=1,nper(i)
      xx(j)=x(i,j)
    5 yy(j)=y(i,j)
      call MoveAbs2(xx(1), yy(1))
c     call Polylineabs2(xx,yy,nper(i))
	write(0,'(2a)') ' tailles des traits et espaces du tirete',
     *		' (typiquement: 1 a 10):'
	read(*,*) dplein, dvid
	dplein=dplein*.001
	dvid=dvid*.001
	call tirete(xx,yy,nper(i),dplein,dvid)
        call CloseRetainSeg(iseg)
    4 continue
      
      write(*,'(2a,$)') ' entrer 0 pour sortir, 1 pour un autre dessin,'
     *     ,' 2 pour superposer a la meme echelle:'
      read(*,*) ibid
      if (ibid.eq.0) goto 100
      go to 10
*
  999 call findes
      end 
c-----------------------------------------------------------
	subroutine tirete(x,y,n,dplein,dvid)
	common/fenetre/xmin,xmax,ymin,ymax
	dimension x(1),y(1),dt(2)

c	programme de trace (x,y) en tirete regulier
c	dplein fixe la longueur du tirete plein
c	dvid celle du vide
c	L'unite est telle que la diagonale vaut 1

	ipos=1
	dcur=0.
	dt(1)=dplein
	dt(2)=dvid
	xint=x(1)
	yint=y(1)
CJJ	call InqWindow(xmin,xmax,ymin,ymax)
	xlen=xmax-xmin
	ylen=ymax-ymin

	dtrait=dt(2-ipos)
	call MoveAbs2(x(1),y(1))

	do 1 i=2,n

    2	continue
	d=sqrt(((x(i)-xint)/xlen/2.)**2+((y(i)-yint)/ylen/2.)**2)
	if((dcur+d).le.dtrait) then
		dcur=dcur+d
		xint=x(i)
		yint=y(i)
		if(ipos.eq.1) then
			call LineAbs2(xint,yint)
			else
			call MoveAbs2(xint,yint)
			endif
		goto 1
		endif

c	cas ou on doit changer l'etat de la plume

	xint=xint+(dtrait-dcur)*(x(i)-xint)/d
	yint=yint+(dtrait-dcur)*(y(i)-yint)/d
	if(ipos.eq.1) then
		call LineAbs2(xint,yint)
		else
		call MoveAbs2(xint,yint)
		endif
	ipos=1-ipos
	dcur=0.
	dtrait=dt(2-ipos)
	goto 2

    1	continue
	return
	end
*-----------------------------------------------------------
       subroutine gradue(tmin,tmax,cmin,cmax,isens)

        character*80 texte,carte

        call MoveAbs2(tmin,cmin)
        call LineAbs2(tmax,cmin)
        call LineAbs2(tmax,cmax)
        call LineAbs2(tmin,cmax)
        call LineAbs2(tmin,cmin)

        na=log10(tmax-tmin)
        if((tmax-tmin).lt.1.) na=na-1
        dx=10.**na
	xnx=(int(tmin/dx)+1)*dx
	ncar2=1
	if(dx.ge.(tmax-tmin)/2) then
		dx=dx/2
		if(na.le.0) na=na-1
	endif
	if(dx.ge.(tmax-tmin)/2) then
		dx=dx/2
                ncar2=2
		if(na.le.0) na=na-1
		if(na.eq.1) na=-1
	endif
	if(dx.le.(tmax-tmin)/6) dx=dx*2
    1  	if((xnx-dx).ge.tmin) then
		xnx=xnx-dx
		goto 1
	endif
        ncar1=log10(max(abs(tmin),abs(tmax)))
        ncar=abs(na)+3
	if(na*ncar1.lt.0) ncar=ncar+ncar1
        aa=min(abs(tmin),abs(tmax))
        write(texte,*)'(f',ncar,'.',-min(0,na),')'
        if(ncar.gt.7.and.aa/dx.lt.10.) then
	    ncar=6+ncar2
            write(texte,*) '(e',ncar,'.',ncar2,')'
        endif

	dtest=dx
   11	if(dtest.le.1.) then
	dtest=dtest*10.
	goto  11
	endif
   33	if(dtest.gt.10.) then
	dtest=dtest*0.1
	goto 33
	endif
	igrad=dtest+.5

c	igrad= nb de graduations intermediaires
	isup=(xnx-tmin)*igrad/dx
	xnx=xnx-(isup*dx)/igrad
	ni=-isup
c        call InqTextExtent2('a',ddx,ddy)
         ddy=ddx*1.3*(cmax-cmin)/(tmax-tmin)
   21   write(carte,texte) xnx
        call MoveAbs2(xnx,cmin)
        if(2*ni-igrad*int(2*ni/igrad).eq.0) then
c	    graduation principale
c        call InqTextExtent2(carte(1:ncar),ddx,bid)
        	call MoveAbs2(xnx,cmin)
        	call LineAbs2(xnx,cmin-ddy)
       	call MoveAbs2(xnx,cmax)
        	call LineAbs2(xnx,cmax+ddy*.7)
		if(isens.eq.1) then
		   call MoveAbs2(xnx-ddx/2,cmin-2.*ddy)
		   else
		   call MoveAbs2(xnx-ddx/2,cmax+ddy)
		   endif
        if(ni-igrad*int(ni/igrad).eq.0) call Text(carte(1:ncar))
		else
c	    graduation intermediaire
        	call MoveAbs2(xnx,cmin)
        	call LineAbs2(xnx,cmin-ddy*.5)
       		call MoveAbs2(xnx,cmax)
        	call LineAbs2(xnx,cmax+ddy*.5*.7)
		endif
	ni=ni+1
        xnx=xnx+dx/igrad
        if(xnx.lt.tmax) goto 21

        na=log10(cmax-cmin)
        if((cmax-cmin).lt.1.) na=na-1
        dx=10.**na
	ncar2=1
	xnx=(int(cmin/dx)+1)*dx
	if(dx.ge.(cmax-cmin)/2) then
		dx=dx/2
		if(na.le.0) na=na-1
	endif
	if(dx.ge.(cmax-cmin)/2) then
		dx=dx/2
		ncar2=2
		if(na.le.0) na=na-1
		if(na.eq.1) na=-1
	endif
	if(dx.le.(cmax-cmin)/6) dx=dx*2
    2  	if((xnx-dx).ge.cmin) then
		xnx=xnx-dx
		goto 2
	endif
        ncar1=log10(max(abs(cmin),abs(cmax)))
        ncar=abs(na)+3
	if(na*ncar1.lt.0) ncar=ncar+ncar1
        aa=min(abs(cmin),abs(cmax))
        write(texte,*)'(f',ncar,'.',-min(0,na),')'
        if(ncar.gt.7.and.aa/dx.lt.10.) then
	    ncar=6+ncar2
            write(texte,*) '(e',ncar,'.',ncar2,')'
        endif

	dtest=dx
  111	if(dtest.le.1.) then
	dtest=dtest*10.
	goto 111
	endif
  333	if(dtest.gt.10.) then
	dtest=dtest*0.1
	goto 333
	endif
	igrad=dtest+.5

	isup=(xnx-cmin)*igrad/dx
	xnx=xnx-(isup*dx)/igrad
	ni=-isup
   22   write(carte,texte) isens*xnx
c       call InqTextExtent2(carte(1:ncar),ddx,ddy)
        if(2*ni-igrad*int(2*ni/igrad).eq.0) then
c	    graduation principale
        	call MoveAbs2(tmin,xnx)
        	call LineAbs2(tmin-1.5*ddx/ncar,xnx)
        	call MoveAbs2(tmin-ddx*(ncar+1)/ncar,xnx)
        if(ni-igrad*int(ni/igrad).eq.0) call Text(carte(1:ncar))
        	call MoveAbs2(tmax,xnx)
        	call LineAbs2(tmax+1.5*ddx/ncar*.7,xnx)
		else
c	    graduation intermediaire
        	call MoveAbs2(tmin,xnx)
        	call LineAbs2(tmin-.75*ddx/ncar,xnx)
        	call MoveAbs2(tmax,xnx)
        	call LineAbs2(tmax+.75*ddx/ncar*.7,xnx)
		endif
	ni=ni+1
        xnx=xnx+dx/igrad
        if(xnx.lt.cmax) goto 22

        return
        end
