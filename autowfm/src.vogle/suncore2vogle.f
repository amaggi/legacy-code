c Compilation sous Linux 07/02/2002
c tout les call linewidth() ont ete mis en commentaires
c la routine n'etant pas identifiee dans les librairies X
c de linux
c Bug fixes F.Tilmann
c   - setfont (23.3.99)
c-----------------------------------------------------------------------

c     Routines de conversion d'ordres Suncore en Vogle et Postscript.
 
c-----------------------------------------------------------------------
c     Reprend toutes les routines Suncore de "convsuncore.h", plus
c     quelques complementaires (point, rectangle, circle) pour une 
c     conversion en ordres Vogle et PostScript (limitation au 2D).
c     Un fichier PostScript est cree a chaque appel de "initdes".
c     Le fichier PostScript n'est clos qu'a l'appel de "findes".
c     Pour plus de details sur la syntaxe et l'utilisation des rou-
c     tines, se reporter aux commentaires en entete de la routine.
c    ! Non standard f77 calls: adapt comments cSUN and cHP below !
c     Tested on Sun and HP.  Hugues Dufumier, Fevrier 98.
 
c-----------------------------------------------------------------------
 
c-----------------------------------------------------------------------
      subroutine initdes(X0,X1,Y0,Y1,IVIEW,INPUT,IDIM,ICOUL)
c-----------------------------------------------------------------------
*********** Procedure simplifiee pour dessiner sous VOGLE **************
*           prevoit l'utilisation de la souris comme "locator"
*           possibilite de dessiner en couleur, en 2D et 3D
*           selection table de gris ou de couleur (incluses ou defaut)
*           definit un "viewport" standard pour ronds ronds
*           initialise un segment permanent (dessin direct)
*           permet une sortie postcript (argument.ps)
*
* X0,X1,Y0,Y1 : fenetre 2D dans l'espace utilisateur.
* IVIEW   definit la forme du viewport et sa taille (de |IVIEW|=1 a 9)
*         <0: diverses tailles utiles orthonormees (ronds ronds), 
*         >0: diverses tailles fixes rapport 0.75 (occupation maxi).
*         (9 = plein ecran en landscape, 7 = tient dans format A4)
* INPUT=0 pas d'effet ici: souris toujours initialisee.
* IDIM=0  dessin en deux dimensions, trois dimensions sinon (incomplet).
* ICOUL=  0: dessin noir et blanc sur ecran monochrome (2 couleurs),
*         1: palette de gris pour ecran couleur (128 teintes de gris),
*         2: palette arc en ciel pour ecran couleur (58 couleurs),
*         3: palette feux rouge-orange-vert ecran couleur (128 couleurs),
*         4: palette tomo rouge-bleu pour ecran couleur (128 couleurs),
*         5: palette mer bleu-vert pour ecran couleur (128 couleurs),
*         autres: palette Vogle par defaut (8 couleurs de base).

      character fich*80,nomps*70,cday*24,device*10,open*1
      common/psdim/PSXDIM,PSYDIM
      common/line/ilin,lx,ly

c     adapt this parameter to your screen size (IVIEW=9 -> full screen)
      parameter(SCREENSIZE=1140)

c     palettes predefinies: dans tous les cas 0=noir, ncol-1=blanc
      parameter(NCOU=58)
      common/rgb/icol,ncol,red(256),gre(256),blu(256)
      integer iredc(NCOU),igrec(NCOU),ibluc(NCOU)
      integer iredv(8),igrev(8),ibluv(8),getdepth

c     table "arc en ciel" (degrade 50 couleurs + 8 speciales) (ICOUL=2)
      data iredc/0,255,255,255,255,255,255,255,255,255,255,250,230,200,
     &150,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,150,170,190,210,200,180,
     &160,140,120,100,80,70,60,50,40,30,0,90,150,210,255,0,200,0,255/
      data igrec/0,0,70,100,120,140,160,180,200,220,240,250,255,255,255
     &,255,255,255,255,255,255,255,250,250,240,220,200,180,150,130,100,0
     &,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,90,150,210,250,160,0,0,255/
      data ibluc/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,120,140,180,200,220,
     &240,250,250,255,255,255,255,255,255,255,255,250,250,250,250,240,
     &230,210,190,170,150,130,130,120,110,100,90,80,0,90,150,210,180,40,
     &0,150,255/

c     table Vogle 8 couleurs par defaut (ICOUL<0 ou >5)
      data iredv/0,255,0,255,0,255,0,255/
      data igrev/0,0,255,255,0,0,255,255/
      data ibluv/0,0,0,0,255,255,255,255/

cJJL2002
	data iopen/0/

      icol=ICOUL
      iwin=abs(IVIEW)
      ilin=0
      lx=0
      ly=0

c     au premier appel de initdes
      if(iopen.eq.0) then

cSUN:
       call getarg(0,nomps)
cHP:   call igetarg(1,nomps,70)
c      initialisation Vogle pour recherche device
       print*,'% Vogle: Set the environment variable VDEVICE in'
       print*,'% your .cshrc file to the desired output device.'
cJJL2002       call vinit(device)
       call vinit("X11")
       call vgetdev(device)
       INDEXMAX=2**(getdepth()-1)
       call vexit
c     aux appels suivants de initdes
      else
       call findes
      endif

c     initialisation du fichier postscript sous unite 66
cSUN:
      call fdate(cday)
cHP:  call date(cday)
      write(open,'(i1)') mod(iopen,10)
      iopen=iopen+1
      fich=nomps(1:lnblnk(nomps))//'.postscript'
      print*,'% Output PostScript File: ',fich(:lnblnk(fich))
      open(66,file=fich)
      write(66,'(a)')'%!PS-SuncoreToVogle+PS'
      write(66,'(a)')'%statusdict /manualfeed true put'
c     gsave initial permettant la concatenation de fichiers PS
      write(66,'(a)')'gsave newpath'
      write(66,'(a)')'/s {stroke} def'
CJJ      write(66,'(a)')'/n {newpath} def'
      write(66,'(a)')'/m {moveto} def'
      write(66,'(a)')'/l {lineto} def'
      write(66,'(a)')'/rm {rmoveto} def'
      write(66,'(a)')'/rl {rlineto} def'
      write(66,'(a)')'/gs {gsave} def'
      write(66,'(a)')'/gr {grestore} def'
      write(66,'(a)')'/cp {closepath} def'
      write(66,'(a)')'/mx {matrix} def'
CJJ      write(66,'(a)')'/sg {setgray} def'
      write(66,'(a)')'/sc {setrgbcolor} def'
      write(66,'(a)')'/sh {show} def'
      write(66,'(a)')'/f {fill} def'
      write(66,'(a)')'/centr {stringwidth 2 div neg exch 2 div neg exch
     &rm} def'
      write(66,'(a)')'/right {stringwidth neg exch neg exch rm} def'
      write(66,'(a)')'/p {/yy exch def /xx exch def s newpath xx yy
     & m xx yy 1 0 360 arc} def'
      write(66,'(a)')'/r {/ys exch def /yi exch def /xs exch def /xi exc
     &h def s newpath xi yi m xs yi l xs ys l xi ys l cp} def'
      write(66,'(a)')'/txt {/hy exch def /hx exch def /ang exch def /dy
     &exch def /dx exch def currentfont dx dy mx translate ang mx rotate
     & mx concatmatrix hx hy mx scale exch mx concatmatrix makefont setf
     &ont} def'

c     echelle postscript suivant option IVIEW
c     (de 1 a 9 = plein ecran en landscape, 7 tient dans un A4)
c     landscape
      if(X1-X0.ge.Y1-Y0) then
        PSXDIM=7260.*float(iwin)/9.
        PSYDIM=5900.*float(iwin)/9.
        xscl=0.141
        yscl=xscl*PSXDIM/PSYDIM
        write(66,'(a)')'580 20 translate 90 rotate'
      else
c     portrait
        PSXDIM=5900.*float(iwin)/9.
        PSYDIM=7260.*float(iwin)/9.
        yscl=0.141
        xscl=yscl*PSYDIM/PSXDIM
        write(66,'(a)')'20 20 translate'
      endif
      write(66,'(2f7.4,a)') xscl,yscl,' scale'

c     viewport orthonorme ou occupation maximale papier
      if(IVIEW.lt.0) then
c     viewport pour fenetre orthonormee (ronds ronds)
        wplarg=(X1-X0)/amax1((X1-X0),(Y1-Y0))
        wphaut=(Y1-Y0)/amax1((X1-X0),(Y1-Y0))
      else
c     fenetre non orthonormee (occupation PS max mais ellipses possibles)
        if(X1-X0.ge.Y1-Y0) then
          wplarg=1.0
          wphaut=0.7
         else
          wplarg=0.7
          wphaut=1.0
         endif
      endif

c     le viewport maxi en Postscript A4 est de 1 par sqrt(2)/2 = 0.7
      if(X1-X0.ge.Y1-Y0.and.wphaut.gt.0.7) then
        wplarg=wplarg*0.7/wphaut
        wphaut=0.7
      endif
      if(X1-X0.lt.Y1-Y0.and.wplarg.gt.0.7) then
        wphaut=wphaut*0.7/wplarg
        wplarg=0.7
      endif

c     creation fenetre Vogle avec memes proportions que viewport
c     coin haut gauche, taille fenetre suivant option IVIEW
      call prefposition(0,0)
      FAC=SCREENSIZE*float(iwin)/9.
      call prefsize(nint(FAC*wplarg),nint(FAC*wphaut))
      call vinit(device)

c     definit une fenetre 2D ou 3D sans perspective
      if(IDIM.eq.0) then
        call setviewport2(0.0,wplarg,0.0,wphaut)
        call setwindow(X0,X1,Y0,Y1)
      else
        call setviewport3(0.0,wplarg,0.0,wphaut,0.0,1.0)
        call setwindow3(X0,X1,Y0,Y1,-100.,100.)
      endif

c     choix table de gris ou de couleur
      if(INDEXMAX.eq.1.and.icol.ne.0) then
        print*,'% Monochrome device: Grey/Color option ICOUL set to 0 !'
        icol=0
        endif
      if(icol.eq.0) then
c       noir et blanc uniquement pour ecran monochrome
        ncol=2
      else if(icol.eq.1) then
c       table de gris pour ecran couleur (128 nuances)
        ncol=min0(128,INDEXMAX)
        do 51 i=1,ncol
        red(i)=float(i-1)/float(ncol-1)
        gre(i)=float(i-1)/float(ncol-1)
  51    blu(i)=float(i-1)/float(ncol-1)
      else if(icol.eq.2) then
c       table de couleur "arc en ciel" (50+8 couleurs)
        ncol=min0(NCOU,INDEXMAX)
        do 52 i=1,ncol
        red(i)=float(iredc(i))/255.
        gre(i)=float(igrec(i))/255.
  52    blu(i)=float(ibluc(i))/255.
      else if(icol.eq.3) then
c       table de couleur feux rouge-orange-vert (128 couleurs)
        ncol=min0(128,INDEXMAX)
        n2=ncol/2
        do 53 i=1,ncol
        red(i)=1.
        gre(i)=1.
        if(i.le.n2+2) gre(i)=float(i-2)/float(n2+1)
        if(i.gt.n2-2) red(i)=float(ncol-1-i)/float(n2+1)
  53    blu(i)=0.4-0.4*abs(n2-i)/float(n2)
      else if(icol.eq.4) then
c       table de couleur tomo rouge-bleu (128 couleurs)
        ncol=min0(128,INDEXMAX)
        n2=ncol/2
        do 54 i=2,ncol-1
        red(i)=1.
        blu(i)=1.
        if(i.le.n2+2) blu(i)=float(i-2)/float(n2+1)
        if(i.gt.n2-2) red(i)=float(ncol-1-i)/float(n2+1)
  54    gre(i)=1.-0.8*abs(n2-i)/float(n2)
      else if(icol.eq.5) then
c       table de couleur mer bleu-vert (128 couleurs)
        ncol=min0(128,INDEXMAX)
        n2=ncol/2
        do 55 i=1,ncol
        blu(i)=1.
        gre(i)=1.
        if(i.le.n2+2) gre(i)=float(i-2)/float(n2+1)
        if(i.gt.n2-2) blu(i)=float(ncol-1-i)/float(n2+1)
  55    red(i)=0.6-0.6*abs(n2-i)/float(n2)
      else
c       charge la table Vogle 8 couleurs par defaut
        ncol=8
        do 59 i=1,ncol
        red(i)=float(iredv(i))/255.   
        gre(i)=float(igrev(i))/255.   
  59    blu(i)=float(ibluv(i))/255. 
      endif
c     dans tous les cas, indice 0 = noir, indice ncol-1 = blanc
      red(1)=0.
      gre(1)=0.
      blu(1)=0.
      red(ncol)=1.
      gre(ncol)=1.
      blu(ncol)=1.

c     chargement table de gris ou de couleur (sinon table par defaut)
      print*,'% Vogle output Device: ',device,ncol,' colors.'
      if(icol.ge.0.and.icol.le.5) then
        do 60 i=0,ncol-1
        ir=nint(255.*red(i+1))
        ig=nint(255.*gre(i+1))
        ib=nint(255.*blu(i+1))
  60    call mapcolor(i,ir,ig,ib)
      endif

c     couleurs de fond et de dessin par defaut: noir sur blanc
      call color(ncol-1)
      call clear
      call settextindex(0)
      call setfillindex(0)
      call setlineindex(0)

c     dimension de base d'un "dot" pour les linestyles
      dashlen=abs(X1-X0)/6000.
      call setdash(dashlen)
      call linestyle('')
      write(66,'(a)') '5 setlinewidth'

c     font par defaut
      call font("times.r")
      call ycentertext()
      write(66,'(a)')'/Times-Roman findfont 45 scalefont setfont'
      labely=iconvy(Y1)+18
      labelx=iconvx(X1)
      write(66,*) '     2',labely,' m ('//fich(:lnblnk(fich))//') sh'
      write(66,*) labelx,labely,' m (suncore2vogle - '//cday//') right'
      write(66,*) '(suncore2vogle - '//cday//') sh'

c     clipping
      call clipping(on)
      write(66,*) iconvx(X0)-1,iconvy(Y0)-1,' m'
      write(66,*) iconvx(X1)+1,iconvy(Y0)-1,' l'
      write(66,*) iconvx(X1)+1,iconvy(Y1)+1,' l'
      write(66,*) iconvx(X0)-1,iconvy(Y1)+1,' l'
      write(66,*) iconvx(X0)-1,iconvy(Y0)-1,' l  cp clip s'

c     gsave memoire de configuration pour chaque newframe
      write(66,'(a)')'gs'
c     ouvre un segment permanent
CJJ      call Pcreateretseg(1)

      return
      end

c-----------------------------------------------------------------------
      subroutine findes
c-----------------------------------------------------------------------
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup

      if(nseg(1000).eq.-1) call closetempseg()
c     retrieves initial configuration to allow concatenation of PS files
      write(66,'(a)')'showpage gr gr'
      close(66)
      call vexit
      return
      end

c-----------------------------------------------------------------------

c     Routines de conversion d'ordres graphiques Suncore en Vogle et PS

c-----------------------------------------------------------------------
      subroutine setviewport2(x1,x2,y1,y2)
      common/ps/iseg,nseg(1000),xa,xb,ya,yb,xinf,xsup,yinf,ysup
c     le viewport en PS et Suncore va de 0 a 1, en Vogle de -1 a 1
      sc=2./amin1(x2-x1,y2-y1)
      call viewport(sc*x1-1.,sc*x2-1.,sc*y1-1.,sc*y2-1.)
      xa=x1
      xb=x2
      ya=y1
      yb=y2
      return
      end
c-----------------------------------------------------------------------
      subroutine setviewport3(x1,x2,y1,y2,z1,z2)
c     a le meme effet que setviewport2
      call setviewport2(x1,x2,y1,y2)
      return
      end
c-----------------------------------------------------------------------
      subroutine setwindow(xi,xs,yi,ys)
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      call ortho2(xi,xs,yi,ys)
      xinf=xi
      xsup=xs
      yinf=yi
      ysup=ys
      return
      end
c-----------------------------------------------------------------------
      subroutine setwindow3(xi,xs,yi,ys,zi,zs)
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      call ortho(xi,xs,yi,ys,zi,zs)
      xinf=xi
      xsup=xs
      yinf=yi
      ysup=ys
      return
      end
c-----------------------------------------------------------------------
      subroutine createretainseg(n)
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      common/line/ilin,lx,ly

      if(n.le.0) return
      iseg=iseg+1
      nseg(iseg)=-n
      write(66,'(a,i3.3,a,$)') '/seg',n,' { '
cVogl l'ordre equivalent en vogle serait le suivant, mais il est inutile
	call makeobj(n)
c     repart du point precedent pour eviter les pertes de memoire PS
      write(66,*) 'newpath ',lx,ly,' m'
      ilin=0
      return
      end
c-----------------------------------------------------------------------
      subroutine closeretainseg(nn)
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
c     nseg < 0 ==> segment ouvert
c     nseg > 0 ==> segment deja clos
c     nn=0  ==> clore le segment en cours
c     nseg = -9999  ==> segment deja clos et deleted

      if(iseg.eq.0) return
      n=nn
c     fermeture du segment courant: retrouve son numero
      if(n.eq.0) then
        if(nseg(iseg).ge.0.or.nseg(iseg).eq.-9999) return
        n=-nseg(iseg)
        endif
c     fermeture du segment n
      do 1 i=1,iseg
      if(n.ne.-nseg(i)) goto 1
        write(66,'(a,i3.3)') ' stroke } def seg', n
cVogl l'ordre equivalent en vogle serait le suivant, mais il ne fait que
cVogl retarder l'affichage de l'objet, qui intervient apres sa fermeture
	call closeobj
	call callobj(n)
      write(66,'(a)') 's'
      nseg(i)=n
      return
    1 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine delretainsegment(n)
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      common/rgb/icol,ncol,red(256),gre(256),blu(256)

      if(n.le.0) return
c     effacement du segment n (apres fermeture si encore ouvert)
      do 1 i=1,iseg
      if(nseg(i).eq.-9999.or.n.ne.abs(nseg(i))) goto 1
      if(n.eq.-nseg(i)) call closeretainseg(n)
        write(66,'(a,i3.3,a)') '1 setgray seg',n,' 0 setgray'
cVogl l'ordre equivalent en vogle serait le suivant, mais il ne provoque
cVogl nullement l'effacement de l'objet de l'ecran -> sans interet aucun
	call color(ncol-1)
	call callobj(n)
	call delobj(n)
	call vflush
	call color(0)
CJJ      call vflush()
      nseg(i)=-9999
      return
    1 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine delallretainsegs()
c     en Postscript et Vogle, equivaut a un NewFrame (donc redondant)
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      common/rgb/icol,ncol,red(256),gre(256),blu(256)

      if(iseg.lt.1) return
      do i=1,iseg
      if(nseg(i).ne.-9999) call delretainsegment(abs(nseg(i)))
      enddo
c     nettoie ecran et recharge font et parametres courants
      write(66,*) 'showpage gr gs newpath'
      call setfont(-1)
      call color(ncol-1)
      call clear
      call color(0)
      iseg=0
      return
      end
c-----------------------------------------------------------------------
      subroutine createtempseg()
c     le segment 1000 represente le segment temporaire
c     segments temporaires sans interet en Vogle, a part le flush
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      common/line/ilin,lx,ly
      nseg(1000)=-1
c     repart du point precedent pour eviter les pertes de memoire PS
      write(66,*) 'newpath ',lx,ly,' m'
      ilin=0
      return
      end
c-----------------------------------------------------------------------
      subroutine closetempseg()
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
c     nseg < 0 ==> segment ouvert
c     nseg > 0 ==> segment deja clos
c     nseg = -9999  ==> segment deja clos et deleted

      if(nseg(1000).ne.-1) return
CJJ      call vflush()
      write(66,*) 's'
      nseg(1000)=1
      return
      end
c-----------------------------------------------------------------------
      subroutine newframe()
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      common/rgb/icol,ncol,red(256),gre(256),blu(256)

      if(nseg(1000).eq.-1) call closetempseg()
c     nettoie ecran et recharge font et parametres courants
      write(66,*) 'showpage gr gs newpath'
      call setfont(-1)
      call color(ncol-1)
      call clear
      call color(0)
      return
      end
c-----------------------------------------------------------------------
      subroutine lineabs2(x,y)
      common/line/ilin,lx,ly

CJJ      call setlineindex(-1)
      call draw2(x,y)
      write(66,*) iconvx(x),iconvy(y),' l'
      lx=iconvx(x)
      ly=iconvy(y)
      call Testmaxlin
      return
      end
c-----------------------------------------------------------------------
      subroutine linerel2(dx,dy)
      common/line/ilin,lx,ly

CJJ      call setlineindex(-1)
      call rdraw2(dx,dy)
      write(66,*) iconrx(dx),iconry(dy),' rl'
      lx=lx+iconrx(dx)
      ly=ly+iconry(dy)
      call Testmaxlin
      return
      end
c-----------------------------------------------------------------------
      subroutine moveabs2(x,y)
      common/line/ilin,lx,ly

      write(66,*) iconvx(x),iconvy(y),' m'
      call move2(x,y)
      lx=iconvx(x)
      ly=iconvy(y)
      call Testmaxlin
      return
      end
c-----------------------------------------------------------------------
      subroutine moverel2(dx,dy)
      common/line/ilin,lx,ly

      write(66,*) iconrx(dx),iconry(dy),' rm'
      call rmove2(dx,dy)
      lx=lx+iconrx(dx)
      ly=ly+iconry(dy)
      call Testmaxlin
      return
      end
c-----------------------------------------------------------------------
      subroutine punkt(x,y)
c     subroutine maison (non Suncore) pour tracer un point

      write(66,*) iconvx(x),iconvy(y),' p'
      call point2(x,y)
      return
      end
c-----------------------------------------------------------------------
      subroutine rectangle(xi,xs,yi,ys,ifill)
c     subroutine maison (non Suncore) pour tracer un rectangle
c     ifill=0: tour seul, 1: remplissage seul, 2: remplissage+tour
c     point courant non modifie
      real xi,xs,yi,ys

      if(ifill.eq.0) then
       call setlineindex(-1)
       call polyfill(0)
       call rect(xi,yi,xs,ys)
       write(66,'(4i5,$)') iconvx(xi),iconvx(xs),iconvy(yi),iconvy(ys)
       write(66,'(a)') ' r s newpath'
       endif
      if(ifill.eq.1) then
       call setfillindex(-1)
       call rect(xi,yi,xs,ys)
       write(66,'(4i5,$)') iconvx(xi),iconvx(xs),iconvy(yi),iconvy(ys)
       write(66,'(a)') ' r f s newpath'
       endif
      if(ifill.eq.2) then
       call setfillindex(-1)
       call rect(xi,yi,xs,ys)
       write(66,'(4i5,$)') iconvx(xi),iconvx(xs),iconvy(yi),iconvy(ys)
       write(66,'(a)') ' r gs f s newpath gr'
       call setlineindex(-1)
       call polyfill(0)
       call rect(xi,yi,xs,ys)
       write(66,'(a)') 's newpath'
       endif
cVogl Permet affichage synchronise. A enlever si objet ouvert (makeobj)
cVogl call vflush()

      return
      end
c-----------------------------------------------------------------------
      subroutine cercle(x0,y0,r,ifill)
c     subroutine maison (non Suncore) pour tracer un cercle
c     ifill=0: tour seul, 1: remplissage seul, 2: remplissage+tour
c     point courant non modifie
      parameter(NARC=50)
      real x,y,x0,y0,r

      call circleprecision(NARC)
      x1=iconvx(x0+r)
      y1=iconvy(y0)
c     retrace le cercle, car la fonction arc en PS donnerait une ellipse
      if(ifill.eq.0) then
       call setlineindex(-1)
       call polyfill(0)
       call circle(x0,y0,r)
       write(66,*) 's newpath ',x1,y1,' m'
       do i=1,NARC
       x=x0+r*cos(float(i)*6.2832/float(NARC))
       y=y0+r*sin(float(i)*6.2832/float(NARC))
       write(66,*) iconvx(x),iconvy(y),' l'
       enddo
       write(66,'(a)') 'cp s newpath'
       endif
      if(ifill.eq.1) then
       call setfillindex(-1)
       call circle(x0,y0,r)
       write(66,*) 's newpath ',x1,y1,' m'
       do i=1,NARC
       x=x0+r*cos(float(i)*6.2832/float(NARC))
       y=y0+r*sin(float(i)*6.2832/float(NARC))
       write(66,*) iconvx(x),iconvy(y),' l'
       enddo
       write(66,'(a)') 'cp f s newpath'
       endif
      if(ifill.eq.2) then
       call setfillindex(-1)
       call circle(x0,y0,r)
       write(66,*) 's newpath ',x1,y1,' m'
       do i=1,NARC
       x=x0+r*cos(float(i)*6.2832/float(NARC))
       y=y0+r*sin(float(i)*6.2832/float(NARC))
       write(66,*) iconvx(x),iconvy(y),' l'
       enddo
       write(66,'(a)') 'cp gs f s newpath gr'
       call setlineindex(-1)
       call polyfill(0)
       call circle(x0,y0,r)
       write(66,'(a)') 's newpath'
       endif
cVogl Permet affichage synchronise. A enlever si objet ouvert (makeobj)
CJJ      call vflush()

      return
      end
c-----------------------------------------------------------------------
      subroutine polylineabs2(x,y,n)
c     courbe f(x,y) demarrant au point courant
      real x(*),y(*)
      common/line/ilin,lx,ly

c     retrouve couleur de trait et position courantes
CJJ      call setlineindex(-1)
      write(66,*) 's newpath ',lx,ly,' m'

      do i=1,n
      write(66,*) iconvx(x(i)),iconvy(y(i)),' l'
      call Testmaxlin
      call draw2(x(i),y(i))
      enddo
      lx=iconvx(x(n))
      ly=iconvy(y(n))
      write(66,*) 's newpath ',lx,ly,' m'
cVogl Permet affichage synchronise. A enlever si objet ouvert (makeobj)
CJJ      call vflush()
      ilin=0
      return
      end
c-----------------------------------------------------------------------
      subroutine polylinerel2(dx,dy,n)
c     courbe f(dx,dy) demarrant au point courant
      real dx(*),dy(*)
      common/line/ilin,lx,ly

c     retrouve couleur de trait et position courantes
CJJ      call setlineindex(-1)
      write(66,*) 's newpath ',lx,ly,' m'

      do i=1,n
      write(66,*) iconrx(dx(i)),iconry(dy(i)),' rl'
      call Testmaxlin
      call rdraw2(dx(i),dy(i))
      lx=lx+iconrx(dx(i))
      ly=ly+iconry(dy(i))
      enddo
      write(66,*) 's newpath ',lx,ly,' m'
cVogl Permet affichage synchronise. A enlever si objet ouvert (makeobj)
CJJ      call vflush()
      ilin=0
      return
      end
c-----------------------------------------------------------------------
      subroutine polygonabs2(x,y,n)
c     Polygone plein. Augmenter la dimension de points si n > 1000
c     Remplit d'abord, dessine le contour (continu en Vogle) ensuite.
      real x(*),y(*),points(2,128)
      common/line/ilin,lx,ly

      if(n.gt.128) print*,'% ! Vogle cannot fill polygons with more than
     & 128 points -> truncation to 128 !'
      m=min(n,128)
      do i=1,m
      points(1,i)=x(i)
      points(2,i)=y(i)
      enddo
c     retrouve couleur de fill courante
      call setfillindex(-1)
      call poly2(m,points)
      write(66,*) 's newpath ',iconvx(x(1)),iconvy(y(1)),' m'
c     dessin sans stroke intermediaire pour fill, car dim.max 128
      do i=2,m
      write(66,*) iconvx(x(i)),iconvy(y(i)),' l'
      enddo
c     trace du tour, en trait continu comme Suncore
      write(66,*) 'cp gs f s newpath gr [] 0 setdash'
      call setlineindex(-1)
      call polyfill(0)
      call poly2(m,points)
      lx=iconvx(x(1))
      ly=iconvy(y(1))
      call setlinestyle(-2)
cVogl Permet affichage synchronise. A enlever si objet ouvert (makeobj)
CJJ      call vflush()
      return
      end
c-----------------------------------------------------------------------
      subroutine polygonrel2(dx,dy,n)
c     Polygone plein. Augmenter la dimension de points si n > 1000
c     Remplit d'abord, dessine le contour (continu en Vogle) ensuite.
      real dx(*),dy(*),points(2,128)
      common/line/ilin,lx,ly

      if(n.gt.128) print*,'% ! Vogle cannot fill polygons with more than
     & 128 points -> truncation to 128 !'
      m=min(n,128)
      call getgp2(x0,y0)
      points(1,1)=x0
      points(2,1)=y0
      do i=2,m
      points(1,i)=points(1,i-1)+dx(i)
      points(2,i)=points(2,i-1)+dy(i)
      enddo
c     retrouve couleur de fill et position courantes
      call setfillindex(-1)
      call poly2(m,points)
      write(66,*) 's newpath ',lx,ly,' m'
c     dessin sans stroke intermediaire pour fill, car dim.max 128
      do i=1,m
      write(66,*) iconrx(dx(i)),iconry(dy(i)),' rl'
      enddo
c     trace du tour, en trait continu comme Suncore
      write(66,*) 'cp gs f s newpath gr [] 0 setdash'
      call setlineindex(-1)
      call polyfill(0)
      call poly2(m,points)
      call setlinestyle(-2)
cVogl Permet affichage synchronise. A enlever si objet ouvert (makeobj)
CJJ      call vflush()
      return
      end
c-----------------------------------------------------------------------
      subroutine Testmaxlin
c     inclut un stroke toutes les 500 commandes pour eviter surcharge PS
      parameter(MAXLIN=500)
      common/line/ilin,lx,ly
      ilin=ilin+1
      if(ilin.gt.MAXLIN) then
        write(66,*) 'gs s newpath gr'
        ilin=0
        endif
      return
      end
c-----------------------------------------------------------------------
      subroutine setcharprecision(n)
c     STRING=0 CHARACTER=1
      common/txt/ichpr,ijust,ifont,chx,chy,angle,modif
      if(n+1.ne.n0) modif=1
      ichpr=n
      n0=n+1
      return
      end
c-----------------------------------------------------------------------
      subroutine setfont(ni)
c     recherche les polices les plus ressemblantes suivant charprecision
c     STRING: times ou futura etroit, taille fixe
c     CHARACTER: 0-5 = roman, greek, script, oldenglish, stick, symbols
c     si n<0, recharge la font precedente
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      common/txt/ichpr,ijust,ifont,chx,chy,angle,modif
      save n0

      write(*,'(a)')'SETFONT'
      n=ni
      modif=2
      write(66,*)'N0',n0
      if(n.lt.0) n=n0
      ifont=n
      if(ichpr.eq.1) then
       if(n.eq.0) call font("times.r")
       if(n.eq.0) write(66,'(a,$)')'/Bookman-Light findfont'
       if(n.eq.1) call font("greek")
       if(n.eq.1) write(66,'(a,$)')'/Symbol findfont'
       if(n.eq.2) call font("cursive")
       if(n.eq.2) write(66,'(a,$)')'/ZapfChancery-MediumItalic 
     & findfont 1.2 scalefont'
       if(n.eq.3) call font("gothic.eng")
       if(n.eq.3) write(66,'(a,$)')'/Bookman-DemiItalic findfont'
       if(n.eq.4) call font("futura.l")
       if(n.eq.4) write(66,'(a,$)')'/Helvetica findfont'
       if(n.eq.5) call font("symbolic")
       if(n.eq.5) write(66,'(a,$)')'/ZapfDingbats findfont'
      else
       if(n.eq.0) then
        call font("times.rb")
        xmsize=0.018*(xsup-xinf)
        ymsize=0.019*(ysup-yinf)
        write(66,'(a,$)')'/NewCenturySchlbk-Bold findfont  60 scalefont'
        endif
       if(n.eq.1) then
        call font("futura.l")
        xmsize=0.013*(xsup-xinf)
        ymsize=0.013*(ysup-yinf)
        write(66,'(a,$)')'/Helvetica findfont 75 scalefont'
        endif
       if(n.eq.2) then
        call font("futura.l")
        xmsize=0.010*(xsup-xinf)
        ymsize=0.009*(ysup-yinf)
        write(66,'(a,$)')'/Helvetica findfont 50 scalefont'
        endif
       if(n.eq.3) then
        call font("futura.m")
        xmsize=0.013*(xsup-xinf)
        ymsize=0.013*(ysup-yinf)
        write(66,'(a,$)')'/Helvetica-Bold findfont 85 scalefont'
        endif
       if(n.eq.4) then
        call font("times.r")
        xmsize=0.013*(xsup-xinf)
        ymsize=0.014*(ysup-yinf)
        write(66,'(a,$)')'/Palatino-Roman findfont 90 scalefont'
        endif
       if(n.eq.5) then
        call font("times.rb")
        xmsize=0.015*(xsup-xinf)
        ymsize=0.014*(ysup-yinf)
        write(66,'(a,$)')'/Palatino-Bold findfont 90 scalefont'
        endif
      call textsize(xmsize,ymsize)
      endif
      write(66,'(a)')' setfont s newpath'
      n0=n
      return
      end
c-----------------------------------------------------------------------
      subroutine setcharsize(hx,hy)
c     Sans effet si charprecision=STRING.
c     Fonts Suncore environ 2 fois plus grandes que celles de Vogle.
      common/txt/ichpr,ijust,ifont,chx,chy,angle,modif
      if(ichpr.eq.0) return
      if(hx.ne.hx0.or.hy.ne.hy0) modif=1
      call textsize(2.*hx,2.*hy)
      chx=2.*float(iconrx(hx))
      chy=2.*float(iconry(hy))
      if(chx.eq.0.) chx=2.
      if(chy.eq.0.) chy=2.
      hx0=hx
      hy0=hy
      return
      end
c-----------------------------------------------------------------------
      subroutine setcharjust(n)
c     off=0 left=1 center=2 right=3
c     Cette fonction n'etant pas operationnelle en Suncore, elle est
c     programmee a l'affichage du texte en fonction de son extension.
c     En Vogle, la justification verticale equivalente est centree.
      common/txt/ichpr,ijust,ifont,chx,chy,angle,modif
      if(n.eq.0) call centertext(off)
      if(n.eq.1) call leftjustify()
      if(n.eq.2) call xcentertext()
      if(n.eq.3) call rightjustify()
      call ycentertext()
      ijust=n
      return
      end
c-----------------------------------------------------------------------
      subroutine setcharpath2(dx,dy)
c     Sans effet si charprecision=STRING.
      common/txt/ichpr,ijust,ifont,chx,chy,angle,modif
      if(ichpr.eq.0) return
      if(dx.ne.dx0.or.dy.ne.dy0) modif=1
      if(dx.eq.0.) angle=90.
      if(dx.ne.0.) angle=atan2(dy,dx)*180./3.141592654
      call textang(angle)
      dx0=dx
      dy0=dy
      return
      end
c-----------------------------------------------------------------------
      subroutine setcharup2(dx,dy)
c     Sans effet en Postscript et Vogle:
c     caracteres perpendiculaires a la chaine
      return
      end
c-----------------------------------------------------------------------
      subroutine text(carte)
      character*(*) carte
      character*100 cardeb,carfin,carps
      common/txt/ichpr,ijust,ifont,chx,chy,angle,modif
      common/line/ilin,lx,ly

c     verification des paires de parentheses
      nc=lnblnk(carte)
      npar=0
      ndeb=0
      nfin=0
      do 1 i=1,nc
      if(carte(i:i).eq.'(') npar=npar+1
      if(carte(i:i).eq.')') npar=npar-1
    1 continue
      if(npar.gt.0) then
      nfin=npar
      do 2 i=1,npar
    2 carfin(i:i)=')'
      endif
      if(npar.lt.0) then
      ndeb=-npar
      do 3 i=1,ndeb
    3 cardeb(i:i)='('
      endif
      carps=' ('//cardeb(:ndeb)//carte(:nc)//carfin(:nfin)//')'
      nps=lnblnk(carps)

c     retrouve caracteristiques du texte PS si modifiees
      if(modif.eq.1) call setfont(ifont)
      call settextindex(-1)
c     positionnement & justification
      if(ichpr.eq.0) then
        write(66,*) lx,ly-40,' m'
      else
        if(modif.ge.1) write(66,*) lx,ly,' m',0,nint(-chy/2.),angle,
     &    int(chx),int(chy),' txt'
        if(modif.eq.0.and.ijust.ge.2) write(66,*) lx,ly,' m'
        if(ijust.eq.2) write(66,'(2a)') carps(:nps),' centr'
        if(ijust.eq.3) write(66,'(2a)') carps(:nps),' right'
      endif
      write(66,'(2a)') carps(:nps),' sh s newpath'

c     parametres de trait affectent texte Vogle: necessite reinitialiser
      call linestyle('1')
c      call linewidth(0)
      call drawstr(carte)
      call setlinestyle(-1)
      call setlinewidth(-1.)
cVogl Permet affichage synchronise. A enlever si objet ouvert (makeobj)
CJJ      call vflush()
      modif=0
      return
      end
c-----------------------------------------------------------------------
      subroutine setmarkersymbol(nn)
c     les marqueurs sont les codes Ascii n de 32 a 127, traces par texte
c     NB: il existe une font "markers" Vogle, mais differents de Suncore.
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      character*1 asci(96),mrq
      common/marqueur/mrq
      data asci/' ','!','#','#','$','%','&',' ','[',']','*','+',' ',
     &'-','.','/','0','1','2','3','4','5','6','7','8','9',':',';','<',
     &'=','>','?','@','A','B','C','D','E','F','G','H','I','J','K','L',
     &'M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','[','/',
     &']','^','_',' ','a','b','c','d','e','f','g','h','i','j','k','l',
     &'m','n','o','p','q','r','s','t','u','v','w','x','y','z','{','|',
     &'}','~','^'/
      n=nn
      if(n.lt.32.or.n.gt.127) n=32
      call font("futura.m")
      xmsize=0.015*(xsup-xinf)
      ymsize=0.015*(ysup-yinf)
      call textsize(xmsize,ymsize)
      call centertext(on)
      write(66,*) 's newpath /Helvetica findfont 75 scalefont setfont'
      mrq=asci(n-31)
      return
      end
c-----------------------------------------------------------------------
      subroutine markerabs2(x,y)
c     dessin du marqueur defini avant par mrq , taille fixe 75 points
      character*1 mrq
      common/marqueur/mrq
      write(66,*) iconvx(x)-37,iconvy(y)-37,' m'
      write(66,'(3a)')' (',mrq,') sh'
c     parametres de trait affectent texte Vogle: necessite reinitialiser
      call linestyle('1')
c      call linewidth(0)
      call move2(x,y)
      call drawchar(mrq)
      call setlinestyle(-1)
      call setlinewidth(-1.)
      return
      end
c-----------------------------------------------------------------------
      subroutine setlinestyle(istyle)
c     En Suncore SOLID=0, DOTTED=1, DASHED=2, DOTDASHED=3
c     En Vogle et Postscript, on reproduit des patterns similaires.
c     Si n=-1, on reprend le linestyle precedent en Vogle seulement
c     Si n=-2, on reprend le linestyle precedent en Vogle et Postscript
      common/line/ilin,lx,ly
      n=istyle
      if(n.lt.0) n=n0
      if(istyle.ne.-1) then
      if(n.eq.0) write(66,'(a,$)') 's newpath [] 0 setdash'
      if(n.eq.1) write(66,'(a,$)') 's newpath [20] 0 setdash'
      if(n.eq.2) write(66,'(a,$)') 's newpath [80 50] 0 setdash'
      if(n.eq.3) write(66,'(a,$)') 's newpath [20 40 90 40] 0 setdash'
      write(66,*) lx,ly,' m'
      endif
      if(n.eq.0) call linestyle('')
      if(n.eq.1) call linestyle('100')
      if(n.eq.2) call linestyle('11110000')
      if(n.eq.3) call linestyle('100011111000')
      ilin=0
      n0=n
      return
      end
c-----------------------------------------------------------------------
      subroutine setlinewidth(width)
c     Conversion PS et Vogle equivalentes pour fenetres de tailles 
c     similaires et facteurs respectifs de 50 et 7, adaptables.
c     Si width=-1. recharge la largeur de trait precedente en Vogle seul
c     Si width=-2. recharge la largeur de trait precedente en Vogle & PS
      common/line/ilin,lx,ly
      parameter(psfac=50.,vofac=7.) 
      if(width.ge.0.) rwidth=width
      if(width.lt.0.) rwidth=width0
      if(width.ne.-1.) then
       factor=rwidth*psfac
       write(66,*) 's newpath',factor,' setlinewidth',lx,ly,' m'
       endif
      iw=nint(rwidth*vofac)
c      call linewidth(iw)
      width0=rwidth
      ilin=0
      return
      end
c-----------------------------------------------------------------------
      subroutine setfillindex(index)
c     Si index<0 on recharge la couleur de remplissage precedente
c     on limite les s n ici pour rapidite et permettre fill+contour
      common/curcol/indcur
      if(indx.ne.indcur) write(66,'(a,$)') 's newpath '
      if(index.ge.0) indx=index
      if(index.lt.0) indx=indxf
      call setcolor(indx)
      call polyfill(1)
      indxf=indx
      end
c-----------------------------------------------------------------------
      subroutine settextindex(index)
c     Si index<0 on recharge la couleur de texte precedente
      common/curcol/indcur
      if(index.ge.0) indx=index
      if(index.lt.0) indx=indxt
      if(indx.ne.indcur) write(66,'(a,$)') 's newpath '
      call setcolor(indx)
      indxt=indx
      end
c-----------------------------------------------------------------------
      subroutine setlineindex(index)
c     Si index<0 on recharge la couleur de trait precedente
      if(index.ge.0.and.index.ne.indxl) write(66,'(a,$)') 's newpath '
      if(index.ge.0) indx=index
      if(index.lt.0) indx=indxl
      call setcolor(indx)
      indxl=indx
      end
c-----------------------------------------------------------------------
      subroutine setcolor(index)
c     Routine maison pour definir la couleur ou le gris en PS et Vogle.
c     En Suncore et Vogle, index varie de 0 a ncol-1,
c     en PS les indices des tableaux RGB vont de 1 a ncol.
      common/rgb/icol,ncol,red(256),gre(256),blu(256)
      common/line/ilin,lx,ly
      common/curcol/indcur

      indx=iabs(index)
      if(indx.gt.ncol-1) then
        if(icol.ne.0) write(*,'(2(a,i4),a)')'% ! index ',indx,
     &  ' > ncol (Max.',ncol,' colors) !'
        indx=0
        endif
      call color(indx)

c     si changement de couleur, recharge couleur et position courante
      if(indx.eq.indcur) return
      j=indx+1
      if(icol.le.1) write(66,'(f5.2,a,$)') red(j),' setgray'
      if(icol.ge.2) write(66,'(3f5.2,a,$)') red(j),gre(j),blu(j),' sc'
      write(66,*) lx,ly,' m'
      indcur=indx
      ilin=0

      return
      end
c-----------------------------------------------------------------------
c Modif eric 6/02/2002 pour compilation sous Linux
c (routine jamais appellee par le code)
c     subroutine inqtextextent2(carte,dx,dy)
c     En PostScript, retourne dx,dy dans les variables nommees dx,dy.
c     character*(*) carte
c     call getfontsize(dx1,dy)
c     dx=strlength(carte)
c     write(66,*) '('//carte(:lnblnk(carte))//')',
c    &' stringwidth /dy exch def /dx exch def'
c     return
c     end
c-----------------------------------------------------------------------
      subroutine inqcurrpos2(x,y)
c     En Suncore et Vogle, retourne la position World courante dans x,y.
c     En PostScript, se repositionne au point courant.
      common/line/ilin,lx,ly
      call getgp2(x,y)
      write(66,*) 's newpath ',lx,ly,' m'
      return
      end
c-----------------------------------------------------------------------
      subroutine awtbuttongetloc2(msec,loc,ibut,x,y)
c     cette fonction n'a pas d'effet en PostScript
c   ! En Suncore renvoie des coordonnees NDC, en Vogle des coord. World.
c     En sortie, ibut est le numero du bouton presse, x,y la position du
c     locator (souris). En Vogle, pas de prise en compte du temps msec:
c     on le remplace par un comptage +/- equivalent de 0 a msec/2.
c     Numero de locator loc egalement non pris en compte.
      integer locator
      do j=0,nint(msec/2.)
      ibut=locator(x,y)
      if(ibut.ne.0) goto 1
      enddo
c     en Vogle  : gauche, milieu, droit = 1, 2, 4 + combinaisons
c     en Suncore: gauche, milieu, droit = 1, 2, 3 seulement
c     renvoie 99 si plusieurs boutons appuyes simultanement
   1  if(ibut.eq.3.or.ibut.ge.5) ibut=99
      if(ibut.eq.4) ibut=3
      return
      end
c-----------------------------------------------------------------------
      subroutine mapndctoworld2(x,y,xw,yw)
c     cette fonction n'a pas d'effet en PostScript
c   ! On ne fait pas de conversion ici, puisque cette fonction est en
c     general utilisee apres un AwtButtonGetLoc2 pour convertir les
c     coordonnees Suncore NDC en World. Vogle a renvoye directement
c     les coordonnees World a l'appel de awtbuttongetloc2.
      xw=x
      yw=y
      return
      end
c-----------------------------------------------------------------------
      subroutine awaitanybutton(msec,ibut)
c     cette fonction n'a pas d'effet en PostScript
c     En Vogle, meme effet que awtbuttongetloc2(msec,loc,ibut,x,y),
c     sauf que les coordonnees x,y ne sont pas retournees.
      call awtbuttongetloc2(msec,loc,ibut,x,y)
      return
      end
c-----------------------------------------------------------------------
      function iconvx(x)
      common/psdim/PSXDIM,PSYDIM
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      iconvx=(x-xinf)*PSXDIM*(x2-x1)/(xsup-xinf)+PSXDIM*x1
      return
      end
c-----------------------------------------------------------------------
      function iconvy(y)
      common/psdim/PSXDIM,PSYDIM
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      iconvy=(y-yinf)*PSYDIM*(y2-y1)/(ysup-yinf)+PSYDIM*y1
      return
      end
c-----------------------------------------------------------------------
      function iconrx(x)
      common/psdim/PSXDIM,PSYDIM
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      iconrx=x*PSXDIM*(x2-x1)/(xsup-xinf)
      return
      end
c-----------------------------------------------------------------------
      function iconry(y)
      common/psdim/PSXDIM,PSYDIM
      common/ps/iseg,nseg(1000),x1,x2,y1,y2,xinf,xsup,yinf,ysup
      iconry=y*PSYDIM*(y2-y1)/(ysup-yinf)
      return
      end
c-----------------------------------------------------------------------
c     function lnblnk (string)
c     fonction rajoutee car non standard Fortran 77 (inconnue sur HP)
c     character*(*) string
c      do i=len(string),1,-1
c     if(string(i:i).ne.' ') goto 1
c     enddo
c  1  lnblnk=i
c     return
c     end
