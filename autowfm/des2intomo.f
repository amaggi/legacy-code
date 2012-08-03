c derniere modif 24 Juin 2000
c     Version amelioree de readdesaus2.
c     Doit aussi permettre de faire un test synthetique de 
c     l'algorithme de regionalisation.

      character bcd*100,sta*4,sto*4,texte1*80,bidc*20
      character*4 stat(400)
      character texte2*80
      character*254 entrc(60000),car
      character*254 entrend
      character*254 coco,adrentr,adrres,entr
      real      late(60000),lone(60000),prof(40), lats, lons
      dimension p(80),er(80)
      dimension somVs(40),somQ(40),Vs(40),Q(40),erQ(40),erVs(40)
      integer compt,comptra,nsta,nbtra
      real uplat,downlat,eastlon,westlon,evlo
       
      yerr=0
      nmo=5
      in=10
      inc=11
      in2=12
      so=15
      sox=16
      write(*,*)'-------------------------------------'
      write(*,*)'ATTENTION aux dimensions de tableaux:'
      write(*,*)'Pas plus de 400 stations et 60000 trajets!'
      write(*,*)'Pas plus de 40 couches'
      write(*,*)'253 caracteres au plus pour les fic des'
      write(*,*)'-------------------------------------'


      write (*,*)'Enter the coordinates that delineate the region'
      write (*,*)'under study (uplat,downlat,westlon,eastlon):'
      write (*,*)'For Australia it was 20 -66 89 189 '
      write (*,*)'For Africa it is 45 -30 0 110 '
      read(*,*)uplat,downlat,westlon,eastlon
      if(eastlon.lt.0) eastlon=eastlon+360
      if(westlon.gt.eastlon) westlon=westlon-360

c     lecture dans le fichier lstcoord des coordonnees
c     des epicentres et des stations.
      write(*,*)'Voulez vous  creer des fichiers intomodes'
      write(*,*)'synthetiques ou toutes les donnees ont la meme' 
      write(*,*)'vitesse?(1=oui,0=non)' 
      read(*,*)itest
      if (itest.eq.1)then
         write(*,*)'Entrez cette vitesse :' 
         read(*,*)ctest
      endif
      write(*,*)'Deux parametres sont traite par ce pgm:'
      write(*,*)'Vs et Xi ou Vs et Q. Entrez 1 si le 2eme'
      write(*,*)'parametre est Q.'
      read(*,*) iderparaeqQ
      open(so,status='new',file='intomodesVs')
      if (iderparaeqQ.eq.1) then
        open(sox,status='new',file='intomodeslQ')
      else
        open(sox,status='new',file='intomodesXi')
      endif
      write(*,*)'Voulez-vous eliminer les trajets pour lesquels'
      write(*,*)'le modele de facteur de qualite est peu credible?' 
      write(*,*)'(1=oui,0=non)'
      read(*,*) iq 
      write(*,*)'Nom du fichier coord'
      read(*,*)coco
      open(inc,status='old',file=coco)
      read (inc,*) nsta,nbtra
      write(so,*)nsta,nbtra           
      write(sox,*)nsta,nbtra           
      do 5 i=1,nsta
        read(inc,2003)stat(i),lats,lons
c       write(*,*)i, nsta, stat(i),lats,lons
        write(so,2003)stat(i),lats,lons
        write(sox,2003)stat(i),lats,lons
5     continue
2003    format(a4,2x,2f9.4)
        write(*,2003)stat(nsta),lats,lons
      do 7 j=1,nbtra
        read(inc,'(a)')entrc(j)
        read(inc,*)late(j),lone(j)
      if(lone(j).lt.0.)then 
c          write(*,*)'LONE(j) EVLO',lone(j),evlo 
           evlo=lone(j)+360 
c          write(*,*)'EVLO ',evlo
      else
           evlo=lone(j)
      endif
c petit beug trouve le 6 aout 2001 coorige
c     if(evlo.gt.eastlon) evlo=lone(j)-360
      if(evlo.gt.eastlon) evlo=evlo-360
c     write(*,*)uplat, downlat, westlon, eastlon
c     write(*,*)late(j), evlo
      if(late(j).gt.uplat.or.late(j).lt.downlat
     *.or.evlo.lt.westlon)then
c    *.or.evlo.lt.westlon.or.evlo.gt.eastlon)then
          write(*,*)'Agrandir la fenetre de regio pour ce trajet:'
          write(*,'(a,f,f)')entrc(j),late(j),evlo
          write(*,*)'WESTLON',westlon,' EASTLON',eastlon
      endif
7     continue
      close(inc)
      write(*,*)'Lecture du fichier coord terminee'
      do 9 j=1,40 
      somVs(j)=0.
      somQ(j)=0.
9     continue
      comptra=0
      do 10 compt=1, nbtra
c         write(*,*)'entrez le nom du fichier des:'
          read(*,'(a)') entr
      if(itest.eq.1) then
c         write(*,*)'On prend 1 couche pour ce test synthe'
          ncou=2
          ncou2= 2*ncou
      endif
      if(itest.ne.1) then
c	comment out the whole of this mes - only want to read the des file
c        i2=lnblnk(entr) 
c        i1=index(entr,'/19')
c        ir1=index(entr,'/20')
c        if(i1.eq.0) then
c           i1=ir1
c           if(ir1.eq.0.)stop ' PROBLEME'
c        elseif (ir1.ne.0.and.ir1.lt.i1) then
c          i1=ir1
c        endif

c        entrend=entr(i1+1:i2)
c        entrend=entrend(:lnblnk(entrend))
c        write(*,'(a)')'FIN',entrend
c         if (iderparaeqQ.eq.1.and.iq.eq.1) then
c          adrres=entr(:lnblnk(entr))//"/res."//entrend  
c          write(*,'(a)')adrres 
c          open(in2,status='old',file=adrres)
c           read(in2,*) bidc,bidc,fact
c           if((fact.ge.2.).or.(fact.lt.0.5)) goto 199
c          close(in2)
c         endif
c        comptra=comptra+1
c        adrentr=entr(:lnblnk(entr))//'/des.'//entrend
c        write(*,'(a)')'ADRENTR',adrentr

         open(in,status='old',file=entr)
         read(in,'(a)') bcd  
c        write(*,'(a)') bcd  
         read(in,*)bid 
         read(in,*) ncou,bido
c        write(*,*) ncou,bido
         ncou2= 2*ncou

c index cherche la premiere occurence d'une sous-chaine de caracteres
c dans une chaine donnee. Sous chaine pas trouvee=>resultat=0
c sous chaine trouvee=>resultat=rang du 1er caractere de la sous-chaine
c trouve.
  112    read(in,'(a)',end=199) texte1
c        if(texte1(:16).ne.'  Modele inverse') goto 112
         if(index(texte1,'Modele inverse').eq.0) goto 112
         read(in,900) (p(j), j=1,10)
c        write(*,900) (p(j), j=1,10)
         read(in,900) (p(j), j=11,20)
         read(in,900) (p(j), j=21,30)
         read(in,900) (p(j), j=31,40)
         read(in,900) (p(j), j=41,50)
         read(in,900) (p(j), j=51,60)
         read(in,950) (p(j), j=61,68)

  114    read(in,'(a)',end=199) texte2
c        if(texte2(:16).ne.'  Erreurs-modele') goto 114
         if(index(texte2,'Erreurs-modele').eq.0) goto 114
         read(in,900) (er(j), j=1,10)
c        write(*,900) (er(j), j=1,10)
         read(in,900) (er(j), j=11,20)
         read(in,900) (er(j), j=21,30)
         read(in,900) (er(j), j=31,40)
         read(in,900) (er(j), j=41,50)
         read(in,900) (er(j), j=51,60)
         read(in,950) (er(j), j=61,68)
900      format(10f8.4)
950      format(8f8.4)
         close(in)

         k=1
         do 115 j=1,ncou2
         if (j.le.ncou) then
            Vs(j)=p(j)
            erVs(j)=er(j)
         endif
         if (j.gt.ncou.and.j.le.ncou2) then
          if(iderparaeqQ.eq.1) then
c           le pgm Montagner est fait pour regionaliser 1/Q
c           avec en input Q. En wfm on sort log10Q
c           on repasse en Q et en er(Q)(en supposant que la matrice
c           de covariance cov(Q) est diagonale).
            Q(k)=10**p(j)
            erQ(k)=log(10.)*Q(k)*er(j)
            if(iq.eq.1.and.((Q(k).gt.3000).or.(Q(k).lt.20))) goto 199
            k=k+1
          else
c           dans ce cas le parametre Q(k) contient en fait Xi  
            Q(k)=p(j)
            erQ(k)=er(j)               
c           TEST: on va mettre 1/Xi dans Q(k)
c           Q(k)=1/p(j)
c           erQ(k)=er(j)/(p(j)**2)               
            k=k+1
          endif
         endif
115      continue
         k=1
         do 116 j=1,ncou2
         if (j.le.ncou) then
            somVs(j)=somVs(j)+(1/Vs(j))
         endif
         if (j.gt.ncou.and.j.le.ncou2) then
            somQ(k)=somQ(k)+(1/Q(k))
            k=k+1
         endif
116      continue
      else
         k=1
         do 120 j=1,ncou2
         if (j.le.ncou) then
            Vs(j)=ctest
            erVs(j)=0.1  
         endif
         if (j.gt.ncou.and.j.le.ncou2) then
            Q(k)=200.    
            erQ(k)=20.
            k=k+1
         endif
120      continue

       endif
c index cherche la premiere occurence d'une sous-chaine de caracteres
c dans une chaine donnee. Sous chaine pas trouvee=>resultat=0
c sous chaine trouvee=>resultat=rang du 1er caractere de la sous-chaine
c trouve.
      itest2=0
      itrouve=0
      do 121 ii=1,nsta 
      sto=stat(ii)(:lnblnk(stat(ii)))
      itest2=index(entr,sto)
      if (itest2.ne.0) then
           itrouve=itrouve+1
           sta=stat(ii)
      endif 
121   continue
      if (itrouve.ne.1) then
          write(*,*)'ITROUVE',itrouve
          write(*,'(a)')entr
          write(*,'(a)')sto
          stop'Station redondante ou non trouvee'
      endif
      
      comptra=comptra+1
c     write(*,*)'COMPTRA',comptra
      if (comptra.lt.10) then                           
        write(so,1000) sta,".z0",comptra,entr(1:lnblnk(entr))
        write(sox,1000) sta,".z0",comptra,entr(1:lnblnk(entr))
      elseif (comptra.lt.100) then                           
        write(so,1001) sta,".z",comptra,entr(1:lnblnk(entr))
        write(sox,1001) sta,".z",comptra,entr(1:lnblnk(entr))
      elseif (comptra.lt.1000) then
        write(so,1002) sta,".z",comptra,entr(1:lnblnk(entr))
        write(sox,1002) sta,".z",comptra,entr(1:lnblnk(entr))
      elseif (comptra.lt.10000) then
        write(so,1003) sta,".z",comptra,entr(1:lnblnk(entr))
        write(sox,1003) sta,".z",comptra,entr(1:lnblnk(entr))
      elseif (comptra.lt.100000) then
        write(so,1004) sta,".z",comptra,entr(1:lnblnk(entr))
        write(sox,1004) sta,".z",comptra,entr(1:lnblnk(entr))
      endif
1000  format(a,a,i1,5x,a)
1001  format(a,a,i2,5x,a)
1002  format(a,a,i3,5x,a)
1003  format(a,a,i4,5x,a)
1004  format(a,a,i5,5x,a)

c     write(*,*) "*****", entr(:lnblnk(entr)), "*****"
      do 122 j=1,nbtra
       car=entrc(j)
       car=car(:lnblnk(car)-4)
c      write(*,*) car(:lnblnk(car))
       if (index(entr(:lnblnk(entr)),car(:lnblnk(car))).ne.0) then
c       if (car(:lnblnk(car)-4).eq.entrend(:lnblnk(entrend))) then
c         write(*,*)'Found path'           
          write(so,*)late(j),lone(j)
          write(sox,*)late(j),lone(j)
c         have written the path lat and lon - now exit the loop
          goto 1222
       endif
122   continue

1222  write(so,2000) (Vs(j),j=1,ncou)
      write(so,2000) (erVs(j),j=1,ncou)
      write(sox,2010) (Q(j),j=1,ncou)
      write(sox,2010) (erQ(j),j=1,ncou)
2000  format(33f8.4)
2010  format(33f14.4)
199   continue
10    continue
      write(*,*)'Nombre de trajets retenus:',comptra
      write(*,*)'moyenne des vitesse ... '
      do 123 j=1,ncou 
      prof(1)=40
      prof(2)=50
      if(j.gt.2.and.j.lt.29)prof(j)=prof(j-1)+25.
      if(j.ge.29)prof(j)=prof(j-1)+50.
      write(*,2005),prof(j),1/(somVs(j)/comptra)
123   continue
      write(*,*)'et du deuxieme parametre traite (calculee en 1/X)'
       do 124 j=1,ncou 
      prof(1)=40 
      prof(2)=50 
      if(j.gt.2.and.j.lt.29)prof(j)=prof(j-1)+25.
      if(j.ge.29)prof(j)=prof(j-1)+50.
      write(*,2005),prof(j),1/(somQ(j)/comptra)
124   continue
2005  format(f7.2,f9.3,' .05 0. 0. 0.005')
      close(so)
      close(sox)
      end
