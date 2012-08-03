c derniere modif 20/11/03     
c     Version de des2intomoGLOB.f
c     qui genere un fichier intomodesVS.V3
c     avec un format simplifie(coordonnee de la station
c     donnee pour chaque evenement)
c     Cette version lit les corrdonnees
c     des epicentres et stations dans  un fichier de type nomcoo
c     'sans liste de station au debut'
c------------------------------------------------------
c     Doit aussi permettre de faire un test synthetique de 
c     l'algorithme de regionalisation.
c     Version pour inversion globale (21/10/03)

      character bcd*100,texte1*80,bidc*20
      character texte2*80,coord*60
      character*41 entrc(120000),car
      character*120 entrend,entrstart
      character*120 coco,adrentr,adrres,entr
      real      late(120000),lone(120000),prof(40)
      real      latsn(120000),lonsn(120000)
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
      write(*,*)'Pas plus de 2000 stations et 120000 trajets!'
      write(*,*)'Pas plus de 40 couches'
      write(*,*)'40 caracteres au plus pour les fic des'
      write(*,*)'-------------------------------------'


c------------------------------------------------------------------
c     Set up boundaries of region - the whole world
c------------------------------------------------------------------
      uplat=90.
      downlat=-90.
      westlon=-180.
      eastlon=+180.
      if(eastlon.lt.0) eastlon=eastlon+360
      if(westlon.gt.eastlon) westlon=westlon-360

c------------------------------------------------------------------
c     Get the input information from the user
c------------------------------------------------------------------
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
      write(*,*)'parametre est Q, autrement cest Xi.'
      read(*,*) iderparaeqQ
      open(so,status='new',file='intomodesVs.V3')
      if (iderparaeqQ.eq.1) then
        open(sox,status='new',file='intomodeslQ.V3')
      else
        open(sox,status='new',file='intomodesXi.V3')
      endif
      write(*,*)'Voulez-vous eliminer les trajets pour lesquels'
      write(*,*)'le modele de facteur de qualite est peu credible?' 
      write(*,*)'(1=oui,0=non)'
      read(*,*) iq 

      write(*,*)'Nom du fichier coord (nomcoo + ntraj)'
      read(*,*)coco

c------------------------------------------------------------------
c     Read the coordinate file
c------------------------------------------------------------------
      open(inc,status='old',file=coco)
      read(inc,*)nbtra
      write(so,*)nbtra
      write(sox,*)nbtra
      do 7 j=1,nbtra
        read(inc,'(a41,a60)') entrc(j),coord
        read(coord,*)late(j),lone(j),latsn(j),lonsn(j)
        if(lone(j).lt.0.)then 
          evlo=lone(j)+360 
        else
          evlo=lone(j)
        endif
        if(evlo.gt.eastlon) evlo=evlo-360
        if(late(j).gt.uplat.or.late(j).lt.downlat
     *     .or.evlo.lt.westlon)then
          write(*,*)'Agrandir la fenetre de regio pour ce trajet:'
          write(*,'(a,f,f)')entrc(j),late(j),evlo
          write(*,*)'WESTLON',westlon,' EASTLON',eastlon
        endif
7     continue
      close(inc)
      write(*,*)'Lecture du fichier coord terminee'
c------------------------------------------------------------------



c------------------------------------------------------------------
c     Initialise sums of Vs and Q for final averages
c------------------------------------------------------------------
      do 9 j=1,40 
        somVs(j)=0.
        somQ(j)=0.
9     continue

c------------------------------------------------------------------
c     START LOOP OVER DES FILES
c------------------------------------------------------------------
      comptra=0
      do 10 compt=1, nbtra
        read(*,'(a)') entr
        write(*,*) entr
c------------------------------------------------------------------
c       if doing synthetic then ....
c------------------------------------------------------------------
        if(itest.eq.1) then
          ncou=2
          ncou2= 2*ncou
        endif
c------------------------------------------------------------------
c       if doing real data then ....
c------------------------------------------------------------------
        if(itest.ne.1) then
c------------------------------------------------------------------
c         find the date parte of the filename
c------------------------------------------------------------------
          i2=lnblnk(entr) 
          i1=index(entr,'/19')
          ir1=index(entr,'/20')
          if(i1.eq.0) then
            i1=ir1
            if(ir1.eq.0.)stop ' PROBLEME'
          elseif (ir1.ne.0.and.ir1.lt.i1) then
            i1=ir1
          endif
          entrend=entr(i1+1:i2)
          entrstart=entr(:i1)
          entrend=entrend(:lnblnk(entrend))
c------------------------------------------------------------------
c         If we are throwing out paths on the basis of Q models,
c         we have to read the Qmodels
c------------------------------------------------------------------
          if (iderparaeqQ.eq.1.and.iq.eq.1) then
            adrres=entr(:lnblnk(entr))//"/res."//entrend  
            write(*,'(a)')adrres 
            open(in2,status='old',file=adrres)
            read(in2,*) bidc,bidc,fact
            if((fact.ge.2.).or.(fact.lt.0.5)) goto 199
            close(in2)
          endif
c------------------------------------------------------------------
c         Read appropriate des file
c------------------------------------------------------------------
          adrentr=entrstart(:lnblnk(entrstart))//'/des.'//entrend

          open(in,status='old',file=adrentr)
          read(in,'(a)') bcd  
          read(in,*)bid 
          read(in,*) ncou,bido
          ncou2= 2*ncou

c index cherche la premiere occurence d'une sous-chaine de caracteres
c dans une chaine donnee. Sous chaine pas trouvee=>resultat=0
c sous chaine trouvee=>resultat=rang du 1er caractere de la sous-chaine
c trouve.
  112    read(in,'(a)',end=199) texte1
         if(index(texte1,'Modele inverse').eq.0) goto 112
         read(in,900) (p(j), j=1,10)
         read(in,900) (p(j), j=11,20)
         read(in,900) (p(j), j=21,30)
         read(in,900) (p(j), j=31,40)
         read(in,900) (p(j), j=41,50)
         read(in,900) (p(j), j=51,60)
         read(in,950) (p(j), j=61,68)

  114    read(in,'(a)',end=199) texte2
         if(index(texte2,'Erreurs-modele').eq.0) goto 114
         read(in,900) (er(j), j=1,10)
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

c-Ecriture dans le fichier de sortie
      
      comptra=comptra+1
      if (comptra.lt.10) then                           
        write(so,1000) "z0",comptra,entrend
        write(sox,1000) "z0",comptra,entrend
      elseif (comptra.lt.100) then                           
        write(so,1001) "z",comptra,entrend
        write(sox,1001) "z",comptra,entrend
      elseif (comptra.lt.1000) then
        write(so,1002) "z",comptra,entrend
        write(sox,1002) "z",comptra,entrend
      elseif (comptra.lt.10000) then
        write(so,1003) "z",comptra,entrend
        write(sox,1003) "z",comptra,entrend
      elseif (comptra.lt.100000) then
        write(so,1004) "z",comptra,entrend
        write(sox,1004) "z",comptra,entrend
      elseif (comptra.lt.1000000) then
        write(so,1005) "z",comptra,entrend
        write(sox,1005) "z",comptra,entrend
      endif
1000  format(a,i1,5x,a)
1001  format(a,i2,5x,a)
1002  format(a,i3,5x,a)
1003  format(a,i4,5x,a)
1004  format(a,i5,5x,a)
1005  format(a,i6,5x,a)

      write(so,*)late(compt),lone(compt),latsn(compt),lonsn(compt)
      write(sox,*)late(compt),lone(compt),latsn(compt),lonsn(compt)

      write(so,2000) (Vs(j),j=1,ncou)
      write(so,2000) (erVs(j),j=1,ncou)
      write(sox,2010) (Q(j),j=1,ncou)
      write(sox,2010) (erQ(j),j=1,ncou)
2000  format(34f8.4)
2010  format(34f14.4)
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
