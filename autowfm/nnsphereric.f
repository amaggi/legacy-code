c------------------------------------------------------------------------
c
c	A simple interface to NNsphere initialization routines 
c	with all information passed between routines by common blocks.
c
c	This version is for locating a point on the surface of
c	a sphere only. It calculates the convex hull of input nodes
c	and uses 2-D NN routines for calculating the natural neighbour
c	lists.
c
c	Also performs debug checking routines if required
c
c	All points read in Latitude Longitude pairs.
c
c	INPUT PARAMETERS (idebug,iwrite,mode_del)
c	
c	If idebug = 1       then error checking is performed on the 
c		            tetrahedralization and the results are 
c		            written to standard out.
c
c	If iwrite = 1       then triangles vertices and neighbours
c		            are written to standard out.
c
c	If mode_del = 1     then Delaunay is read in from file `del.in'
c	            = 0     then Delaunay is calculated internally.
c
c
c	INPUT FILE: del.in (only read in if mode_del = 1)
c
c		  
c	INPUT FILE: del.in 
c
c	     INPUT               			Explanation 
c            list of triangles faces on convex hull 
c       	  format: (output format of `qhull')
c		  nt					:number of faces
c		  n1,n2,n3				:nodes at vertices
c		  ...
c		  
c
c	Comments:
c
c						M. Sambridge, Feb. 1999.
c  Slightly modified by eric
c------------------------------------------------------------------------
c
c       SUBROUTINE NNSPHERERIC_INIT
        SUBROUTINE NNSPHERERICINIT
     &            (idebug,iwrite,mode_del,latnode,lonnode,nperic)
c  
        include		'nn.param'
c       implicit real*8(a-h,o-z)

        common/nnsetuparrays/points(3,np_max),
     &                       vertices(4,nt_max),
     &                       neighbour(3,nt_max),
     &                       coords(2,np_max),
     &                       centres(2,nt_max),
     &                       nnn(np_max+1),
     &                       nnlist(nwork3d),
     &                       ntlist(nwork3d)

        common/nnparam/np,nt

c       real latnode(nperic),lonnode(nperic)
        real*8 latnode(nperic),lonnode(nperic)
        
c       real*8          alat,alon,alat2,alon2
	real*8		points
        real*8          coords
        real*8          centres
c       real*8          depth,rad
c       real*8		degtorad,radtodeg
	integer		vertices
	integer		neighbour
	integer		nnn
	integer		nnlist
	integer		ntlist

	logical		nnwrite

c					local parameters
        logical         consistent
        logical         debug,writeout,usefindnode,extendtetra
        integer		cyc(4)
        data		cyc/2,3,4,1/
c       data            degtorad/0.017453292/
c       data            radtodeg/57.29578/
c       data            rad/6371.0/
        data            rad/1.d0/
        data            pi/3.141592653589793116/
        data            pihalf/1.570796326794896558/

        degtorad=pi/180.
        radtodeg=180./pi
c
c     write(*,*)'NNSPHERE PARA',idebug,iwrite,mode_del,nperic
      if(mode_del.eq.0)then
        lu_in = 0
      else if(mode_del.eq.1)then
        lu_in = 50
        open(lu_in,file='del.in',status='old')
      end if
 
c     open(51,file='nodes.in',status='old')
c     read(51,*)
c     read(51,*)np
      np=nperic
      nd = 2
 
      if(nd.gt.nd_max)then
         write(*,*)' Too many dimensions in input file'
         stop
      endif
c     else if(np.gt.np_max)then
c        write(*,*)' Too many points in input file'
c        stop
c     end if
c     write(*,*)' reading data...'
 
c                                               read in data rbox format
c     do i=1,np
c	 read(51,*)(coords(j,i),j=1,2)
c     end do
c     close(51)
c
c						convert input points to 
c						Cartesian reference frame
c						Note that this assumes
c						lat NOT colat
c                                               BUT in forcbrut.c it is a colatitude
c                                               and a longitude that have been passed in
c                                               radians to this rsubroutine.    
      depth = 0.
      do i=1,np
c eric : latnode(it's a colat) and lonnode are already in radian in forcbrut 
	  alat = latnode(i)
	  alon = lonnode(i)
          coords(1,i)=90.-(alat*radtodeg)
          coords(2,i)=alon*radtodeg
          alat2=alat
          alon2=alon
c         alat2 = (90.-alat)*degtorad
c         alon2 = alon*degtorad
c         points(1,i) = (rad-depth)*(dsin(alat2)*dcos(alon2))
c         points(2,i) = (rad-depth)*(dsin(alat2)*dsin(alon2))
c         points(3,i) = (rad-depth)*(dcos(alat2))
          points(1,i) = (rad-depth)*dble(sin(alat2)*cos(alon2))
          points(2,i) = (rad-depth)*dble(sin(alat2)*sin(alon2))
          points(3,i) = (rad-depth)*dble(cos(alat2))
      end do
 
      lud = 6
 
c                                       set debug mode on for nn routines
      nnwrite = .false.
      debug = .false.
      writeout = .true.
      writeout = .false.
      usefindnode = .true.
      extendtetra = .false.
      if(idebug.eq.1)debug = .true.
      if(iwrite.eq.1)writeout = .true.
c     if(icall_find_node.eq.1)usefindnode = .true.
c     if(iextend_outside_hull.eq.1)extendtetra = .true.
      if(debug)write(*,*)' Total number of points read in:',np
c
c                                       Perform setup of nn interpolation
c                                       (i.e. calculate Delaunay triangles,
c                                       build neighbour)
 
      if(debug)write(*,*)' calling nn3d_setup'
c     write(*,*)'Input nnsurface',np,nt_max,np_max,nnpn_max,nwork3d,
c    & lu_in,nt,
c    & loc
c     do i=1,np
c      write(*,*)(points(ji,i),ji=1,3)
c     end do
c     call nnsurface_setup
      call nnsurfacesetup
     &     (np,nt_max,np_max,nnpn_max,nwork3d,
     &      points,centres,lu_in,nt,vertices,
     &      neighbour,loc,nnn,nnlist,ntlist,iwrite)
 

      if(mode_del.ne.0)close(lu_in)

      if(debug)write(*,*)' finished nnsurface_setup'
 
      if(debug)write(*,fmt='(/"  Number of points = ",i7,
     &              " Triangles on surface = ",i7/)')np,nt
      
c					abort checking
      if(.not.debug)return
c                                       check calculation of
c                                       neighbour array

c     call check_neighbour
      call checkneighbour
     &     (neighbour,2,vertices,nt,.false.,consistent)

      if(consistent)then
         write(*,*)' neighbour matrix consistent'
         write(*,*)
      end if

c
c					if extendtetra is set then
c					neighbour array has negative
c					entries
c
c                                       write out vertices and neighbours
c
      if(writeout)then

         open(7,file='VerticeNeighbourArrays')
         write(7,*)' Vertices and neighbour arrays'
         do 10 i=1,nt
            write(7,*)i,' :',vertices(1,i),vertices(2,i),
     &                vertices(3,i),' n:',
     &                neighbour(1,i),neighbour(2,i),
     &                neighbour(3,i)
 10      continue
         close(7)
       
         open(9,file='NpointCoords')
         write(9,*)'npoint, coords'
         do i=1,np
            write(9,*)i,' :',points(1,i),points(2,i),
     &                points(3,i),' c:',
     &                coords(1,i),coords(2,i)
         end do
         close(9)

c                                       special write out of nnlist
         open(8,file='listvoro')
         write(8,*)'indice of each voronoi cell, number of 
     & vortices, lat,lon vortices'
c note that ntlist(nnn(i+1)-1)=0 : this is why we stop at nnn(i+1)-2
         do i=1,np
          write(8,300)i,nnn(i+1)-1-nnn(i),
     & (centres(1,ntlist(k)),centres(2,ntlist(k)),k=nnn(i),nnn(i+1)-2)
         end do

  300    format(1x,i4,1x,i4,50(f8.2,1x,f8.2,2x))
c 300    format(1x,i4,':',50(i4,1x))

         close(8)
      end if
c
 100  format(i2,' : ',25(i2,1x))
 101  format(i3,' : ',20(i3,1x))
 102  format(i4,' : ',15(i4,1x))
 103  format(i6,' : ',10(i6,1x))
 200  format(1x,//' Error - maximum number of faces  
     &              on convex hull exceeded    
     &              increase size of parameter nf_max and recompile'/)
c     do 1040 jj=1,10
c     do 1040 ii=1,4
c1040  write(*,*)'vertices',ii,jj,' : ',vertices(ii,jj)
      return
c     end
      END
c------------------------------------------------------------------------
c
c	nnsurface_setup - Performs all setup procedures for 
c			  natural neighbour location on a sphere.
c
c	Input:
c	        np			number of nodes
c		points(3,np)	        array of node co-ordinates
c	        nt_max			maximum number of triangles
c		np_max			maximum number of nodes 
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                       for nd=2, set to ~20 in calling program)
c		nmax			maximum sum of the number of neighbours 
c					per node (set to ~10*np_max in 
c					calling program)
c	        nd			number of dimensions 
c	        lu_in			logical unit of Delaunay input file
c
c	Output:
c	        nt			number of triangles
c	        vertices(4,nt)		array of spherical triangle vertices 
c	        neighbour(3,nt)		array of neighbouring s-triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j.
c		loc			an initial guess tetrahedron for point 
c					location routine `tetloc' used by nn3D
c					(set somewhere near the centre)
c
c	Comments:
c		 In calling program arrays should be dimensioned:
c
c		 real*8		points(3,np_max)
c		 integer	vertices(4,nt_max)
c		 integer	neighbour(3,nt_max)
c
c		 If lu_in > 0 the the Deluanay tessellation is read in from
c		 logical unit `lu_in' instead of being calculated by 
c		 routine qhullf. This can be useful if qhullf fails because
c		 of precision errors and the Deluanay may be determined
c		 externally to this program using a double precision version
c		 or another algorithm, e.g. Fortune's sweepline method.
c
c		 It is assumed that the read in format has one triangle 
c		 line represented as by four numbered from ZERO. 
c
c		 Three other arrays are calculated as a `by product' of
c		 of the routine build_nn. These must be dimensioned
c		 in the calling program in the following way:
c
c		 integer nnn(np_max+1)  : number of neighbours per node
c		 integer nnlist(nmax) : nn list for nodes
c		 integer ntlist(nmax)   : work array
c
c		 Note that array vertices has to be dimensioned (4,nt)
c		 so that correct
c
c		 Calls are made to: qhullf, build_nn.
c		 
c					    M. Sambridge, RSES.
c
c------------------------------------------------------------------------
c
c       Subroutine nnsurface_setup
	Subroutine nnsurfacesetup
     &             (np,nt_max,np_max,nnpn_max,nmax,
     &              points,centres,lu_in,nt,vertices,
     &              neighbour,loc,nnn,nnlist,ntlist,iwrite)

	real*8		points(3,*)
c       real*8          centres(3,*)
c    Modif eric
        real*8          centres(2,*)
	integer		vertices(4,*)
	integer		neighbour(3,*)
	integer		nnn(*)
	integer		nnlist(*)
	integer		ntlist(*)
	logical		nnwrite

	common/nnswitches/nnwrite,lud

        if(lu_in.le.0)then
c                                       calculate Convex hull 
 
           call qhullf(np,3,3,nt_max,1,points,nt,vertices)

        else
c                                       read in Delaunay vertices
 
           nt = 0
           read(lu_in,*)

  1        read(lu_in,*,end=3,err=100)
     &     vertices(1,nt+1),vertices(2,nt+1),vertices(3,nt+1)

           nt = nt + 1
           if(nt.ge.nt_max) go to 200
           go to 1
  3        continue

        end if

        if(nt.ge.nt_max) go to 200
 
c					adjust array vertices to
c					range from nodes 1 to np
	do 5 i=1,nt
           do 5 j=1,3
	   vertices(j,i) = vertices(j,i) + 1
 5      continue

c                                       Build neighbour matrix
c					for nodes on a sphere
c					Also builds nnn and nnlist
c					used by find_node_sph

c       call build_nn
        call buildnn
     &       (np,vertices,nt,np_max,nmax,
     &        neighbour,nnn,nnlist,ntlist)

c                                       calculate latitude (not colat)
c                                       and longitude of circumcentres

c       call ccentres_sphere
        call ccentressphere
     &       (points,vertices,nt,neighbour,centres,iwrite,50)


c					set initial guess for 
c					triangle location procedure
        loc = nt/2

	return

  100   write(*,*)
     &  'Error in nnsurface_setup:',
     &  ' read error in Delaunay input file'
        stop

  200   write(*,*) 'Error in nnsurface_setup:',
     &             ' too many triangles on sphere'
        write(*,*) 'Remedy: increase size of parameter nt_max'
        write(*,*) '        in calling program.'
        stop 

	end
c------------------------------------------------------------------------
c
c       build_nn - Builds neighbour array for Delaunay triangulation in 2-D.
c
c       Input:
c               vertices(4,nt)          array of triangle vertices
c               nt                      number of triangles
c               np_max                  maximum number of nodes
c               nmax                    maximum total number of neighbours 
c                                       per node (should be set to 
c                                       3*nt_max + np_max in calling program)
c
c       Output:
c               neighbour(3,nt)         array of neighbouring triangles
c               nnlist(namx)            array of neighbouring nodes
c               ntlist(namx)            array of triangles attached 
c                                       to each node
c
c       Comments:
c                Assumes input list of vertices in anticlockwise sequence
c                and produces an anticlockwise list of neighbour triangles.
c                The value of neighbour(i,j) is the index of the neighbouring
c                triangle opposite node i in triangle j.
c
c                Three temporary work arrays are used and must be dimensioned
c                in the calling program in the following way:
c
c                integer nnn(np_max+1)  : number of neighbours per node
c                integer nnlist(nmax)   : natural neighbours per node
c                integer ntlist(nmax)   : triangles attached to each node
c                                         ON A SPHERE ONLY !
c
c                The value of nmax should be set to (3*nt_max + np_max)
c                in the calling program.
c
c                No calls to other routines.
c
c                                       M. Sambridge, RSES, April 1994.
c                                       (using ideas by J.Braun)
c
c------------------------------------------------------------------------
c
c       Subroutine build_nn
        Subroutine buildnn
     &             (np,vertices,nt,np_max,nmax,
     &              neighbour,nnn,nnlist,ntlist)
c
        integer         vertices(4,*)
        integer         neighbour(3,*)
        integer         nnn(*)
        integer         nnlist(*)
        integer         ntlist(*)
        logical         nnwrite

        common/nnswitches/nnwrite,lud

        if(nnwrite)write(*,*)' Building neighbour v ...'
c
c                                       initialize neighbour list
        do 5 i = 1,3
           do 5 j = 1,nt
              neighbour(i,j) = 0
 5      continue
c                                       initialize work arrays
        do i = 1,nmax
           nnlist(i) = 0
           ntlist(i) = 0
        end do

        do 7 i = 1,np
           nnn(i) = 0
 7      continue

        do 10 it = 1,nt
           i1 = vertices(1,it)
           i2 = vertices(2,it)
           i3 = vertices(3,it)
           nnn(i1) = nnn(i1) + 1
           nnn(i2) = nnn(i2) + 1
           nnn(i3) = nnn(i3) + 1
 10     continue

c                                       turn nnn into a running sum
        itemp = nnn(1)+1
        nnn(1) = 1
        do 20 j = 2,np+1
           itemp2  = itemp 
           itemp   = itemp + nnn(j)+1
           nnn(j) = itemp2 + 1
 20     continue
c       write(*,*)' size of array =',nnn(np+1)-1
c       write(*,*)' 3nt+np        =',3*nt+np

        if(nnn(np+1).ge.nmax)then
           write(*,*)'Error: array sizes too small in subroutine '
     &               ,'build_nn'
           write(*,*)'       maximum number of neighbours for all nodes'
           write(*,*)'       is too small: current value =',nmax
           write(*,*)'       Increase size of parameter nmax'
           write(*,*)'       to at least',nnn(np+1)
           write(*,*)'       This will be satisfied if nmax is set'
           write(*,*)'       to 3*nt_max+np_max in calling program' 
           stop
        end if

        do 25 it = 1,nt
           i1 = vertices(1,it) 
           i2 = vertices(2,it) 
           i3 = vertices(3,it) 
c                                               compare neighbours i1 i2
c                                               (remove go to ?)
           j1 = nnn(i1)
           j2 = nnn(i1+1) - 1 
           jt = 0
           do 30 j = j1,j2
              if(nnlist(j).eq.0)then
                 nnlist(j) = i2
                 ntlist(j) = it
c                                               if we have recorded connection
c                                               then jump out of loop
                 go to 31
              else if(nnlist(j).eq.i2.and.ntlist(j).ne.it)then
                 jt = ntlist(j)
                 go to 31
              end if
  30       continue
  31       continue
c                                               if neighbours are found then
c                                               skip second loop 
           if(jt.eq.0)then
              j1 = nnn(i2)
              j2 = nnn(i2+1) - 1 
              do 32 j = j1,j2
                 if(nnlist(j).eq.0)then
                    nnlist(j) = i1
                    ntlist(j) = it
c                                               if we have inserted connection
c                                               then jump out of loop
                    go to 33
                 end if
  32          continue
           end if
  33       continue

           if(jt.ne.0)then
c                                               found neighbours it,jt with
c                                               common nodes i1 and i2
              neighbour(3,it) = jt
              k1 = vertices(1,jt)
              k2 = vertices(2,jt)
              k3 = vertices(3,jt)
              if(k1.ne.i1.and.k1.ne.i2)then 
                 neighbour(1,jt) = it
              else if(k2.ne.i1.and.k2.ne.i2)then 
                 neighbour(2,jt) = it
              else
                 neighbour(3,jt) = it
              end if
           end if
c                                               compare neighbours i1 i3
           jt = 0
           j1 = nnn(i1)
           j2 = nnn(i1+1) - 1 
           do 130 j = j1,j2
              if(nnlist(j).eq.0)then
                 nnlist(j) = i3
                 ntlist(j) = it
c                                               if we have recorded connection
c                                               then jump out of loop
                 go to 131
              else if(nnlist(j).eq.i3.and.ntlist(j).ne.it)then
                 jt = ntlist(j)
                 go to 131
              end if
  130      continue
  131      continue
c                                               if neighbours are found then
c                                               skip second loop 
           if(jt.eq.0)then
              j1 = nnn(i3)
              j2 = nnn(i3+1) - 1 
              do 132 j = j1,j2
                 if(nnlist(j).eq.0)then
                    nnlist(j) = i1
                    ntlist(j) = it
c                                               if we have inserted connection
c                                               then jump out of loop
                    go to 133
                 end if
  132         continue
              end if
  133      continue
           if(jt.ne.0)then
c                                               found neighbours it,jt with
c                                               common nodes i1 and i3
             neighbour(2,it) = jt
             k1 = vertices(1,jt) 
             k2 = vertices(2,jt)
             k3 = vertices(3,jt)
             if(k1.ne.i1.and.k1.ne.i3)then 
                neighbour(1,jt) = it
             else if(k2.ne.i1.and.k2.ne.i3)then 
                neighbour(2,jt) = it
             else
                neighbour(3,jt) = it
             end if
           end if
c                                               compare neighbours i2 i3
           jt = 0
           j1 = nnn(i2)
           j2 = nnn(i2+1) - 1 
           do 230 j = j1,j2
              if(nnlist(j).eq.0)then
                 nnlist(j) = i3
                 ntlist(j) = it
c                                               if we have recorded connection
c                                               then jump out of loop
                 go to 231
              else if(nnlist(j).eq.i3.and.ntlist(j).ne.it)then
                 jt = ntlist(j)
                 go to 231
              end if
  230      continue
  231      continue
c                                               if neighbours are found then
c                                               skip second loop 
           if(jt.eq.0)then
              j1 = nnn(i3)
              j2 = nnn(i3+1) - 1 
              do 232 j = j1,j2
                 if(nnlist(j).eq.0)then
                    nnlist(j) = i2
                    ntlist(j) = it
c                                               if we have inserted connection
c                                               then jump out of loop
                    go to 233
                 end if
  232         continue
              end if
  233      continue
           if(jt.ne.0)then
c                                               found neighbours it,jt with
c                                               common nodes i2 and i3
             neighbour(1,it) = jt
             k1 = vertices(1,jt) 
             k2 = vertices(2,jt)
             k3 = vertices(3,jt)
             if(k1.ne.i2.and.k1.ne.i3)then 
                neighbour(1,jt) = it
             else if(k2.ne.i2.and.k2.ne.i3)then 
                neighbour(2,jt) = it
             else
                neighbour(3,jt) = it
             end if
           end if

 25     continue

c                                               build ntlist as the
c                                               indices of triangles
c                                               attached to each node
c
c                                               Note this assumes that the
c                                               Number of triangles attached
c                                               to each node is equal to
c                                               the number of neighbours
c                                               (BUT we're leaving an extra
c                                               component in ntlist that is zero)
c                                               and hence only works for
c                                               nodes on a sphere 
c                                               (since there is no boundary).
c                                               It will not work for
c                                               nodes in a plane.
        do i=np+1,2,-1
           nnn(i) = nnn(i-1)
        end do
        do i = 1,nmax
           ntlist(i) = 0
        end do
        do i=1,nt
           i1 = vertices(1,i)
           i2 = vertices(2,i)
           i3 = vertices(3,i)
           ntlist(nnn(i1+1)) = i
           ntlist(nnn(i2+1)) = i
           ntlist(nnn(i3+1)) = i
           nnn(i1+1) = nnn(i1+1) + 1
           nnn(i2+1) = nnn(i2+1) + 1
           nnn(i3+1) = nnn(i3+1) + 1
        end do
        nnn(1) = 1
        do i=2,np+1
           nnn(i) = nnn(i) + 1
        end do

        if(nnwrite)write(*,*)' built neighbour nnlist and ntlist'

        return
        end
c
c------------------------------------------------------------------------
c
c       Subroutine ccentres_sphere - calculates the circum centres
c                                    of spherical triangles
c
c       Input:
c
c       Output:
c
c
c                                       M. Sambridge, RSES, Sept. 1998.
c
c------------------------------------------------------------------------
c
c       Subroutine ccentres_sphere
        Subroutine ccentressphere
     &             (points,vertices,nt,neighbour,centres,iwrite,lu)

        real*8          points(3,*)
        real*8          centres(2,*)
        real*8          p(2,3),a(3),b(3),an(3)
        real*8          d1,d3,d5
        real*8          acroblen,adotb,blen,alen,dan
        real*8          c(2),xc,yc,zc,rc
        real*8          pi,rad2deg
        integer         vertices(4,*)
        integer         neighbour(3,*)

        pi = 3.141592653589793d0
        rad2deg = 180.d0/pi

        do i=1,nt

           a(1) = points(1,vertices(2,i))-points(1,vertices(1,i))
           a(2) = points(2,vertices(2,i))-points(2,vertices(1,i))
           a(3) = points(3,vertices(2,i))-points(3,vertices(1,i))
           b(1) = points(1,vertices(3,i))-points(1,vertices(1,i))
           b(2) = points(2,vertices(3,i))-points(2,vertices(1,i))
           b(3) = points(3,vertices(3,i))-points(3,vertices(1,i))

           alen  = dsqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
           blen  = dsqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))
           adotb = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
           d1 = a(1)*b(2) - a(2)*b(1)
           d1 = d1*d1
           d3 = a(2)*b(3) - a(3)*b(2)
           d3 = d3*d3
           d5 = a(3)*b(1) - a(1)*b(3)
           d5 = d5*d5
           acroblen = dsqrt(d1+d3+d5)
           p(1,1) = 0.d0
           p(2,1) = 0.d0
           p(1,2) = adotb/blen
           p(2,2) = acroblen/blen
           p(1,3) = blen
           p(2,3) = 0.d0
           an(1)  = a(1) - p(1,2)*b(1)/blen
           an(2)  = a(2) - p(1,2)*b(2)/blen
           an(3)  = a(3) - p(1,2)*b(3)/blen
           dan    = dsqrt(an(1)*an(1)+an(2)*an(2)+an(3)*an(3))
           an(1)  = an(1)/dan
           an(2)  = an(2)/dan
           an(3)  = an(3)/dan

           call circum(p,c)

           xc = points(1,vertices(1,i))
     &          + c(1)*b(1)/blen + c(2)*an(1)
           yc = points(2,vertices(1,i))
     &          + c(1)*b(2)/blen + c(2)*an(2)
           zc = points(3,vertices(1,i))
     &          + c(1)*b(3)/blen + c(2)*an(3)
c
c                                               check distance
c          d1 = (xc-points(1,vertices(1,i)))
c          d2 = (yc-points(2,vertices(1,i)))
c          d3 = (zc-points(3,vertices(1,i)))
c           da = d1*d1 + d2*d2 +d3*d3
c          d1 = (xc-points(1,vertices(2,i)))
c          d2 = (yc-points(2,vertices(2,i)))
c          d3 = (zc-points(3,vertices(2,i)))
c           db = d1*d1 + d2*d2 +d3*d3
c          d1 = (xc-points(1,vertices(3,i)))
c          d2 = (yc-points(2,vertices(3,i)))
c          d3 = (zc-points(3,vertices(3,i)))
c           dc = d1*d1 + d2*d2 +d3*d3
c           write(*,*)i,':',da,db,dc

c                                               convert point to
c                                               spherical co-ords
           rc = dsqrt(xc*xc+yc*yc+zc*zc)
           centres(1,i) = 90.d0 - rad2deg*dacos(zc/rc)
           centres(2,i) = rad2deg*datan2(yc,xc)
           if(centres(2,i).lt.0.d0)then
              centres(2,i) = 360.d0+centres(2,i)
           end if

        end do

c                                               write out edges of
c                                               voronoi cells
c
        if(iwrite.eq.1)then

c          open(lu,file='voronoi_edges.out',status='unknown')
           open(lu,file='voronoi.xy',status='unknown')
           write(lu,'(1a)')'>'

           do i=1,nt
              do j=1,3
                k1 = neighbour(j,i)
c                                               lat1,lon1 lat2,lon2

c               write(lu,*)centres(1,i),centres(2,i),
c    &                    centres(1,k1),centres(2,k1)
                if(k1.lt.i)then
                  write(lu,*)centres(1,i),centres(2,i)
                  write(lu,*)centres(1,k1),centres(2,k1)
                  write(lu,'(1a)')'>'
                end if

              end do
           end do
           close(lu)

        end if


        return
        end
c
c------------------------------------------------------------------------
c
c       Circum - calculates circumcentre of a triangle
c
c
c       Input:
c               p(2,3)                  array of node points
c
c       Output:
c               centre(2)              centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of
c                                       circumcircle about input Delaunay
c                                       triangle
c
c       Comments:
c
c                No calls to other routines.
c
c                                       M. Sambridge, RSES, Sept. 1999.
c
c------------------------------------------------------------------------
c
        Subroutine circum(p,centre)
c
        real*8          p(2,3)
        real*8          centre(2)
        real*8          x1,x2,x3,y1,y2,y3,x,y
        real*8          dx2m1,dx2p1,dy2m1,dy2p1
        real*8          dx3m1,dx3p1,dy3m1,dy3p1
        real*8          denom
c                                               Find centre of
c                                               Circumcircle
           x1 = p(1,1)
           x2 = p(1,2)
           x3 = p(1,3)
           y1 = p(2,1)
           y2 = p(2,2)
           y3 = p(2,3)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1
           x = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1*0.5d0
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0*dy2m1)/
     &         (denom)

           y = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1)*0.5d0)/
     &         (denom)

           centre(1) = x
           centre(2) = y
           x1 = x - x1
           y1 = y - y1
c          centre(3) = x1*x1 + y1*y1

        return
        end

c------------------------------------------------------------------------
c
c	find_node_sph - locates the nearest node to input point 
c			on a sphere
c
c	Input:
c		xd(2)  			co-ordinates of input points	
c					to be passed to misfit routine
c		points(2,np)		array of node points	
c	        node			first guess of node nearest xd(2)
c		nnn			Integer pointer array used for nnlist
c		nnlist			Integer array containing list of natural
c				        neighbours to each node 
c				        (built by build_nn)
c
c	Output:
c	        node			index of nearest node (or Voronoi cell
c				        containing input point.
c	        walk			number of cells encountered in 
c					directed walk 
c		xlon			Longitude of nearest node
c		xlat			Latitude of nearest node
c
c	Comments:
c		 No calls to other routines.
c
c					M. Sambridge, RSES, Feb. 1998.
c
c------------------------------------------------------------------------
c
c       Subroutine find_node_sph (xd,nodein,walk,xlon,xlat)
c       Subroutine find_node_sph (intrala,xd,nodein,walk,xlon,xlat,idim)
        Subroutine findnodesph (intrala,xd,nodein,walk,xlon,xlat,idim)

        include         'nn.param'

        common/nnsetuparrays/points(3,np_max),
     &                       vertices(4,nt_max),
     &                       neighbour(3,nt_max),
     &                       coords(2,np_max),
     &                       centres(2,nt_max),
     &                       nnn(np_max+1),
     &                       nnlist(nwork3d),
     &                       ntlist(nwork3d)

	common/walkpath/path(np_max)

        real*8          points
        real*8          coords
        real*8          centres
        integer         vertices
        integer         neighbour
        integer         nnn
        integer         nnlist
        integer         ntlist

c       real*8          xd(2)
        real*8          xd(idim)
        real*8          amisfit
        real*8          amisfit_min
        real*8          xlat
        real*8          xlon
        integer         walk
        integer         path
c ##############################MODIF ERIC############################
      if(intrala.eq.0) then
c         write(1000,*)' Vertices and neighbour arrays'
c          nt= 9122  
c          np= 4568 
c         do 10 i=1,nt
c            write(1000,*)i,' :',vertices(1,i),vertices(2,i),
c     &                vertices(3,i),' n:',
c     &                neighbour(1,i),neighbour(2,i),
c     &                neighbour(3,i)
c 10      continue
c                                       special write out of nnlist
c         write(1001,*)'nnlist'
c         do i=1,np
c            write(1001,300)i,(nnlist(k),k=nnn(i),nnn(i+1)-1)
c           write(1001,300)i,(ntlist(k),k=nnn(i),nnn(i+1)-1)
c         end do
c         write(1002,*)'npoint, coords'
c         do i=1,np
c            write(1002,*)i,' :',points(1,i),points(2,i),
c     &                points(3,i),' c:',
c     &                coords(1,i),coords(2,i)
c         end do
c
      endif
  300    format(1x,i4,':',50(i4,1x))
c#######################################################################
        walk = 1
	path(walk) = nodein

c       write(*,*)'  INSIDE FINDNODE :  XD(1)', xd(1),' XD(2)',xd(2)
c    & ,' NODE',nodein
c                                       Use a walk through natural neighbour
c                                       list to find nearest node
        noden = nodein
        if(noden.le.0)noden = 1

        call misfit(xd,noden,amisfit)

c       write(*,*)'inside find_node_sph starting node=',noden,
c    &           ' d=',amisfit
c       write(*,*)'xd1,2,nodein,walk,xlon,xlat,idim',xd(1)
c    &  ,xd(2),nodein,walk,xlon,xlat,idim

	amisfit_min = amisfit

 20     continue
        do i=nnn(noden),nnn(noden+1)-1
           n = nnlist(i)
           if(n.ne.0)then
	      call misfit(xd,n,amisfit)
	      if(amisfit.lt.amisfit_min)then
                 noden=n
		 amisfit_min = amisfit
c                write(*,*)'find_node_sph : newnode=',
c    &           noden,' misfit ',amisfit_min
                 walk = walk + 1
	         path(walk) = n
                 go to 20
              end if
           end if
        end do
c                                       found optimum node
        nodein = noden
	xlat = coords(1,nodein)
	xlon = coords(2,nodein)

c       write(*,*)'findnode : found node=',nodein,' d ',amisfit_min

        return
        end
c
c------------------------------------------------------------------------
c
c	check_neighbour - checks whether the array neighbour is
c			  consistent.
c
c	Input:
c	        neighbour(3,nt)		array of neighbouring tetrahedra.	
c					Neighbour(i,j) is the tetrahedron
c					opposite node i in tetrahedron j.
c	        nd			number of dimensions
c	        vertices(4,nt)		array of tetrahedron vertices	
c	        nt			number of tetrahedra
c
c	Output:
c               consistent		logical set to true if the
c					test is passed
c
c	Comments:
c		 This is a utility routine used to test whether the
c		 calculated array neighbour is internally consistent.
c		 A 2-D and 3-D version is contained with the one
c	  	 routine. The routine is not called by any nn routines
c		 but can be placed after the nn setup routine in the
c		 calling program to check the consistency of the
c		 neighbour matrix.
c
c		 Currently works for 2-D or 3-D cases only.
c
c		 If special option ignoreneg = .true. then
c		 negative entries in the 3D neighbour array
c		 are treated as zeros. This is so that the extend
c		 tetrahedra option can be used.
c
c		 Note that this is a special version of
c		 check_neighbour with array vertices dimensioned (4,*)
c		 so that it will work with internally generated
c		 convex hull vertices.
c	
c	         Calls no other routines.
c		 
c					M. Sambridge, RSES, May 1995.
c
c------------------------------------------------------------------------
c
c       Subroutine check_neighbour
	Subroutine checkneighbour
     &             (neighbour,nd,vertices,nt,ignoreneg,consistent)

	integer		neighbour(nd+1,*)
c	integer		vertices(nd+1,*)
	integer		vertices(4,*)
	logical		consistent
	logical		ignoreneg

        consistent = .true.

c						perform 2-D test
        if(nd.eq.2)then

	do 10 i=1,nt
	   do 20 j=1,3
              k = neighbour(j,i)
              j1 = mod(j,3)+1
              j2 = mod(j1,3)+1
              i1 = vertices(j1,i)
              i2 = vertices(j2,i)
              if(k.eq.0)go to 20
	      do 30 m=1,3
                 m1 = mod(m,3)+1
                 m2 = mod(m1,3)+1
                 k1 = vertices(m1,k) 
                 k2 = vertices(m2,k) 
                 if(neighbour(m,k).eq.i)then
                    if(k1.eq.i1.and.k2.eq.i2.or.
     &                 k1.eq.i2.and.k2.eq.i1)go to 20
                 else if(m.eq.3)then
                 write(*,*)' 2D neighbour matrix inconsistent'
                 write(*,*)' triangle',i,' has neighbour',k,' index',j
                 write(*,*)' triangle',k,' has no neighbour',i
                 consistent = .false.
                 end if
 30           continue
              consistent = .false.
              write(*,*)
     &        ' neighbour found but vertices inconsistent'
              write(*,*)
     &        ' neighbours triangle',i,' and',k
              write(*,*)
     &        ' vertices of triangle',i,' :',i1,i2,vertices(j,i)
              write(*,*)
     &        ' vertices of triangle',k,' :',k1,k2,vertices(m,k)
c             write(*,*)i,':',vertices(1,i),vertices(2,i),vertices(3,i)
c             write(*,*)k,':',vertices(1,k),vertices(2,k),vertices(3,k)
 20        continue
 10     continue

c						perform 3-D test
	else if(nd.eq.3)then

	do 110 i=1,nt
	   do 120 j=1,4
              k = neighbour(j,i)
              j1 = mod(j,4)+1
              j2 = mod(j1,4)+1
              j3 = mod(j2,4)+1
              i1 = vertices(j1,i)
              i2 = vertices(j2,i)
              i3 = vertices(j3,i)
              if(k.eq.0)go to 120
              if(ignoreneg.and.k.lt.0)go to 120
	      do 130 m=1,4
                 m1 = mod(m,4)+1
                 m2 = mod(m1,4)+1
                 m3 = mod(m2,4)+1
                 k1 = vertices(m1,k) 
                 k2 = vertices(m2,k) 
                 k3 = vertices(m3,k) 
                 if(neighbour(m,k).eq.i)then
                    if((k1.eq.i1.and.k2.eq.i2.and.k3.eq.i3).or.
     &                 (k1.eq.i1.and.k2.eq.i3.and.k3.eq.i2).or.
     &                 (k1.eq.i2.and.k2.eq.i1.and.k3.eq.i3).or.
     &                 (k1.eq.i2.and.k2.eq.i3.and.k3.eq.i1).or.
     &                 (k1.eq.i3.and.k2.eq.i1.and.k3.eq.i2).or.
     &                 (k1.eq.i3.and.k2.eq.i2.and.k3.eq.i1))go to 120
                 else if(m.eq.4)then
                 write(*,*)' 3D neighbour matrix inconsistent'
                 write(*,*)' triangle',i,' has neighbour',k,' index',j
                 write(*,*)' triangle',k,' has no neighbour',i
                 consistent = .false.
                 end if
 130          continue
              consistent = .false.
              write(*,*)
     &        ' neighbour found but vertices inconsistent'
              write(*,*)
     &        ' neighbours triangle',i,' and',k
              write(*,*)
     &        ' vertices of triangle',i,' :',i1,i2,i3,vertices(j,i)
              write(*,*)
     &        ' vertices of triangle',k,' :',k1,k2,k3,vertices(m,k)
 120       continue
 110    continue

	end if

	return
	end
c------------------------------------------------------------------------
c
c	find_node_brute - locates the nearest node to input point 
c			  on a sphere using a simple exhaustive search.
c
c	Input:
c		xd(2)  			co-ordinates of input points	
c	        node			first guess of node nearest xd(2)
c		nnn			Integer pointer array used for nnlist
c		nnlist			Integer array containing list of natural
c				        neighbours to each node 
c				        (built by build_nn)
c
c	Output:
c	        node			index of nearest node (or Voronoi cell)
c				        containing input point.
c
c	Comments:
c		 No calls to other routines.
c
c					M. Sambridge, RSES, Feb. 1998.
c
c------------------------------------------------------------------------
c
c       Subroutine find_node_brute (xd,node)
        Subroutine findnodebrute (xd,node)

        include         'nn.param'

        common/nnsetuparrays/points(3,np_max),
     &                       vertices(4,nt_max),
     &                       neighbour(3,nt_max),
     &                       coords(2,np_max),
     &                       centres(2,nt_max),
     &                       nnn(np_max+1),
     &                       nnlist(nwork3d),
     &                       ntlist(nwork3d)

        common/nnparam/np,nt

        real*8          points
        real*8          coords
        real*8          centres
        integer         vertices
        integer         neighbour
        integer         nnn
        integer         nnlist
        integer         ntlist

        real*8          xd(2)
        real*8          amisfit
        real*8          amisfit_min
c                                       Use a walk through natural neighbour
c                                       list to find nearest node

        call misfit(xd,1,amisfit_min)
	node = 1

	do i=2,np

           call misfit(xd,i,amisfit)

	   if(amisfit.lt.amisfit_min)then
	      amisfit_min = amisfit
              node = i
	   end if

	end do

        return
        end
c
c------------------------------------------------------------------------
c
c	Subroutine misfit_L2 - defines misfit criteria for minimization
c
c			    Here we use the distance from the input
c			    point to the input node.
c
c	Input:
c
c	xd			: current point
c	node			: current node
c
c	Output:
c
c	amisfit			: misfit value 
c
c------------------------------------------------------------------------
c
c       Subroutine misfit_L2(x,node,amisfit)
	Subroutine misfitL2(x,node,amisfit)

        include         'nn.param'

        common/nnsetuparrays/points(3,np_max),
     &                       vertices(4,nt_max),
     &                       neighbour(3,nt_max),
     &                       coords(2,np_max),
     &                       centres(2,nt_max),
     &                       nnn(np_max+1),
     &                       nnwork(nwork3d),
     &                       ntlist(nwork3d)
c
        data            degtorad/0.017453292/
c       data            rad/6371.0/
        data            rad/1.0/

        real*8 		x(2)
        real*8 		xd(3)
	real*8		points
	real*8		coords
        real*8          centres
	real*8		amisfit
	integer		vertices
	integer		neighbour
	integer		nnn
	integer		nnwork
	integer		ntlist
c
        depth = 0.0
        alat = x(1)
        alon = x(2)
        alat2 = (90-alat)*degtorad
        alon2 = alon*degtorad
        xd(1) = (rad-depth)*dble(sin(alat2)*cos(alon2))
        xd(2) = (rad-depth)*dble(sin(alat2)*sin(alon2))
        xd(3) = (rad-depth)*dble(cos(alat2))

        a = points(1,node)-xd(1)
        b = points(2,node)-xd(2)
        c = points(3,node)-xd(3)

        amisfit = a*a + b*b + c*c
c
	return
	end
c------------------------------------------------------------------------
c
c       Subroutine misfit - defines misfit criteria for minimization
c
c                           Here we use the geodesic distance from the input
c                           point to the input node on the surface of
c                           the sphere. Input point is in degrees. Output
c                           distance is in radians.
c       Input:
c
c       xd                      : current point xd(1) = lat, xd(2) = long
c       node                    : current node
c
c       Output:
c
c       amisfit                 : geodesic distance on surface of sphere
c                                 in radians
c
c                                               M. Sambridge, RSES, Feb. 1999
c
c------------------------------------------------------------------------
        Subroutine misfit(x,node,amisfit)

        include         'nn.param'
        
        implicit real*8(a-h,o-z)

        common/nnsetuparrays/points(3,np_max),
     &                       vertices(4,nt_max),
     &                       neighbour(3,nt_max),
     &                       coords(2,np_max),
     &                       centres(2,nt_max),
     &                       nnn(np_max+1),
     &                       nnwork(nwork3d),
     &                       ntlist(nwork3d)
c
c       data            degtorad/0.017453292/
c       data            rad/6371.0/
        data            rad/1.0/
        data            pi/3.141592653589793116/
        data            pihalf/1.570796326794896558/

        real*8          x(2)
c       real*8          xd(3)
        real*8          points
        real*8          coords
        real*8          centres    
        real*8          alat1,alat2,alon1,alon2
        real*8          amisfit
        integer         vertices
        integer         neighbour
        integer         nnn
        integer         nnwork
        integer         ntlist
c
        depth = 0.0
        alat1 = x(1)
        alon1 = x(2)
        alat2 = coords(1,node)
        alon2 = coords(2,node)
c       write(*,*)'X1,Y1,X2,Y2',alat1,alon1,alat2,alon2

c to radians:
        degtorad = pi/180.d0
        x1=alon1*degtorad
        x2=alon2*degtorad
        y1=alat1*degtorad
        y2=alat2*degtorad

c                                       change lat to colat
        thet1 = pihalf - y1
        thet2 = pihalf - y2

        carc = dcos(thet1)*dcos(thet2)
     &         + dsin(thet1)*dsin(thet2)*dcos(x2-x1)

        amisfit = dacos(carc)
c
        return
        end
c------------------------------------------------------------------------
c
c       Tlocsph - locates the triangle on a sphere containing input
c                 point x,y
c
c       Input:
c               x,y                     co-ordinates of input points
c                                       lat (not colat) and long (degrees)
c               points(3,np)            array of node points (x,y,z)
c               vertices(4,nt)          array of triangle vertices
c               neighbour(3,nt)         array of neighbouring triangles
c                                       Neighbour(i,j) is the triangle
c                                       opposite node i in triangle j,
c                                       stored counterclockwise about j.
c               rad                     radius of sphere
c               xo                      centre of sphere in Cartesian
c                                       co-ordinates
c               iloc                     first guess of triangle containing
c                                       (x, y).
c
c       Output:
c               iloc                     index of triangle containing
c                                       input point.
c
c       Comments:
c                If (x,y) is outside convex hull iloc is a `nearby' triangle
c                on the hull and out is set to true.
c
c                No calls to other routines.
c
c                                       M. Sambridge, RSES, Nov 1999.
c
c------------------------------------------------------------------------
c
        SUBROUTINE TLOCSPH3 (iloc,x,y)
c
        include		'nn.param'
        common/nnsetuparrays/points(3,np_max),
     &                       vertices(4,nt_max),
     &                       neighbour(3,nt_max),
     &                       coords(2,np_max),
     &                       centres(2,nt_max),
     &                       nnn(np_max+1),
     &                       nnlist(nwork3d),
     &                       ntlist(nwork3d)
      
        real*8          x
        real*8          y
        real*8          xo(3)
        real*8          points
        real*8          coords 
        real*8          centres
        integer         vertices
        integer         neighbour
        integer         nnn
        integer         nnlist
        integer         ntlist
        logical         out
        integer         c1(3)
        data            c1/2,3,1/
        data            degtorad/0.017453292/
c
c                                               convert input points to
c                                               Cartesian reference frame
c                                               Note that this assumes
c                                               lat NOT colat
        xo(1) = 0.0
        xo(2) = 0.0
        xo(3) = 0.0
        rad = 1.0
        out = .false.
        depth = 0.0
        write(*,*)
        write(*,*)'tloc,x,y',iloc,x,y
c       iloc=2
c       write(*,*)'BJ tloc,x,y',iloc,x,y
        alat2 = (90-x)*degtorad
        alon2 = y*degtorad
        xv = xo(1) + (rad-depth)*sin(alat2)*cos(alon2)
        yv = xo(2) + (rad-depth)*sin(alat2)*sin(alon2)
        zv = xo(3) + (rad-depth)*cos(alat2)

 10     continue
        if(out)then
           write(*,*)' Error in subroutine tlocsph'
           write(*,*)' can not find triangle containing input point'
           stop
        end if
c                                       point is outside convex hull
           i1 = vertices(1,iloc)
           i2 = vertices(2,iloc)
           i3 = vertices(3,iloc)
           write(*,*) i1,i2,i3,iloc
           write(*,*)coords(1,i1),coords(2,i1) 
           write(*,*)coords(1,i2),coords(2,i2) 
           write(*,*)coords(1,i3),coords(2,i3) 

        do 20 i=1,3
           j = c1(i)
           k = c1(j)
           i1 = vertices(i,iloc)
           i2 = vertices(j,iloc)
           i3 = vertices(k,iloc)
           x1 = points(1,i1)
           y1 = points(2,i1)
           z1 = points(3,i1)
           x2 = points(1,i2)
           y2 = points(2,i2)
           z2 = points(3,i2)
           x3 = points(1,i3)
           y3 = points(2,i3)
           z3 = points(3,i3)
           x4 = sngl(xo(1))
           y4 = sngl(xo(2))
           z4 = sngl(xo(3))
           ax = (y2-y4)*(z3-z4)-(z2-z4)*(y3-y4)
           ay = (z2-z4)*(x3-x4)-(x2-x4)*(z3-z4)
           az = (x2-x4)*(y3-y4)-(y2-y4)*(x3-x4)
           prod = (x2-x1)*ax+(y2-y1)*ay
     &            +(z2-z1)*az
           if(prod.lt.0.0)then
              ax = -ax
              ay = -ay
              az = -az
           end if

           b=ax*x4+ay*y4+az*z4

           prod = ax*xv + ay*yv + az*zv

           if(prod.gt.b)then
              if(neighbour(i,iloc).eq.0)then
                 out = .true.
              else
                 iloc = neighbour(i,iloc)
              end if
              go to 10
           end if
 20     continue

        write(*,*)'tloc',iloc
        return
        end
c
c------------------------------------------------------------------------
c
c	find_sphtriangle - A front end routine to tlocsph which locates
c		           the spherical triangle containing the 
c			   input point alat alon. Returns triangle
c			   index loc.
c
c                                       M. Sambridge, RSES, Nov 1999.
c
c------------------------------------------------------------------------
c
c	Subroutine find_sphtriangle(alat,alon,loc)
 	Subroutine tlocsph(loc,alat,alon)

        include         'nn.param'

        common/nnsetuparrays/points(3,np_max),
     &                       vertices(4,nt_max),
     &                       neighbour(3,nt_max),
     &                       coords(2,np_max),
     &                       centres(2,nt_max),
     &                       nnn(np_max+1),
     &                       nnlist(nwork3d),
     &                       ntlist(nwork3d)

	real*8		points
        real*8          coords
        real*8          centres
	integer		vertices
	integer		neighbour
	integer		nnn
	integer		nnlist
	integer		ntlist
	integer		loc

c       loc = 1
	write(*,*)' alat ',alat
	write(*,*)' alon ',alon
	write(*,*)' loc ',loc

        call tlocsph2
     &       (alat,alon,points,vertices,neighbour,loc)

	return
	end
c
c
c------------------------------------------------------------------------
c
c	Tlocsph - locates the triangle on a sphere containing input
c		  point x,y
c
c	Input:
c		x,y			co-ordinates of input points	
c					lat (not colat) and long
c		points(3,np)		array of node points (x,y,z)	
c	        vertices(4,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c	        rad			radius of sphere 
c	        xo			centre of sphere in Cartesian 
c					co-ordinates 
c	        loc			first guess of triangle containing
c					(x, y).
c
c	Output:
c	        loc			index of triangle containing 
c					input point.
c
c	Comments:
c		 If (x,y) is outside convex hull loc is a `nearby' triangle
c		 on the hull and out is set to true.
c
c		 No calls to other routines.
c
c                                       M. Sambridge, RSES, Nov 1999.
c
c------------------------------------------------------------------------
c
	Subroutine tlocsph2(x,y,points,vertices,neighbour,loc)
c
	real*8		points(3,*)
	real*8		xo(3)
	integer		vertices(4,*)
	integer		neighbour(3,*)
        logical		out
        integer		c1(3)
        data		c1/2,3,1/
        data            degtorad/0.017453292/
c
c                                               convert input points to
c                                               Cartesian reference frame
c                                               Note that this assumes
c                                               lat NOT colat
	xo(1) = 0.0
	xo(2) = 0.0
	xo(3) = 0.0
	rad = 1.0
	out = .false.
	depth = 0.0
	xv = rad*cos(x)
	alat2 = (90-x)*degtorad
        alon2 = y*degtorad
        xv = xo(1) + (rad-depth)*sin(alat2)*cos(alon2)
        yv = xo(2) + (rad-depth)*sin(alat2)*sin(alon2)
        zv = xo(3) + (rad-depth)*cos(alat2)

 10     continue
	if(out)then
	   write(*,*)' Error in subroutine tlocsph'
	   write(*,*)' can not find triangle containing input point'
	   stop
	end if
c					point is outside convex hull
        do 20 i=1,3
	   j = c1(i)
	   k = c1(j)
           i1 = vertices(i,loc)
           i2 = vertices(j,loc)
           i3 = vertices(k,loc)
           x1 = points(1,i1)
           y1 = points(2,i1)
           z1 = points(3,i1)
           x2 = points(1,i2)
           y2 = points(2,i2)
           z2 = points(3,i2)
           x3 = points(1,i3)
           y3 = points(2,i3)
           z3 = points(3,i3)
           x4 = sngl(xo(1))
           y4 = sngl(xo(2))
           z4 = sngl(xo(3))
           ax = (y2-y4)*(z3-z4)-(z2-z4)*(y3-y4)
           ay = (z2-z4)*(x3-x4)-(x2-x4)*(z3-z4)
           az = (x2-x4)*(y3-y4)-(y2-y4)*(x3-x4)
           prod = (x2-x1)*ax+(y2-y1)*ay
     &            +(z2-z1)*az
           if(prod.lt.0.0)then
              ax = -ax
              ay = -ay
              az = -az
           end if

           b=ax*x4+ay*y4+az*z4

	   prod = ax*xv + ay*yv + az*zv

	   if(prod.gt.b)then
	      if(neighbour(i,loc).eq.0)then
                 out = .true.
	      else
	         loc = neighbour(i,loc)
	      end if
	      go to 10
	   end if
 20     continue
	
	return
	end
