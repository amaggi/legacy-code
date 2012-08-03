c     $Id: cluster.f,v 1.1 2005/02/01 01:17:03 alessia Exp $
c----------------------------------------------------------------------
      subroutine cluster(lat,lon,npts,c_rad,itype,n_clust,c_index)
c----------------------------------------------------------------------
c
c   Clusters points defined by lat, lon into clusters of radius c_rad
c   using distaz to calculate distances
c
c   Input via call:
c   lat, lon  =  coordinates of points (arrays of length npts)
c   npts      =  number of points
c   c_rad     =  radius
c   itype     =  radius type 1=degree 2=km
c
c   Output via call:
c   n_clust   =  number of clusters
c   c_index   =  integer index of cluster for each point (array of length npts)

c   dkm -=- epicentral distance in kilometers
c
c   Note: Single precision version.
c
c----------------------------------------------------------------------
      real lat,lon,c_rad
      integer npts,itype,n_clust,c_index

      real azm,bzm,ddg,dkm

      dimension lat(*),lon(*),c_index(*)

c     initialise indexes
      n_clust=0
      do 10 i=1,npts
        c_index(i)=0
10    continue

c     loop over points
      do 20 i=1,npts
c       if this point is not already assigned to a cluster then use it as root
c       for a new cluster
        if (c_index(i).eq.0) then
          n_clust=n_clust+1
          c_index(i)=n_clust
          do 21 j=i,npts
c           if this point is not already assigned to a cluster then test it 
            if (c_index(j).eq.0) then
              call distaz(lat(i),lon(i),lat(j),lon(j),azm,bzm,ddg,dkm)
              if((itype.eq.1 .and. ddg.le.c_rad) .or. 
     *           (itype.eq.1 .and. dkm.le.c_rad) ) c_index(j)=n_clust
            endif
21        continue
        endif
20    continue

      return
      end subroutine
c----------------------------------------------------------------------
