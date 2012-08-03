c     $Id: read_txt.f,v 1.1 2005/02/01 01:17:03 alessia Exp $
c----------------------------------------------------------------------
      subroutine readtxt(ifile,itype,npts,nhead,lat,lon,dep,mb,ename)
c----------------------------------------------------------------------
c
c   Clusters points defined by lat, lon into clusters of radius c_rad
c   using distaz to calculate distances
c
c   Input via call:
c   ifile     =  file index for .txt file
c   itype     =  1=jdy 2=wideformat
c   npts      =  number of events in .txt file
c   nhead     =  number of header lines in .txt file
c
c   Output via call:
c   lat,lon   =  lat/lon arrays length npts
c   dep       =  depth array (length npts)
c   ename     =  event name (length npts)
c----------------------------------------------------------------------
      real lat,lon,mb
      integer npts,nhead,itype,ifile
      integer dep
      character*14 ename

      dimension lat(*),lon(*),dep(*),ename(*),mb(*)

      do 10 i=1,nhead
        read(ifile,*)
10    continue

c     loop over points
      do 20 i=1,npts
        read(ifile,1001) i1,i2,i3,i4,i5,i6,f1,lat(i),lon(i),dep(i),
     *                   mb(i), ename(i)
20    continue

1001  format(i4,1x,'(',i3,')',1x,4(i2,1x),f4.1,2x,f7.3,2x,f8.4,1x,i4,
     *       2x,f3.1,1x,a14)

      return
      end subroutine
c----------------------------------------------------------------------
