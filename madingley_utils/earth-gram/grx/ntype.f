c
c*******************************************************************************
c
c    Subroutine ntype
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine ntype(xtype,ytype)
c
c    Subroutine ntype defines the type of plot scales (linear or log).
c
c    Inputs  - xtype  = CHARACTER*3 string which defines whether the
c                       horizontal axis is linear or logarithmic.
c                       = 'LIN' - linear
c                       = 'LOG' - logarithmic
c              ytype  = CHARACTER*3 string which defines whether the
c                       vertical axis is linear or logarithmic.
c                       = 'LIN' - linear
c                       = 'LOG' - logarithmic
c
      character*(*) xtype, ytype
c
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,iitran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
c
      character*1 c
c
      ixtype = 0
      ilen = len(xtype)
      if (ilen .gt. 0) then
        c = xtype(1:1)
        if (c .eq. 'l' .or. c .eq. 'L') then
          if (ilen .gt. 1) then
            c = xtype(2:2)
            if (c .eq. 'o' .or. c .eq. 'O') then
              if (ilen .gt. 2) then
                c = xtype(3:3)
                if (c .eq. 'g' .or. c .eq. 'G') then
                  ixtype = 1
                end if
              end if
            end if
          end if
        end if
      end if
c
      iytype = 0
      ilen = len(ytype)
      if (ilen .gt. 0) then
        c = ytype(1:1)
        if (c .eq. 'l' .or. c .eq. 'L') then
          if (ilen .gt. 1) then
            c = ytype(2:2)
            if (c .eq. 'o' .or. c .eq. 'O') then
              if (ilen .gt. 2) then
                c = ytype(3:3)
                if (c .eq. 'g' .or. c .eq. 'G') then
                  iytype = 1
                end if
              end if
            end if
          end if
        end if
      end if
c
      return
      end
      subroutine gettype (xtype, ytype)
      character*(*) xtype, ytype
      common /pdim/ xdim,ydim,xlow,ylow,rxdim,rydim,rxlow,rylow,
     1              xbl,ybl,xbh,ybh,xbm,ybm,iitran,tangle,
     2              ca,sa,cellht,cellwd,ixtype,iytype
      if (ixtype .eq. 1) then
        xtype = 'LOG'
      else
        xtype = 'LIN'
      end if
      if (iytype .eq. 1) then
        ytype = 'LOG'
      else
        ytype = 'LIN'
      end if
      return
      end
