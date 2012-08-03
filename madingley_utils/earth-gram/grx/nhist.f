      subroutine nhistbin (ndata, data, iwt, wts, nbin, bin, binw,
     +                     nbmax, nbout, bout, bwout, yout)
c
      integer              ndata
      real*4                      data(ndata)
      integer                           iwt
      real*4                                 wts(ndata)
      integer                                     nbin
      real*4                                          bin(*)
      real*4                                               binw(*)
      integer              nbmax
      integer                     nbout
      real*4                            bout(nbmax)
      real*4                                  bwout(nbmax)
      real*4                                         yout(nbmax)
c
c    nhistbin will compute weighted histogram sums. nhistplot
c    can be called to make the plot.
c
c    Inputs -	ndata	= The number of input data points.
c		data(ndata)
c			= The data values.
c		iwt	= Weighting flag. If set (1), then the
c			  wts() array is used to determine data
c			  weights. If not set (0), then all data
c			  weights are one by default.
c		wts(ndata)
c			= The data weights. This is only used
c			  when iwt = 1.
c		nbin	= The number of bins to hold histogram
c			  sums.
c		bin(nbin)
c			= The center points of each bin.
c		binw(nbin)
c			= The widths of each bin.
c			Note: if nbin is 0, then bin(1) and
c			binw(1) specify a starting center point
c			and a constant width to define as many
c			bins as are required by the data.
c		nbmax	= Maximum dimensions of bout() and yout().
c
c    Outputs -	nbout	= The output number of binned points. This
c			  will be the minimum of nbin and nbmax
c			  when nbin > 0 or equal to the minimum
c			  of nbmax and the number of bins required
c			  when nbin = 0.
c		bout(nbout)
c			= The bin center values.
c		bwout(nbout)
c			= The bin widths.
c		yout(nbout)
c			= The histogram weighted sums.
c
      nbout = 0
      if (ndata .lt. 1) return
      do 10  i = 1, nbmax
        yout(i) = 0.0
        if (nbin .lt. 1) then
          bout(i) = bin(1) + (i-1)*binw(1)
          bwout(i) = binw(1)
        else
          if (i .gt. nbin) go to 11
          bout(i) = bin(i)
          bwout(i) = binw(i)
        end if
   10 continue
   11 continue
      binwi = 1.0/binw(1)
      do 100  i = 1, ndata
        if (nbin .lt. 1) then
          x = (data(i)-bin(1))*binwi + 1.5
          j = x
        else
          do 120  j = 1, nbin
            if (data(i) .ge. bin(j)-0.5*binw(j) .and.
     +          data(i) .lt. bin(j)+0.5*binw(j)) go to 121
  120     continue
          go to 100
  121     continue
        end if
        if (j .gt. nbmax) go to 100
        if (j .lt. 1) go to 100
        if (j .gt. nbout) nbout = j
        if (iwt .eq. 1) then
          yout(j) = yout(j) + wts(i)
	else
          yout(j) = yout(j) + 1.0
	end if
  100 continue
c
      return
      end
c
      subroutine nhistplot (nhist, bhist, bwhist, yhist, 
     +                      itype, ifill, thick)
c
      integer               nhist
      real*4                       bhist(nhist)
      real*4                              bwhist(nhist)
      real*4                                      yhist(nhist)
      integer               itype
      integer                      ifill
      real*4                              thick
c
c    nihstplot will plot histogram bars.
c
c    Inputs -	nhist	= No. of histogram bins to plot.
c		bhist(nhist)
c			= The bin center values.
c		bwhist(nhist)
c			= The bin widths.
c		yhist(nhist)
c			= The histogram sums (heights).
c		itype	= Plot type.
c			  = 0 - Plot histogram bars vertically.
c			  = 1 - Plot histogram bars horizontally.
c		ifill	= Fill type.
c			  = 0 - No fill.
c			  = 1 - Fill bars with foreground color.
c			  = 2 - Fill bars with background color.
c		thick	= Line thickness.
c
      do 100  i = 1, nhist
        if (yhist(i) .eq. 0.0) go to 100
        if (bwhist(i) .eq. 0.0) go to 100
        if (itype .eq. 0) then
          x1 = bhist(i) - 0.5*bwhist(i)
          x2 = bhist(i) + 0.5*bwhist(i)
          y1 = 0.0
          y2 = yhist(i)
        else
          y1 = bhist(i) - 0.5*bwhist(i)
          y2 = bhist(i) + 0.5*bwhist(i)
          x1 = 0.0
          x2 = yhist(i)
        end if
        if (x2 .lt. x1) then
          x = x1
          x1 = x2
          x2 = x
        end if
        if (y2 .lt. y1) then
          y = y1
          y1 = y2
          y2 = y
        end if
        if (ifill .eq. 1) then
          call getbac (bh, bl, bs)
          call getfor (fh, fl, fs)
          call setbac (fh, fl, fs)
          call clrrgn (x1, x2, y1, y2)
          call setbac (bh, bl, bs)
        else if (ifill .eq. 2) then
          call clrrgn (x1, x2, y1, y2)
        end if
        call box (x1, x2, y1, y2, thick, 0, 0)
  100 continue
c
      return
      end
