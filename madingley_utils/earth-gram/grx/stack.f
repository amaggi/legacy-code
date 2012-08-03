      subroutine stack(xdim,ydim,xlow,ylow,xmin,xmax,pos,
     1  vdim,ymin,ymax)
c
c    routine stack sets dimensions and scales etc.
c
      if (pos .ge. 1.) go to 100
      if (pos .le. 0.)  go to 100
      yydim = 2.*vdim
      yylow = ydim*(1.-pos)
      yylow = yylow - vdim + ylow
      ylw2 = yylow
      if (yylow .lt. ylow) ylw2 = ylow
      if (yylow .gt. ylow+ydim)  go to 100
      ydm2 = yydim - ylw2 + yylow
      if (ylw2 + ydm2 .gt. ylow+ydim) ydm2 = ylow+ydim - ylw2
      call setdim(xdim,ydm2,xlow,ylw2)
      yymax = ymax*(ylw2+ydm2-yylow-yydim/2.)*2./yydim
      yymin = ymin*(yylow+yydim/2.-ylw2)*2./yydim
      call setscl(xmin,xmax,yymin,yymax)
  100 return
      end
