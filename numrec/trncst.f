      SUBROUTINE trncst(x,y,iorder, ncity,n,de)
      INTEGER ncity,iorder(ncity),n(6)
      REAL de,x(ncity),y(ncity)
      INTEGER ii,j
      REAL xx(6),yy(6),alen,x1,x2,y1,y2
      alen(x1,x2,y1,y2)=sqrt((x2-x1)**2+(y2-y1)**2)
      n(4)=1+mod(n(3),ncity)
      n(5)=1+mod((n(1)+ncity-2),ncity)
      n(6)=1+mod(n(2),ncity)
      do 11 j=1,6
        ii=iorder(n(j))
        xx(j)=x(ii)
        yy(j)=y(ii)
11    continue
      de=-alen(xx(2),xx(6),yy(2),yy(6))-alen(xx(1),xx(5),yy(1),
     *yy(5))-alen(xx(3),xx(4),yy(3),yy(4))+alen(xx(1),xx(3),yy(1),
     *yy(3))+alen(xx(2),xx(4),yy(2),yy(4))+alen(xx(5),xx(6),yy(5),yy(6))
      return
      END
