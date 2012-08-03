      subroutine ngrid (xsthick, xshue, xslight, xssat, ixsstyl,
     +                  xnthick, xnhue, xnlight, xnsat, ixnstyl,
     +                  ysthick, yshue, yslight, yssat, iysstyl,
     +                  ynthick, ynhue, ynlight, ynsat, iynstyl)
c
      real*4            xsthick, xshue, xslight, xssat
      integer                                           ixsstyl
      real*4            xnthick, xnhue, xnlight, xnsat
      integer                                           ixnstyl
      real*4            ysthick, yshue, yslight, yssat
      integer                                           iysstyl
      real*4            ynthick, ynhue, ynlight, ynsat
      integer                                           iynstyl
c
c    ngrid specifies grid drawing prior to calls to axis.
c
c    Inputs -	xsthick	= Thickness of X-axis grid lines corresponding
c			  to small tic marks. If < 0.0, then no
c			  grid lines are drawn for small tic marks.
c		xshue	= Color hue of X-axis small tic mark grid lines.
c		xslight	= Color lightness of X-axis small tic mark grid lines.
c		xssat	= Color saturation of X-axis small tic mark grid lines.
c		ixsstyl	= Line style of X-axis small tic mark grid lines.
c		xnthick, xnhue, xnlight, xnsat, ixnstyl
c			= Same for X-axis grid lines corresponding
c			  to large (or numbered) tic marks.
c		ysthick, yshue, yslight, yssat, iysstyl,
c		ynthick, ynhue, ynlight, ynsat, iynstyl
c			= Same for Y-axis grid lines.
c
      common /ngrcom/ xsth, xsh, xsl, xss, ixss,
     +                xnth, xnh, xnl, xns, ixns,
     +                ysth, ysh, ysl, yss, iyss,
     +                ynth, ynh, ynl, yns, iyns
c
      xsth = xsthick
      xsh = xshue
      xsl = xslight
      xss = xssat
      ixss = ixsstyl
      xnth = xnthick
      xnh = xnhue
      xnl = xnlight
      xns = xnsat
      ixns = ixnstyl
      ysth = ysthick
      ysh = yshue
      ysl = yslight
      yss = yssat
      iyss = iysstyl
      ynth = ynthick
      ynh = ynhue
      ynl = ynlight
      yns = ynsat
      iyns = iynstyl
c
      return
      end
