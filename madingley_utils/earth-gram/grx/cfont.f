c
c*******************************************************************************
c
c    Subroutine cfont
c
c    Author - Danny Harvey
c
c*******************************************************************************
c
      subroutine cfont(jfont)
c
c    Subroutine CFONT changes the font used for text labelling.
c
c    input  - jfont  = font number
c
c                      Hershey vector font family:
c                      = 1 - KRoman
c                      = 2 - KGreek
c                      = 3 - SRoman
c                      = 4 - SGreek
c                      = 5 - SScript
c                      = 6 - PIRoman
c                      = 7 - PIGreek
c                      = 8 - PIItalics
c                      = 9 - PNRoman
c                      = 10 - PNGreek
c                      = 11 - PNItalics
c                      = 12 - DRoman
c                      = 13 - CScript
c                      = 14 - CCyrilic
c                      = 15 - TRoman
c                      = 16 - TItalics
c                      = 17 - Gothic German
c                      = 18 - Gothic English
c                      = 19 - Gothic Italian
c                      = 20 - Map symbols
c                      = 21 - Math, Astronomical, and Astrological symbols
c                      These are the original Hershey vector fonts for both
c                      X-windows and PostScript.
c
c                      Hershey-PostScript font family:
c                      112 - DRoman, Heletica
c                      115 - TRoman, Times-Roman
c                      116 - TItalics, Times-Italic
c                      130 - TRoman, Palatino-Roman
c                      131 - TItalics, Palatino-Italic
c                      132 - TRoman, Palatino-Bold
c                      These use the original Hershey vector fonts for the
c                      screen, but PostScript fonts for hard copy output.
c
c    This subroutine allows for more than one font to be stored in core
c    memory at a time depending on how the parameter variable NPMAX
c    is set, however this will only affect the response time and will be
c    otherwise transparent to the user.
c
c    parameters  -  NPMAX  = maximum number of words in the character stroke
c                            array, IPOINT.  This must be set large enough
c                            to accomodate AT LEAST the largest single font
c                            that will be used (e.g. font 18 ~ 9000).
c                            An error exit will occur if this number is set
c                            too small.
c                   LUFNT  = FORTRAN logical unit number for the binary
c                            font file
c
      parameter  (NPMAX = 10000)
      parameter  (LUFNT = 9)
c
      common /fonts/ ifoff,ipoint(NPMAX)
      character*40 path
      integer*2 ipoint
      integer*2 ifoff
c
      integer*2      jjfont
c
      integer*2  itot(21),ifont(21),nfonts
      integer*2  jfoff(21)
      integer*2 ifnt, nrd
c
      data  (itot(i),i=1,21) /
     1       1943,2043,2443,2770,3919,3823,3957,4026,4407,4624,
     2       4710,4805,5659,4960,7572,7437,8938,8645,7153,3036,
     3       4227  /
      data  nfonts  / 0 /
      data  ifoff   / -1 /
      data  jfoff(1) / 1 /
c
c    look to see if font is already loaded
c
      jjfont = jfont
      llfont = jfont
      if (llfont .lt. 100) llfont = 0
      if (llfont .eq. 0 .or. llfont .eq. 112 .or. llfont .eq. 115
     1                  .or. llfont .eq. 116 .or. llfont .eq. 130
     1                  .or. llfont .eq. 131 .or. llfont .eq. 132) then
      else
	write (6, '(/a)') 'cfont: Illegal font.'
	llfont = 0
      end if
      call hdstlf (llfont)
      if (jjfont .eq. 112) jjfont = 12
      if (jjfont .eq. 115) jjfont = 15
      if (jjfont .eq. 116) jjfont = 16
      if (jjfont .eq. 130) jjfont = 15
      if (jjfont .eq. 131) jjfont = 16
      if (jjfont .eq. 132) jjfont = 15
      if (jjfont .lt. 1)  jjfont = 1
      if (jjfont .gt. 21)  jjfont = 1
      if (nfonts .eq. 0)  go to 11
      do 10  i = 1,nfonts
      if (jjfont .eq. ifont(i))  go to 20
   10 continue
   11 continue
c
c    font not loaded - first determine where font can fit
c
      jtot = itot(jjfont) + 191
      if (jtot .gt. NPMAX)  go to 900
c
c    now shift other fonts down deleting the ones that wont fit
c
      ktot = jtot
      if (nfonts .eq. 0)  go to 105
      do 100  i = 1,nfonts
      ktot = ktot + itot(ifont(i)) + 191
      if (ktot .gt. NPMAX)  go to 110
  100 continue
  105 i = nfonts + 1
  110 nfonts = i
      if (nfonts .eq. 1)  go to 121
      nfm1 = nfonts - 1
      do 120  i = 1,nfm1,1
      j = nfonts - i + 1
      jm1 = j - 1
      jfoff(j) = jfoff(jm1) + jtot
      ifont(j) = ifont(jm1)
      j1 = jfoff(j) - 1
      j2 = jfoff(jm1) - 1
      jjtot = itot(ifont(j)) + 191
      do 125  k = 1,jjtot
  125 ipoint(j1+k) = ipoint(j2+k)
  120 continue
  121 continue
c
c    read in the font data from the binary file into the first slot
c
      call datafile (path, 'GRX_FONTPATH', 'fonts.bin')
      if (path .eq. ' ') then
997   	   write(0,*) 'cfont:  Error opening fonts.bin'
           stop
      endif

        open(unit=LUFNT,file=path,access='sequential',
     2       form='unformatted',status='old',err=997)

998   rewind LUFNT
      j = jjfont - 1
      if (j .eq. 0)  go to 150
c
c    position file to correct record
c
      do 155  i = 1,j
      read (LUFNT,end=910,err=920) nrd,ifnt
      do 156  jj = 1, nrd
  156 read (LUFNT,end=910,err=920)
  155 continue
  150 continue
c
c    read file
c
      read (LUFNT,end=910,err=920) nrd,ifont(1),(ipoint(i),i=1,191)
      j1 = 192
      do 157  jj = 1, nrd
      j2 = j1 + 4000 - 1
      if (j2 .gt. jtot) j2 = jtot
      read (LUFNT,end=910,err=920) (ipoint(i),i=j1,j2)
  157 j1 = j2 + 1
      if (ifont(1) .ne. jjfont)  go to 930
      jfoff(1) = 1
      ifoff = 1
c
c    close file and exit
c
      close(LUFNT)
      return
c
c    font already loaded
c
   20 ifoff = jfoff(i)
      return
c
c    error exits
c
c    font buffer array too small for single font
c
  900 write (6,901) jjfont,jtot,NPMAX
  901 format (' CFONT: font no. ',i2,' requires ',i5,' words',
     1  ' in font buffer'/'        buffer set at ',i5,' too small',
     2  ' - run stopped')
      stop
c
c    premature end of file on input font file
c
  910 write (6,911)
  911 format (' CFONT: unexpected EOF on input font file',
     1 ' - run stopped')
      stop
c
c    read error on input font file
c
  920 write (6,921)
  921 format (' CFONT: read error on input font file',
     1  ' - run stopped')
      stop
c
c    font no. on input file doesnt match requested font no.
c
  930 write (6,931) jjfont,ifont(1)
  931 format (' CFONT: requested font no. ',i2,' doesnt match',
     1  ' font no. ',i2,' on input file'/'        - run stopped')
      stop
c
c    end of subroutine CFONT
c
      end
