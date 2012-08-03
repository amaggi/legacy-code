c
c*******************************************************************************
c
c    Subroutine datafile
c
c    Author - Xuelin Qiu
c
c*******************************************************************************
c
c $Id: datafile.f,v 1.1.1.1 2002/07/12 11:15:19 maggi Exp $
c $Log: datafile.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:19  maggi
c
c
c Revision 1.1  2002/05/23 10:28:28  maggi
c Initial revision
c
c
      subroutine datafile(path, pname, fname)
c
c    Subroutine DATAFILE gets the full path of the binary font data file.
c
c    input  - pname  = path name for the binary font data file
c             fname  = file name for the binary font data file
c
c    output - path   = full path for the binary font data file
c

      character*40 path, pname, fname

      pname = '/dione/Codes/SWcodes/grx'
      path  = '/dione/Codes/SWcodes/grx/fonts.bin'

      return
      end
