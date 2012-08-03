c
c*******************************************************************************
c
c    Subroutine datafile
c
c    Author - Xuelin Qiu
c
c*******************************************************************************
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
