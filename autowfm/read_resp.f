	program read_resp
c
c	writes the pole/zero files for each station in the form required 
c	by the automated waveform inversion code.  If the RESP file is in 
c	Hz the units are changed to radiands.
c
	character station*4,network*2,stream*2,channel*3,units*3,type*1
	character start*8,end*8
	character outfile*14
	real*8 a0,frequency,zreal(25),zimag(25),preal(40),pimag(40)
	
	read(5,'(a2)') network
	read(5,'(a)') station
	read(5,'(a2)') stream
	read(5,'(a)') channel

c
c	Make output file name and open input file
c
	nchar = lnblnk(station)
	if(nchar .lt. 4) read(station(3:3),'(a)') station(4:4)
	outfile = '??.????.??.???'
	write(outfile(1:2),'(a)') network
	write(outfile(4:7),'(a)') station
	if(stream .eq. '') then
	   write(outfile(9:10),'(a2)') '..'
	else
	   write(outfile(9:10),'(a2)') stream
	endif
	if(channel .eq. 'LHZ') write(outfile(12:14),'(a)') 'lhz'
c	write(outfile(16:30),'(a)') date
	open(2,file=outfile,status='unknown')
c
c	open input files
c
	open(10,file='start',status='old')
	open(11,file='end',status='old')
	open(12,file='a0',status='old')
	open(13,file='sensitivity',status='old')
	open(14,file='units',status='old')
	open(15,file='frequency',status='old')
	open(16,file='type',status='old')
	open(17,file='nzeros',status='old')
	open(18,file='zeros',status='old')
	open(19,file='npoles',status='old')
	open(20,file='poles',status='old')
c
c	read in parameters
c
10	read(10,'(a)',end=99) start
	read(11,'(a)') end
	read(12,*) a0
	read(13,*) sensitivity
	read(14,'(a3)') units
c  freq. needed if to calculate Hz--> rad	
	read(15,*) frequency   
	read(16,'(a1)') type
	read(17,*) nzeros
	do 20 ii=1,nzeros
20	   read(18,*) zreal(ii),zimag(ii)
	read(19,*) npoles
	do 30 jj=1,npoles
30	   read(20,*) preal(jj),pimag(jj)

	write(2,'(a,1x,a,1x,a8,''-> '',a8)') station, channel,
     &		start, end
	write(2,'(''-2'')')
c
c	If the input polezero file is in Hz convert to radians
c
	if(type .eq. 'B') then
	   write(0,*) type,a0,nzeros,npoles,frequency
	   call convert_hz_rad(a0,nzeros,zreal,
     &		zimag,npoles,preal,pimag,frequency)
	endif
c
c	If the units are m/s add a zero to go to displacement
c
	if(units .eq. 'M/S') then
	   do 35 ii=1,nzeros
	      zreal(nzeros+2-ii) = zreal(nzeros+1-ii)
35	      zimag(nzeros+2-ii) = zimag(nzeros+1-ii)
	   nzeros = nzeros + 1
	   zreal(1) = 0.
	   zimag(1) = 0.
	endif
c
c	Output results
c
	write(2,'(''ZEROS '',i2)') nzeros
	do 40 ii=1,nzeros
40	   write(2,*) zreal(ii),zimag(ii)
	write(2,'(''POLES '',i2)') npoles
	do 50 jj=1,npoles
50	   write(2,*) preal(jj),pimag(jj)
	constant = a0 * sensitivity / 1e6
	write(2,'(''CONSTANT '',e12.6,'' muM'')') constant
c
c	read extra line from the files with units of Hz
c
	if(network .eq. 'II') then
	   read(12,*) a0
	   read(14,'(a3)') units
	   read(15,*) frequency
	   read(16,'(a1)') type
	   read(17,*) nzeros
	   read(19,*) npoles
	endif

	go to 10
99	stop
	end
c ---------------------------------------------------------------------
	subroutine convert_hz_rad(a0,noz,zr,zi,nop,pr,pi,freq)
c
c	all poles and zeros must be specified explicitly; generates new 
c	p-z file and Ao for rad/sec output; input poles and zeros are 
c	assumed to be in Hz.
c
	integer nop,noz
	real*8 a,b,c,d,pr,pi,zr,zi,tpi,prtem,pitem
	real*8 bot,a0,zrtem,zitem,freq
	dimension pr(40),pi(40),zr(25),zi(25)
c      2 pi	
	tpi=6.28319
c  mul zeroes and poles wirh 2pi
	do 10 ii=1,noz
	   zr(ii) = zr(ii) * tpi
10	   zi(ii) = zi(ii) * tpi
	do 20 ii=1,nop
	   pr(ii) = pr(ii) * tpi
20	   pi(ii) = pi(ii) * tpi

c	work out numerator  (zeroes) = (w-z1) * (w-z2) * (w-z3)..., w=2*pi*f

	zrtem=zr(1)
	zitem=zi(1)
	zrtem=-zrtem
	zitem=((tpi*freq)-zitem)

	do 40 j=2,noz
	   c=-zr(j)
	   d=(tpi*freq)-zi(j)
	   a=zrtem
	   b=zitem
	   zrtem=((a*c)-(b*d))
	   zitem=((a*d)+(b*c))
40	   continue
     
c	work out denominator  (poles) (w-p1) * (w-p2) * (w-p3)..., w=2*pi*f

	prtem=pr(1)
	pitem=pi(1)
	prtem=-prtem
	pitem=((tpi*freq)-pitem)

	do 50 j=2,nop
	   c=-pr(j)
	   d=(tpi*freq)-pi(j)
	   a=prtem
	   b=pitem
	   prtem=((a*c)-(b*d))
	   pitem=((a*d)+(b*c))
50	   continue

	a=zrtem
	b=zitem
	c=prtem
	d=pitem

c	take modulus of numerator and denominator

	top=sqrt((a*a)+(b*b))
	bot=sqrt((c*c)+(d*d))
      
	write(*,*) 'Will convert Hz --> rad'
	write(*,*) 'top=',top,' bot=',bot

c	calculate normalization

	a0=bot/top

	return
	end
