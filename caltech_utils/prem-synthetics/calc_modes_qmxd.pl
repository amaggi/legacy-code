#!/usr/bin/perl
use vars qw($opt_m $opt_e $opt_S $opt_d $opt_c $opt_n $opt_s $opt_O $opt_C);

use Time::Local;
use POSIX;
use Math::Trig;
use File::Basename;
use Getopt::Std;
use lib '/opt/seismo-util/lib/perl/';
use CMT_TOOLS;

sub Usage {
  print STDERR <<END;
------------------------------------------------------------------------------------
  USAGE: calc_modes_qmxd.pl -m CMTfile -n ns -s s/s [-t Tmin] [-T Tmax]  
         -O output_dir -c comps -C [-S station_file] |[-d datafiles]
  FUNCTION:  Wrapper for QmXD, a mode summation code.
	 Writes the RECORDHEADERS file needed for the mode summation and runs QmXD.
  OPTIONS:
<--EVENT INFO-->
      -m CMTfile -- CMTSOLUTION file in Harvard format (multiple sources ok)
	     The t=0 of synthetics is taken to be the centroid time of the 
             first source. [Default: ./CMTSOLUTION]
<--STATION INFO-->
      -S station_file -- each entry looks like:
	     AAK   II   42.6390    74.4940 1645.0    0.0
	     name  netw   lat        lon    altitude burial-depth
      -d  -- Use the stations of the given datafiles and
	     the components given by the -c option.
      -c comps -- e.g. LHZ/LHE/LHN or LHZ/LHT
<--SYNTHETICS INFO-->
      -n ns  -- total number of samples in synthetic
      -s s/s -- samples per second
      -t Tmin -- min period to consider (default 8s)
      -T Tmax -- max period to consider (default 8000s)
<--OUTPUT--> Writes into ascii files with names like AAK.II.LHZ.qmxd.
      -O output_dir -- where to put synthetics.
      -C -- calculate modes, if not specified the script only 
	    writes RECORDHEADERS file
  *** Vala, June 28, 2002 last update: March 24, 2005 ***
END
  exit 1;
}

#Checking command line:
@ARGV>0 or Usage();
if(!getopts('m:e:S:dc:n:s:t:T:O:C')){ die "Check input arguments\n";}
if(!$opt_m and !(-e "CMTSOLUTION")){die "Specify moment tensor file, -m (default ./CMTSOLUTION )\n";}
if(!$opt_n){die "Specify total number of samples\n";}
if(!$opt_s){die "Specify number of samples per second\n";}
if(!$opt_c){die "Specify components (-c)\n";}
if((!$opt_S and !$opt_d)|($opt_S and $opt_d)){die "Specify one of -S and -d \n";}
if($opt_O){$output_dir = $opt_O;}
else{$output_dir = ".";}
if($opt_t > 8){$Tmin = $opt_t;}
elsif($opt_t && $opt_t < 8){die "Minimum Tmin is 8 s\n";}
else{$Tmin = "8";}
if($opt_T && $opt_T > $opt_t){$Tmax = $opt_T;}
elsif($opt_T && $opt_T < $opt_t){die "Tmax should be larger than Tmin\n";}
else{$Tmax = "8000";}
if($Tmin !~/\./){$Tmin = "$Tmin.";}
if($Tmax !~/\./){$Tmax = "$Tmax.";}
print "Period range of synthetics: $Tmin - $Tmax seconds\n";
#other settings:
@data_files = @ARGV;
$saclst = "/opt/seismo-util/bin/saclst";
$rec_file = "RECORDHEADERS";
$spsec = $opt_s;
$nsp = $opt_n;

#EVENT INFO

if($opt_m){$cmt_file = $opt_m;}
else{$cmt_file = "./CMTSOLUTION";}
print "Using event info from file: $cmt_file\n";

#-evaluate start time of synthetics:

($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$evid,$tshift,$hdur,$evla,$evlo,$evdp) = &get_cmt($cmt_file);
($year,$jday,$hr,$min,$sec,$msec) = &tdiff($oyear,$omonth,$oday,$ohr,$omin,$osec,$omsec,$tshift);

print "tshift and hdur of first source:$tshift $hdur \n";
$dt=1/$opt_s;
if((90-abs($evla)) < 0.001){die "Will not work for events located at the poles, move earthquake 0.001 degree of pole and it will work \n";}
if(2*$dt > $Tmin){die "Nyquist violated: You have to have more than two samples per wavelength (2*dt<Tmin) dt: $dt  Tmin:$Tmin\n";}

##STATION INFO and writing to RECORDHEADERS-file
open(REC,">$rec_file");
$format = "%-8s %-5s %-8s %8.4f %9.4f %6.1f %6.1f %5.1f %6.1f %12.4f %7.0f %4.0f %03.0f %02.0f %02.0f %06.3f\n";
if($opt_c){
  @comps = split(/\//,$opt_c);
  if($opt_S){
    die "Can\'t find station file: $opt_S\n" unless -e $opt_S;
    print "Using station-file $opt_S and components @comps \n";
    open(STAT,"$opt_S");
    @stations = <STAT>;
    close(STAT);
  }
  elsif($opt_d){
    print "Using stations for files given on command line and component(s) @comps \n";
    @stations = `$saclst kstnm knetwk stla stlo stel f @data_files | awk \'{print \$2, \$3, \$4, \$5, \$6}\'`;
  }
  else{
    die "Ups, something went wrong on the command line, stopped";
  }
  foreach $station (@stations){
    chomp($station);
    ($stat,$netw,$stla,$stlo,$stel,$stbur) = split(/ +/,$station);
    if(90-(abs($stla)) < 0.001 ){
      die "Will not work for station $stat $netw since it is located at the pole move it 0.001 degree of pole and it will work \n"
    }
    $stbur = "0" unless $bur !~ /^\s*$/;
    $rec_string = "";
    foreach $comp (@comps){
      if($comp =~/Z$/){
	$stazi = 0;
	$stdip = -90;
      }
      elsif($comp =~/E$/){
	$stazi = 90;
	$stdip = 0;
      }
      elsif($comp =~/N$/){
	$stazi = 0;
	$stdip = 0;
      }
      elsif($comp =~/R$/){
	($gcarc,$baz)=&gcarc_backaz($stla,$stlo,$evla,$evlo);
	$stazi = $baz-180;
	$stazi += 360 unless $stazi>0;
	$stdip = 0;
      }
      elsif($comp =~/T$/){
	($gcarc,$baz)=&gcarc_backaz($stla,$stlo,$evla,$evlo);
	$stazi = $baz-90;
	$stazi += 360 unless $stazi>0;
	$stdip = 0;
      }
      else{
	die "components have to end with Z E N R or T \n";
      }
      print "$stat, $netw, $comp, $stla, $stlo,$stel,$stbur,$stazi,$stdip,$spsec,$nsp,$year,$jday, $hr, $min\n";
      $rec_string = sprintf("$format",$stat, $netw, $comp, $stla, $stlo,$stel,$stbur,$stazi,$stdip,$spsec,$nsp,$year,$jday, $hr, $min, "$sec.$msec");
      if($records !~ /$rec_string/){
	print REC "$rec_string";
	$records = "$records $rec_string";
      }
    }
  }
  close(REC);
}
elsif($opt_D){
  print "Using stations and components for files given on command line\n";
  @stations = `$saclst kstnm knetwk stla stlo stel kcmpnm cmpaz cmpinc f @data_files | awk \'{print \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9}\'`;
  foreach $station (@stations){
    ($stat,$netw,$stla,$stlo,$stel,$comp, $stazi,$stdip) = split(/ +/,$station);
    $stbur = "0";
    $stdip = $stdip-90;
    printf REC ("$format",$stat, $netw, $comp, $stla, $stlo,$stel,$stbur,$stazi,$stdip,$spsec,$nsp,$year,$jday, $hr, $min, "$sec.$msec");
  }
  close(REC);
}
##Running QmXD


if($opt_C){
  print "Calculating mode-synthetics\n";
  print "Refer to $rec_file to see what records are being calculated\n";
  $command = "/opt/Modes/QmXD/bin/QmXDF $cmt_file RECORDHEADERS ${output_dir} /opt/Modes/QmXD/catalogues/seigsml_prem_an.750.8s /opt/Modes/QmXD/U4L8CHEBYPART /opt/Modes/QmXD/s+p.run11.06 9 $Tmin $Tmax 1 0 0 0 0 1";
  print "Running command: $command\n";
  `$command`;
}
else{
  print "Only writing record headers file, use -C to calculate the modes.\n";
}



##SUBROUTINES

sub gcarc_backaz {
local($stla,$stlo,$evla,$evlo)=@_;
#using names from formula in Lay & Wallace, page 134.

$A=deg2rad($stlo-$evlo);
$b=deg2rad(90-$evla);
$c=deg2rad(90-$stla);
$a=acos(cos($b)*cos($c)+sin($b)*sin($c)*cos($A)); #gcarc
$B=acos((cos($b)-cos($a)*cos($c))/(sin($a)*sin($c))); #backaz
$gcarc = rad2deg($a);
$backaz = rad2deg($B);
if(($stlo-$evlo<180 & $stlo-$evlo>0)|($stlo+360-$evlo<180 & $stlo+360-$evlo>0)){
  $backaz=360-$backaz;
}
return ($gcarc,$backaz);
}


