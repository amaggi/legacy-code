#!/usr/bin/env python
# encoding: utf-8

import sys, os, getopt, AM_geopath, AM_subs, request_io, datetime

class Request_Pars(object):
  TYPE_CMT=0
  TYPE_DATETIME=1
  def __init__(self,opts):
    # go throught the option list setting the parameters
    temp_dict={}
    for (option,parameter) in opts:
      temp_dict[option]=parameter

    try:
      if temp_dict.has_key('--Tcmt'):
        self.type=Request_Pars.TYPE_CMT
        self.cmt_file=temp_dict['--Tcmt']
      elif temp_dict.has_key('--Tdatetime'):
        self.type=Request_Pars.TYPE_DATETIME
        self.datetime_string=temp_dict['--Tdatetime']
      else:
        raise getopt.GetoptError('No reference time')
      self.label=temp_dict['--label']
      self.datacenter=temp_dict['--datacenter']
      self.sta_file=temp_dict['--sta']
      self.cmp=temp_dict['--cmp']
      self.sec_before=int(float(temp_dict['-b']))
      self.sec_after =int(float(temp_dict['-a']))
      self.username =temp_dict['-u']

      self.print_pars()

    except KeyError:
      raise getopt.GetoptError('Incorrect command line')
    
  def print_pars(self): 
    print "==========================================================="
    print "Request data center : %s"%(self.datacenter)
    print "Request label : %s"%(self.label)
    if self.type == Request_Pars.TYPE_CMT:
      print "Reference time from cmtsolution file %s"%(self.cmt_file)
    elif self.type == Request_Pars.TYPE_DATETIME:
      print "Reference time from string: %s"%(self.datetime_string)
    print "Time span before reference time : %d"%(self.sec_before)
    print "Time span after  reference time : %d"%(self.sec_after)
    print "Stations from file %s"%(self.sta_file)
    print "==========================================================="


  def set_times(self):
    if self.type == Request_Pars.TYPE_CMT:
      ref_time = AM_geopath.CMTEvent(self.cmt_file).ref_time
    else:
      ref_time = AM_subs.string_to_datetime(self.datetime_string)


    self.start_time=ref_time + datetime.timedelta(0,-1*self.sec_before)
    self.end_time=ref_time + datetime.timedelta(0,self.sec_after)

  def set_stations(self):
    try:
      self.stations=AM_geopath.StationList([])
      self.stations.read_from_STATIONS_file(self.sta_file)
      # check various station selection parameters and revise self.stations
    except:
      self.stations=None
      raise

  def message_body(self):
    lines=[]
    for sta in self.stations.sta_list:
      lines.append("%-5s %-2s %4d %2d %2d %2d %2d %4.1f %4d %2d %2d %2d %2d %4.1f %1d %s\n"%\
                   (sta.staname,sta.net,self.start_time.year,self.start_time.month,\
                    self.start_time.day,self.start_time.hour,self.start_time.minute,\
                    float(self.start_time.second),self.end_time.year,self.end_time.month,\
                    self.end_time.day, self.end_time.hour, self.end_time.minute, \
                    float(self.end_time.second),1,self.cmp))
    return ''.join(lines)

if __name__== '__main__':

  try:
    (opts,args_proper) = getopt.getopt(sys.argv[1:],'b:a:u:',['Tcmt=','Tdatetime=','label=','datacenter=','cmp=','sta='])
    request_parameters=Request_Pars(opts)
    request_parameters.set_times()
    request_parameters.set_stations()
    msg_body=request_parameters.message_body()
    message=request_io.construct_breqfast_email(request_parameters.username,request_parameters.label,request_paramters.datacenter,msg_body)
    print message
    print "\nSending email..."
    request_io.send_breqfast_email(request_parameters.username,request_parameters.datacenter,message)
  

  except getopt.GetoptError:
    print "Usage: data_request.py --T[cmt CMTSOLUTION |datetime yyyymmddhhmmss] --label a_request_label --datacenter iris|orfeus --sta STATIONS_file --cmp [B|L]H[E|N|Z|?] -b n_seconds_before -a n_seconds_after -u username"
