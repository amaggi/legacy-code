#!/usr/bin/env python
# encoding: utf-8

import sys, os, email.Message, email.Utils, smtplib

class BREQfast_dbase(object):
  headers={}
  emails={}
  LABEL_FIELD=7
  def __init__(self):

    # Set up breqfast headers (tedious at best!)
    h=[]
    h.append(".NAME Alessia Maggi\n")
    h.append(".INST Ecole et Obervatoire des Sciences de la Terre\n")
    h.append(".MAIL 5 rue Rene Descartes, 67084 Strasbourg Cedex, FRANCE\n")
    h.append(".EMAIL alessia@sismo.u-strasbg.fr\n")
    h.append(".PHONE +33-390-245028\n")
    h.append(".FAX +33-390-240125\n")
    h.append(".MEDIA FTP\n")
    h.append(".LABEL XXLABELXX\n")
    h.append(".END\n")

    self.headers["alessia"]=h
    self.emails["alessia"]= email.Utils.formataddr(("Alessia Maggi","alessia@sismo.u-strasbg.fr"))

    h=[]
    h.append(".NAME Luis Rivera\n")
    h.append(".INST Ecole et Obervatoire des Sciences de la Terre\n")
    h.append(".MAIL 5 rue Rene Descartes, 67084 Strasbourg Cedex, FRANCE\n")
    h.append(".EMAIL luis@sismo.u-strasbg.fr\n")
    h.append(".PHONE +33-390-240047\n")
    h.append(".FAX +33-390-240125\n")
    h.append(".MEDIA FTP\n")
    h.append(".LABEL XXLABELXX\n")
    h.append(".END\n")

    self.headers["luis"]=h
    self.emails["luis"]= email.Utils.formataddr(("Luis Rivera","luis@sismo.u-strasbg.fr"))


    h=[]
    h.append(".NAME Joanne Buckenmeyer\n")
    h.append(".INST Ecole et Obervatoire des Sciences de la Terre\n")
    h.append(".MAIL 5 rue Rene Descartes, 67084 Strasbourg Cedex, FRANCE\n")
    h.append(".EMAIL joanne.buckenmeyer@hotmail.fr\n")
    h.append(".PHONE +33-390-240091\n")
    h.append(".FAX +33-390-240125\n")
    h.append(".MEDIA FTP\n")
    h.append(".LABEL XXLABELXX\n")
    h.append(".END\n")

    self.headers["joanne"]=h
    self.emails["joanne"]= email.Utils.formataddr(("Joanne Buckenmeyer","joanne.buckenmeyer@hotmail.fr"))



    self.emails["iris"]= email.Utils.formataddr(("BREQfast IRIS","breq_fast@iris.washington.edu"))
    self.emails["orfeus"]= email.Utils.formataddr(("BREQfast ORFEUS","breq_fast@knmi.nl"))


dbase=BREQfast_dbase()
#############################################

# Make email

def construct_breqfast_email(username,label,datacenter,request_lines=''):
  try:
    # get the relevant breqfast header
    header=dbase.headers[username]
    # overwrite the label field
    header[BREQfast_dbase.LABEL_FIELD]=".LABEL %s\n"%(label)
    # create the header field
    breqfast_string=''.join(header)
    # create the full email...
    m=email.Message.Message()
    m.add_header("From",dbase.emails[username])
    m.add_header("To", dbase.emails[datacenter]) 
    m.add_header("Bcc",dbase.emails[username])
    m.add_header("Subject", "Data Request")
    m.set_payload(breqfast_string+request_lines)

    return m.as_string()

  except KeyError:
    print "Username %s unkown.  Insert relevant parameters in request_io.py."%(username)
    raise

def send_breqfast_email(username,datacenter,msg_string):
  s=smtplib.SMTP(host="130.79.8.200",port=25)
  s.sendmail(dbase.emails[username],dbase.emails[datacenter],msg_string)
  s.quit()

if __name__== '__main__':
  username="alessia"
  datacenter="orfeus"
#  username="jjl"
  label="testing"
  lines="This is a \n set of lines \n to send\n\n"
  msg_string=construct_breqfast_email(username,label,datacenter,lines)
#  send_breqfast_email(username,msg_string)
  print msg_string
