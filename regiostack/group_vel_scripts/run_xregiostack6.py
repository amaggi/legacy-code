#!/usr/bin/env python
# encoding: utf-8

from os import system
from dsp2intomodes import period2key

periods = [20, 50, 100, 150]

for per in periods:
    key = period2key(per)
    fname = "input_%s.txt" % key
    system("xregiostack6 < %s" % fname)

