#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 29 11:54:39 2018

@author: anin
"""

total = 0


def get_sec(time_str):
    h, m, s = time_str.split(':')
    return float(h) * 3600 + float(m) * 60 + float(s)

counter = 0
with open('time_info.txt', 'r') as inp, open('total_time.txt', 'w') as outp:
   for line in inp:
       try:
           num = get_sec(line.strip())
           total += num
           counter += 1
           
       except ValueError:
           print('{} is not a number!'.format(line))
           
   outp.write('#Generations %s: %s seconds'%(counter,total))

print('{} Generations took: {} seconds'.format(counter,total))