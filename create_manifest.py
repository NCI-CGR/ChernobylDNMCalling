#!/usr/bin/python

import csv
import sys
import re
import os
import errno
from collections import defaultdict

triosdir = os.getenv('TRIOSDIR')
bamdir = os.getenv('BAMDIR')

d = defaultdict(dict)
f= open("combined_manifest.txt","r")
fl =f.readlines()[1:]
f.close()
for x in fl:
	cols=x.rstrip().split('\t')
	if(bool(re.search('c', cols[3]))):
		d[cols[2]][cols[1]] = cols[0]

	if(bool(re.search('mo', cols[3]))):
		d[cols[2]]['mom'] = cols[0]
	
	if(bool(re.search('fa', cols[3]))):
                d[cols[2]]['dad'] = cols[0]

for key in (sorted(d.keys())):
	for key2 in (sorted(d[key].keys())):
		if(bool(re.search('c', key2))):
			filename = triosdir + "/" + key2 + "/manifest.txt"
			os.makedirs(os.path.dirname(filename), exist_ok=True)
			f = open(filename, "w")
			wrt_str = "outId\t" + key2 + "\n" + "child\t" + d[key][key2] + "\n" + "mom\t" + d[key]['mom'] + "\n" + "dad\t" + d[key]['dad'] + "\n" + "bamDir\t" + bamdir + "\n" + "dataDir\t" + triosdir + "/" + key2 + "\n"
			f.write(wrt_str)
			f.close()