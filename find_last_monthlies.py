from __future__ import division, print_function

import sys

from sciamachy import get_darkstateid, SdmfDataMissingError
from vardark import read_ch8_darks

n_skipped = 0
petje = 1.0
for orbit in range(51567,53000):
	stateid = get_darkstateid(petje, orbit)
	#print(petje, stateid)
	try:
		jds_, readouts_, noise_, tdet_, pet_, coadd_ = read_ch8_darks((orbit,orbit), stateid)
	except SdmfDataMissingError as e:
		#print("skip orbit", orbit)
		n_skipped += 1
	if jds_.size > 6:
		print(orbit, jds_.size)

		

