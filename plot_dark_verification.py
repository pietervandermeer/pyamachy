from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import h5py

from scia_dark_functions import scia_dark_fun2
from sciamachy import petcorr

pixel = 621
orbit1 = 23213
orbit2 = 23628

plot_x = (
23627.0359079,
23627.0508471,
23627.1287297,
23627.1436695,
23627.2141677,
23627.2291075,
23628.0359149,
23628.0508547,
23628.0703416,
23628.0852814,
23628.1083425,
23628.1232817,
23628.1501586,
23628.165099,
23628.1845859,
23628.1994412,
23628.2189274,
23628.2338679,
23628.4304718,
23628.445411,
23628.4648985,
23628.4798384,
23628.4993252,
23628.514265,
23628.5337518,
23628.5486917,
23628.5681785,
23628.5831177,
23628.6026045,
23628.6175444,
23628.6577049,
23628.6726441,
23628.6921316,
23628.7070708,
23628.7265582,
23628.7414974,
23628.7609849,
23628.7759247,
23628.7954116,
23628.8103514,
23628.8298382,
23628.844778,
23628.8642649,
23628.8792047,
23628.8986909,
23628.9136314,
23628.9331176,
23628.948058,
23628.9675449,
23629.0359057,
23629.0508462,
23629.0739073,
23629.0888465,
23629.1232971,
23629.1382369,
23629.1577244,
23629.1726643,
23629.1956402,
23629.2105807,
23629.2300669,
23629.2450067,
23629.4305034,
23629.4454439,
23629.4649307,
23629.4798705,
23629.4993567,
23629.5142972,
23629.533784,
23629.5487232,
23629.5682107,
23629.5831505,
23629.6026374,
23629.6175772,
23629.6576872,
23629.6726271,
23629.6921133,
23629.7070538,
23629.7265399,
23629.7414798,
23629.7609666,
23629.7759071,
23629.7953939,
23629.8103337,
23629.8298199,
23629.8447604,
23629.8642466,
23629.8791871,
23629.8986732,
23629.9136131,
23629.9331006,
23629.9480404,
23630.0358868,
23630.0508272,
23630.0851542,
23630.1000947,
23630.1268007,
23630.1417405,
23630.1612273,
23630.1761665,
23630.2029595,
23630.2178993,
23630.2373855,
23631.0577685,
23631.072709,
23631.0921952,
23631.107135,
23631.1777176,
23631.1926574)

plot_y = (
5155.36027714,
5156.281324,
5152.10056723,
5151.10209589,
5147.25793281,
5149.23048671,
5155.56168693,
5150.83877849,
5153.03037651,
5152.37561238,
5151.1228498,
5154.78168736,
5155.15691191,
5152.81023293,
5153.83308252,
5146.29484031,
5149.79217638,
5147.73085007,
5128.37723154,
5127.14902238,
5127.14921845,
5127.21656476,
5124.71470205,
5124.71064452,
5127.4493777,
5128.41862337,
5125.6933972,
5127.14902238,
5128.65294786,
5124.10863634,
5129.38036962,
5131.125213,
5134.09492307,
5135.93638408,
5135.24765191,
5133.69769675,
5137.88064495,
5139.24498189,
5140.3630695,
5143.18691189,
5145.93117026,
5146.4955097,
5146.53148876,
5148.19875236,
5149.0921307,
5148.5668094,
5150.41987097,
5155.28287141,
5147.40947901,
5152.95655884,
5152.17787962,
5150.79678104,
5150.94156037,
5148.56465215,
5146.83224274,
5149.69391578,
5149.57211412,
5146.08320532,
5143.2573909,
5140.61483174,
5139.44760904,
5124.39156644,
5125.78251275,
5121.9599832,
5118.26475204,
5115.24013131,
5119.73502243,
5118.87626243,
5118.43213968,
5118.72373851,
5118.03080089,
5121.20714078,
5123.37643777,
5120.85711794,
5120.77067228,
5123.94083871,
5122.90853548,
5128.27897094,
5131.29553727,
5128.90666553,
5133.80145751,
5133.69357006,
5133.60078811,
5135.27453914,
5134.66971971,
5139.361398,
5136.6744559,
5139.78866051,
5140.98542181,
5142.54582366,
5141.85368402,
5141.96945924,
5143.89072308,
5144.12532617,
5144.15795598,
5142.5453348,
5141.25167583,
5142.09460707,
5139.079552,
5135.34933453,
5134.10099328,
5136.12661988,
5142.19531196,
5142.08665629,
5140.03846731,
5141.78614163,
5136.52748403,
5132.86467404,
)

plot_yerr = (
2.66881458748,
3.08130067384,
2.45260464733,
2.83742312152,
2.07002520193,
2.82631972414,
2.7279602857,
2.85443442202,
1.87396075225,
2.92194454066,
2.10671466907,
2.59903238396,
1.98404323699,
3.03944205346,
2.11394299895,
2.82154770807,
2.27235733115,
2.68679035227,
2.3015072029,
2.17882778613,
2.37912450403,
2.5868817623,
1.98150632761,
2.91322022081,
2.28550376362,
2.41319691479,
2.20290239513,
3.0957369306,
1.99274518467,
2.97123505987,
2.65006728544,
2.72960012569,
1.98072195762,
2.69148422084,
2.48805215553,
3.18803127757,
2.13706711416,
3.21176517483,
1.89583109286,
2.15811296334,
2.25532509253,
2.77463554943,
2.22430094411,
2.8016437552,
2.63419131225,
2.59172475638,
3.00831073296,
2.65824811349,
2.42323946769,
2.70571029188,
2.94762423895,
2.28660823104,
2.6100908525,
2.03201333379,
2.78898911802,
2.14819867533,
2.63185481144,
2.51290517878,
2.68254553308,
2.11353362751,
2.99172093678,
2.61593143673,
2.48120507437,
2.08093859081,
2.76054701319,
2.45663248031,
2.51712759162,
2.50449551991,
2.76020789606,
2.16803683885,
2.39793879488,
2.08940529911,
2.53059955127,
2.41383036964,
2.56915924787,
2.09730771933,
2.92262779357,
2.19103348797,
2.93450788627,
2.41901017145,
2.70172584489,
2.34449263558,
2.77067559819,
2.59721063046,
2.72087030922,
2.1275205246,
2.75369058585,
2.72926096497,
3.05818598199,
2.3274802091,
2.73295783932,
2.39013361182,
3.03494367313,
2.63762024557,
2.70388027823,
2.68267712158,
2.44606474907,
2.67289302498,
2.78773589024,
1.87861839662,
2.99060847009,
2.3872735019,
2.43660813036,
2.66882765974,
1.96199601928,
3.21694107248,
2.30920362454,
2.7210538483
)

xnew = (
23627.0,
23627.02,
23627.04,
23627.06,
23627.08,
23627.1,
23627.12,
23627.14,
23627.16,
23627.18,
23627.2,
23627.22,
23627.24,
23627.26,
23627.28,
23627.3,
23627.32,
23627.34,
23627.36,
23627.38,
23627.4,
23627.42,
23627.44,
23627.46,
23627.48,
23627.5,
23627.52,
23627.54,
23627.56,
23627.58,
23627.6,
23627.62,
23627.64,
23627.66,
23627.68,
23627.7,
23627.72,
23627.74,
23627.76,
23627.78,
23627.8,
23627.82,
23627.84,
23627.86,
23627.88,
23627.9,
23627.92,
23627.94,
23627.96,
23627.98,
23628.0,
23628.02,
23628.04,
23628.06,
23628.08,
23628.1,
23628.12,
23628.14,
23628.16,
23628.18,
23628.2,
23628.22,
23628.24,
23628.26,
23628.28,
23628.3,
23628.32,
23628.34,
23628.36,
23628.38,
23628.4,
23628.42,
23628.44,
23628.46,
23628.48,
23628.5,
23628.52,
23628.54,
23628.56,
23628.58,
23628.6,
23628.62,
23628.64,
23628.66,
23628.68,
23628.7,
23628.72,
23628.74,
23628.76,
23628.78,
23628.8,
23628.82,
23628.84,
23628.86,
23628.88,
23628.9,
23628.92,
23628.94,
23628.96,
23628.98,
23629.0,
23629.02,
23629.04,
23629.06,
23629.08,
23629.1,
23629.12,
23629.14,
23629.16,
23629.18,
23629.2,
23629.22,
23629.24,
23629.26,
23629.28,
23629.3,
23629.32,
23629.34,
23629.36,
23629.38,
23629.4,
23629.42,
23629.44,
23629.46,
23629.48,
23629.5,
23629.52,
23629.54,
23629.56,
23629.58,
23629.6,
23629.62,
23629.64,
23629.66,
23629.68,
23629.7,
23629.72,
23629.74,
23629.76,
23629.78,
23629.8,
23629.82,
23629.84,
23629.86,
23629.88,
23629.9,
23629.92,
23629.94,
23629.96,
23629.98,
23630.0,
23630.02,
23630.04,
23630.06,
23630.08,
23630.1,
23630.12,
23630.14,
23630.16,
23630.18,
23630.2,
23630.22,
23630.24,
23630.26,
23630.28,
23630.3,
23630.32,
23630.34,
23630.36,
23630.38,
23630.4,
23630.42,
23630.44,
23630.46,
23630.48,
23630.5,
23630.52,
23630.54,
23630.56,
23630.58,
23630.6,
23630.62,
23630.64,
23630.66,
23630.68,
23630.7,
23630.72,
23630.74,
23630.76,
23630.78,
23630.8,
23630.82,
23630.84,
23630.86,
23630.88,
23630.9,
23630.92,
23630.94,
23630.96,
23630.98,
23631.0
)

trend_y = (
5153.38585077,
5153.38585077,
5153.38585077,
5153.38585077,
5153.38585077,
5153.38585077,
5153.38585077,
5153.36981134,
5153.3537719,
5153.33773247,
5153.32169304,
5153.30565361,
5153.28961418,
5153.27357474,
5153.25753531,
5153.24149588,
5153.22545645,
5153.20941702,
5153.19337758,
5153.17733815,
5153.16129872,
5153.14525929,
5153.12921986,
5153.11318043,
5153.09714099,
5153.08110156,
5153.06506213,
5153.0490227,
5153.03298327,
5153.01694383,
5153.0009044,
5152.98486497,
5152.96882554,
5152.95278611,
5152.93674667,
5152.92070724,
5152.90466781,
5152.88862838,
5152.87258895,
5152.85654952,
5152.84051008,
5152.82447065,
5152.80843122,
5152.79239179,
5152.77635236,
5152.76031292,
5152.74427349,
5152.72823406,
5152.71219463,
5152.6961552,
5152.68011576,
5152.66407633,
5152.6480369,
5152.63199747,
5152.61595804,
5152.59991861,
5152.58387917,
5152.49245631,
5152.40103344,
5152.30961058,
5152.21818771,
5152.12676485,
5152.03534198,
5151.94391911,
5151.85249625,
5151.76107338,
5151.66965052,
5151.57822765,
5151.48680479,
5151.39538192,
5151.30395906,
5151.21253619,
5151.12111332,
5151.02969046,
5150.93826759,
5150.84684473,
5150.75542186,
5150.663999,
5150.57257613,
5150.48115327,
5150.3897304,
5150.29830753,
5150.20688467,
5150.1154618,
5150.02403894,
5149.93261607,
5149.84119321,
5149.74977034,
5149.65834748,
5149.56692461,
5149.47550175,
5149.38407888,
5149.29265601,
5149.20123315,
5149.10981028,
5149.01838742,
5148.92696455,
5148.83554169,
5148.74411882,
5148.65269596,
5148.56127309,
5148.46985022,
5148.37842736,
5148.28700449,
5148.19558163,
5148.10415876,
5148.0127359,
5147.87016532,
5147.72759474,
5147.58502416,
5147.44245359,
5147.29988301,
5147.15731243,
5147.01474186,
5146.87217128,
5146.7296007,
5146.58703012,
5146.44445955,
5146.30188897,
5146.15931839,
5146.01674782,
5145.87417724,
5145.73160666,
5145.58903608,
5145.44646551,
5145.30389493,
5145.16132435,
5145.01875378,
5144.8761832,
5144.73361262,
5144.59104204,
5144.44847147,
5144.30590089,
5144.16333031,
5144.02075974,
5143.87818916,
5143.73561858,
5143.593048,
5143.45047743,
5143.30790685,
5143.16533627,
5143.0227657,
5142.88019512,
5142.73762454,
5142.59505396,
5142.45248339,
5142.30991281,
5142.16734223,
5142.02477166,
5141.88220108,
5141.7396305,
5141.59705992,
5141.45448935,
5141.31191877,
5141.16934819,
5141.02677762,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
5140.88420704,
)

plot_x_ = np.array(plot_x)
plot_y_ = np.array(plot_y)
plot_yerr_ = np.array(plot_yerr)

f = h5py.File("interpolated_monthlies_long.h5")
ds_amps = f["amps"]
ds_amp2 = f["amp2"]
ds_phases = f["phases"]
ds_orbits = f["orbits"]
orbits = ds_orbits[:]
#idx = np.where(orbits == orbit1)[0]
idx = orbits == orbit1
amp1 = ds_amps[idx,pixel]
amp2 = ds_amp2[idx]
phase1 = ds_phases[idx,0]
phase2 = ds_phases[idx,1]

#p = (0, plot_y_.mean(), amp1, 0, phase1, amp2, phase2)
p = (0, 0, amp1, 0, phase1, amp2, phase2)

#x_axis = np.linspace(orbit2-1,orbit2+2,100)
x_axis = np.array(xnew)
l = x_axis.size
pets = np.ones(l) # - petcorr # data points are normalized to BU/s, so set PET = 1.0 s

x = x_axis, pets

wave = scia_dark_fun2(p, x) 
wave = wave + trend_y -11 # add interpolated dc offset (daily+seasonal variation)

plt.cla()
plt.ticklabel_format(useOffset=False)
plt.title("Monthly orbit "+str(orbit2)+" predicted by vardark data of monthly orbit "+str(orbit1)+", pixel "+str(pixel)+", ch8\n\n")
plt.xlabel("Orbit number")
plt.ylabel("Dark current (BU/s)")
plt.errorbar(plot_x, plot_y, yerr=plot_yerr, fmt='bo')
plt.plot(x_axis, wave, 'g-')
plt.show()

