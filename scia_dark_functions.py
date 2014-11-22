import numpy as np
from numpy import cos, pi

def scia_dark_fun1(p, x):
    ao = p[0]
    dc = p[1]
    amp1 = p[2] 
    trend = p[3] 
    phase_shift1 = p[4]
    #amp2 = p[5]
    #phase_shift2 = p[6]
    #print(ao,dc,amp1,trend,phase_shift1)
    # Extract exposure information and orbit phase from x
    orbit_phase, pet, coadd = x
    #print(orbit_phase)
    n_x = orbit_phase.size # nr of datapoints

    dark = np.zeros(n_x)
    dark += dc
#    print(dark)
    #print(amp1 * cos(2*pi*(orbit_phase + phase_shift1)))
    dark += amp1 * cos(2*pi*(orbit_phase + phase_shift1))
#    print(dark)
    #print(amp1 * amp2 * cos(4*pi*(orbit_phase + phase_shift2)))
    #dark += amp1 * amp2 * cos(4*pi*(orbit_phase + phase_shift2))
    dark += trend * orbit_phase
    return dark*pet + ao

# no coadding
def scia_dark_fun2(p, x):
    ao = p[0]
    dc = p[1]
    amp1 = p[2] 
    trend = p[3]
    phase_shift1 = p[4]
    amp2 = p[5]
    phase_shift2 = p[6]
    #print(ao,dc,amp1,trend,phase_shift1,amp2,phase_shift2)
    # Extract exposure information and orbit phase from x
    orbit_phase, pet = x
    n_x = orbit_phase.size # nr of datapoints

    dark = np.zeros(n_x)
    dark += dc
    dark += amp1 * cos(2*pi*(orbit_phase + phase_shift1))
    dark += amp1 * amp2 * cos(4*pi*(orbit_phase + phase_shift2))
    dark += trend * orbit_phase
    return dark*pet + ao

# normalized (no ao, trend or pet).. only the wave itself
def scia_dark_fun2n(p, x):
    ao = p[0]
    dc = p[1]
    amp1 = p[2] 
    trend = p[3]
    phase_shift1 = p[4]
    amp2 = p[5]
    phase_shift2 = p[6]
    #print(ao,dc,amp1,trend,phase_shift1,amp2,phase_shift2)
    # Extract exposure information and orbit phase from x
    orbit_phase = x
    n_x = orbit_phase.size # nr of datapoints

    dark = np.zeros(n_x)
    #dark += dc
    dark += amp1 * cos(2*pi*(orbit_phase + phase_shift1))
    dark += amp1 * amp2 * cos(4*pi*(orbit_phase + phase_shift2))
    #dark += trend * orbit_phase
    return dark

# normalized (no ao, trend or pet).. but with dc
def scia_dark_fun2m(p, x):
    ao = p[0]
    dc = p[1]
    amp1 = p[2] 
    trend = p[3]
    phase_shift1 = p[4]
    amp2 = p[5]
    phase_shift2 = p[6]
    #print(ao,dc,amp1,trend,phase_shift1,amp2,phase_shift2)
    # Extract exposure information and orbit phase from x
    orbit_phase = x
    n_x = orbit_phase.size # nr of datapoints

    dark = np.zeros(n_x)
    dark += dc
    dark += amp1 * cos(2*pi*(orbit_phase + phase_shift1))
    dark += amp1 * amp2 * cos(4*pi*(orbit_phase + phase_shift2))
    dark += trend * orbit_phase
    return dark

