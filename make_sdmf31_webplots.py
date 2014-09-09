"""
dumps evolution plots for all channels of dark current, transmission and so on.
Author : Pieter van der Meer, SRON, 2011
"""

# TODO: if plots seem odd, try '--median' option on each!

import sciamachy_module
import envisat
from plot_darkcurrent31_evo import DarkcurrentEvo
from plot_transmission31_evo import TransmissionEvo

# make channel average dark current plots
def make_darkcurrent_evoplots():

    print 'initialising dark current object'
    dc = DarkcurrentEvo(args=['--noscreen','--filter','--data=dc'])
    yranges = ((0,.2),(0,.15),(0,1.2),(0,1.2),(0,.8),(0,2),(0,3000),(0,1500),(2500,6000))

    print 'creating channel average plots of dark current'
    i=0
    for chan in sciamachy_module.channels:

        print 'channel '+chan

        #
        # full orbit range 
        #

        # add --median for median instead of goodpix channel average
        dc.set_args(['--channel='+chan,'--data=dc'])
        dc.set_yrange(yranges[i])
        dc.save_image('chan'+chan+'_darkcurrent_avg')

        #
        # last months
        #

        # add --median for median instead of goodpix channel average
        dc.set_args(['--channel='+chan, '--last-orbits=1300','--data=dc'])
        dc.set_yrange(yranges[i])
        dc.save_image('chan'+chan+'_darkcurrent_avg_last')

        i+=1

# make channel average analog offset plots
def make_analogoffset_evoplots():

    print 'initialising analog offset object'
    dc = DarkcurrentEvo(args=['--noscreen','--filter','--data=ao'])
    yranges = ((1940,2020),(2840,2900),(3540,3580),(2400,2420),(2590,2640),(3900,4200),(3600,4400),(3000,4000),(3100,3800))

    print 'creating channel average plots of dark current'
    i = 0
    for chan in sciamachy_module.channels:

        print 'channel '+chan

        #
        # full orbit range 
        #

        # add --median for median instead of goodpix channel average
        dc.set_args(['--channel='+chan,'--data=ao'])
        dc.set_yrange(yranges[i])
        dc.save_image('chan'+chan+'_analogoffset_avg')

        #
        # last months
        #

        # add --median for median instead of goodpix channel average
        dc.set_args(['--channel='+chan, '--last-orbits=1300','--data=ao'])
        dc.set_yrange(yranges[i])
        dc.save_image('chan'+chan+'_analogoffset_avg_last')

        i+=1

# make channel average transmission plots
def make_transmission_evoplots():

    print 'initialising transmission object'
    dc = TransmissionEvo(args=['--noscreen'])

    print 'creating channel average plots of dark current'
    for chan in sciamachy_module.channels:

        print 'channel '+chan

        #
        # full orbit range 
        #

        # add --median for median instead of goodpix channel average
        dc.set_args(['--legend', '--channel='+chan])
        dc.save_image('chan'+chan+'_transmission_avg')

        #
        # last month
        #

        # add --median for median instead of goodpix channel average
        dc.set_args(['--legend', '--channel='+chan, '--last-orbits=1300'])
        dc.save_image('chan'+chan+'_transmission_avg_last')

if __name__ == '__main__':
    make_darkcurrent_evoplots()
    make_analogoffset_evoplots()
    make_transmission_evoplots()

