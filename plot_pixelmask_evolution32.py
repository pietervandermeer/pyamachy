#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Plots evolution of SDMF 3.2 (smooth) pixel mask (channel 8 only).
"""

from __future__ import division, print_function

import ConfigParser
import sys
import argparse
import logging
import h5py
import numpy
import wx

from viewers import GUIViewer, DumpingViewer
import envisat
import sciamachy
import distinct_colours

#- functions -------------------------------------------------------------------

class PlotPixelmaskEvolution():
    """
    class that plots evolution of pixel masks (SDMF3.1) for channel 8.
    uses Viewer class to do the actual visualisation. uses config file and handles
    command-line arguments.
    """

    def __init__(self):

        # Indicate data isn't loaded yet.
        self.loaded = False
        
        #
        # Get arguments from command line, and set defaults
        #
        
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('-o', '--output', dest='output_fname', type=str)
        parser.add_argument('--config', dest='config_file', type=file, 
                            default='default3.2.cfg')
        parser.add_argument('-v', '--verbose', dest='verbose', 
                            action='store_true')
        parser.add_argument('-V', '--version', action='version', 
                            version='%(prog)s 0.1')
        parser.add_argument('--noscreen', action='store_true')
        parser.add_argument('--last-orbits', dest='last_orbits', default='0', 
                            type=int, 
                            help='displays only last 3 months of data')
        # parameters used especially for the GUI
        parser.add_argument('-l', '--legend', action='store_true', 
                            dest='legend', help='displays legend')
        parser.add_argument('-f', '--filter', action='store_true', 
                            dest='filter', help='filters out bad orbits')
        self.args = parser.parse_args()

        #
        # Parse config file, exit if unsuccessful
        #
        
        try:
            self.cfg = self.get_config(self.args.config_file)
        except ConfigParser.NoOptionError, ex:
            msg = "There was a missing option in the configuration file '"
            logging.exception(msg+self.args.config_fname+"'!")
            raise

        #
        # initialise orbit filter to filter out contaminated and useless orbits
        #

        #ideal_orbits = numpy.arange(46000)+5000
        #c = numpy.in1d(ideal_orbits, self.orbits)
        #print 'number of orbits present (ideally 46000) = ', c.sum()

        self.orbfilter = sciamachy.orbitfilter()
        fname = self.cfg['db_dir']+self.cfg['pixelmask_fname']
        fmask = h5py.File(fname, 'r')
        self.orbits = (fmask['orbits'])[:]
        fmask.close()

        self.mask = self.orbfilter.get_monthly_orbit_filter(self.orbits)
        self.mask &= self.orbfilter.get_quality_orbit_filter(self.orbits, margin=100)

        #
        # get 8 distinct colour-blind friendly colours for scatter plot
        #

        self.cols = distinct_colours.get_distinct(8)

        #
        # Pass control options on to GUIViewer
        #
        
        if self.args.noscreen:
            self.view = DumpingViewer(self.args.output_fname, self)
        else:
            control_params = [{'name':'legend','action_descr':'Show legend',
                              'handler':self.OnToggleLegend,'type':'bool',
                              'value':self.args.legend},
                              {'name':'filter',
                              'action_descr':'Filter bad orbits',
                              'handler':self.OnToggleFilter,
                              'value':self.args.filter,'type':'bool'}]
            self.view = GUIViewer(self, title="Pixel quality evolution", 
                                  image_basename=self.args.output_fname,
                                  params=control_params)

    
    # load configuration from file
    def get_config(self, config_file):
        parser=ConfigParser.SafeConfigParser()
        parser.readfp(config_file)
        dict = {}
        dict['db_dir'] = parser.get('Global','masterdirectory')
        dict['pixelmask_fname'] = parser.get('Global','smoothmask_file')
        dict['monthlies_fname'] = parser.get('Global','monthlies_file')
        return dict
    
    # set parameter
    def set_param(self, name, value):
        if name == 'channel':
            #print value
            self.chan_idx = (numpy.where(numpy.array(self.channels)==value))[0]
            #print numpy.where(self.channels == value)
            self.args.channel = value
            self.load()
        else:
            print("unknown parameter "+name)

    # loads data and store it for use in 'plot' method.
    def load(self):
        
        #
        # Would be nice to have data class, but in this simple case that
        # seems like overkill. 
        #
        
        fname = self.cfg['db_dir']+self.cfg['pixelmask_fname']
        fmask = h5py.File(fname, 'r')
        
        chanmasks_inv = (fmask['invalid'])[:,:]
        chanmasks_sat = (fmask['saturation'])[:,:]
        chanmasks_err = (fmask['darkError'])[:,:]
        chanmasks_res = (fmask['darkResidual'])[:,:]
        chanmasks_wls = (fmask['wlsResponse'])[:,:]
        chanmasks_sun = (fmask['sunResponse'])[:,:]
        chanmasks_com = (fmask['combined'])[:,:]
        chanmasks_comf = (fmask['combinedFlag'])[:,:]
        
        self.totalbad_inv = chanmasks_inv.sum(axis=1) / 1024
        self.totalbad_sat = chanmasks_sat.sum(axis=1) / 1024
        self.totalbad_err = chanmasks_err.sum(axis=1) / 1024
        self.totalbad_res = chanmasks_res.sum(axis=1) / 1024
        self.totalbad_wls = chanmasks_wls.sum(axis=1) / 1024
        self.totalbad_sun = chanmasks_sun.sum(axis=1) / 1024
        self.totalbad_com = chanmasks_com.sum(axis=1) / 1024

        self.totalbad_comf = chanmasks_comf.sum(axis=1) / 1024
        
        fmask.close()
        
        self.loaded = True

    # plot method. this should be about the same as in your stand-alone script.
    # use input object 'fig' as a subplot (axis) object to plot with.
    def display(self, fig):
        
        #
        # load data
        #

        if not self.loaded:
            self.load()
        
        if self.args.filter:
            orbits_ = self.orbits[self.mask]
        else:
            orbits_ = self.orbits
        
        #
        # plot data
        #
        
        fig.cla()
        fig.set_title("Evolution of SDMF 3.2 pixel quality in channel 8\n\n")
        fig.set_xlabel("Orbit number")
        fig.set_ylabel("Pixel quality (channel mean)")
        if self.args.last_orbits > 0:
            ma = max(orbits_)
            fig.set_xlim(ma-self.args.last_orbits,ma)
        if self.args.filter:
            fig.plot(orbits_,self.totalbad_comf[self.mask],color=self.cols[0],
                        label='Bad pixel count (combined flag)') #,marker='.',ls='none',markersize=1
#            fig.plot(orbits_,self.totalbad_com[self.mask],color=self.cols[1],
#                        label='Total pixel quality (combined)')
            fig.plot(orbits_,self.totalbad_wls[self.mask],color=self.cols[2],
                        label='Total pixel quality (WLS response)')
            fig.plot(orbits_,self.totalbad_sun[self.mask],color=self.cols[3],
                        label='Total pixel quality (Sun response)')
            fig.plot(orbits_,self.totalbad_inv[self.mask],color=self.cols[4],
                        label='Total pixel quality (invalid)')
#            fig.plot(orbits_,self.totalbad_err[self.mask],color=self.cols[5],
#                        label='Total pixel quality (dark error)')
            fig.plot(orbits_,self.totalbad_res[self.mask],color=self.cols[6],
                        label='Total pixel quality (dark residual)')
            fig.plot(orbits_,self.totalbad_sat[self.mask],color=self.cols[7],
                        label='Total pixel quality (saturated)')
        else:
            print(orbits_.size,self.totalbad_comf.size)
            fig.plot(orbits_,self.totalbad_comf,color=self.cols[0],
                        label='Bad pixel count (combined flag)')
#            fig.plot(orbits_,self.totalbad_com,color=self.cols[1],
#                        label='Total pixel quality (combined)')
            fig.plot(orbits_,self.totalbad_wls,color=self.cols[2],
                        label='Total pixel quality (WLS response)')
            fig.plot(orbits_,self.totalbad_sun,color=self.cols[3],
                        label='Total pixel quality (Sun response)')
            fig.plot(orbits_,self.totalbad_inv,color=self.cols[4],
                        label='Total pixel quality (invalid)')
#            fig.plot(orbits_,self.totalbad_err,color=self.cols[5],
#                        label='Total pixel quality (dark error)')
            fig.plot(orbits_,self.totalbad_res,color=self.cols[6],
                        label='Total pixel quality (dark residual)')
            fig.plot(orbits_,self.totalbad_sat,color=self.cols[7],
                        label='Total pixel quality (saturated)')
        if not hasattr(self, 'ax2'):
            self.ax2 = fig.twiny()
        if not hasattr(self, 'ay2'):
            self.ay2 = fig.twinx()
        self.ax2.set_xlabel("Date")
        self.ay2.set_ylabel("Bad pixel count")
        self.ay2.set_ylim((0,1024))
        dates = envisat.convert_orbit_to_jd(self.orbits)
        self.ax2.plot_date(dates,self.totalbad_com,visible=False)
        x1, x2 = fig.get_xlim()
        self.ax2.set_xlim(list(envisat.convert_orbit_to_jd((x1,x2))))
        fig.set_ylim((0,1.0))
        self.ax2.grid(True)
        #fig.grid(True,which='minor')
        if self.args.legend:
            fig.legend(loc='lower right', scatterpoints=10)
            #fig.legend()

    # execute this when Show legend check box is clicked
    def OnToggleLegend(self, event):
        self.args.legend = not self.args.legend
        self.view.refresh_plot()

    # execute this when Filter orbits check box is clicked
    def OnToggleFilter(self, event):
        self.args.filter = not self.args.filter
        self.view.refresh_plot()

#- main code -------------------------------------------------------------------

if __name__ == '__main__':

    # Every wxWidgets application must have a class derived from wx.App
    class TestApp(wx.App):

        # wxWindows calls this method to initialize the application
        def OnInit(self):
            return True

    trans = PlotPixelmaskEvolution()
    app = TestApp(0)
    app.MainLoop()
