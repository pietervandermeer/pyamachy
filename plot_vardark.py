# -*- coding: iso-8859-1 -*-
#
# COPYRIGHT (c) 2014 SRON (pieter.van.der.meer@sron.nl)
#
#   This is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License, version 2, as
#   published by the Free Software Foundation.
#
#   The software is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, 
#   Boston, MA  02111-1307, USA.
#

"""
plots computed vardark product and overplots state executions.
"""

from __future__ import print_function
from __future__ import division

import ConfigParser
import sys
import argparse
import logging
import h5py
import numpy
import wx
import matplotlib

from viewers import GUIViewer, DumpingViewer
import envisat
import distinct_colours
from sciamachy_module import NonlinCorrector, read_extracted_states, petcorr, orbitfilter
from vardark_module import extract_dark_states, fit_monthly, fit_eclipse_orbit, compute_trend
from scia_dark_functions import scia_dark_fun1, scia_dark_fun2

# Used to guarantee to use at least Wx2.8
import wxversion
#wxversion.ensureMinimal('2.8')
# uncomment the following to use wx rather than wxagg
#matplotlib.use('WX')
#from matplotlib.backends.backend_wx import FigureCanvasWx as FigureCanvas
# comment out the following to use wx rather than wxagg
matplotlib.use('WXAgg')

#- functions -------------------------------------------------------------------

class VarDarkPlotter():

    def __init__(self):

        #
        # Set defaults
        #

        # Indicate data isn't loaded yet.
        self.loaded = False
        
        #
        # Get arguments from command line, and set defaults
        #
        
        parser = argparse.ArgumentParser(description=
                 'Displays SCIAMACHY dark state measurements corrected by vardark product per pixel per orbit')
        parser.add_argument('-o', '--output', dest='output_fname', type=str)
        parser.add_argument('--config', dest='config_file', type=file, 
                            default='default3.1.cfg')
        parser.add_argument('-v', '--verbose', dest='verbose', 
                            action='store_true')
        parser.add_argument('-V', '--version', action='version', 
                            version='%(prog)s 0.1')
        parser.add_argument('--noscreen', action='store_true')
        # parameters used especially for the GUI
        parser.add_argument('-l', '--legend', action='store_true', 
                            dest='legend', help='displays legend')
        parser.add_argument('-O', '--orbit', dest='orbit', type=int, help='sets orbit number', default=24044)
        parser.add_argument('-P', '--pixnr', dest='pixnr', type=int, help='sets pixel number', default=620)
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
        # extract orbit list and turn into strings (for the gui)
        #

        fname = self.cfg['db_dir']+self.cfg['extract_fname']
        fextract = h5py.File(fname, 'r')
        group = fextract['State_08']
        orbits = (numpy.unique(((group['orbitList'])[:]))).tolist()
        self.orbits = map(str, orbits)
        fextract.close()

        #
        # get 6 distinct colour-blind friendly colours for scatter plot
        #

        self.cols = distinct_colours.get_distinct(6)

        #
        # instantiate helper objects
        #

        self.nlc = NonlinCorrector()
        self.ofilt = orbitfilter()

        #
        # Pass control options on to GUIViewer
        #
        
        if self.args.noscreen:
            self.view = DumpingViewer(self.args.output_fname, self)
        else:
            control_params = [{'name':'orbit','action_descr':'Choose Orbit',
                              'handler':self.OnOrbit,'type':'int'},
                              {'name':'legend','action_descr':'Show legend',
                              'handler':self.OnToggleLegend,'type':'bool',
                              'value':self.args.legend},
                             ]
            self.view = GUIViewer(self, title="Vardark-corrected darks", 
                                  image_basename=self.args.output_fname,
                                  params=control_params)
    
    # load configuration from file
    def get_config(self, config_file):
        parser=ConfigParser.SafeConfigParser()
        parser.readfp(config_file)
        dict = {}
        dict['db_dir'] = parser.get('Global','masterdirectory')
        dict['extract_fname'] = parser.get('Global','extract_file')
        dict['monthlies_fname'] = parser.get('Global','monthlies_file')
        return dict
    
    # set parameter
    def set_param(self, name, value):
        if name == 'pixel':
            self.args.pixel = value
        elif name == 'orbit':
            self.args.orbit = value
            self.load()
            self.view.refresh_plot()
        else:
            print("unknown parameter "+name)

    # loads data and store it for use in 'plot' method.
    def load(self):

        #
        # compute vardark
        #

        normal_orbit = self.args.orbit
        monthly_orbit = self.ofilt.get_closest_monthly(self.args.orbit)
        use_short_states = False
        use_long_states = True
        pixnr = 590

        channel_phase, channel_phase2, aos, lcs, amps, channel_amp2, trends = fit_monthly(monthly_orbit, shortFlag=use_short_states, longFlag=use_long_states)
        self.aos = aos
        self.lcs = lcs
        self.amps = amps
        self.channel_amp2 = channel_amp2
        self.trends = trends
        self.channel_phase = channel_phase
        self.channel_phase2 = channel_phase2
        print('channel_phase=', channel_phase)
        print('channel_phase2=', channel_phase)
        print('aos=', aos)
        print('lc=', lcs)
        print('amp=', amps)
        print('trend=', trends)

        #plt.cla()
        #plt.scatter(numpy.arange(1024), aos, c='b')
        #plt.scatter(numpy.arange(1024), lcs, c='g')
        #plt.show()

        print('ao=', aos[pixnr], 'lc=', lcs[pixnr], 'amp=', amps[pixnr], 'trend=', trends[pixnr])

        # fit constant part of lc and trend
        x, lcs_fit, trends_fit, readouts, sigmas = fit_eclipse_orbit(normal_orbit, aos, lcs, amps, channel_phase, shortFlag=use_short_states, longFlag=use_long_states)
        #trends_fit = numpy.zeros(n_pix) # just to illustrate difference
        self.lcs_fit = lcs_fit
        self.trends_fit = trends_fit
        self.readouts = readouts
        self.sigmas = sigmas 

        self.readout_phases, readout_pets, readout_coadd = x

        # directly compute constant part of lc and trend for averaged eclipse data points
        trends_lin, lcs_lin = compute_trend(normal_orbit, aos, amps, channel_phase, shortFlag=use_short_states, longFlag=use_long_states)
        self.trends_lin = trends_lin
        self.lcs_lin = lcs_lin

        self.loaded = True

    # plot method. this should be about the same as in your stand-alone script.
    # use input object 'fig' as a subplot (axis) object to plot with.
    def display(self, fig):
        
        #
        # load data
        #

        if not self.loaded:
            self.load()

        #
        # plot data
        #

        # print(trends_ecl.size, trends_ecl)
        # plt.cla()
        # plt.scatter(numpy.arange(n_pix), trends_ecl, c='b')
        # plt.scatter(numpy.arange(n_pix), lcs_ecl, c='g')
        # plt.show()

        pixnr = self.args.pixnr

        pfit = self.aos[pixnr], self.lcs_fit[pixnr], self.amps[pixnr], self.trends_fit[pixnr], self.channel_phase, self.channel_amp2, self.channel_phase2
        plin = self.aos[pixnr], self.lcs_lin[pixnr], self.amps[pixnr], self.trends_lin[pixnr], self.channel_phase, self.channel_amp2, self.channel_phase2
        #plin = aos[pixnr], lcs_lin[pixnr], amps[pixnr], 0, channel_phase

        pts_per_orbit = 50
        n_orbits = 2
        total_pts = n_orbits*pts_per_orbit
        orbphase = numpy.arange(total_pts)/float(pts_per_orbit)
        coadd = numpy.ones(total_pts)
        pets = numpy.array([1/16., 1/8., 1/2., 1]) - petcorr
        cols = ['b','g','r','y','k','m','#ff00ff','#ffff00']
        
        fig.cla()
        fig.set_title("Vardark correction of dark states, orbit "+str(self.args.orbit)+", pix "+str(self.args.pixnr)+"\n")
        fig.set_xlabel("Orbit phase")
        fig.set_ylabel("Dark signal (BU)")
        fig.set_xlim([0,2])
        #fig.set_ylim([cons-amp, cons+amp])
        #print(numpy.sqrt(self.noise[:,pixnr]))
        #print(numpy.sqrt(self.simunoise[pixnr]))
        for i in range(len(pets)):
            x = orbphase, numpy.zeros(total_pts)+pets[i] #, coadd
            fig.plot(orbphase, scia_dark_fun2(pfit, x), c=cols[i], label="fit orbvar "+str(pets[i]))
            fig.plot(orbphase, scia_dark_fun2(plin, x), c=cols[i], marker='+', label="lin orbvar "+str(pets[i]))
        fig.errorbar(self.readout_phases, self.readouts[:,pixnr], yerr=self.sigmas[:,pixnr], ls='none', marker='o', label="dark states")

        if self.args.legend:
            fig.legend(loc='upper right', scatterpoints=10)
            #fig.legend()

    # execute this when Show legend check box is clicked
    def OnToggleLegend(self, event):
        self.args.legend = not self.args.legend
        self.view.refresh_plot()

    # execute this when Options->Choose orbit is selected
    def OnOrbit(self, event):
        # Create a list of choices
        choices = self.orbits
        # Create the dialog
        dialog = wx.SingleChoiceDialog(None, 'Pick Orbit...', 'Orbit nr', choices)
        # Show the dialog
        ret = dialog.ShowModal()
        if ret == wx.ID_OK:
            orbit = dialog.GetStringSelection()
            self.set_param('orbit', int(orbit))
        else:
            print("orbit not accepted")
        # Destroy the dialog
        dialog.Destroy()

    # skip through pixels and orbits using cursor keys (first focus the canvas)
    def OnKey(self, event):
        key = event.GetKeyCode()
        print(key)
        if key == 314:
            self.args.pixnr -= 1
            if self.args.pixnr < 0:
                self.args.pixnr = 0
            self.view.refresh_plot()
        elif key == 316:
            self.args.pixnr += 1
            if self.args.pixnr > 1023:
                self.args.pixnr = 1023
            self.view.refresh_plot()
        if key == 315:
            self.args.orbit += 1
            if self.args.orbit < 0:
                self.args.orbit = 0
            self.loaded = False
            self.view.refresh_plot()
        elif key == 317:
            self.args.orbit -= 1
            if self.args.orbit > 60000:
                self.args.orbit = 60000
            self.loaded = False
            self.view.refresh_plot()

#- main ------------------------------------------------------------------------

#parser = argparse.ArgumentParser()
#parser.add_argument("pixnr")
#args = parser.parse_args()
#print(args.pixnr)
#check_eclipse_calib(int(args.pixnr))

if __name__ == '__main__':

    # Every wxWidgets application must have a class derived from wx.App
    class TestApp(wx.App):

        # wxWindows calls this method to initialize the application
        def OnInit(self):
            return True

    plotter = VarDarkPlotter()
    app = TestApp(0)
    app.MainLoop()
