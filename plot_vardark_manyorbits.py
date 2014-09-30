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

class MVarDarkPlotter():

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
        parser.add_argument('-O', '--orbit', dest='orbit', type=int, help='sets orbit number', default=24200)
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
        self.old_monthly = -1
        self.orb_window = 20

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

        if self.old_monthly != monthly_orbit:
            (lst) = fit_monthly(monthly_orbit, shortFlag=use_short_states, longFlag=use_long_states)
            channel_phase, channel_phase2, aos, lcs_fit, amps, channel_amp2, trends_fit = lst
            self.aos = aos
            self.lcs_fit = lcs_fit
            self.amps = amps
            self.channel_amp2 = channel_amp2
            self.trends_fit = trends_fit
            self.channel_phase = channel_phase
            self.channel_phase2 = channel_phase2
            self.old_monthly = monthly_orbit
        else:
            aos = self.aos
            lcs_fit = self.lcs_fit
            amps = self.amps
            channel_amp2 = self.channel_amp2
            trends_fit = self.trends_fit
            channel_phase = self.channel_phase
            channel_phase2 = self.channel_phase2
        print('channel_phase=', channel_phase)
        print('channel_phase2=', channel_phase2)
        print('aos=', aos)
        print('lc=', lcs_fit)
        print('amp=', amps)
        print('trend=', trends_fit)

        print('ao=', aos[pixnr], 'lc=', lcs_fit[pixnr], 'amp=', amps[pixnr], 'trend=', trends_fit[pixnr])

        self.trends_lin = list()
        self.lcs_lin = list()
        self.readouts_list = list()
        self.sigmas_list = list()
        self.phases_list = list()
        for i_orb in range(self.orb_window):
            print("window orbit", i_orb)
            orb = normal_orbit+i_orb
            (lst) = extract_dark_states(orb, shortFlag=use_short_states, longFlag=use_long_states)
            n_exec, self.polar_phases, readout_pet, readout_coadd, self.readouts, self.sigmas, self.readout_phases = lst
            self.phases_list.append(self.readout_phases + float(i_orb))
            self.readouts_list.append(self.readouts)
            self.sigmas_list.append(self.sigmas)

            # directly compute constant part of lc and trend for averaged eclipse data points
            trends, lcs = compute_trend(orb, aos, amps, channel_amp2, channel_phase, channel_phase2, shortFlag=use_short_states, longFlag=use_long_states)
            self.trends_lin.append(trends)
            self.lcs_lin.append(lcs)

        dump_pixel()

        self.loaded = True

    def dump_pixel(self):
        pixnr = self.args.pixnr
        for i_orb in range(self.orb_window):
            phi, jds, readouts, sigmas, tdet = extract_two_dark_states_(orbit, stateid1, stateid2)
            # TODO!
            print(phi, rp)

    # plot method. this should be about the same as in your stand-alone script.
    # use input object 'fig' as a subplot (axis) object to plot with.
    def display(self, fig):
        
        #
        # load and prepare data
        #

        if not self.loaded:
            self.load()

        pixnr = self.args.pixnr

        pfit = self.aos[pixnr], self.lcs_fit[pixnr], self.amps[pixnr], self.trends_fit[pixnr], self.channel_phase, self.channel_amp2, self.channel_phase2

        #
        # plot data
        #

        pts_per_orbit = 50
        n_orbits = 1
        total_pts = n_orbits*pts_per_orbit+1
        orbphase = .12 + (numpy.arange(total_pts)/float(pts_per_orbit))
        coadd = numpy.ones(total_pts)
        pets = numpy.array([1/2., 1]) - petcorr
        cols = ['b','g','r','y','k','m','#ff00ff','#ffff00']
        
        fig.cla()
        fig.set_title("Vardark correction of dark states, orbit "+str(self.args.orbit)+", pix "+str(self.args.pixnr)+"\n")
        fig.set_xlabel("Orbit phase")
        fig.set_ylabel("Dark signal (BU)")
        fig.set_xlim([0,self.orb_window+1])

        for i_orb in range(self.orb_window):

            for i_pet in range(len(pets)):

                lc = (self.lcs_lin[i_orb])[pixnr]
                trend = (self.trends_lin[i_orb])[pixnr]
                plin = self.aos[pixnr], lc, self.amps[pixnr], trend, self.channel_phase, self.channel_amp2, self.channel_phase2
                ph = orbphase+float(i_orb)
                x = orbphase, numpy.zeros(total_pts)+pets[i_pet] #, coadd

                if self.readout_phases.size < 20: 
                    # not a montly calib. orbit, linear computed trend, please
                    fig.plot(ph, scia_dark_fun2(plin, x), c=cols[i_pet], marker='+', label="lin orbvar "+str(pets[i_pet]))
                else:
                    fig.plot(ph, scia_dark_fun2(pfit, x), c=cols[i_pet], label="fit orbvar "+str(pets[i_pet]))

            ph = (self.phases_list[i_orb])
            rd = (self.readouts_list[i_orb])[:,pixnr]
            si = (self.sigmas_list[i_orb])[:,pixnr]
            fig.errorbar(ph, rd, yerr=si, ls='none', marker='o', label="dark states")

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

    plotter = MVarDarkPlotter()
    app = TestApp(0)
    app.MainLoop()
