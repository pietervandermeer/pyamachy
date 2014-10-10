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
plots computed vardark product and overplots state executions. also has a mode to show residuals ('R' key)
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
from sciamachy_module import NonlinCorrector, petcorr, orbitfilter
from vardark_module import fit_monthly, fit_eclipse_orbit, compute_trend, get_darkstateid, trending_phase, AllDarks
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
        parser.add_argument('-O', '--orbit', dest='orbit', type=int, help='sets orbit number', default=35466)
        parser.add_argument('-P', '--pixnr', dest='pixnr', type=int, help='sets pixel number', default=597)
        self.args = parser.parse_args()
        self.residual_mode = False

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
        self.orb_window = 14

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
        # load initial orbit window
        #

        orbit_range = [self.args.orbit, self.args.orbit+self.orb_window]
        alldarks.lump(orbit_range)
        alldarks.finalize()

        #
        # initialize GUI kdb state
        #

        self.shift = False

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

    def dump_pixel(self):
        pixnr = self.args.pixnr
        orbit = self.args.orbit
        for i_orb in range(self.orb_window):
            orb = orbit + i_orb
            phi, jds, readouts, sigmas, tdet, coadd = extract_two_dark_states_(orb, 8, 63) # TODO: correct state ids for entire mission
            # TODO!
            phi += i_orb
            for i in range(phi.size):
                print(phi[i], readouts[i][pixnr]-self.aos[pixnr], sigmas[i][pixnr], tdet[i])

    # loads data and store it for use in 'plot' method.
    def load(self):

        #
        # compute vardark
        #

        normal_orbit = self.args.orbit
        self.monthly_orbit = self.ofilt.get_closest_monthly(self.args.orbit)
        use_short_states = False
        use_long_states = True
        pixnr = self.args.pixnr

        if self.old_monthly != self.monthly_orbit:
            new_monthly_range = [self.monthly_orbit-self.orb_window/2, self.monthly_orbit+self.orb_window/2]
            alldarks.get_range(new_monthly_range) # autolump new range into buffers
            (lst) = fit_monthly(alldarks, self.monthly_orbit, shortFlag=use_short_states, longFlag=use_long_states, verbose=self.args.verbose)
            channel_phase, channel_phase2, aos, lcs_fit, amps, channel_amp2, trends_fit = lst
            self.aos = aos
            self.lcs_fit = lcs_fit
            self.amps = amps
            self.channel_amp2 = channel_amp2
            self.trends_fit = trends_fit
            self.channel_phase = channel_phase
            self.channel_phase2 = channel_phase2
            self.old_monthly = self.monthly_orbit
        else:
            aos = self.aos
            lcs_fit = self.lcs_fit
            amps = self.amps
            channel_amp2 = self.channel_amp2
            trends_fit = self.trends_fit
            channel_phase = self.channel_phase
            channel_phase2 = self.channel_phase2
        if self.args.verbose:
            print('channel_phase=', channel_phase)
            print('channel_phase2=', channel_phase2)
            print('aos=', aos)
            print('lc=', lcs_fit)
            print('amp=', amps)
            print('amp2=', channel_amp2)
            print('trend=', trends_fit)
            print('ao=', aos[pixnr], 'lc=', lcs_fit[pixnr], 'amp=', amps[pixnr], 'trend=', trends_fit[pixnr])

        self.trends_lin = list()
        self.lcs_lin = list()
        self.readouts_list = list()
        self.sigmas_list = list()
        self.phases_list = list()
        self.rd_pets = list()
        for i_orb in range(self.orb_window):
            #print("window orbit", i_orb)
            orb = normal_orbit+i_orb

            #(lst) = extract_two_dark_states_(orb, stateid1=s1, stateid2=s2)
            #self.readout_phases, jds, self.readouts, self.sigmas, tdet, readout_pet, readout_coadd = lst
            (lst) = alldarks.get_range([orb,orb+.999])
            n_exec, dummy, readout_pet, readout_coadd, self.readouts, self.sigmas, self.readout_phases = lst

            self.phases_list.append(self.readout_phases)
            self.readouts_list.append(self.readouts)
            self.sigmas_list.append(self.sigmas)
            self.rd_pets.append(readout_pet)

            # directly compute constant part of lc and trend for averaged eclipse data points
            trends, lcs = compute_trend(alldarks, orb, aos, amps, channel_amp2, channel_phase, channel_phase2, shortFlag=use_short_states, longFlag=use_long_states)
            self.trends_lin.append(trends)
            self.lcs_lin.append(lcs)

        #self.dump_pixel()

        self.loaded = True

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
        orbphase = (numpy.arange(total_pts)/float(pts_per_orbit))
        coadd = numpy.ones(total_pts)
        pets = numpy.array([1/2., 1]) - petcorr
        cols = ['b','g','r','y','k','m','#ff00ff','#ffff00']

        fig.cla()
        fig.set_xlim([0+self.args.orbit, self.orb_window+self.args.orbit])
        fig.set_xlabel("Orbit")
        fig.ticklabel_format(useOffset=False)

        if self.residual_mode:

            fig.set_title("Residuals of vardark correction, pix "+str(self.args.pixnr)+"\n")
            fig.set_ylabel("Dark residual (BU)")
            fig.plot([0,100000], [0,0], label="zero")

            for i_orb in range(self.orb_window):

                # dark states
                ph_st = (self.phases_list[i_orb])
                rd = (self.readouts_list[i_orb])[:,pixnr]
                si = (self.sigmas_list[i_orb])[:,pixnr]
                pe = (self.rd_pets[i_orb])

                # vardark model
                lc = (self.lcs_lin[i_orb])[pixnr]
                trend = (self.trends_lin[i_orb])[pixnr]
                plin = self.aos[pixnr], lc, self.amps[pixnr], trend, self.channel_phase, self.channel_amp2, self.channel_phase2
                x_model = ph_st + trending_phase -float(i_orb), pe
                #print(self.channel_phase, ph_st, rd, pe, scia_dark_fun2(plin, x_model))
                x_plot = ph_st #+ self.args.orbit
                fig.errorbar(x_plot, rd-scia_dark_fun2(plin, x_model), yerr=si, ls='none', marker='o', c=cols[0]) #marker='+', 

        else:

            fig.set_title("Vardark correction of dark states, pix "+str(self.args.pixnr)+"\n")
            fig.set_ylabel("Dark signal (BU)")

            for i_orb in range(self.orb_window):

                # dark states
                ph_st = (self.phases_list[i_orb])
                rd = (self.readouts_list[i_orb])[:,pixnr]
                si = (self.sigmas_list[i_orb])[:,pixnr]
                #+self.args.orbit
                fig.errorbar(ph_st, rd, yerr=si, ls='none', marker='o', label="dark states")

                # vardark model
                for i_pet in range(len(pets)):

                    lc = (self.lcs_lin[i_orb])[pixnr]
                    trend = (self.trends_lin[i_orb])[pixnr]
                    plin = self.aos[pixnr], lc, self.amps[pixnr], trend, self.channel_phase, self.channel_amp2, self.channel_phase2
                    ph = orbphase+float(i_orb)  + trending_phase
                    x_model = orbphase+ trending_phase, numpy.zeros(total_pts)+pets[i_pet] #, coadd

                    x_plot = ph + self.args.orbit  
                    fig.plot(x_plot, scia_dark_fun2(plin, x_model), c=cols[i_pet], linewidth=2.0, label="lin orbvar "+str(pets[i_pet])) #marker='+', 
                    #fig.plot(ph, scia_dark_fun2(pfit, x_model), c=cols[i_pet], label="fit orbvar "+str(pets[i_pet]))

        if self.args.legend:
            fig.legend(loc='upper right', scatterpoints=10)

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
        if key == 314: #<- : pixnr down
            self.args.pixnr -= 1
            if self.args.pixnr < 0:
                self.args.pixnr = 0
            self.view.refresh_plot()
        elif key == 316: #-> : pixnr up
            self.args.pixnr += 1
            if self.args.pixnr > 1023:
                self.args.pixnr = 1023
            self.view.refresh_plot()
        elif key == 315: #arrow up : orbit up
            self.args.orbit += 1
            if self.args.orbit < 0:
                self.args.orbit = 0
            self.loaded = False
            self.view.refresh_plot()
        elif key == 78: # 'n'ext monthly
            print("curr monthly=",self.monthly_orbit)
            new_monthly = self.ofilt.get_next_monthly(self.monthly_orbit)
            self.args.orbit = new_monthly - (self.orb_window/2)
            print("ORBIT=", self.args.orbit)
            self.loaded = False
            self.view.refresh_plot()
        elif key == 80: # 'p'rev monthly
            new_monthly = self.ofilt.get_previous_monthly(self.monthly_orbit)
            self.args.orbit = new_monthly - (self.orb_window/2)
            print("ORBIT=", self.args.orbit)
            self.loaded = False
            self.view.refresh_plot()
        elif key == 317: #arrow down : orbit down
            self.args.orbit -= 1
            if self.args.orbit > 60000:
                self.args.orbit = 60000
            self.loaded = False
            self.view.refresh_plot()
        elif key == 82: #R(esiduals) (toggle)
            self.residual_mode = not self.residual_mode
            self.view.refresh_plot()
        elif key == 306: 
            self.shift = not self.shift

#- main ------------------------------------------------------------------------

# initialize without any dark data, just specify PET list
alldarks = AllDarks([1.0, 0.5])

#print("import darks state def 1")
#alldarks = AllDarks([5750,43362])
#alldarks = AllDarks([5750,43362])
#alldarks = AllDarks([23000,25000])
#print("import darks state def 2")
#alldarks.lump([43362,53200])
#alldarks.lump([43362,43400])
#print("finalize")
#alldarks.finalize()

if __name__ == '__main__':

    # Every wxWidgets application must have a class derived from wx.App
    class TestApp(wx.App):

        # wxWindows calls this method to initialize the application
        def OnInit(self):
            return True

    plotter = MVarDarkPlotter()
    app = TestApp(0)
    app.MainLoop()
