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

from __future__ import print_function
from __future__ import division

import ConfigParser
import sys
import argparse
import logging
import h5py
import numpy
import wx
import ctypes as ct 
import matplotlib

from viewers import GUIViewer, DumpingViewer
import envisat
import distinct_colours
from sciamachy_module import NonlinCorrector, read_extracted_states
from simudark_module import simudark_orbvar_function, read_simudark, test_simudark_orbvar_function

# Used to guarantee to use at least Wx2.8
import wxversion
#wxversion.ensureMinimal('2.8')
# uncomment the following to use wx rather than wxagg
#matplotlib.use('WX')
#from matplotlib.backends.backend_wx import FigureCanvasWx as FigureCanvas
# comment out the following to use wx rather than wxagg
matplotlib.use('WXAgg')

#- functions -------------------------------------------------------------------

# class that plots evolution of pixel masks (SDMF3.1) per detector channel.
# uses Viewer class to do the actual visualisation. uses config file and handles
# command-line arguments.
class PlotSimudarkCorrection():

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
                 'Displays SCIAMACHY dark state measurements corrected by simudark product per pixel per orbit')
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
        self.args = parser.parse_args()

        # TODO configurable
        self.args.orbit = 24000
        self.args.pixnr = 620

        #
        # Parse config file, exit if unsuccessful
        #
        
        try:
            self.cfg = self.get_config(self.args.config_file)
        except ConfigParser.NoOptionError, ex:
            msg = "There was a missing option in the configuration file '"
            logging.exception(msg+self.args.config_fname+"'!")
            raise

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
        # instantiate non-linearity corrector
        #

        self.nlc = NonlinCorrector()

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
            self.view = GUIViewer(self, title="Dark-corrected darks", 
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
        # load dark state executions (entire ch8) and correct non-linearity
        #
        
        fname = self.cfg['db_dir']+self.cfg['extract_fname']
        orbrange = [self.args.orbit,self.args.orbit]
        states = read_extracted_states(orbrange, 8, fname, readoutMean=True, readoutNoise=True)
        state_mtbl = states['mtbl']
        print(states['readoutMean'].shape)
        readouts = states['readoutMean']
        noise = states['readoutNoise']
        print(readouts.shape)
        for idx_exec in range(readouts.shape[0]):
            readouts[idx_exec,:] = self.nlc.correct(readouts[idx_exec,:])
        pet = states['pet'][0]
        coadd = states['coadd'][0]
        self.readouts = readouts[:,7*1024:8*1024] #/ coadd (was already done?)
        self.noise = noise[:,7*1024:8*1024] / numpy.sqrt(coadd)
        self.state_phases = state_mtbl['orbitPhase'][:]

        #
        # load simudark
        #

        petcorr = 1.18125e-3
        pet -= petcorr
        print("pet=",pet)
        simudark = read_simudark(self.args.orbit, ao=True, lc=True, amp1=True, sig_ao=True, sig_lc=True, sig_amp1=True)
        simudark_mtbl = simudark['mtbl']
        d = {}
        d['phases'] = self.state_phases
        d['pet'] = numpy.array(pet) * numpy.ones(1024) # TODO: get pet from state data!
        d['amp1'] = simudark['amp1'][:]
        d['amp2'] = simudark_mtbl['AMP2'] #[:]
        d['phase1'] = simudark_mtbl['PHASE1'] #[:]
        d['phase2'] = simudark_mtbl['PHASE2'] #[:]
        funk = numpy.array(simudark_orbvar_function(d))
        # add ao and lc offsets to orbital variation
        print("funk.shape=", funk.shape)
        darklevel = simudark['ao']+simudark['lc']*pet
        self.simunoise = simudark['sig_ao'] + simudark['sig_lc']*pet*pet + simudark['sig_amp1']
        print("darklevel.shape=", darklevel.shape)
        for i in range(funk.shape[0]):
            funk[i,:] += darklevel
        self.simudark = funk

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
        
        fig.cla()
        fig.set_title("Simudark correction of dark states, orbit "+str(self.args.orbit)+", pix "+str(self.args.pixnr)+"\n")
        fig.set_xlabel("Orbit phase")
        fig.set_ylabel("Dark signal (BU)")
        pixnr = self.args.pixnr
        print(numpy.sqrt(self.noise[:,pixnr]))
        print(numpy.sqrt(self.simunoise[pixnr]))
        fig.errorbar(self.state_phases, self.readouts[:,pixnr], yerr=numpy.sqrt(self.noise[:,pixnr]),
            ls='none', color=self.cols[0], marker='o', label='State execution')
        fig.errorbar(self.state_phases, self.simudark[:,pixnr], yerr=numpy.sqrt(self.simunoise[pixnr]), 
            ls='none', color=self.cols[1], marker='o', label='Simudark orbvar')
        if self.args.legend:
            fig.legend(loc='upper left', scatterpoints=10)
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

    plotter = PlotSimudarkCorrection()
    app = TestApp(0)
    app.MainLoop()
