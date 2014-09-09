# -*- coding: iso-8859-1 -*-
import ConfigParser
import sys
import argparse
import logging
import h5py
import numpy
import wx

from viewers import GUIViewer, DumpingViewer
import envisat
import sciamachy_module

#- functions -------------------------------------------------------------------

# class that plots transmission per pixel for each detector channel.
# uses Viewer class to do the actual visualisation. uses config file and handles
# command-line arguments.
class PlotWLS():

    def __init__(self):

        #
        # Set defaults
        #
        self.channels = ['1','2','3','4','5','6','6+','7','8']
        self.pixranges = [numpy.arange(1024), \
                          1024+numpy.arange(1024), \
                          1024*2+numpy.arange(1024), \
                          1024*3+numpy.arange(1024), \
                          1024*4+numpy.arange(1024), \
                          1024*5+numpy.arange(795), \
                          1024*5+795+numpy.arange(229), \
                          1024*6+numpy.arange(1024), \
                          1024*7+numpy.arange(1024), \
                         ]

        # Indicate data isn't loaded yet.
        self.loaded = False
        
        #
        # Get arguments from command line, and set defaults
        #
        
        parser = argparse.ArgumentParser(description=
                 'Displays evolution of SCIAMACHY channel pixel evolution')
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
        parser.add_argument('-c', '--channel', default='1', 
                            help='sets detector channel', choices=self.channels)
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
        # initialise correctors for memory and non-linearity
        #

        self.mc = sciamachy_module.MemCorrector()
        self.nlc = sciamachy_module.NonlinCorrector()

        #
        # set default channel
        #

        self.chan_idx = (numpy.where(numpy.array(self.channels) 
                         == self.args.channel))[0]

        #
        # load WLS spectra
        #
        
        fname = self.cfg['db_dir']+self.cfg['extract_fname']
        indexref = 560 #orb 50717
        index49k = 541 #orb 49007
        index47k = 519 #orb 46982
        fex   = h5py.File(fname, 'r')
        group = fex['State_61'] # WLS state
        pixrange = self.pixranges[self.chan_idx] # pixels in sel. channel
        specref = (group['readoutMean'])[indexref,:]
        spec49k = (group['readoutMean'])[index49k,:]
        spec47k = (group['readoutMean'])[index47k,:]
        fex.close()
        
        #
        # load analog offset spectra
        #
        
        fname =  self.cfg['db_dir']+self.cfg['dark_fname']
        fdark = h5py.File(fname, 'r')
        group = fdark['DarkFit/analogOffset']
        aoref = group[50717-1,:]
        ao49k = group[49007-1,:]
        ao47k = group[46982-1,:]
        fdark.close()
        
        #
        # correct for memory and non-linearity
        #
        
        #print "specref.shape=",specref.shape
        specref = self.mc.correct(specref)
        specref = self.nlc.correct(specref)
        spec47k = self.mc.correct(spec47k)
        spec47k = self.nlc.correct(spec47k)
        spec49k = self.mc.correct(spec49k)
        spec49k = self.nlc.correct(spec49k)
        aoref = self.mc.correct(aoref)
        aoref = self.nlc.correct(aoref)
        ao47k = self.mc.correct(ao47k)
        ao47k = self.nlc.correct(ao47k)
        ao49k = self.mc.correct(ao49k)
        ao49k = self.nlc.correct(ao49k)
        
        #
        # compute ao-corrected relative spectra
        #
        
        specref -= aoref
        spec47k -= ao47k
        spec49k -= ao49k
        self.quo47k = spec47k/specref
        self.quo49k = spec49k/specref

        #
        # Pass control options on to GUIViewer
        #
        
        if self.args.noscreen:
            self.view = DumpingViewer(self.args.output_fname, self)
        else:
            control_params = [{'name':'channel','action_descr':'Choose channel',
                              'handler':self.OnChannel,'type':'str'},
                              {'name':'legend','action_descr':'Show legend',
                              'handler':self.OnToggleLegend,'type':'bool',
                              'value':self.args.legend}]
            self.view = GUIViewer(self, title="WLS transmission spectra", 
                                  image_basename=self.args.output_fname,
                                  params=control_params)

    
    # load configuration from file
    def get_config(self, config_file):
        parser=ConfigParser.SafeConfigParser()
        parser.readfp(config_file)
        dict = {}
        dict['db_dir'] = parser.get('Global','masterdirectory')
        dict['extract_fname'] = parser.get('Global','extract_file')
        dict['dark_fname'] = parser.get('Global','dark_file')
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
            print "unknown parameter "+name

    # loads data and store it for use in 'plot' method.
    def load(self):
        # dummy.. do everything at the initialisation as it's really not much
        # data
        self.loaded = True

    # plot method. this should be about the same as in your stand-alone script.
    # use input object 'fig' as a subplot (axis) object to plot with.
    def display(self, fig):
        
        #
        # load data
        #

        if not self.loaded:
            self.load()
        pixrange = self.pixranges[self.chan_idx] # pixels in sel. channel
        
        #
        # plot data
        #
        
        fig.cla()
        fig.set_title("WLS transmission in channel "+self.args.channel+"\n\n\n")
        fig.set_xlabel("Pixel number")
        fig.set_ylabel("WLS spec/ref")
        fig.scatter(pixrange,self.quo47k[pixrange],color='r',s=1,label='47k/ref')
        fig.scatter(pixrange,self.quo49k[pixrange],color='g',s=1,label='49k/ref')
        x1, x2 = fig.get_xlim()
        #fig.set_ylim((0,65535))
        fig.grid(True)
        #fig.grid(True,which='minor')
        if self.args.legend:
            fig.legend(loc='upper left', scatterpoints=10)
            #fig.legend()

    # execute this when Show legend check box is clicked
    def OnToggleLegend(self, event):
        self.args.legend = not self.args.legend
        self.view.refresh_plot()

    # execute this when Options->Choose channel is selected
    def OnChannel(self, event):
        # Create a list of choices
        choices = self.channels
        # Create the dialog
        dialog = wx.SingleChoiceDialog(None, 'Pick Channel...', 
                                       'Channel Picker', choices)
        # Show the dialog
        ret = dialog.ShowModal()
        if ret == wx.ID_OK:
            chan = dialog.GetStringSelection()
            self.set_param('channel', chan)
            self.view.refresh_plot()
        # Destroy the dialog
        dialog.Destroy()

#- main code -------------------------------------------------------------------

if __name__ == '__main__':

    # Every wxWidgets application must have a class derived from wx.App
    class TestApp(wx.App):

        # wxWindows calls this method to initialize the application
        def OnInit(self):
            return True

    trans = PlotWLS()
    app = TestApp(0)
    app.MainLoop()
