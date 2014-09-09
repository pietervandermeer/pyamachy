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

"""
orbit viewer for SDMF3.1 products.
class that plots various types of SDMF data per pixel for each detector channel.
plots a single orbit of data. chooses daily/weekly calibration orbit closest to 
the user-specified uses Viewer class to do the actual visualisation. uses config 
file and handles command-line arguments.

author: Pieter van der Meer, SRON, 2011
"""
class OrbitViewer():

    def __init__(self):

        # Indicate data isn't loaded yet.
        self.loaded = False
        self.datatypes = ['ao','dc','suntrans','wlstrans','ppg','mask']        
        self.datanames = ['Analog offset','Dark current','Sun-over-ESM transmission','WLS transmission','Pixel-to-pixel gain','Dead pixel mask']
        self.units = ['(BU)','(BU/s)','[%]','[%]','','']

        #
        # Get arguments from command line, and set defaults
        #
        
        parser = argparse.ArgumentParser(description=
                 'Displays SDMF3.1 products per channel')
        parser.add_argument('-o', '--output', dest='output_fname', type=str)
        parser.add_argument('--config', dest='config_file', type=file, 
                            default='default3.1.cfg')
        parser.add_argument('-v', '--verbose', dest='verbose', 
                            action='store_true')
        parser.add_argument('-V', '--version', action='version', 
                            version='%(prog)s 0.1')
        parser.add_argument('--noscreen', action='store_true')
        parser.add_argument('--orbit', default='19820', 
                            help='sets orbit number', type=int)
        # parameters used especially for the plot
        parser.add_argument('--data', action='store', type=str, default='ao',
                            help='data type to display', choices=self.datatypes)
        parser.add_argument('-l', '--legend', action='store_true', 
                            dest='legend', help='displays legend')
        parser.add_argument('-c', '--channel', default='1', 
                            help='sets detector channel', choices=sciamachy_module.channels)
        self.args = parser.parse_args()
        data = self.args.data

        #
        # Parse config file, exit if unsuccessful
        #
        
        try:
            self.cfg = self.get_config(self.args.config_file)
        except ConfigParser.NoOptionError, ex:
            msg = "There was a missing option in the configuration file '"
            logging.exception(msg+str(self.args.config_file)+"'!")
            raise

        #
        # set default channel
        #

        self.chan_idx = (numpy.where(numpy.array(sciamachy_module.channels) 
                         == self.args.channel))[0]


        #
        # set plot parameters for user-specified data type
        #

        idx = (numpy.where(numpy.array(self.datatypes) == data))[0]
        self.dataname = self.datanames[idx]
        self.dataunit = self.units[idx]

        #
        # load spectra of user-specified type
        #

        if data == 'ppg':
            fname = self.cfg['db_dir']+self.cfg['ppg_fname']
            fdat = h5py.File(fname, 'r')
            self.orbits = fdat['orbitList'][:]
            dset = fdat['pixelGain']
        elif data == 'dc':
            fname = self.cfg['db_dir']+self.cfg['dark_fname']
            fdat = h5py.File(fname, 'r')
            n_orbits = len(fdat['DarkFit/metaTable'])
            self.orbits = numpy.arange(n_orbits, dtype=int)
            dset = fdat['DarkFit/darkCurrent']
        elif data == 'ao':
            fname = self.cfg['db_dir']+self.cfg['dark_fname']
            fdat = h5py.File(fname, 'r')
            n_orbits = len(fdat['DarkFit/metaTable'])
            self.orbits = numpy.arange(n_orbits, dtype=int)
            dset = fdat['DarkFit/analogOffset']
        elif data == 'suntrans':
            fname = self.cfg['db_dir']+self.cfg['transmission_fname']
            fdat = h5py.File(fname, 'r')
            self.orbits = fdat['Transmission/orbitList'][:]
            dset = fdat['Transmission/transmission']
        elif data == 'wlstrans':
            fname = self.cfg['db_dir']+self.cfg['transmission_fname']
            fdat = h5py.File(fname, 'r')
            self.orbits = fdat['WLStransmission/orbitList'][:]
            dset = fdat['WLStransmission/transmission']
        elif data == 'mask':
            fname = self.cfg['db_dir']+self.cfg['pixelmask_fname']
            fdat = h5py.File(fname, 'r')
            self.orbits = fdat['orbitalMask/orbitList'][:]
            criteria = sciamachy_module.mask_criteria
            dset0 = fdat['orbitalMask/'+criteria[0]]
            dset1 = fdat['orbitalMask/'+criteria[1]]
            dset2 = fdat['orbitalMask/'+criteria[2]]
            dset3 = fdat['orbitalMask/'+criteria[3]]
            dset4 = fdat['orbitalMask/'+criteria[4]]
            dset5 = fdat['orbitalMask/'+criteria[5]]
        else:
            raise ValueError("unknown data type!")

        # find closest index
        dist = numpy.abs(self.orbits - self.args.orbit)
        indices = numpy.argsort(dist)
#        idx = numpy.where(self.orbits == self.args.orbit)
        idx = indices[0]
        self.orbit = self.orbits[idx]
#        print "idx=",idx
        if data == 'dc' or data == 'ao':
            self.dat = dset[idx,:]
        elif data == 'mask':
            self.dat1 = dset0[:,idx].astype(float)
            self.dat2 = dset1[:,idx].astype(float)
            self.dat3 = dset2[:,idx].astype(float)
            self.dat4 = dset3[:,idx].astype(float)
            self.dat5 = dset4[:,idx].astype(float)
            self.dat6 = dset5[:,idx].astype(float)
            idx_ = numpy.where(self.dat1 == 0)
            self.dat1[idx_] = numpy.nan
            idx_ = numpy.where(self.dat2 == 0)
            self.dat2[idx_] = numpy.nan
            idx_ = numpy.where(self.dat3 == 0)
            self.dat3[idx_] = numpy.nan
            idx_ = numpy.where(self.dat4 == 0)
            self.dat4[idx_] = numpy.nan
            idx_ = numpy.where(self.dat5 == 0)
            self.dat5[idx_] = numpy.nan
            idx_ = numpy.where(self.dat6 == 0)
            self.dat6[idx_] = numpy.nan
        else:
            self.dat = dset[:,idx]
        fdat.close()

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
            self.view = GUIViewer(self, title="PPG spectra", 
                                  image_basename=self.args.output_fname,
                                  params=control_params)

    
    # load configuration from file
    def get_config(self, config_file):
        parser=ConfigParser.SafeConfigParser()
        parser.readfp(config_file)
        dict = {}
        dict['db_dir'] = parser.get('Global','masterdirectory')
        dict['ppg_fname'] = parser.get('Global','ppg_file')
        dict['dark_fname'] = parser.get('Global','dark_file')
        dict['pixelmask_fname'] = parser.get('Global','pixelmask_official')
        dict['transmission_fname'] = parser.get('Global','transmission_file')
        return dict
    
    # set parameter
    def set_param(self, name, value):
        if name == 'channel':
            #print value
            self.chan_idx = (numpy.where(numpy.array(sciamachy_module.channels)==value))[0]
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
        pixrange = sciamachy_module.pixranges[self.chan_idx] # pixels in sel. channel
        
        #
        # plot data
        #

        ylabel = self.args.data.upper()+" "+self.dataunit

        fig.cla()
        fig.set_title(self.dataname+" in channel "+self.args.channel
                      +", orbit "+str(self.orbit)+"\n")
        fig.set_xlabel("Pixel number")
        fig.set_ylabel(ylabel)
        if self.args.data == 'mask':
            criteria = sciamachy_module.mask_criteria
            fig.set_yticks((6,5,4,3,2,1))
            labels = fig.set_yticklabels(criteria)
            fig.scatter(pixrange,self.dat1[pixrange]+5,color='k',s=1,label=criteria[0])
            fig.scatter(pixrange,self.dat2[pixrange]+4,color='r',s=1,label=criteria[1])
            fig.scatter(pixrange,self.dat3[pixrange]+3,color='g',s=1,label=criteria[2])
            fig.scatter(pixrange,self.dat4[pixrange]+2,color='b',s=1,label=criteria[3])
            fig.scatter(pixrange,self.dat5[pixrange]+1,color='m',s=1,label=criteria[4])
            fig.scatter(pixrange,self.dat6[pixrange]+0,color='y',s=1,label=criteria[5])
        else:
            fig.scatter(pixrange,self.dat[pixrange],color='r',s=1,label=ylabel)
        #x1, x2 = fig.get_xlim()
#        fig.set_xlim((0,1024))
        fig.set_xlim((pixrange.min(),pixrange.max()))
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
        choices = sciamachy_module.channels
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

    view = OrbitViewer()
    app = TestApp(0)
    app.MainLoop()
