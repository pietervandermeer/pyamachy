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
class that plot's evolution of SCIAMACHY's dark current per detector channel.
Channel averages are done with good pixels (using SDMF3.0 pixelmask!). Channel
median is an alternative, but probably somewhat biased. That's why I chose a
good pixel average.
uses Viewer class to do the actual visualisation. uses config file and 
handles command-line arguments.
This is data from SDMF 3.1. Aim: plot dark current directly without 
intermediate database product.

Note: when doing full orbit range (>50k orbits) RAM usage is extensive:
4GB RAM is recommended!

Author : Pieter van der Meer, SRON, 2011
"""
class DarkcurrentEvo():

    # NOTE: '--filter' and '--data' settings are *fixed* after object 
    # initialisation! set_args() won't change that!
    def __init__(self, pixelmask=None, args=None):

        # Indicate data isn't loaded yet.
        self.loaded = False
        self.yrange = None
        
        #
        # Get arguments from command line, and set defaults
        #

        parser = argparse.ArgumentParser(description=
               'Displays evolution of SCIAMACHY channel dark current')
        parser.add_argument('-o', '--output', dest='output_fname', type=str)
        parser.add_argument('--config', dest='config_file', type=file, 
                            default='default3.1.cfg')
        parser.add_argument('-v', '--verbose', dest='verbose', 
                            action='store_true')
        parser.add_argument('-V', '--version', action='version', 
                            version='%(prog)s 0.99')
        parser.add_argument('--noscreen', action='store_true')
        # parameters specific to this type of plot
        parser.add_argument('--last-orbits', dest='last_orbits', default='0',
                            type=int, 
                            help='displays only last 3 months of data')
        parser.add_argument('-f', '--filter', action='store_true', 
                            dest='filter', help='filters out bad orbits')
        parser.add_argument('--median', action='store_true', dest='median', 
                            help='use channel median instead of goodpix average')
        parser.add_argument('-l', '--legend', action='store_true', 
                            dest='legend', help='displays legend')
        parser.add_argument('-c', '--channel', default='1', 
                            help='sets detector channel', 
                            choices=sciamachy_module.channels)
        parser.add_argument('--data', default='dc', 
                            help='sets data type (dark current or analog offset)', 
                            choices=['ao','dc'])
        parser.add_argument('-p', '--pixel', default='all', type=str,
           help='sets detector pixel ("all", [0..1023], [0..794] for ch6, [0..229] for ch 6+)')
        self.parser = parser
        if args==None:
            self.set_args(sys.argv[1:])
        else:
            self.set_args(args)

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
        # load all dark current orbits
        #

        print 'load dark orbit list'
        fdark = h5py.File(self.cfg['db_dir']+
                           self.cfg['dark_fname'], 'r')
        gfit = fdark["DarkFit"]
        shape = (gfit["metaTable"]).shape
        n_orbits = shape[0]
        n_orbits = 10000 # less orbits, good for testing!
        self.darkorbits = numpy.arange(n_orbits)

        #
        # load dpm data for selected orbits
        #

        print 'load pixel mask data'
        if pixelmask==None:
            pixelmask = sciamachy_module.pixelmask()
        idx = numpy.in1d(pixelmask.orbits, self.darkorbits)
        self.maskorbits = pixelmask.orbits[idx]
        self.msk = pixelmask.msk[idx,:]

        #
        # load dark current data for good orbits
        # 

        print "load dark current data"
        revidx = numpy.in1d(self.darkorbits, self.maskorbits)
        if self.args.data == 'ao':
            self.dc = (gfit["analogOffset"])[revidx,:]
        elif self.args.data == 'dc':
            self.dc = (gfit["darkCurrent"])[revidx,:]
        self.orbits = self.darkorbits[revidx]

        fdark.close()

        #
        # filter out contaminated and useless orbits
        #

        if self.args.filter:
            print 'apply orbit quality filter'
            orbfilter = sciamachy_module.orbitfilter()
            orbmask = orbfilter.get_quality_orbit_filter(self.orbits)
            self.orbits = self.orbits[orbmask]
            self.dc = self.dc[orbmask,:]
            self.msk = self.msk[orbmask,:]

        #
        # pass control options on to GUIViewer
        #
        
        if self.args.noscreen:
            self.view = DumpingViewer(self.args.output_fname, self)
        else:
            control_params = [{'name':'channel','action_descr':'Choose channel',
                              'handler':self.OnChannel,'type':'str'},
                              {'name':'pixel','action_descr':'Choose pixel(s)',
                              'handler':self.OnPixel,'type':'str'},
                              {'name':'legend','action_descr':'Show legend',
                              'handler':self.OnToggleLegend,'type':'bool',
                              'value':self.args.legend}]
            self.view = GUIViewer(self, title="Dark current", 
                                  image_basename=self.args.output_fname,
                                  params=control_params)

    # parse and check and set arguments specified in array 'args'
    def set_args(self, args):
        self.loaded = False
        self.args = self.parser.parse_args(args)
        # check pixel range for active channel
        if self.args.pixel != 'all':
            try:
                pix = int(self.args.pixel)
            except:
                logging.exception("wrong pixel nr (not 'all', or integer)")
                raise
            if pix < 0:
                logging.exception("negative pixel number not allowed!")
                raise ValueError("negative pixel number not allowed!")
            if self.args.channel == '6+':
                if pix > 229:
                    logging.exception("wrong pixel nr for ch 6+ (>229)!")
                    raise ValueError("wrong pixel nr for ch 6+ (>229)!")
            elif self.args.channel == '6':
                if pix > 794:
                    logging.exception("wrong pixel nr for ch 6+ (>794)!")
                    raise ValueError("wrong pixel nr for ch 6+ (>794)!")
            else:
                if pix > 1023:
                    logging.exception("wrong pixel nr (>1023)!")
                    raise ValueError("wrong pixel nr (>1023)!")
    
    # load configuration from file
    def get_config(self, config_file):
        parser=ConfigParser.SafeConfigParser()
        parser.readfp(config_file)
        dict = {}
        dict['db_dir'] = parser.get('Global','masterdirectory')
        dict['pixelmask_fname'] = parser.get('Global','pixelmask_official')
        dict['dark_fname'] = parser.get('Global','dark_file')
        return dict
    
    # loads data and store it for use in 'plot' method.
    def load(self):

        #
        # load dark formatted to missing mask data
        #

        # get data per pixel / channel
        #print 'slice dark data'
        chan_idx = self.get_chan_idx()
        pixrange = sciamachy_module.pixranges[chan_idx] # pixels in sel. channel
        if self.args.pixel == 'all':
            if self.args.median:
                # if you like the median (smooth, but probably biased)
                self.dc_ = numpy.median(self.dc[:,pixrange], axis=1)
            else:
                # if you like the good pixel mean (accurate, but some variability with time)
                chmsk = self.msk[:,pixrange]
                chunk = numpy.ma.array(self.dc[:,pixrange], mask=chmsk)
                self.dc_ = chunk.mean(1)
        else:
            chanpix = int(self.args.pixel)
            pix = pixrange[chanpix]
            self.dc_ = (self.dc)[:,pix]

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
        # prepare titles
        #

        if self.args.data == 'ao':
            label = 'Analog offset'
            unit = '(BU)'
        elif self.args.data == 'dc':
            label = 'Dark current'
            unit = '(BU/s)'

        if self.args.pixel == 'all':
            title = label+" of channel "+self.args.channel+"\n\n\n"
        else:
            title = label+" of channel "+self.args.channel+\
                    " pixel "+self.args.pixel+ "\n\n\n"

        #
        # plot data
        #

        fig.cla()
        fig.set_title(title)
        fig.set_xlabel("Orbit number")
        fig.set_ylabel(label+" "+unit)
        if self.args.last_orbits > 0:
            ma = max(self.orbits)
            fig.set_xlim(ma-self.args.last_orbits,ma)
        fig.scatter(self.orbits,self.dc_,color='r',s=1,label=label)
        if self.yrange != None:
            fig.set_ylim(self.yrange)
        if not hasattr(self, 'ax2'):
            self.ax2 = fig.twiny()
        self.ax2.set_xlabel("Date")
        dates = envisat.convert_orbit_to_jd(self.orbits)
        self.ax2.plot_date(dates,self.dc_,visible=False)
        x1, x2 = fig.get_xlim()
        self.ax2.set_xlim(list(envisat.convert_orbit_to_jd((x1,x2))))
        self.ax2.grid(True)
        #fig.grid(True,which='minor')

    # dump PNG image and thumbnail
    def save_image(self,basename):
        self.view.show()
        self.view.save_image(basename)

    # get index of channel name
    def get_chan_idx(self):
        chan = self.args.channel
        return (numpy.where(numpy.array(sciamachy_module.channels)==chan))[0]

    def set_yrange(self, yrange):
        self.yrange = yrange

    #
    # -- GUI functionality, don't call these manually --
    #

    # set parameter (from GUI)
    def set_param(self, name, value):
        if name == 'channel':
            self.args.channel = value
            self.load()
        elif name == 'pixel':
            self.args.pixel = value
            self.load()
        else:
            print "unknown parameter "+name

    # execute this when Show legend check box is clicked
    def OnToggleLegend(self, event):
        self.args.legend = not self.args.legend
        self.view.refresh_plot()

    # execute this when Options->Choose channel is selected
    def OnChannel(self, event):
        # Create the dialog
        dialog = wx.SingleChoiceDialog(None, 'Pick Channel...', 
                                       'Channel Picker', 
                                       sciamachy_module.channels)
        # Show the dialog
        ret = dialog.ShowModal()
        if ret == wx.ID_OK:
            chan = dialog.GetStringSelection()
            self.set_param('channel', chan)
            self.view.refresh_plot()
        # Destroy the dialog
        dialog.Destroy()

    # execute this when Options->Choose pixel is selected
    def OnPixel(self, event):
        # Create a list of choices
        if self.args.channel == '6+':
            n_pix = 229
        else:
            n_pix = 1024
        lst = [''] * (n_pix+1)
        lst[0] = 'all'
        for i in range(n_pix):
            lst[i+1] = str(i)
        # Create the dialog
        dialog = wx.SingleChoiceDialog(None, 'Pick pixel...', 
                                       'Pixel Picker', lst)
        # Show the dialog
        ret = dialog.ShowModal()
        if ret == wx.ID_OK:
            pix = dialog.GetStringSelection()
            self.set_param('pixel', pix)
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

    dc = DarkcurrentEvo()
    app = TestApp(0)
    app.MainLoop()

