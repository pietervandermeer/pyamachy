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
class that plot's evolution of SCIAMACHY's transmission per detector channel.
uses Viewer class to do the actual visualisation. uses config file and 
handles command-line arguments.
This is data from SDMF 3.1. Aim: plot transmission directly without 
intermediate database product.
Author : Pieter van der Meer, SRON, 2011
"""
class TransmissionEvo():

    def __init__(self, pixelmask=None, args=None):

        #
        # Set defaults
        #

        # Indicate data isn't loaded yet.
        self.loaded = False
        
        #
        # Get arguments from command line, and set defaults
        #
        
        parser = argparse.ArgumentParser(description='Displays evolution of SCIAMACHY channel transmission')
        parser.add_argument('-o', '--output', dest='output_fname', type=str)
        parser.add_argument('--config', dest='config_file', type=file, 
                            default='default3.1.cfg')
        parser.add_argument('-v', '--verbose', dest='verbose', 
                            action='store_true')
        parser.add_argument('-V', '--version', action='version', 
                            version='%(prog)s 0.99')
        parser.add_argument('--noscreen', action='store_true')
        # parameters used especially for the plot
        parser.add_argument('--last-orbits', dest='last_orbits', default='0', type=int, help='displays only last 3 months of data')
        parser.add_argument('-l', '--legend', action='store_true', dest='legend', help='displays legend')
        parser.add_argument('-c', '--channel', default='1', help='sets detector channel', choices=sciamachy_module.channels)
        parser.add_argument('-p', '--pixel', default='all', type=str,
           help='sets detector pixel ("all", [0..1023], [0..794] for ch6, [0..229] for ch 6+)')
        # parse
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
        # load all transmission data
        #

        print 'load transmission data'
        ftrans = h5py.File(self.cfg['db_dir']+self.cfg['transmission_fname'], 'r')
        gwls = ftrans["WLStransmission"]
        gsun = ftrans["Transmission"]
        self.wlsorbits = (gwls["orbitList"])[:]
        self.sunorbits = (gsun["orbitList"])[:]
        self.wlstrans = (gwls["transmission"])[:]
        self.suntrans = (gsun["transmission"])[:]

        #
        # load dpm data for selected orbits
        #

#        self.allorbits = numpy.union1d(self.sunorbits, self.wlsorbits)

        print 'load pixel mask data'
# SDMF 3.1 pixelmask: not intended for this use without SNR figures!
#        fmask = h5py.File(self.cfg['db_dir']+self.cfg['pixelmask_fname'], 'r')
# SDMF 3.0 pixelmask: very useful indeed
        fmask = h5py.File('/SCIA/share/SDMF/3.0/sdmf_pixelmask.h5', 'r')
        gmask = fmask['orbitalMask']
        orbitlist = (gmask['orbitList'])[:]
#        idx = numpy.in1d(orbitlist, self.allorbits)
        wlsidx = numpy.in1d(orbitlist, self.wlsorbits)
        sunidx = numpy.in1d(orbitlist, self.sunorbits)
#        self.combined = (gmask['combined'])[:,idx]
        self.wlsmsk = (gmask['combined'])[:,wlsidx]
        self.sunmsk = (gmask['combined'])[:,sunidx]
        fmask.close()

        #
        # reload transmission formatted to missing mask data?
        #

        print 'reload transmission data'
        revsunidx = numpy.in1d(self.sunorbits, orbitlist)
        revwlsidx = numpy.in1d(self.wlsorbits, orbitlist)
        self.wlsorbits = (gwls["orbitList"])[revwlsidx]
        self.sunorbits = (gsun["orbitList"])[revsunidx]
        self.wlstrans = (gwls["transmission"])[:,revwlsidx]
        self.suntrans = (gsun["transmission"])[:,revsunidx]
        
        ftrans.close()

        #
        # Pass control options on to GUIViewer
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
            self.view = GUIViewer(self, title="Transmission", 
                                  image_basename=self.args.output_fname,
                                  params=control_params)
    
    # set arguments
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
        dict['transmission_fname'] = parser.get('Global','transmission_file')
        dict['pixelmask_fname'] = parser.get('Global','pixelmask_official')
        return dict
    

    # loads data and store it for use in 'plot' method.
    def load(self):
        self.loaded = True

    # plot method. this should be about the same as in your stand-alone script.
    # use input object 'fig' as a subplot (axis) object to plot with.
    def display(self, fig):
        
        #
        # load data
        #

        if not self.loaded:
            self.load()
        
        # get data per pixel / channel
        chan_idx = self.get_chan_idx()
        pixrange = sciamachy_module.pixranges[chan_idx] # pixels in sel. channel
        if self.args.pixel == 'all':
            chmsk = self.sunmsk[pixrange,:]
            chunk = numpy.ma.array(self.suntrans[pixrange,:], mask=chmsk)
            suntrans = chunk.mean(0)
            chmsk = self.wlsmsk[pixrange,:]
            chunk = numpy.ma.array(self.wlstrans[pixrange,:], mask=chmsk)
            wlstrans = chunk.mean(0)
        else:
            chanpix = int(self.args.pixel)
            pix = pixrange[chanpix]
            suntrans = self.suntrans[pix,:]
            wlstrans = self.wlstrans[pix,:]

        #
        # plot data
        #

        if self.args.pixel == 'all':
            title = "Transmission of channel "+self.args.channel+"\n\n\n"
        else:
            title = "Transmission of channel "+self.args.channel+\
                    " pixel "+self.args.pixel+ "\n\n\n"

        fig.cla()
        fig.set_title(title)
        fig.set_xlabel("Orbit number")
        fig.set_ylabel("Transmission [%]")
        if self.args.last_orbits > 0:
            ma = max(self.sunorbits)
            fig.set_xlim(ma-self.args.last_orbits,ma)
        fig.scatter(self.sunorbits,suntrans*100,color='r',s=1,label='Sun-over-ESM')
        fig.scatter(self.wlsorbits,wlstrans*100,color='b',s=1,label='WLS')
        if not hasattr(self, 'ax2'):
            self.ax2 = fig.twiny()
        self.ax2.set_xlabel("Date")
        dates = envisat.convert_orbit_to_jd(self.sunorbits)
        self.ax2.plot_date(dates,suntrans,visible=False)
        x1, x2 = fig.get_xlim()
        self.ax2.set_xlim(list(envisat.convert_orbit_to_jd((x1,x2))))
        fig.set_ylim((0,120))
        self.ax2.grid(True)
        #fig.grid(True,which='minor')
        if self.args.legend:
            fig.legend(loc='lower left', scatterpoints=10)
            #fig.legend()

    # dump PNG image and thumbnail
    def save_image(self,basename):
        self.view.show()
        self.view.save_image(basename)

    # get index of channel name
    def get_chan_idx(self):
        chan = self.args.channel
        return (numpy.where(numpy.array(sciamachy_module.channels)==chan))[0]

    #
    # -- GUI functionality, don't call these manually --
    #

    # set parameter
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
                                       'Channel Picker', sciamachy_module.channels)
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

    trans = TransmissionEvo()
    app = TestApp(0)
    app.MainLoop()
