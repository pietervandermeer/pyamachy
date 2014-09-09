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

#- functions -------------------------------------------------------------------

# class that plot's evolution of SCIAMACHY's transmission per detector channel.
# uses Viewer class to do the actual visualisation. uses config file and handles
# command-line arguments.
class PlotTransmission():

    def __init__(self):

        #
        # Set defaults
        #
        self.channels = ['1','2','3','4','5','6','6+','7','8']

        # Indicate data isn't loaded yet.
        self.loaded = False
        
        #
        # Get arguments from command line, and set defaults
        #
        
        parser = argparse.ArgumentParser(description='Displays evolution of SCIAMACHY channel transmission')
        parser.add_argument('-o', '--output', dest='output_fname', type=str)
        parser.add_argument('--config', dest='config_file', type=file, 
                            default='default.cfg')
        parser.add_argument('-v', '--verbose', dest='verbose', 
                            action='store_true')
        parser.add_argument('-V', '--version', action='version', 
                            version='%(prog)s 0.1')
        parser.add_argument('--noscreen', action='store_true')
        parser.add_argument('--last-orbits', dest='last_orbits', default='0', type=int, help='displays only last 3 months of data')
        # parameters used especially for the GUI
        parser.add_argument('-l', '--legend', action='store_true', dest='legend', help='displays legend')
        parser.add_argument('-c', '--channel', default='1', help='sets detector channel', choices=self.channels)
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
        # Pass control options on to GUIViewer
        #
        
        if self.args.noscreen:
            self.view = DumpingViewer(self.args.output_fname, self)
        else:
            control_params = [{'name':'channel','action_descr':'Choose channel',
                              'handler':self.OnChannel,'type':'str'},
                              {'name':'legend','action_descr':'Show legend',
                              'handler':self.OnToggleLegend,'type':'bool','value':self.args.legend}]
            self.view = GUIViewer(self, title="Transmission", 
                                  image_basename=self.args.output_fname,
                                  params=control_params)
    
    # load configuration from file
    def get_config(self, config_file):
        parser=ConfigParser.SafeConfigParser()
        parser.readfp(config_file)
        dict = {}
        dict['db_dir'] = parser.get('Global','masterdirectory')
        dict['chanev_fname'] = parser.get('Global','chanev_file')
        return dict
    
    # set parameter
    def set_param(self, name, value):
        if name == 'channel':
            self.args.channel = value
            self.load()
        else:
            print "unknown parameter "+name

    # loads data and store it for use in 'plot' method.
    def load(self):
        
        #
        # Would be nice to have data class, but in this simple case that
        # seems like overkill. 
        #
        
        fchan = h5py.File(self.cfg['db_dir']+self.cfg['chanev_fname'], 'r')
        
        chanev = fchan['Channel'+self.args.channel]
        
        self.orbits   = chanev["ORBIT"]
        self.suntrans = chanev["TRANSMISSION_GP"]
        self.wlstrans = chanev["WLSTRANSMISSION_GP"]
        idx = numpy.where(self.suntrans > 0.0)
        self.suntrans = self.suntrans[idx]
        self.sunorbits = self.orbits[idx]
        idx = numpy.where(self.wlstrans > 0.0)
        self.wlsorbits = self.orbits[idx]
        self.wlstrans = self.wlstrans[idx]
        
        fchan.close()
        
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
        fig.set_title("Transmission of channel "+self.args.channel+"\n\n\n")
        fig.set_xlabel("Orbit number")
        fig.set_ylabel("Transmission [%]")
        if self.args.last_orbits > 0:
            ma = max(self.orbits)
            fig.set_xlim(ma-self.args.last_orbits,ma)
        fig.scatter(self.sunorbits,self.suntrans*100,color='r',s=1,label='Sun-over-ESM')
        fig.scatter(self.wlsorbits,self.wlstrans*100,color='b',s=1,label='WLS')
        if not hasattr(self, 'ax2'):
            self.ax2 = fig.twiny()
        self.ax2.set_xlabel("Date")
        dates = envisat.convert_orbit_to_jd(self.sunorbits)
        self.ax2.plot_date(dates,self.suntrans*100,visible=False)
        x1, x2 = fig.get_xlim()
        self.ax2.set_xlim(list(envisat.convert_orbit_to_jd((x1,x2))))
        fig.set_ylim((0,120))
        self.ax2.grid(True)
        #fig.grid(True,which='minor')
        if self.args.legend:
            fig.legend(loc='lower left')
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

    trans = PlotTransmission()
    app = TestApp(0)
    app.MainLoop()
