# -*- coding: iso-8859-1 -*-
import ConfigParser
import sys
import argparse
from logging import log
import h5py
import numpy
import wx

from viewers import GUIViewer, DumpingViewer
from envisat import convert_orbit_to_jd

#- functions -------------------------------------------------------------------

class SMRPlotter():

    def __init__(self):
        
        #
        # Indicate data isn't loaded yet.
        #
        
        self.loaded = False
        
        #
        # Get arguments from command line, and set defaults
        #
        
        parser = argparse.ArgumentParser(description=
                                         'Displays evolution of the SCIAMACHY Sun Mean Reference per pixel')
        parser.add_argument('-p', '--pixel', default='0', type=int)
        parser.add_argument('-o', '--output', dest='output_fname', type=str)
        parser.add_argument('--config', dest='config_file', type=file, 
                            default='default.cfg')
        parser.add_argument('-v', '--verbose', dest='verbose', 
                            action='store_true')
        parser.add_argument('-V', '--version', action='version', 
                            version='%(prog)s 0.1')
        parser.add_argument('--noscreen', action='store_true')
        parser.add_argument('--last-orbits', dest='last_orbits', default='0', type=int)
        self.args = parser.parse_args()
        
        #
        # Parse config file, exit if unsuccessful
        #
        
        try:
            self.cfg = self.get_config(self.args.config_file)
        except ConfigParser.NoOptionError, ex:
            msg = "There was a missing option in the configuration file '"
            log.exception(msg+self.args.config_fname+"'!")
            raise
        
        if self.args.noscreen:
            self.view = DumpingViewer(self.args.output_fname, self)
        else:
            control_params = [{'name':'pixel','action_descr':'Choose pixel',
                              'handler':self.OnPixel,'type':'int'}]
            # pass on GUI options: pixel selection
            self.view = GUIViewer(self, title="Sun Mean Reference", 
                                  image_basename=self.args.output_fname,
                                  params=control_params)

    # load configuration from file
    def get_config(self, config_file):
        parser=ConfigParser.SafeConfigParser()
        parser.readfp(config_file)
        dict = {}
        dict['db_dir'] = parser.get('Global','masterdirectory')
        dict['smr_fname'] = parser.get('Global','smr_file')
        return dict
    
    # set parameter
    def set_param(self, name, value):
        if name == 'pixel':
            self.args.pixel = value
            self.load()
        else:
            print "unknown parameter "+name

    # loads data and store it for use in 'plot' method.
    def load(self):
        
        #
        # Would be nice to have data class, but in this simple case that
        # seems like overkill. 
        #
        
        fsmr = h5py.File(self.cfg['db_dir']+self.cfg['smr_fname'], 'r')
        
        if not int(self.args.pixel) in range(8192):
            raise ValueError("Wrong pixel: "+self.args.pixel)
        
        orbits = fsmr["orbitList"]
        smr = fsmr["SMR"]
        n_orbits = orbits.shape[0]

        self.orbits = orbits[0:n_orbits]
        self.smr = smr[0:8192,0:n_orbits]
        
        fsmr.close()
        
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

        pix = int(self.args.pixel)
        fig.cla()
        fig.set_title("Sun Mean Reference of pixel "+str(self.args.pixel)+"\n\n\n")
        fig.set_xlabel("Orbit")
        fig.set_ylabel("Sun Mean Reference [BU]")

        # dual x-axis plot
        # multiple instances mess up the axis (overdraws!) make sure it's a singleton!
        if not hasattr(self, 'ax2'):
            self.ax2 = fig.twiny()
        self.ax2.set_xlabel("Date")
        fig.scatter(self.orbits,self.smr[pix],color='r',s=20)#,edgecolors='black' # nice, but costly for the eps and printer!
        if self.args.last_orbits > 0:
            ma = max(self.orbits)
            fig.set_xlim(ma-self.args.last_orbits,ma)
        self.ax2.plot_date(convert_orbit_to_jd(self.orbits),self.smr[pix],visible=False)
        x1, x2 = fig.get_xlim()
        self.ax2.set_xlim(list(convert_orbit_to_jd((x1,x2))))

        y1, y2 = fig.get_ylim()
        fig.set_ylim(y1,y2)

    # execute this when Options->Choose pixel is selected
    def OnPixel(self, event):
        # Create a list of choices
        choices = [str(x) for x in range(8192)]
        # Create the dialog
        dialog = wx.SingleChoiceDialog(None, 'Pick Pixel...', 
                                       'Pixel Picker', choices)
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

    trans = SMRPlotter()
    app = TestApp(0)
    app.MainLoop()
