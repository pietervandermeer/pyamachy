"""
channel average transmission (Sun-over-ESM and WLS) evolution plot.

this calculates channel averages from scratch, and is not perfect at this.
however, it is a nice illustration of numpy and h5py power.

Pieter van der Meer, SRON - Netherlands Institute for Space Research, 2011

"""

import numpy as np
import numpy.lib.arraysetops as aso
import h5py
import wx
# toggle between direct pylab and plotpanel. i don't like plotpanel stuff
# since it looks like you need to make all the zoom functionality yourself. 
# that's pretty harsh.

import matplotlib
import pylab as p

from wxPlotPanel import PlotPanel
#from wxxPlotPanel import PlotPanel
import ConfigParser

#-- classes --------------------------------------------------------------------

# pretty useless plot panel using embedded plot in wx. you need to add
# zoom functionality yourself. go figure.
if __name__ == '__main__':
    class PyamachyPlotPanel (PlotPanel):
        import numpy as np
        """Plots several lines in distinct colors."""
        def __init__( self, parent, point_lists, clr_list, **kwargs ):
            self.parent = parent
            self.point_lists = point_lists
            self.clr_list = clr_list

            # initiate plotter
            PlotPanel.__init__( self, parent, **kwargs )
            self.SetColor( (255,255,255) )

        def draw( self ):
            """Draw data."""
            if not hasattr( self, 'subplot' ):
                self.subplot = self.figure.add_subplot( 111 )

            for i, pt_list in enumerate( self.point_lists ):
                plot_pts = np.array( pt_list )
                clr = [float( c )/255. for c in self.clr_list[i]]
                self.subplot.scatter( plot_pts[:,0], plot_pts[:,1], color=clr )

# Create a new frame class, derived from the wxPython Frame.
class MyFrame(wx.Frame):

    def __init__(self, parent, id, title):
        # First, call the base class' __init__ method to create the frame
        wx.Frame.__init__(self, parent, id, title)

        # Associate some events with methods of this class
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_MOVE, self.OnMove)

        # Add a panel and some controls to display the size and position
        panel = wx.Panel(self, -1)
        label1 = wx.StaticText(panel, -1, "Size:")
        label2 = wx.StaticText(panel, -1, "Pos:")
        self.sizeCtrl = wx.TextCtrl(panel, -1, "", style=wx.TE_READONLY)
        self.posCtrl = wx.TextCtrl(panel, -1, "", style=wx.TE_READONLY)
        self.panel = panel

        # Use some sizers for layout of the widgets
        sizer = wx.FlexGridSizer(2, 2, 5, 5)
        sizer.Add(label1)
        sizer.Add(self.sizeCtrl)
        sizer.Add(label2)
        sizer.Add(self.posCtrl)

        border = wx.BoxSizer()
        border.Add(sizer, 0, wx.ALL, 15)
        panel.SetSizerAndFit(border)
        self.Fit()


    # This method is called by the System when the window is resized,
    # because of the association above.
    def OnSize(self, event):
        size = event.GetSize()
        self.sizeCtrl.SetValue("%s, %s" % (size.width, size.height))

        # tell the event system to continue looking for an event handler,
        # so the default handler will get called.
        event.Skip()

    # This method is called by the System when the window is moved,
    # because of the association above.
    def OnMove(self, event):
        pos = event.GetPosition()
        self.posCtrl.SetValue("%s, %s" % (pos.x, pos.y))

# Every wxWidgets application must have a class derived from wx.App
class MyApp(wx.App):

    # wxWindows calls this method to initialize the application
    def OnInit(self):

        # Create an instance of our customized Frame class
        frame = MyFrame(None, -1, "This is a test")
        frame.Show(True)

        # Tell wxWindows that this is our main window
        self.SetTopWindow(frame)

        # Return a success flag
        return True

#-- functions ------------------------------------------------------------------

# this is the mean() function with the twist that it discards invalid values 
# note: assumes float arrays
# usage: 
#  a = [[1,2],[3,np.nan]]
#  print mean_nan(a, axis=0)
def mean_nan(a, axis=None):
    import numpy as np
    c=np.ma.fix_invalid(a)
    return c.mean(axis=axis)

# good pixel mean
def mean_gp(a, mask, axis=None):
    import numpy as np
    c=np.ma.array(a,mask=mask)
    return c.mean(axis=axis)

def get_config(config_fname):
    global db_dir, transmission_fname, mask_fname
    #config_fname = 'default.cfg' # TODO: or use commandline config?
    parser=ConfigParser.SafeConfigParser()
    parser.read([config_fname])
    
    db_dir             = parser.get('Global','masterdirectory')
    transmission_fname = parser.get('Global','transmission_file')
    mask_fname         = parser.get('Global','pixelmask_file')
#    bla_fname          = parser.get('Global','bla_file') # test! hehe

#-- globals --------------------------------------------------------------------

# defaults globals. will be overridden by config file.
db_dir             = "/SCIA/share/SDMF/3.0/"
transmission_fname = "sdmf_transmission.h5"
mask_fname         = "sdmf_pixelmask.h5"
channel            = 8

#-- main code ------------------------------------------------------------------

def main():
    config_fname = 'default.cfg'
    # parse config file
    try:
        get_config(config_fname)
    except ConfigParser.NoOptionError, ex:
        print "There was missing option in the configuration file '"+\
            config_fname+"'!"
        print ex
        return
    
    app = MyApp(0)     # Create an instance of the application class
    app.MainLoop()     # Tell it to start processing events
    
    ch_start = (channel-1)*1024
    ch_end   = channel*1024
    ft = h5py.File(db_dir+transmission_fname, 'r')
    fmask = h5py.File(db_dir+mask_fname, 'r')
    
    #open pixel mask database 
    smoothgroup = fmask["smoothMask"]
    mask_orblist = smoothgroup["orbitList"]
    # get actual pixelmask dataset (not into memory because it's gigantic!)
    masks = smoothgroup["combined"]
    
    #
    # plot sun transmission
    #
    
    sungroup = ft["Transmission"]
    orblist = sungroup["orbitList"]
    trans = sungroup["transmission"]
    
    #
    # find matches between sun data and smoothmask data
    #
    
    print "matching Sun-over-ESM orbits with smoothmask orbits..."
    idx  = np.in1d(mask_orblist, orblist, assume_unique=True)
    ridx = np.in1d(orblist, mask_orblist, assume_unique=True)
    print "done."
    
    print "slicing mask data..."
    ch1masks = masks[ch_start:ch_end, idx]
    print "done."
    
    # good pixel channel mean
    mn = mean_gp(trans[ch_start:ch_end,ridx], ch1masks, axis=0) 
    orblist_ = orblist[ridx]
    
    p.scatter(orblist_,mn,c='r')
    #x0 = orblist_
    #y0 = mn
    
    #
    # plot wls transmission
    #
    
    wlsgroup = ft["WLStransmission"]
    orblist = wlsgroup["orbitList"]
    trans = wlsgroup["transmission"]
    
    #
    # find matches between WLS data and smoothmask data
    #
    
    print "matching WLS orbits with smoothmask orbits..."
    idx  = np.in1d(mask_orblist, orblist, assume_unique=True)
    ridx = np.in1d(orblist, mask_orblist, assume_unique=True)
    print "done."
    
    print "slicing mask data..."
    ch1masks = masks[ch_start:ch_end, idx]
    print "done."
    
    # good pixel channel mean
    mn = mean_gp(trans[ch_start:ch_end,ridx], ch1masks, axis=0) 
    orblist_ = orblist[ridx]
    
    p.scatter(orblist_,mn,c='b')
    #x1 = orblist_
    #y1 = mn
    
    #
    # finish the plot
    #
    
    p.title('Channel evolution of transmission in channel '+str(channel))
    p.show()
    
    #
    # close data files
    #
    
    ft.close()
    fmask.close()
    
    #
    # done processing, plot it (in case we're using PlotPanel stuff)
    #
    
    #points = [[(xi,yi) for xi,yi in zip( x0, y0 )],
            #[(xi,yi) for xi,yi in zip( x1, y1 )]]
    #clrs = [[225,200,160], [219,112,147]]
    #app = wx.PySimpleApp( 0 )
    #frame = wx.Frame( None, wx.ID_ANY, 'Channel evolution of transmission in channel',size=(300,300))
    #panel = PyamachyPlotPanel( frame, points, clrs )
    #frame.Show()
    #app.MainLoop()

if __name__ == "__main__":
    main()
