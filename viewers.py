# -*- coding: iso-8859-1 -*-
"""
Package that features two graph viewers (1 interactive based on WxPython, 
1 non-interactive that writes straight to PNG)
Both Viewers have the same interface (using duck typing).
"""

#!/usr/bin/env python

from __future__ import print_function

__version__ = '1.0'
__all__ = [
    'GUIViewer',
]

# Used to guarantee to use at least Wx2.8
import wxversion
#wxversion.ensureMinimal('2.8')
from numpy import arange, sin, pi
import matplotlib

# uncomment the following to use wx rather than wxagg
#matplotlib.use('WX')
#from matplotlib.backends.backend_wx import FigureCanvasWx as FigureCanvas
# comment out the following to use wx rather than wxagg
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
#from wx import Menu, MenuBar
import wx

#- classes ---------------------------------------------------------------------

#
# Classes don't inherit from a superclass (yet), we just use 'Duck typing'
# TODO: inheritance is better since they share the same save_image method.
#

# TODO: accept some kind of list of GUI options

# Every wxWidgets application must have a class derived from wx.App
class DummyApp(wx.App):
    def OnInit(self):
        return True

# image dimensions in inches
image_dims = (20,12)

class DumpingViewer():
    
    def __init__(self, image_basename, customer):
        self.customer = customer
        self.figure = Figure()
        self.canvas = FigureCanvasAgg(self.figure)
        self.axes = self.figure.add_subplot(111)
        # Call back customer to display it's data on our figure.
        self.show()
        if image_basename != None:
            self.save_image(image_basename)

    def show(self):
        self.customer.display(self.axes)

    def save_image(self, basename):
        # set figure size in inches
        #f = self.axes.gcf()
        self.figure.set_size_inches(image_dims)
        # save as a thumbnail (size in inches x 10 DPI : 10 times 10x6 = 100x60)
        self.figure.savefig(basename+'.thumb.png', dpi=10)
        # save 100 DPI full-size (1000x600)
        self.figure.savefig(basename+'.png', dpi=100)
        # save as EPS
        self.figure.savefig(basename+'.eps')

class GUIViewer(wx.Frame):
    
    def __init__(self, customer, title=None, image_basename=None, params=None):
        
        #
        # Set defaults
        #
        
        self.gridEnable = False
        
        self.customer = customer
        
        #
        # If customer is derived from wx.App we can use it as an application.
        # Otherwise, we have to roll our own. 
        #
        
        is_sub = issubclass(customer.__class__, wx.App)
        if not is_sub:
            app = DummyApp(0)
            print("GUIViewer: Initialised wx dummy app")
            app.MainLoop()
        
        wx.Frame.__init__(self,None,-1,title,size=(550,350))
        favicon = wx.Icon('falcbee.gif', wx.BITMAP_TYPE_ANY, 16, 16)
        wx.Frame.SetIcon(self, favicon)
        self.figure = Figure()
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.axes = self.figure.add_subplot(111)
        # Call back customer to display it's data on our figure.
        self.customer.display(self.axes)
        #self.add_annotation()
        
        if image_basename != None:
            self.save_image(image_basename)
        
        self.SetBackgroundColour(wx.NamedColor("WHITE"))
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)

        #
        # Add menu bar
        #
        
        menu = wx.Menu()
        m_about = menu.Append(wx.ID_ABOUT, "&About", 
                              "More information about this program")
        menu.AppendSeparator()
        m_exit = menu.Append(wx.ID_EXIT, "E&xit", "Terminate the program")

        menuBar = wx.MenuBar()
        menuBar.Append(menu, "&File");

        #
        # Add check box bar
        #

        checkSizer = wx.BoxSizer(wx.HORIZONTAL)

        #
        # Add model controls
        #
        
        if params != None:
            opt_menu = wx.Menu()
            for dict in params:
                if dict['type'] != 'bool':
                    m_item = opt_menu.Append(wx.NewId(), dict['name'], 
                                            dict['action_descr'])
                    self.Bind(wx.EVT_MENU, dict['handler'], m_item)
                else:
                    cbox = wx.CheckBox(self, label=dict['action_descr'])
                    cbox.Bind(wx.EVT_CHECKBOX, dict['handler'])
                    cbox.SetValue(dict['value'])
                    checkSizer.Add(cbox, 0, wx.ALL, 5)
            menuBar.Append(opt_menu, "&Options");

        self.SetMenuBar(menuBar)
        
        self.Bind(wx.EVT_MENU, self.OnAbout, m_about)
        self.Bind(wx.EVT_MENU, self.OnClose, m_exit)

        self.canvas.Bind(wx.EVT_KEY_DOWN, self.onKeyEvent)
        
        #
        # Add toolbar
        #

        self.add_toolbar()  # comment this out for no toolbar

        #
        # Add check boxes
        #

        self.sizer.Add(checkSizer)

        self.SetSizer(self.sizer)
        self.Fit()
        self.Show(True)

    def onToggleGrid(self, event):
        print("GUIViewer.onToggleGrid()")
            
    def OnAbout(self, event):
        dlg = wx.MessageDialog(self, "Really great viewer program",
                               "About Me", wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    def OnClose(self, event):
        self.Close(True)

    def onKeyEvent(self,event=None):
        """ capture , act upon keystroke events"""
        if event == None: return
        #print(str(event))
        #key = event.KeyCode()
        #print(key)
        
    def add_toolbar(self):
        self.toolbar = NavigationToolbar2Wx(self.canvas)
        self.toolbar.Realize()
        if wx.Platform == '__WXMAC__':
            # Mac platform (OSX 10.3, MacPython) does not seem to cope with
            # having a toolbar in a sizer. This work-around gets the buttons
            # back, but at the expense of having the toolbar at the top
            self.SetToolBar(self.toolbar)
        else:
            # On Windows platform, default window size is incorrect, so set
            # toolbar width to figure width.
            tw, th = self.toolbar.GetSizeTuple()
            fw, fh = self.canvas.GetSizeTuple()
            # By adding toolbar in sizer, we are able to put it at the bottom
            # of the frame - so appearance is closer to GTK version.
            # As noted above, doesn't work for Mac.
            self.toolbar.SetSize(wx.Size(fw, th))
            self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        # update the axes menu on the toolbar
        self.toolbar.update()

    def EvtText(self, event):
        frstr = event.GetString()
        print('EvtText: %s\n' % frstr)
        try:
            self.freq = float(frstr)
        except:
            print("wrong value")
            return
        self.axes.cla()
        self.refresh_plot()

    def refresh_plot(self):
        self.customer.display(self.axes)
        self.canvas.draw()
        self.Update()

    def save_image(self, basename):
        print(self.figure)
        # set figure size in inches
        #f = self.axes.gcf()
        self.figure.set_size_inches(image_dims)
        # save as a thumbnail (size in inches x 10 DPI : 10 times 10x6 = 100x60)
        self.figure.savefig(basename+'.thumb.png', dpi=10)
        # save 100 DPI full-size (1000x600)
        self.figure.savefig(basename+'.png', dpi=100)
        # save as EPS
        self.figure.savefig(basename+'.eps')

# experimental validator class.. 
#class Validate_TextCtrl(wx.PyValidator):
    #def __init__ (self):
        #wx.PyValidator.__init__(self)
        #print("init validator")
    #def Clone (self):
        #return Validate_TextCtrl ()
    #def TransferToWindow ( self ):
        #return True
    #def TransferFromWindow ( self ):
        #return True
    #def Validate ( self, Win ):
        #print("Validate")
        #return True

    def add_annotation(self):
        #self.figure.text(48800, 20, unicode('unicode: Institut f\374r Festk\366rperphysik', 'latin-1'))
        self.figure.suptitle(unicode('unicode: Institut f\374r Festk\366rperphysik', 'latin-1'))

#- main code -------------------------------------------------------------------

if __name__ == '__main__':

    # Every wxWidgets application must have a class derived from wx.App
    class TestApp(wx.App):

        # wxWindows calls this method to initialize the application
        def OnInit(self):

            'Create the main window and insert the custom frame'
            frame = GUIViewer(self, title="Example GUIViewer")
            frame.Show(True)

            # Tell wxWindows that this is our main window
            self.SetTopWindow(frame)

            return True

        # Call-back method for displaying your data on the Viewer's canvas
        def display(self, fig):
            import numpy
            fig.plot(numpy.arange(100),numpy.arange(100))

    app = TestApp(0)
    app.MainLoop()
