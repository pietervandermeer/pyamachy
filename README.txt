---------------------------------------
Pyamachy product generators and viewers
---------------------------------------

Author: Pieter van der Meer, SRON - Netherlands Institute for Space Research, 
        2011

-------
Purpose
-------

To have a replacement of the IDL SDMF code, especially the viewers, but also
some of the important end products (like the pixel mask). The implementation
in Python should require less code and should be more readable, whilst not 
requiring any licenses like IDL does (not even for development). Also,
Python libraries are more stable than IDL's, which is a very good reason to
make the conversion. 

However, there is one most important reason for making the leap to Python. This
is the amazing plot functionality. It is much more modern than IDL's, not unlike
Matlab's, so has features likes stable GUI-assisted zoom and image dump and 
anti-aliased plotting with modern readable fonts as standard. Also, since it is
so close to Matlab, it should seem familiar to larger audience. 

---------------------------------------
The generators, viewers and other files
---------------------------------------

- generators:

pixelmask_module.py: generates SDMF 3.1 pixelmask (100% identical values, but
 different structure from SDMF 3.1 IDL database)

- viewers:

plot_pixelmask_evolution: viewer for SDMF3.1 pixelmask evolution.

plot_transmission30_evo.py: SDMF3.0 transmission product viewer as evolution 
over years.

plot_smr_date.py: Sun Mean Reference evolution viewer (per pixel)

plot_transmission31_evo.py: SDMF3.1 transmission product viewer as evolution 
over years. Works per channel average/median or per pixel.

plot_darkcurrent31_evo.py: SDMF3.1 dark current evolution viewer (ao or dc).
Works per channel average/median or per pixel. Since this is a lot of data
(Giga Bytes!), initialisation may take a while.

orbit31_viewer.py: SDMF3.1 product viewer. Per orbit, multiple products like
dark current, PPG, transmission, mask.

wls31_test.py: experimental viewer to investigate transmission gain in UV-VIS
in the recent orbits (>50000).

- config files:

default3.1.cfg: configuration file for SDMF 3.1 code
default.cfg: configuration file for SDMF 3.0 code

- modules:

darklimb_io_module.py: module that contains darklimb database reader function
read_statedark_module.py: module that contains statedark db reader function
viewers.py: module that contains viewer classes GUIViewer and DumpingViewer
sciamachy_module.py: module that contains instrument parameters and low-level
 calibration routines.
envisat.py: module that contains data and methods related to the Envisat 
 platform: currently only conversion of orbit nr -> date.

--------------------------
Usage
--------------------------

Everything is written in Python 2.6.6. All of the above programs may be run by 
using:

> python pyprogram.py [parameters]

Using the --help parameter the program will display all the parameters it
will accept with explanation and valid values.

> python pyprogram.py --help

Viewers are fully usable as off-line PNG file generators, in which case they
don't use a screen. This also dumps a thumbnail version! For instance:

> python aviewer.py --noscreen -o outfile.png

Viewer parameters may be specified on the prompt as well as configured in
their GUI (for instance, channel and pixel number, show legend on/off, etc).

> python plot_transmission30_evo.py --channel=8 --legend

This sets the transmission evolution viewer to detector channel 8 and enables
the legend display.

> python make_sdmf31_webplots.py

Makes all the webplots.. This is currently only evolution plots.

----------------
Generator Design
----------------

A generator generates an SDMF product orbit-by-orbit. It contains a 'calculate'
method for this purpose. The 'write' method writes the computed product to
HDF5. 

In case of a product that uses a sliding window of input data (for instance 50
orbits), 'calculate_first' to ingest the first window of data, followed by 
repeated 'calculate_next' calls to slide the window 1 orbit (add 1 new, remove
1 old) is recommended. Since an object can keep all the orbit data to itself,
it hides the complexity from the user and yet is efficient because it need only
read data once when calculating from orbit to orbit.

Typically, a generator contains a main function which instantiates the object
and may handle commandline arguments or run a unit test. 

For an example, look at pixelmask_module.py.

-------------
Viewer Design
-------------

The viewers have been designed around a generic Viewer class. The Viewer class
is implemented as either a 'DumpingViewer' (dumps PNG files to disk) or a 
'GUIViewer' (interactive display in a customisable window). Both of these 
classes have the same API and use "duck typing". 

Upon initialisation, a Viewer class takes a reference to the 'customer': the
class that owns the Viewer. The customer should have methods 'load' and 
'display'. Where 'load' loads the data applicable to the settings from the 
viewer (the user may have changed settings in the GUI) and 'display' runs the
plot code. 

The plot code is designed to use Python's Matplotlib library. This is done to
facilitate going from a stand-alone quick-and-dirty test plot, to the Viewer
environment, with customisable commandline and GUI, and automated image dumps.
Typically the programmer has only to copy an example 'customer' class and 
replace 'load', 'display' and initialisation methods. 

He can then add his own custom parameters (such as pixel number, channel, orbit 
range, etc). Most of the parametrisation can be done by using the standard 
ArgumentParser class from Python. Custom GUI handlers, however, still need to 
be added, but this is a relatively painless affair. 

Configuration file handling is done with Python's standard ConfigParser.

For an example, look at the combination of plot_pixelmask_evolution.py and 
viewers.py.
