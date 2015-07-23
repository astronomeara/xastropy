"""
#;+ 
#; NAME:
#; ginga
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for Ginga stuff
#;   Jul-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import gzip, os
import subprocess, glob

from astropy import units as u
from astropy.io import ascii, fits
from astropy.table import QTable, Table, Column

from ginga.util import grc

from xastropy.xutils import xdebug as xdb

# def bintab_to_table(fits_fil,exten=1, silent=False):
# def table_to_fits(table, outfil, compress=False, comment=None):

#
def show_fits(fits_file, host='localhost', port=9000, **kwargs):
    ''' Read a binary FITS file into an active Ginga window

    Parameters
    ---------
    fits_file: str
      Filename  
    host: str, optional
      Name of the host; Default='localhost'
    port: int, optional
      Value of the port; Default=9000
    ''' 
    # Check for file
    if not os.path.isfile(fits_file) == 0:
        raise ValueError('File={:s} not found!'.format(fits_file))
    # Connect to ginga RC
    client = grc.RemoteClient(host, port)
    # Connect to ginga widget
    method = getattr(client, 'ginga')
    # Set up args 
    command = 'load_file'
    args = [command,fits_file]
    # 
    res = method(*args, **kwargs)

