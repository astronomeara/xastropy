"""
#;+ 
#; NAME:
#; koavis
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for KOA data visualization, starting w/ HIRES
#;   13-Aug-2015 by JMO.  
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from astropy import units as u
from astropy import constants as const

########################## ##########################
########################## ##########################
class KOA_Vis(object):
    """ Class for visualization of data from the KOA

    inspired by the suite of tools built by Nicolas Lehner and JMO for KODIAQ 
    
    Attributes
    ----------
    wrest: float
      Rest wavelength of line analyzed
    vmnx: tuple (vmin,vmax)
      Velocity range for analysis

    JXP on 11 Dec 2014, JMO Aug 12 2015
    """
