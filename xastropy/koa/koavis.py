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
from linetools.spectra import io as lsi

########################## ##########################
########################## ##########################
class KOA_Vis(object):
    """ Class for visualization of data from the KOA

    inspired by the suite of tools built by Nicolas Lehner and JMO for KODIAQ 
    
    Attributes
    ----------
    qsonam: string
      Name of the quasar
    inst: string
      Instrument
    zem: float
      emission redshift

    JMO Aug 13 2015
    """

    #Initialize
    def __init__(self,qsonam,inst,zem):

        #target info
        self.qsonam = qsonam
        self.inst = inst
        self.zem = zem

        # Data
        self.koa_data = {}
        self.keys = ['ra2000','dec2000','zabs','vlow','vhigh']
        self.key_dtype = ['a11','a12','f4','f4','f4']

        # Init
        for key in self.keys:
            self.koa_data[key] = 0

            
    # Access the data and return the value
    def __getitem__(self,item):
        try:
            return self.koa_data[item]
        except KeyError:
            raise KeyError
        
    ########################## ##########################
    def injest_koa_spec(self, specfil, zem, **kwargs):
        """ Injest a spectrum, put it in format for koavis

        Parameters
        ----------
        specfil:  Name of the spectrum (usually a fits file)
        zem:  Quasar emission redshift

        Returns
        -------
        fill of object keys
        Spectrum1D class

        JMO Aug 13 2015
        """

        #read in fits file
        self.spec = lsi.readspec(specfil,efil=efil)

        #fill
        self.koa_data['zem']=zem

    ########################## ##########################
    def koa_pick_redlines(self, spec, zabs)
        """ Given a spectrum and an absorber redshift, pick
        those lines that are redward of lya emission
        (e.g. for easy kinematics work)

        Parameters
        ----------
        spec:  Spectrum1D class
           Input spectrum
        zem:  Quasar emission redshift

        Returns
        -------
        lines:  float array
           list of lines to use
        fills a koa_data zabs
        
        JMO Aug 13 2015
        """


        #fill
        self.koa_data['zabs']=zabs
