"""
#;+ 
#; NAME:
#; cos_halos
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for COS-Halos analysis
#;        Now includes a COS-Dwarfs sub-class
#;   29-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp, pickle, sys, glob
from astropy.io import fits, ascii
from astropy import units as u 
from astropy.table import QTable, Table, Column
#from astropy import constants as const

from linetools.spectra import io as lsio

from xastropy.galaxy.core import Galaxy
#from xastropy.cgm.core import CGM_Abs, CGM_Abs_Survey
from xastropy.cgm.core import CGMAbsSurvey, CGMSys
from xastropy.xutils import xdebug as xdb
from xastropy.igm.abs_sys.ionclms import IonClms
from xastropy.kinematics.absline import KinAbs

from astropy.utils.misc import isiterable

# Path for xastropy
#xa_path = imp.find_module('xastropy')[1]

#def ion_name(ion):
#def photo_cross(Z, ion, E, datfil=None, silent=False):

# Class for COS_Halos Survey
class COSHalos(CGMAbsSurvey):
    """Inherits CGM Abs Survey

    Attributes:
    -----------
    fits_path: str, optional
      Path to the FITS data files for COS-Halos
    """
    # Initialize with a .dat file
    def __init__(self, tree=None, fits_path=None, kin_init_file=None, cdir=None):

        # Generate with type
        CGMAbsSurvey.__init__(self)
        self.survey = 'COS-Halos'
        self.ref = 'Tumlinson+11; Werk+12; Tumlinson+13; Werk+13'
        if cdir is None:
            self.cdir = os.environ.get('DROPBOX_DIR')+'/COS-Halos/'
        # Summary Tables
        if fits_path is None:
            self.fits_path = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Halos/lowions/FITS')
        else:
            self.fits_path = fits_path
        # Kinematics
        if kin_init_file is None:
            self.kin_init_file = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Halos/Kin/'+
                                                  'coshalo_kin_driver.dat')
        else:
            self.kin_init_file = kin_init_file

    # Load from mega structure
    def load_single(self,inp, skip_ions=False):
        """ Load a single COS-Halos sightline
        Appends to cgm_abs list

        Paramaeters
        ----------
        inp: tuple or str
          if tuple -- (field,gal_id)
            field: str 
              Name of field (e.g. 'J0226+0015')
            gal_id: str 
              Name of galaxy (e.g. '268_22')
        """
        # Parse input
        if isinstance(inp,basestring):
            fil = inp
        elif isinstance(inp,tuple):
            field,gal_id  = inp
            tmp = self.fits_path+'/'+field+'.'+gal_id+'.fits'
            fils = glob.glob(tmp)
            if len(fils) != 1:
                raise IOError('Bad field, gal_id: {:s}'.format(tmp))
            fil = fils[0]
        else:
            raise IOError('Bad input to load_single')

        # Read COS-Halos file
        print('cos_halos: Reading {:s}'.format(fil))
        hdu = fits.open(fil)
        summ = hdu[1].data
        galx = hdu[2].data
        self.cgm_abs.append( CGMSys(galx['ra'][0],
            galx['dec'][0],
            summ['zfinal'][0],
            galx['qsora'][0],
            galx['qsodec'][0], 
            galx['zqso'][0]))
        mm = len(self.cgm_abs)-1
        # COS-Halos naming
        self.cgm_abs[mm].field = galx['field'][0]
        self.cgm_abs[mm].gal_id = galx['galid'][0]
        # Galxy properties
        self.cgm_abs[mm].galaxy.halo_mass = summ['LOGMHALO'][0] 
        self.cgm_abs[mm].galaxy.stellar_mass = summ['LOGMFINAL'][0] 
        self.cgm_abs[mm].galaxy.sfr = (galx['SFR_UPLIM'][0], galx['SFR'][0],
            galx['SFR_FLAG'][0]) # FLAG actually gives method used
        # Ions
        if skip_ions is True:
            return
        self.cgm_abs[mm].abs_sys.ions = IonClms()
        all_Z = []
        all_ion = []
        for jj in range(summ['nion'][0]):
            iont = hdu[3+jj].data
            if jj == 0: # Generate new Table
                dat_tab = Table(iont)
            else:
                try:
                    dat_tab.add_row(Table(iont)[0])
                except:
                    xdb.set_trace()
            all_Z.append(iont['zion'][0][0])
            all_ion.append(iont['zion'][0][1])
        # Add Z,ion
        dat_tab.add_column(Column(all_Z,name='Z'))
        dat_tab.add_column(Column(all_ion,name='ion'))
        # Set
        self.cgm_abs[mm].abs_sys.ions._data = dat_tab
        # NHI
        self.cgm_abs[mm].abs_sys.NHI = self.cgm_abs[mm].abs_sys.ions[(1,1)]['CLM']


    # Load from mega structure
    def load_mega(self,flg=1, data_file=None, cosh_dct=None, test=False, **kwargs):
        """ Load the data for COS-Halos

        Paramaeters
        ----------
        flg: integer (1)
          Flag indicating how to load the data
          0 = IDL mega structure
          1 = FITS files from Dropbox
        data_file: string
          Name of data file
        pckl_fil: string
          Name of file for pickling

        JXP on 30 Nov 2014
        """
        #from xastropy.cgm import core as xcc
        #reload(xcc)

        # IDL save file
        if flg == 0:
            raise ValueError('This will not work.')
            if data_file is None:
                data_file = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Halos/lowions/'+
                                            'coshalos_lowmetals_mega.sav')
        elif flg == 1: # FITS files
            # Loop
            if test is True:
                cos_files = glob.glob(self.fits_path+'/J091*.fits') # For testing
            else:
                cos_files = glob.glob(self.fits_path+'/J*.fits')
            # Setup
            self.nsys = len(cos_files)
            # Read
            for fil in cos_files:
                self.load_single(fil, **kwargs)
        else:
            raise ValueError('cos_halos.load_mega: Not ready for this flag {:d}'.format(flg))

    
    ########################## ##########################
    def load_abskin(self,flg=1,kin_file=None,kin_init_file=None):
        """ Load the absorption-line kinematic data for COS-Halos (or COS-Dwarfs)
        Calculate from scratch if needed

        Paramaeters
        ----------
        flg: int, optional 
          Flag indicating how to load the data
            0 = Load from file
            1 = Generate
        kin_init_file: str
          Name of kinematics driver file
        kin_file: str
          Name of kinematics output file [First made for John Forbes]
        """
    
        if flg == 0: # Load
            if kin_file is None:
                kin_file = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Halos/Kin/'+
                                                  'COS-Halos_kin.fits')
            hdu = fits.open(kin_file)
            # Metals
            metals = Table(hdu[1].data)
            for row in metals:
                mt = np.where( (row['field']==self.field) & 
                    (row['gal_id']==self.gal_id))[0]
                xdb.set_trace()

        elif flg == 1: # Generate
            # Read init file
            if kin_init_file is None:
                kin_init_file = self.kin_init_file
            kin_init = ascii.read(kin_init_file,guess=False)
    
            # Loop to my loop
            fgal = zip(self.field, self.gal_id)
            for qq,cgm_abs in enumerate(self.cgm_abs):
                # Match to kin_init
                mt = np.where( (cgm_abs.field == kin_init['QSO']) &
                               (cgm_abs.gal_id == kin_init['Galaxy']) )[0]
                if len(mt) == 0:
                    print('load_kin: No kinematics for {:s}, {:s}'.format(cgm_abs.field,
                                                                          cgm_abs.gal_id))
                    continue
                mt = mt[0]

                # Metals
                if kin_init['flgL'][mt] > 0:
                    wrest = kin_init['mtl_wr'][mt]*u.AA 
                    if wrest.value <= 1:
                        xdb.set_trace()
                    spec = self.load_bg_cos_spec( qq, wrest )
                    vmnx = (kin_init['L_vmn'][mt]*u.km/u.s, kin_init['L_vmx'][mt]*u.km/u.s)
                    # Process
                    cgm_abs.abs_sys.kin['Metal'] = KinAbs(wrest, vmnx)
                    cgm_abs.abs_sys.kin['Metal'].fill_kin(spec, per=0.07)
                    # Save spec
                    cgm_abs.abs_sys.kin['Metal'].spec = spec
                else:
                    # Fill with zeros (for the keys)
                    cgm_abs.abs_sys.kin['Metal'] = KinAbs(0.*u.AA, (0., 0.))

                # HI
                if kin_init['flgH'][mt] > 0:
                    wrest = kin_init['HI_wrest'][mt]*u.AA 
                    if wrest.value <= 1:
                        xdb.set_trace()
                    spec = self.load_bg_cos_spec( qq, wrest )
                    vmnx = (kin_init['HIvmn'][mt]*u.km/u.s, kin_init['HIvmx'][mt]*u.km/u.s) 
                    # Process
                    cgm_abs.abs_sys.kin['HI'] = KinAbs(wrest, vmnx)
                    cgm_abs.abs_sys.kin['HI'].fill_kin(spec, per=0.07)
                    cgm_abs.abs_sys.kin['HI'].spec = spec
                else:
                    # Fill with zeros (for the keys)
                    cgm_abs.abs_sys.kin['HI'] = KinAbs(0.*u.AA, (0., 0.))


            #tmp = cos_halos.abs_kin('Metal')['Dv']
            #xdb.set_trace()
    # 
    def load_gal_spec(self, idx):
        """ Load the galaxy spectrum

        Parameters
        ----------
        idx: int
          Index of the cgm_abs list

        Returns:
        ----------
        spec: XSpectrum1D 
          Splices the blue and red side for LRIS

        JXP on 12 Oct 2015
        """
        # Init
        cgm_abs = self.cgm_abs[idx]
        # Directories
        galdir = os.environ.get('DROPBOX_DIR')+'/coshaloanalysis/'
        fielddir = 'fields/'+cgm_abs.field+'/'
        sysdir = cgm_abs.gal_id+'/spec1d/'
        sysname = cgm_abs.field+'_'+cgm_abs.gal_id

        # Find files
        lris_files = glob.glob(galdir+fielddir+sysdir+sysname+'*corr.fits')
        if len(lris_files) == 0:
            raise ValueError('No LRIS files!')
        elif len(lris_files) == 2:
            lris_files.sort()
            specb = lsio.readspec(lris_files[0]) 
            specr = lsio.readspec(lris_files[1]) 
            spec = specb.splice(specr)
        else:
            raise ValueError('Not sure what happened')

        # Return
        return spec


    def load_bg_cos_spec(self, inp, wrest):
        """ Load the absorption-line kinematic data for COS-Halos
        Calculate from scratch if needed

        Parameters
        ----------
        idx: int or tuple
          int -- Index of the cgm_abs list
          tuple -- (field,gal_id)
        wrest: Quantity
          Rest wavelength for spectrum of interest
    
        JXP on 11 Dec 2014
        """
        if isinstance(inp,int):
            cgm_abs = self.cgm_abs[inp]
        elif isinstance(inp,tuple):
            cgm_abs = self[inp]
        # Directories
        fielddir = 'Targets/'+cgm_abs.field+'/'
        sysdir = cgm_abs.gal_id+'_z{:5.3f}'.format(cgm_abs.galaxy.z)
        sysname = cgm_abs.field+'_'+sysdir

        # Transition
        templ_fil = os.environ.get('DROPBOX_DIR')+'/COS-Halos/Targets/system_template.lst'
        tab = ascii.read(templ_fil)
        mt = np.argmin(np.abs(tab['col1']-wrest.value))
        if np.abs(tab['col1'][mt]-wrest.value) > 1e-2:
            raise ValueError('get_coshalo_spec: wrest={:g} not found!'.format(wrest))
        trans = tab['col2'][mt]+tab['col3'][mt]

        # Read
        slicedir = self.cdir+fielddir+sysdir+'/fitting/'
        slicename = sysname+'_'+trans+'_slice.fits'
        spec = lsio.readspec(slicedir+slicename, flux_tags=['FNORM'], sig_tags=['ENORM'])
        # Fill velocity
        spec.velo = spec.relative_vel((cgm_abs.galaxy.z+1)*wrest)
    
        #spec.qck_plot()
        return spec


    def __getitem__(self, inp):
        '''Grab CgmAbs Class from the list

        Parameters:
        -----------
        ion: tuple
          tuple:  (field,gal_id)

        Returns:
        ----------
        cgm_abs
        '''
        # Generate lists
        fields = np.array([cgm_abs.field for cgm_abs in self.cgm_abs])
        galids = np.array([cgm_abs.gal_id for cgm_abs in self.cgm_abs])
        #
        mt = np.where( (fields==inp[0]) & (galids == inp[1]))[0]
        if len(mt) != 1:
            warn.warning('CosHalos: CGM not found')
            return None
        else:
            return self.cgm_abs[mt]

class COSDwarfs(COSHalos):
    """Inherits COS Halos Class

    Attributes:
    -----------
    fits_path: str, optional
      Path to the FITS data files for COS-Halos
    """
    # Initialize with a .dat file
    def __init__(self, tree=None, fits_path=None, kin_init_file=None, cdir=None):

        # Generate with type
        CGMAbsSurvey.__init__(self)
        self.survey = 'COS-Dwarfs'
        self.ref = 'Bordoloi+14'
        if cdir is None:
            self.cdir = os.environ.get('DROPBOX_DIR')+'/COS-Dwarfs/'
        if fits_path is None:
            self.fits_path = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Dwarfs/Targets/FITS')
        else:
            self.fits_path = fits_path
        # Kinematics
        if kin_init_file is None:
            self.kin_init_file = os.path.abspath(os.environ.get('DROPBOX_DIR')+'/COS-Dwarfs/Kin/cosdwarfs_kin_driver.dat') 
        else:
            self.kin_init_file = kin_init_file

'''            
########################## ##########################
# Testing
if __name__ == '__main__':

    flg_fig = 0 
    #flg_fig += 1  # Load FITS
    #flg_fig += 2  # NHI plot
    flg_fig += 2**2  # Simple Kinematics

    # Load FITS
    if (flg_fig % 2) == 1:
        cos_halos = COSHalos()
        cos_halos.load_mega()
        print(cos_halos)
    
    # Simple rho vs NHI plot
    if (flg_fig % 2**2) >= 2**1:
        cos_halos = COSHalos()
        cos_halos.load_mega()
        x= cos_halos.rho
        y= cos_halos.NHI
        xdb.xplot(x, y, scatter=True)
    #
    # Simple kinematics
    if (flg_fig % 2**3) >= 2**2:
        cos_halos = COSHalos()
        cos_halos.load_mega()#test=True)
        cos_halos.load_abskin()
        # Plot
        mtl_kin = cos_halos.abs_kin('Metal')
        gd = np.where(mtl_kin['flg'] > 0)
        xdb.xplot(cos_halos.NHI[gd], mtl_kin['Dv'][gd], scatter=True)

        HI_kin = cos_halos.abs_kin('HI')
        gd = np.where(HI_kin['flg'] > 0)
        xdb.xplot(cos_halos.NHI[gd], HI_kin['Dv'][gd], scatter=True)
    print('All done')
'''            
