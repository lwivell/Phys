"""
This file creates a flat layered planetary body surface for induction studies
"""

import numpy as np
import json
import scipy.special
import matplotlib.pyplot as plt
import astropy.constants as const

miti = const.eps0.value
magperm = const.mu0.value

class LayeredSystem(object):

    def __init__(self, json=None, maxlayers=1):
        """
        Default constructor of model and layers.
        Should really have 
        """

        self._maxlayers = None
        self._structure = None
        self._name = None

        if json is None:
            self._maxlayers = maxlayers
        else:
            self.load_json(json)

    def load_json(self, filename):


        with open(filename, 'r') as fh:
            data = json.load(fh)

            if 'maxlayers' not in data:       #Checks required value exists
                raise RuntimeError('Must state the amount of layers in the model')
            
            if 'name' in data:
                self._name = data['name']

            self._maxlayers = data['maxlayers']
            self._structure = np.zeros((self._maxlayers, 4))

            for a in data['layers']:
                self._structure[a['index'],0] = a['index']
                self._structure[a['index'],1] = a['thickness']
                self._structure[a['index'],2] = a['conductivity']
                self._structure[a['index'],3] = a['permittivity']

    def wavenumber(angfreq, elecon, dieperm):
        """
        Calculates the wavenumber "k" 

        Parameters
        -------------------------------------------------------
        angfreq
            Angular frequency of the of the wave (omega)
        elecon
            Electrical conductivity (sigma)
        dieperm
            Dielectric permittivity (epsilon)

        Returns
        -------------------------------------------------------
        wavenumber k

        """
        
        k = (((angfreq**2)*magperm*dieperm)-((1j)*angfreq*magperm*elecon))**(1/2)
        #k = (-(1j)*angfreq*magperm*elecon)**(1/2)       # Low frequency approximation
        return k

    def intrin_imp(angfreq, elecon, dieperm):

        """
        Calculates the intrinsic impedance 

        Parameters
        -------------------------------------------------------
        angfreq
            Angular frequency of the of the wave (omega)
        elecon
            Electrical conductivity (sigma)
        dieperm
            Dielectric permittivity (epsilon)

        Returns
        -------------------------------------------------------
        Intrinsic impedance (eta)

        """

        eta = (((1j)*magperm*angfreq)/(elecon+((1j)*angfreq*dieperm)))**(1/2)
        return eta
    
    def appar_imp(wavenumber, intrinimp, prevapparimp, thick):
        """
        Calculates the apparent impedance of a layer in a conductivity layered model

        Parameters
        -------------------------------------------------------------------
        wavenumber
            Wavenumber as calculated for this current layer of model
        intrinimp
            Intrinsic impedance of the layer being modelled
        prevapparimp
            The apparent impedance of the layer below in the model. If the bottom layer, use intrinsic impedance of underlying halfspace
        thick
            Thickness of the layer of the model being calculated

        Returns
        -------------------------------------------------------------------
        apparimp
            Apparent impedance, Z, of the layer being modelled
        """    
        tanhval = ((1j)*wavenumber*thick)          
        if np.linalg.norm(tanhval) >= 100:                                   #This for loop simply clips the tanh value as it can get too large for python to handle, causes overflow error
            if np.real(tanhval) >= 0 and np.imag(tanhval)>=0:
                tanhval = 10 + 10j
            if np.real(tanhval) >= 0 and np.imag(tanhval)<=0:
                tanhval = 10 -10j
            if np.real(tanhval) <=0 and np.imag(tanhval)>=0:
                tanhval = -10+10j
            if np.real(tanhval) <=0 and np.imag(tanhval)<=0:
                tanhval = -10-10j  

        Z =np.complex128(intrinimp * (prevapparimp + (intrinimp*np.tanh(tanhval)))/(intrinimp + (prevapparimp*np.tanh(tanhval))))
        return Z    
    
    def appar_resis(apparimp, angfreq):
        """
        Calculates the apparent resistivity from outermost layer of a model

        Parameters
        --------------------------------------------------------------------
        apparimp
            Apparent impedance of the top layer
        angfreq
            Angular frequency

        Returns
        --------------------------------------------------------------------
        apparesis
            Apparent resistivity of the body (rho)
        """

        apparesis = (np.linalg.norm(apparimp)**2)/(angfreq*magperm)
        return apparesis
    
    def appar_cond(apparimp, angfreq):
        """
        Calculates the apparent conductivity from outermost layer of a model

        Parameters
        --------------------------------------------------------------------
        apparimp
            Apparent impedance of the top layer
        angfreq
            Angular frequency

        Returns
        --------------------------------------------------------------------
        apparcond
            Apparent conductivity of the body (rho)
        """

        apparcond = 1/((np.linalg.norm(apparimp)**2)/(angfreq*magperm))
        return apparcond