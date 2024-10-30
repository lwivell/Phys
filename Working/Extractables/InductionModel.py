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
        """
        Json file must have layers from innermost with n=0 and outermost with n= maxlayers-1
        """

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
                self._structure[a['index'],1] = a['thickness']*1000
                self._structure[a['index'],2] = a['conductivity']
                self._structure[a['index'],3] = a['permittivity']*miti

    def wavenumber(self, angfreq, elecon, dieperm):
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

    def intrin_imp(self, angfreq, elecon, dieperm):

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
    
    def appar_imp(self, wavenumber, intrinimp, prevapparimp, thick):
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
        for n in range(len(tanhval)):          
            if np.linalg.norm(tanhval[n]) >= 100:                                   #This for loop simply clips the tanh value as it can get too large for python to handle, causes overflow error
                if np.real(tanhval[n]) >= 0 and np.imag(tanhval[n])>=0:
                    tanhval[n] = 10 + 10j
                if np.real(tanhval[n]) >= 0 and np.imag(tanhval[n])<=0:
                    tanhval[n] = 10 -10j
                if np.real(tanhval[n]) <=0 and np.imag(tanhval[n])>=0:
                    tanhval[n] = -10+10j
                if np.real(tanhval[n]) <=0 and np.imag(tanhval[n])<=0:
                    tanhval[n] = -10-10j  
        
        Q = np.complex128(prevapparimp + (intrinimp*np.tanh(tanhval)))/(intrinimp + (prevapparimp*np.tanh(tanhval)))

        Z =np.complex128(intrinimp * Q)
        return Z, Q    
    
    def appar_resis(self, apparimp, angfreq):
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
        apparesis = np.zeros(apparimp.shape)
        for n in range(len(apparimp)):
            apparesis[n] = (np.linalg.norm(apparimp[n])**2)/(angfreq[n]*magperm)
        return apparesis
    
    def appar_cond(self, apparimp, angfreq):
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
            Apparent conductivity of the body (sigma)
        """
        apparcond = np.zeros(apparimp.shape)
        for n in range(len(apparimp)):
            apparcond[n] = 1/((np.linalg.norm(apparimp[n])**2)/(angfreq[n]*magperm))     
        return apparcond
    
    def iterate(self, angfreqs):
        """
        This function aims to complete the layered calculation for the loaded model
        """
        for n in range(self._maxlayers):
            k = self.wavenumber(angfreqs, self._structure[n][2], self._structure[n][3])
            intrin = self.intrin_imp(angfreqs, self._structure[n][2], self._structure[n][3])
            if n == 0:
                apparimp, Q = self.appar_imp(k, intrin, intrin, self._structure[n][1])
                prevapparimp = apparimp
            else:
                apparimp, Q = self.appar_imp(k, intrin, prevapparimp, self._structure[n][1])
                prevapparimp = apparimp

        conds = self.appar_cond(apparimp, angfreqs)

        return conds, Q