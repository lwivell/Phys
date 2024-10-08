import numpy as np

def wavenumber(angfreq, magperm, elecon, dieperm):
    """
    Calculates the wavenumber "k" 

    Parameters
    -------------------------------------------------------
    angfreq
        Angular frequency of the of the wave (omega)
    magperm
        Magnetic Permeability  (mu)
    elecon
        Electrical conductivity (sigma)
    dieperm
        Dielectric permittivity (epsilon)

    Returns
    -------------------------------------------------------
    wavenumber k

    """
    
    k = (((angfreq**2)*magperm*dieperm)-((1j)*angfreq*magperm*elecon))**(1/2)
    return k

def intrin_imp(angfreq, magperm, elecon, dieperm):
    """
    Calculates the wavenumber "k" 

    Parameters
    -------------------------------------------------------
    angfreq
        Angular frequency of the of the wave (omega)
    magperm
        Magnetic Permeability  (mu)
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

    Z = intrinimp * (prevapparimp + (intrinimp*np.tanh((1j)*wavenumber*thick)))/(intrinimp + (prevapparimp*np.tanh((1j)*wavenumber*thick)))
    return Z

def appar_resis(apparimp, angfreq, magperm):
    """
    Calculates the apparent resistivity from outermost layer of a model

    Parameters
    --------------------------------------------------------------------
    apparimp
        Apparent impedance of the top layer
    angfreq
        Angular frequency
    magperm
        Magnetic permeability

    Returns
    --------------------------------------------------------------------
    apparesis
        Apparent resistivity of the body (rho)
    """

    apparesis = (np.linalg.norm(apparimp)**2)/(angfreq*magperm)
    return apparesis
