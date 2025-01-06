from typing import Tuple
import numpy as np

def spherical_unit_vectors(pos_theta, pos_phi):   

    """
    Calculates the spherical unit vectors for magnetic field component transformation

    Parameters
    -------------
    pos_theta
        1D array of the theta component of position
    pos_phi
        1D array of the phi component of position

    Returns
    unitr
        unit vector in i, for the radial component
    unittheta
        unit vector in j, for the theta component
    unitphi
        unit vector in k, for the phi component
    """
    unitr = np.array((np.sin(pos_theta)*np.cos(pos_phi), np.sin(pos_theta)*np.sin(pos_phi), np.cos(pos_theta)))
    unittheta = np.array((np.cos(pos_theta)*np.cos(pos_phi), np.cos(pos_theta)*np.sin(pos_phi), -np.sin(pos_theta)))
    unitphi = np.array((-np.sin(pos_phi), np.cos(pos_phi), np.zeros_like(pos_phi)))
    return unitr, unittheta, unitphi

def cartesian_to_spherical(posx, posy, posz):
    """
    Converts positions in spherical coordinates to cartesian coordinates
    
    Parameters
    ------------
    pos_x
        1D array of the x position
    pos_y
        1D array of the y position
    pos_z
        1D array of the z position  

    Returns
    -------------
    pos_r
        1D array of the radial component of position
    pos_theta
        1D array of the theta component of position
    pos_phi
        1D array of the phi component of position
    
    """
    posr = np.sqrt((posx**2)+(posy**2)+(posz**2))
    postheta = np.arccos(posz/(np.sqrt((posx**2)+(posy**2)+(posz**2))))
    posphi = np.sign(posy)*np.arccos(posx/(np.sqrt((posx**2)+(posy**2)))) 

def spherical_to_cartesian_pos(pos_r, pos_theta, pos_phi):
    """
    Converts positions in spherical coordinates to cartesian coordinates
    
    Parameters
    ------------
    pos_t
        1D array of the times at which positions were measured
    pos_r
        1D array of the radial component of position
    pos_theta
        1D array of the theta component of position
    pos_phi
        1D array of the phi component of position

    Returns
    -------------
    pos_x
        1D array of the x position
    pos_y
        1D array of the y position
    pos_z
        1D array of the z position
    
    """
    pos_x = pos_r * np.cos(pos_phi) * np.sin(pos_theta)  #shifts are negative as it moves measurement point not centre itself
    pos_y = pos_r * np.sin(pos_phi) * np.sin(pos_theta)
    pos_z = pos_r * np.cos(pos_theta)
    return pos_x, pos_y, pos_z

def cartdepth_from_sph(bodrad, normrad):
    """
    Converts normalised spherical radial position, into cartesian depth

    Parameters
    -------------------------------------------------------------------
    bodrad
        Radius of the body being modelled
    normrad
        Normalised radial position

    Returns
    -------------------------------------------------------------------
    cartdepth
        The cartesian depth
    """

    cartdepth = bodrad * ((1-normrad**3)/(2+normrad**3))
    return cartdepth

def cartcon_from_sphcon(sphcon, normrad):
    """
    Converts spherical conductivity to cartesian conductivity

    Parameters
    --------------------------------------------------------------------
    sphcon
        Spherical conductivity as a function of normalised radius
    normrad
        Normalised radial position

    Returns
    --------------------------------------------------------------------
    cartcon
        Conductivity as a function of depth z instead of r'
    
    """

    cartcon = ((2+normrad**3)/(3*normrad))**4 * sphcon
    return cartcon
