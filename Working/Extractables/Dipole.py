from typing import Tuple
import numpy as np

def calc_OTD_components(pos_t, pos_r, pos_theta, pos_phi, dipmag, diptheta, dipphi, x_shift=0.0, y_shift=0.0, z_shift=0.0) -> Tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    """
    Calculates the magnetic field according to any input OTD conditions

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
    dipmag
        magnitude of the dipole moment
    diptheta
        theta tilt of the dipole
    dipphi
        phi tilt of the dipole
    x_shift
        dipole centre shift in the x direction
    y_shift
        dipole centre shift in the y direction
    z_shift
        dipole centre shift in the z direction

    Returns
    --------------
    pos_t
        1D array of the times at which positions were measured
    OTDbr
        1D array of the radial magnetic field components for the offset tilted dipole model
    OTDbtheta
        1D array of the theta magnetic field components for the offset tilted dipole model
    OTDbphi
        1D array of the phi magnetic field components for the offset tilted dipole model
    """

    pos_x, pos_y, pos_z = spherical_to_cartesian_pos(pos_r, pos_theta, pos_phi)  #Converts coordinates and applies dipole centre shifts
    pos_x = pos_x - x_shift       
    pos_y = pos_y - y_shift
    pos_z = pos_z - z_shift


    position = np.array((pos_x, pos_y, pos_z))    #Defining vectors needed for calculation
    OTDfield = (0,0,0)
    DIPmom = np.array((dipmag*np.sin(np.deg2rad(diptheta))*np.cos(np.deg2rad(dipphi)), dipmag*np.sin(np.deg2rad(diptheta))*np.sin(np.deg2rad(dipphi)), dipmag*np.cos(np.deg2rad(diptheta))))

    OTDfield = (3*np.dot(DIPmom, position)*position)/(np.linalg.norm(position, axis=0)**5) - (DIPmom[:,None]/(np.linalg.norm(position, axis=0)**3)) #Calculates OTD field

    unitr, unittheta, unitphi = spherical_unit_vectors(pos_theta, pos_phi) #Calculates unit vectors

    OTDbr = np.array(np.sum(OTDfield* unitr, axis=0))                  #Converts back into spherical coordinates
    OTDbtheta = np.array(np.sum(OTDfield* unittheta, axis=0) )
    OTDbphi = np.array(np.sum(OTDfield* unitphi, axis=0))

    return pos_t, OTDbr, OTDbtheta, OTDbphi

def GeoMag(pos_t, pos_r, pos_theta, pos_phi, DIPtheta, DIPphi, x_shift=0.0, y_shift=0.0, z_shift=0.0):
    """
    Transforms a Geographic coordinate into an OTD centred system

    Parameters
    -----------------------------------------
    pos_t
        1D array of the times at which positions were measured
    pos_r
        1D array of the radial component of position
    pos_theta
        1D array of the theta component of position
    pos_phi
        1D array of the phi component of position
    DIPtheta
        theta tilt of the dipole
    DIPphi
        phi tilt of the dipole
    x_shift
        dipole centre shift in the x direction
    y_shift
        dipole centre shift in the y direction
    z_shift
        dipole centre shift in the z direction

    Returns
    pos_t
        1D array of the times at which positions were measured
    magposr
        1D array of radial positions
    magpostheta
        1D array of colatitude component
    magposphi
        1D array of azimuthal component 
    """

    posx, posy, posz = spherical_to_cartesian_pos(pos_r, pos_theta, pos_phi)
    
    lambdaD = np.deg2rad(DIPphi) 
    phiD = (np.pi/2) - np.deg2rad(DIPtheta)

    omega = lambdaD + (np.pi/2)
    theta = (np.pi/2) - phiD
    phi = -(np.pi/2)
    
    
    a = (np.cos(phi)*np.cos(omega)) - (np.sin(phi)*np.sin(omega)*np.cos(theta))
    b = (np.cos(phi)*np.sin(omega)) + (np.sin(phi)*np.cos(omega)*np.cos(theta))
    c = np.sin(phi)*np.sin(theta)
    d = (-np.sin(phi)*np.cos(omega)) - (np.cos(phi)*np.sin(omega)*np.cos(theta))
    e = (-np.sin(phi)*np.sin(omega)) + (np.cos(phi)*np.cos(omega)*np.cos(theta))
    f = np.cos(phi)*np.sin(theta)
    g = np.sin(omega)*np.sin(theta)
    h = -np.cos(omega)*np.sin(theta)
    i = np.cos(theta)

    RotMat = np.array([[a,b,c],[d,e,f],[g,h,i]])
    geopos = np.vstack((posx, posy, posz))

    RotPos = np.matmul(RotMat, geopos)

    magposx = RotPos[0] - x_shift
    magposy = RotPos[1] - y_shift
    magposz = RotPos[2] - z_shift

    magposr = np.sqrt((magposx**2)+(magposy**2)+(magposz**2))
    magpostheta = np.arccos(magposz/(np.sqrt((magposx**2)+(magposy**2)+(magposz**2))))
    magposphi = np.sign(posy)*np.arccos(posx/(np.sqrt((posx**2)+(posy**2))))
    
    return pos_t, magposr, magpostheta, magposphi, RotPos

def dipLvalues(pos_r, pos_theta):
    """
    Uses radial and colatitude components of coordinate system centred to the dipole and calculates L-values

    Parameters
    -------------------------
    pos_r
        1D array of radial positions
    pos_theta
        1D array of colatitude component

    Returns
    Lval
        1D array of L-values at given radial and colatitude positions
    """

    lat = (np.pi/2) - pos_theta
    Lval = pos_r / (np.cos(lat)**2)

    return Lval