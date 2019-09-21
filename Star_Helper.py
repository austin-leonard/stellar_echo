"""Several functions to help determine stellar parameters mathematically from limited available data. Also contains
the Binary class, which determines the distance between two binary stars given period and mass
"""
import numpy as np
from scipy import constants as c 
import uncertainties


def get_distance(parallax):
    """ If the distance isn't given in Simbad. Input in milliarcseconds. """
    conv = parallax/1000
    distance = 1/conv
    return distance


def pc_to_au(pc):
    """ Given a distance in parsecs, return that distance in AU. """ 
    au = pc*206264.806
    return au


class Star:
    """A  star with a given b-v color and distance. 
    
    Performs the calculation of the radius, luminosity and mass of the star.

    Has some intrinsic error in the calculations -- seems to overestimate the temperature and underestimate the radius,
    
    though from what I've read this often happens when trying to get these values for M dwarf stars."""

    def __init__(self, name, b_mag, v_mag, dist):
        """
        Parameters:

        name: string, name of the star;

        b_mag: float, b-band magnitude;

        v_mag: float, v-band magnitude;

        dist: float, distance to star in parsecs
        """

        self.name = name
        
        self.b_mag = b_mag
        
        self.v_mag = v_mag
        
        self.dist = dist

    def get_temperature(self):
        """ Obtain the temperature given the B-V magnitude."""
        # There is a discrepancy here between the values obtained from this formula
        # and the temperatures that I find for stars on the internet.
        # I still don't know where it comes from. It seems to overestimate --
        # perhaps it's not a good fit for red dwarf stars of very low temperatures?
        # This formula comes from Ballesteros 2012, eq. (14): https://arxiv.org/pdf/1201.1809.pdf
        b_v_mag = self.b_mag - self.v_mag
        temp = 4600*((1/((0.92*b_v_mag)+1.7))+(1/((0.92*b_v_mag)+0.62)))
        print("The temperature of ", self.name, "is ", temp, "K")
        return temp

    def get_absolute_magnitude(self):
        """ Obtain absolute magnitude in the v-band using apparent magnitude in v-band, given in Simbad,
            and distance, which is either directly given or calculated via parallax."""
        abs_mag = self.v_mag + 5 - 5*np.log10(self.dist)
        print("The absolute magnitude in the v-band is ", abs_mag)
        return abs_mag

    def get_luminosity(self, abs_mag):
        """Given the absolute magnitude of the star in the v-band,
        returns the luminosity of the star in the v-band.""" 
        abs_mag_sun = 4.83  # In v-band
        l_star = 10**(0.4*(abs_mag_sun-abs_mag))
        print("The luminosity of ", self.name, "is ", l_star, " lsun")
        return l_star

    def get_radius(self, Temp, l_star):
        """Temperature in K, Luminosity in Lsun.
        Return Radius in Rsun,,,, or is it? . """
        t_sun = 5778 # Kelvin
        radius = ((t_sun/Temp)**2) * np.sqrt(l_star)
        print("The radius of ", self.name, "is ", radius, "rsun")
        return radius
    
    def mlr(self, abs_mag):
        """ Mass-Luminosity Relation. Obtain the mass of a star (a good approximation, anyway)
        given the absolute magnitude in the V-band. Good for values of M_v < 19.
        Comes from Benedict 2016, eq. (11): https://iopscience.iop.org/article/10.3847/0004-6256/152/5/141/pdf """
        c_0 = 0.19226
        c_1 = -0.050737
        c_2 = 0.010137
        c_3 = -0.00075399
        c_4 = -0.000019858
        x_0 = 13.0 
        sub = abs_mag - x_0
        mass = c_0 + c_1 * sub + c_2 * (sub**2) + c_3 * (sub**3) + c_4 * (sub**4)
        print("The mass of ", self.name, "is", mass, "Msun")
        return mass

# Get the parameters of a star
st = Star("gj 4198", 14.8, 13.19, 12.1665)
temperature = st.get_temperature()
mag = st.get_absolute_magnitude()
lum = st.get_luminosity(mag)
rad = st.get_radius(temperature, lum)
m = st.mlr(mag)


class Binary:

    """Given a binary star system, find the separation in AU.
    
    Parameters: 
    
    name: string, system name
    
    period: float, orbital period (days)
    
    mass: float, combined mass of the system (Msun)
    """
    
    def __init__(self, name, period, mass):
        self.name = name
        self.period = period
        self.mass = mass
    
    def semi_major_axis(self):
        """ Period in days, mass in solar masses.
         Used to obtain semi-major axis of binary star orbit. """
        p_yrs = self.period/365
        g = c.g
        constant_term = (g**2)/(4*(c.pi**2))
        a = ((p_yrs**2) * self.mass * constant_term)**(1/3)
        print("The separation of ", self.name, "is", a, "AU")
        return a

# Get the separation of a binary star
pair = Binary("CM Dra", 1.27, 0.4368)
sep = pair.semi_major_axis()
