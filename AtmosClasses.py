#Owen Lehmer 11/1/16
#University of Washington, Seattle, WA
#Department of Earth and Space Sciences

class Layer:
    """
    This class is the basic structure for our model. The atmosphere is 
    comprised of these Layer class objects. They record all relevant info
    for the given layer (i.e. flux up, flux down, pressure, u, du/dr, density,
    T, etc.).

    IMPORTANT: things set to -1 by default may use that -1 as a check to see
               if the variable is set. 
    """
    def __init__(self, T=-1.0, p_top=-1.0, p_bot=-1.0, rho=-1.0, m=-1.0, \
            m_below=-1.0, u=0, du_dr=0, F_uv=0, F_up=0, F_down=0,\
            F_sol=0, h=-1.0, r=-1.0):
        self.T = float(T)          #The temp for the layer (assumed isothermal) [K]
        self.p_top = float(p_top)  #The pressure at the top of the layer [Pa]
        self.p_bot = float(p_bot)  #The pressure at the bottom of the layer [Pa]
        self.rho = float(rho)      #The average density of the layer [kg m-3]
        self.m = float(m)          #The mass of the layer [kg]
        self.m_below = float(m_below) #The total mass below this layer [kg]
        self.u = float(u)          #The radial outflow velocity in the layer [m s-1]
        self.du_dr = float(du_dr)  #The change in u with respect to r (du/dr) [s-1]
        self.F_uv = float(F_uv)    #The downward UV flux from the star [W m-2]
        self.F_sol = float(F_sol)  #The downward visible flux from the star [W m-2]
        self.F_up = float(F_up)    #The upward longwave flux at the layer [W m-2]
        self.F_down = float(F_down) #The downward longwave flux [W m-2]
        self.h = float(h)          #the height of the layer [m]
        self.r = float(r)          #The radial radius of the layer from the surface [m]
                                   #Note about self.r - this is the radius at the
                                   #bottom of the layer



