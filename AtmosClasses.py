class Layer:
    """
    This class is the basic structure for our model. The atmosphere is 
    comprised of these Layer class objects. They record all relevant info
    for the given layer (i.e. flux up, flux down, pressure, u, du/dr, density,
    T, etc.).
    """
    def __init__(self, T=-1.0, p_top=-1.0, p_bot=-1.0, rho=-1.0, u=-1.0, \
            du_dr=-1.0, F_uv=-1.0, F_up=-1.0, F_down=-1.0, r=-1.0):
        self.T = float(T)          #The temp for the layer (assumed isothermal) [K]
        self.p_top = float(p_top)  #The pressure at the top of the layer [Pa]
        self.p_bot = float(p_bot)  #The pressure at the bottom of the layer [Pa]
        self.rho = float(rho)      #The average density of the layer [kg m-3]
        self.u = float(u)          #The radial outflow velocity in the layer [m s-1]
        self.du_dr = float(du_dr)  #The change in u with respect to r (du/dr) [s-1]
        self.F_uv = float(F_uv)    #The downward UV flux from the star [W m-2]
        self.F_up = float(F_up)    #The upward longwave flux at the layer [W m-2]
        self.F_down = float(F_down) #The downward longwave flux [W m-2]
        self.r = float(r)          #The radial radius of the layer from the surface [m]


