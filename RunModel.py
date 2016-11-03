from AtmosClasses import Layer
import matplotlib.pyplot as plt
import numpy as np
import sys
from math import e, pi




################Constants########################
M_Earth = 5.972E24 #mass of Earth in [kg]
R_Earth = 6.371E6 #radius of Earth [m]
GG = 6.67408E-11 #gravitational constant [m3 kg-1 s-2]
R_H2 = 4124.0 #gas constant for H2 [J kg-1 K]



#################End Constants###################


def CalculateTopU(layer):
    #TODO calculate this...
    #TODO this will return cur_loss and u for the top layer
    return 10.0


def UpdateLayerVelocity(layers, ind, cur_loss):
    """
    This function takes the array of layers and the current index to update. 
    It will update the radial outflow velocity, u, at each layer.

    IMPORTANT: this must be called after UpdateR() as the appropriate distance
    for each layer is needed to calculate the new u value.

    Input:
    layers - the array of layers in the atmosphere
    ind - the index of the layer to operate on in layers array
    cur_loss - the current loss rate from the top of the atmosphere. When 
               calling this function on the first layer (ind=0) this function
               will set cur_loss. All subsequent calls should just pass that
               value back to this function.

    Updates layers directly

    Returns:
    cur_loss - the current loss rate [kg s-1] from the top of the atmosphere.
               This is only calculated for the top layer then passed on all 
               subsequent layers for ease of computation.
    """

    if ind==0:
        #this is the top layer, calculate the loss rate from the top layer
        #TODO this will return a tuple with cur_loss
        cur_loss, u = CalculateTopU(layers[ind])
        layers[ind].u = u
    else:
        #not the top layer, calculate u based on cur_loss
        u = cur_loss/(layers[ind].rho*layers[ind].r**2.0)
        layers[ind].u = u

    return cur_loss

def UpdateLayerDuDr(layers, ind):
    """
    This function will update the du/dr value for each layer. The du/dr of each
    layer will be set by looking at the slope of the layer above and below (in
    the case of the top layer the slope will be approximated from just the top
    layer and the one below).

    Input:
    layers - the layers of the amtosphere
    ind - the index of the layer to process
    
    Updates layers directly
    """

    du_dr = 0.0
    if ind == 0:
        #this is the top layer, don't look at a layer above (there isn't one...)
        du_dr = (layers[ind].u-layers[ind+1].u)/(layers[ind].r-layers[ind+1].r)
    elif ind == len(layers)-1:
        #if this is the bottom layer we're in the same situation, just use 2 layers
        du_dr = (layers[indi-1].u-layers[ind].u)/(layers[ind-1].r-layers[ind].r)
    else:
        #this isn't the top or bottom, use the u values form surrounding layers
        du_dr = (layers[indi-1].u-layers[indi+1].u)/(layers[ind-1].r-layers[ind+1].r)

    layers[ind].du_dr = du_dr


def UpdateLayerDistance(layers, ind, R_gas, core_rad):
    """
    This function will calculate the radial distance of the current layer. This
    is done by assuming the layer is approximately isothermal and gravity is
    constant throughout the layer. The first time this function is called g is
    assumed constant throughout the atmosphere. Once r has been calculated 
    gravity can be calculated for each layer (based on previous r).

    IMPORTANT: this function starts from the ground layer first, then works
    up!

    Input:
    layers - the array of layers in the atmosphere
    ind - the index of the current layer in the layers array
    R_gas - the gas constant for the atmosphere
    core_rad - the radius of the planetary core

    Updates layers directly
    """

    #make a temporary index (t_ind) to reverse the direction
    t_ind = len(layers)-ind-1 #reverse the direction, we want to start at the ground

    g = 0.0 #define gravity, g [m s-2]

    if layers[t_ind].r == -1: #not set yet
        g = GG*layers[t_ind].m_below/core_rad**2.0 #gravity!
    else:
        g = GG*layers[t_ind].m_below/layers[t_ind].r**2.0

    #in addition to g at each layer we also have du/dr as an upward acceleration
    #scale height is given by RT/(g-u*du/dr)
    #scale height, the altitude for which pressure decreases by e (~1/3)
    H = R*layers[t_ind].T/(g-layers[t_ind].u*layers[t_ind].du_dr)
    

    p_factor = layers[t_ind].p_bot/layers[t_ind].p_top 

    delta_r = p_factor*H #this is the height of the layer
    layers[t_ind].h = delta_r #record the height of the layer

    if ind = 0:
        #this is the bottom layer
        layers[t_ind].r = core_rad + delta_r/2.0 #the middle of the layer
    else:
        #this isn't the bottom layer
        layers[t_ind].r = layers[t_ind-1].r +layers[t_ind-1].h/2.0 + delta_r/2.0




def UpdateLayerPressure(layers, ind):
    """
    This function will update the pressure bounds of each layer. This is done
    using P=m*(g-u*du/dr) at each layer then summing from top to bottom.

    Inputs:
    layers - the array of layers in the atmosphere
    ind - the index of the current layer to process

    Updates layers directly
    """

    g = GG*layers[ind].m_below/layers[ind].r**2.0
    area = 4.0*pi*layers[ind].r**2.0
    p_layer = layers[ind].m*(g-layers[ind].u*layers[ind].du_dr)

    if ind != 0:
        #this isn't the top layer so add the pressure from the above layers
        p_layer = p_layer + layers[ind-1].p_bot
        layers[ind].p_top = layers[ind-1].p_bot

    layers[ind].p_bot = p_layer

def UpdateTempProfile(layers, T_profile):
    """
    Simply loop over the layers and copy the temperature value to T_profile.
    This will also compute the difference between the old T and the new T.

    Input:
    layers - the atmospheric layers
    T_profile - the temperature at each atmospheric layer

    Updates T_profile directly

    Return:
    diff - the difference between the current T_profile and the old one
    """

    diff = 0.0
    for i in range(0,len(layers)):
        diff += abs(T_profile[i] - layers[i].T)
        T_profile[i] = layers[i].T

    return diff

        
def BalanceAtmosphere(core_mass, core_rad, atmos_mass, R_gas, F_uv, F_sol, F_long,\
        kappa, T_prof_in=[], N=100, iter_lim=200, tol=1.0E-4):
    """
    This is the top level function to balance the atmospheric model. Given the 
    above parameters this will compute the temperature profile, pressure 
    profile, height profile, and the loss rate from the atmosphere. This model uses a 3-stream
    approach to handle UV, visible, and IR flux.

    Inputs:
    core_mass - the mass of the planetary core [kg]
    core_rad - the radius of the planetary core [m]
    atmos_mass - the mass of the atmosphere [kg]
    R_gas - the specific gas constant for the atmosphere [J kg-1 K]
    F_uv - the downward UV flux at the top of the atmosphere [W m-2]
    F_sol - the non-UV non-longwave downward flux at the top of atmos [W m-2]
    F_long - the longwave downward flux at TOA, typically 0 [w m-2]
    kappa - the mass absorption coefficient to use in optical depth [m2 kg-1]
    T_prof_in - optional temperature profile can be passed in
    N - the number of atmospheric layers to use in the simulation, defaults to 100
    iter_lim - the number of iterations after which the model will stop
    tol - the desired tolerance for model accuracy (looks at the temp profile change)

    NOTE: TOA = top of atmosphere, and TOA is at and index of 0 in our layers array

    Output:
    p_profile - the pressure profile of the atmosphere
    r_profile - the radial distance of each layer
    T_profile - the temperature profile of the atmosphere
    mass_flux - the mass flux loss rate from the atmosphere [kg s-1]
    """
    
    if len(T_prof_in>0) and len(T_prof_in) != N:
        #this is a problem, need a temperature for each layer
        print("Temperature profile input had %d layers \
                but the model expects %d"%(len(T_prof_in),N))
        sys.exit(1)

    #initialize our array to hold the layers
    layers = LayersInit(core_mass, core_rad, atmos_mass, T_prof_in)
    #TODO implement LayersInit()!!!!!!

    T_profile = np.zeros(N) #initialize temperature profile
    UpdateTempProfile(layers, T_profile)

    count = 0 #keep track of what iteration we're on
    diff = tol + 1.0 #track our current level of accuracy

    while count<iter_lim and diff>tol:
        """
        This is where all the time will be spent for this algorithm. There are
        lots of seemingly redundant loops below, but they are important. The
        order of the loops also matters. The general steps are:

        1. Calculate the radial distance, r, for each layer
        2. Once r is known, we can calculate velocity (u)
        3. Once u is known for each layer we can find du/dr
        4. With r, u, and du/dr we can calculate the new pressure for the layer
        5. Update the radiative transfer (aka update temperature)

        An aside:
        If you're sufficiently clever you can probably combine some of the loops
        below, but the runtime will still be O(N) so probably not worth the
        effort.
        """

        #1 - first update the distance to each layer
        for i in range(0,N):
            UpdateLayerDistance(layers, i, R_gas, core_rad)

        #2 - next update the radial outflow velocity for each layer
        cur_loss = 0.0 #initialize to 0
        for i in rage(0,N):
            cur_loss = UpdateLayerVelocity(layers, i, cur_loss)

        #3 - update the du/dr for each layer
        for i in range(0,N):
            UpdateLayerDuDr(layers, i)

        #4 - update the pressure for each layer
        for i in range(0,N):
            UpdateLayerPressure(layers, i)

        #5 - calculate the radiative transfer in the layers
        for i in range(0,N):
            UpdateLayerRadTrans() #TODO implement this!

        #Update the temperature profile and calculate the diff
        diff = UpdateTempProfile(layers, T_profile)

        count += 1

    if count == iter_lim:
        #we ended the while loop due to hitting the iter_lim
        print("failed to converge before iteration limit!")

    #get the pressure and distance profiles
    p_profile = np.zeros(N)
    r_profile = np.zeros(N)
    for i in range(0,N):
        p_profile[i] = (layers[i].p_bot+layers[i].p_top)/2.0 #just average it
        r_profile[i] = layers[i].r

    #calculate the mass flux loss rate
    u_top = layers[0].u
    r_top = layers[0].r
    rho_top = layers[0].rho
    mass_flux = 4.0*pi*rho_top*u_top*r_top**2.0 #loss rate [kg s-1]

    return (p_profile, r_profile, T_profile, mass_flux)




        



def RunModel(core_mass=M_Earth, atmos_mass=0.1*M_Earth, R_gas=):

    print("running model")


    #initialize the atmosphere























RunModel()
