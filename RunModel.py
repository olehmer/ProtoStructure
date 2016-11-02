from AtmosClasses import Layer
import matplotlib.pyplot as plt
import sys
from math import e, pi




################Constants########################
M_Earth = 5.972E24 #mass of Earth in [kg]
GG = 6.67408E-11 #gravitational constant [m3 kg-1 s-2]



#################End Constants###################


def CalculateTopU(layer):
    #TODO calculate this...
    #TODO this will return cur_loss and u for the top layer
    return 10.0

def CalculateU(layers, ind, cur_loss):
    """
    This function will calculate the radial outflow velocity (u) for the 
    layer specified by ind. 

    From the equation of mass continuity we know that rho*u*r^2=cur_loss. If
    the current loss rate (cur_loss) is known then, since both radius and 
    density are known, we can calculate u for the current layer.

    Input:
    layers - the array of atmospheric layers
    ind - the index of the current layer (CANNOT BE 0)
    cur_loss - the current loss rate from the top of the atmosphere [kg s-1]

    Updates layers directly
    """

    if ind==0:
        #the comments clearly say not to do this...
        print("Invalid index passed to: CalculateU(layers, ind)")
        sys.exit(1)

    u = cur_loss/(layers[ind].rho*layers[ind].r**2.0)
    layers[ind].u = u




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
        layers[ind].u = CalculateU(layers, ind, cur_loss)

    return cur_loss


def UpdateLayerDistance(layers, ind, R_gas, mass, core_rad):
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
    mass - the mass interior to the current layer
    core_rad - the radius of the planetary core

    Updates layers directly

    Returns:
    mass - the mass of the layer and all interior atmosphere/core
    """

    #make a temporary index (t_ind) to reverse the direction
    t_ind = len(layers)-ind-1 #reverse the direction, we want to start at the ground

    g = 0.0 #define gravity, g [m s-2]

    if layers[t_ind].r == -1: #not set yet
        g = GG*mass/core_rad**2.0 #gravity!
    else:
        g = GG*mass/layers[t_ind].r

    #scale height is given by RT/g
    H = R*layers[t_ind].T/g #scale height, the altitude for which pressure decreases by e (~1/3)

    p_factor = layers[t_ind].p_bot/layers[t_ind].p_top 

    delta_r = p_factor*H #this is the height of the layer
    layers[t_ind].h = delta_r #record the height of the layer

    if ind = 0:
        #this is the bottom layer
        layers[t_ind].r = core_rad + delta_r/2.0 #the middle of the layer
    else:
        #this isn't the bottom layer
        layers[t_ind].r = layers[t_ind-1].r +layers[t_ind-1].h/2.0 + delta_r/2.0

    #calculate the layer mass from mass=(p_bot-p_top)*area/g
    area = 4.0*pi*(layers[t_ind].r-layers[t_ind].h/2.0)**2.0
    layer_mass = (layers[t_ind].p_bot-layers[t_ind].p_top)*area/g

    mass = mass + layer_mass

    return mass






def RunModel(core_mass=M_Earth, atmos_mass=0.1*M_Earth):

    print("running model")

RunModel()
