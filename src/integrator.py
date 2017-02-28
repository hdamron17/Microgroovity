#! /usr/bin/python
""" 
Groovy egg diving numeric integrations for 2017 microgravity project
"""
import math
import numpy as np
import matplotlib.pyplot as plt


max_y_base = 1 - 1/math.sqrt(3) #max y (when multiplied by B in egg_function)
tau = 2 * math.pi #tau revolution

def egg_function(y, width, height):
    """ 
    Creates an egg with the desired parameters
    
    :param y: the desired Y value of the egg function
    :param width: coefficient which increases width of egg (formerly 'A')
    :param height: height of egg (egg is normally y in [0,1] but this changes 1)
    """
    max_y = max_y_base * height
    assert max_y <= y <= height, "y value %f out of range" % y #checks range
    #TODO if y > max_y: do ellipse
    A0 = 1 / (4 * max_y_base * (1-max_y_base) * (1-max_y_base/2))
    return width * math.sqrt(A0 / height * (1 - y/height) * (1 - y/(2*height)))

def groove_function(y, depth, width, height):
    """ 
    Creates a Groove for an egg with the desired parameters
    
    :param y: the Y value of the groove
    :param depth: depth coefficient of groove (#TODO the meaning of this is up to James)
    :return: Returns groove function x value at position y
    """
    #TODO TODO TODO major TODO
    assert max_y_base * height <= y <= height, "y value %f out of range" % y #checks range
    groove = egg_function(y, width, height) - depth
    return groove if groove > 0 else 0

def egg_derivative(y, width, height):
    """ 
    Calculates derivative of egg function at position y
    
    :param y: y position for which to find egg prime
    :param A: A width coefficient for egg function
    :return: Returns derivative of egg function at position y
    """
    assert max_y_base * height <= y <= height, "y value %f out of range" % y #checks range
    A0 = 1 / (4 * max_y_base * (1-max_y_base) * (1-max_y_base/2))
    return width**2 * A0 / (2*height * egg_function(y, width, height)) * (3/2 * y**2 / height**2 - 3 * y / height + 1)
    #The original derivation yields dx/dy but then it is inversed to dy/dx #TODO check validity of this statement
    #TODO if we're going to use a different shape for the groove, this is wrong

def perimeter(y, width, height, groove_angle, n, depth):
    """ 
    Calculates the perimeter of an egg on a cross-section
    
    :param y: the Y value of the desired cross-section
    :param groove_angle: surface angle to center of egg of each groove
    :param n: number of grooves
    :param depth: depth coefficient of groove function
    :return: Returns perimeter of egg intersection with plane at position y
    """
    max_y = max_y_base * height
    #print("%f <= %f <= %f ?" % (max_y, y, height))
    assert max_y <= y <= height, "y value %f out of range" % y #checks range (does not include perimeter below max point)
    w = egg_function(y, width, height) * 2 #egg width
    d = w / 2 - groove_function(y, depth, width, height) #depth = w/2 - groove function
    #separated into parts for typing simplicity
    #derivation by Gracen
    part1 = n * w / 2 * (tau / n - groove_angle)
    part2 = 2*n * math.sqrt(d**2 + 1/8 * w**2 * (1-math.cos(groove_angle)))
    return part1 + part2

def frontal_area(width, height):
    """ 
    Calculates the frontal area of an egg
    
    :param y: the Y value of a desired egg for calculation
    :param width: width coefficient of egg
    :param height: height coefficient of egg
    :return: Returns frontal area of egg (used for drag)
    """
    max_y = max_y_base * height
    radius = egg_function(max_y, height, width)
    return tau / 2 * radius**2 #area of circle

def accel(y, vel, Cd, density, height, width, groove_angle, n, sur_ten, contact_angle, mass, depth):
    """ 
    Calculates acceleration on a specific egg
    
    :param y: y position of egg
    :param vel: y velocity of egg
    :param Cd: drag coefficient
    :param density: density of outer fluid
    :param height: height of egg
    :param width: width of egg
    :param groove_angle: surface angle of groove (from center of egg) - cannot be larger than tau / n
    """
    drag_force = 0.5 * Cd * frontal_area(width, height) * density * vel**2
    if y <= height:
        slope = egg_derivative(y, width, height)
        #print("atan(%s) = %s; cos(atan(%s)) = %s" % (slope, math.atan(slope), slope, math.cos(math.atan(slope)))) #TODO remove
        capillary_force = perimeter(y, width, height, groove_angle, n, depth) * math.cos(contact_angle + math.atan(slope)) * sur_ten
    else:
        capillary_force = 0 #TODO check force after it has gone through if there is any up force
    force = -drag_force + capillary_force #drag pushes up and capillary force pulls down #TODO check signs
    #print("drag=%s, capillary=%s" % (drag_force, capillary_force))
    return force / mass

def v_next(y_next, v_i, a_i, Cd, density, height, width, groove_angle, n, sur_ten, contact_angle, mass, depth, dt):
    u_part = v_i + a_i / 2 * dt
    #TODO = 0 #slope of zero implying it never changes - should be derivative or somethign similar
    if y_next <= height:
        slope = egg_derivative(y_next, width, height)
        z_part = perimeter(y_next, width, height, groove_angle, n, depth) * math.cos(contact_angle + math.atan(slope)) * sur_ten
    else:
        z_part = 0 #TODO check that there are no capillary forces now
    w = egg_function(max_y_base * height, width, height) * 2 #egg width
    A = -dt * Cd * frontal_area(width, height) * density / (4 * mass) # A in quadratic formula Ax^2 + Bx + C
    B = -1
    C = u_part + dt / (2 * mass) * z_part
    v_f = (-B - math.sqrt(B**2 - 4 * A * C)) / (2 * A) #quadratic formula #TODO is it plus or minus?
    return v_f

def integrate(Cd, height, width, groove_angle, n, mass, depth, density, sur_ten, contact_angle, dt=0.01, time=2.2):
    """ 
    Integrates starting with y0 at maximum of egg and v0 = 0
    :param Cd: coefficient of drag (likely approximation of 0.5 for semisphere)
    :param density: water density rho {g/cm^3} #TODO check units on everything
    :param height: height of egg {cm}
    :param width: width of egg {cm}
    :param groove_angle: surface angle of each groove starting from the center (must not exceed tau/n) {rad}
    :param n: integer number of grooves
    :param sur_ten: surface tnsion of water (the Wiki says 71.99 mN/m) {mN/m}
    :TODO continue params
    :return: Returns numpy 2D array with rows: time, y position, y velocity
    """
    data = np.zeros((3, int(2.2 / dt))) #2D array with size 3 rows and row length of time times dt
    data[0,0] = 0 #start time 0
    data[1,0] = max_y_base * height #starting y at max of egg
    data[2,0] = 0 #zero initial velocity
    for i in range(1, data.shape[1]):
        t_i = data[0,i-1] #previous time
        y_i = data[1,i-1] #previous position
        v_i = data[2,i-1] #previous velocity
        a_i = accel(y_i, v_i, Cd, density, height, width, groove_angle, n, sur_ten, contact_angle, mass, depth)

        #leapfrog algorithm (modified to work better)
        data[0,i] = t_i + dt
        data[1,i] = y_f = y_i + v_i * dt + 0.5 * a_i * dt**2
        data[2,i] = v_next(y_f, v_i, a_i, Cd, density, height, width, groove_angle, n, sur_ten, contact_angle, mass, depth, dt)
        #print("yi=%s, vi=%s, ai=%s,    yf=%s, vf=%s" % (y_i, v_i, a_i, y_f, data[2,i]))
    return data

def check_domains(domains):
    """ 
    Checks that all values are within various constraints
    
    :param domains: dictionary of { value:(min,max) }
    :return: returns True if all values are valid, else False
    """
    for key, value in domains.items():
        if value[0] > value[1] or value[0] > key < value[1]:
            return False
        return True

def dive_depth(Cd, height, width, groove_angle, n, mass, depth, density, sur_ten, contact_angle, dt=0.01, time=2.2):
    """ 
    Wrapper over integration method which returns 0 if parameters are invalid else returns the depth of the dive, not including initial position
    
    Note that all parameters are described in integrate() function comments
    :return: returns depth of dive (subtracts starting y) or 0 if parameters are out of bounds
    """
    domains = {
        Cd : (0, 1), #TODO is this a valid bound
        height : (0, 0.1), #TODO find real limit
        width : (0, 0.1), #TODO find real limit
        groove_angle : (0, tau/n),
        n : (0, 20), #TODO decide on actual max n
        mass : (0, 1), #TODO find real limit
        depth : (0, depth / 2 if groove_angle > 0 else 0),
        density : (0, math.inf),
        sur_ten : (0, math.inf),
        contact_angle : (0, tau / 2),
        dt : (0, math.inf),
        time : (0, math.inf)
    }
    if not check_domains(domains):
        return 0
    dive_data = integrate(Cd, height, width, groove_angle, n, mass, depth, density, sur_ten, contact_angle, dt, time)
    y0 = dive_data[1,0] #initial y position
    yf = dive_data[1].max()
    depth = yf - y0
    return depth

def depth_wrapper(params, density, sur_ten, contact_angle, dt=0.01, time=2.2):
    """ 
    Wrapper over the other wrapper which accepts all variable parameters in optimization as a single collection
    
    :param params: tuple containing 7 values which correspond to the first 7 parameters of dive_depth() function
    Note: other parameters are described below integrate() function
    :return: returns dive depth from the dive_depth function
    """
    Cd = params[0]
    height = params[1]
    width = params[2]
    groove_angle = params[3]
    n = params[4]
    mass = params[5]
    depth = params[6]
    return dive_depth(Cd, height, width, groove_angle, n, mass, depth, density, sur_ten, contact_angle, dt, time)

example_params = (0.5, 0.05, 0.05, tau/8, 8, 0.1, 0.001, 1000, 0.0728, tau/20, 0.001, 2.2)

if __name__ == "__main__":
    print(depth_wrapper(example_params[:7], *example_params[7:]))
    print(dive_depth(*example_params))
    data = integrate(*example_params)
    #print("depth: %f" % max(data[1]))
    plt.plot(data[0], data[1], 'r-', label="postion")
    plt.plot(data[0], data[2], 'b-', label="velocity")
    plt.legend()
    plt.show()
