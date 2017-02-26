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
    #TODO figure out what other parameters are required for egg function
    max_y = max_y_base * height
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
    groove = egg_function(y, width, height) - depth
    return groove if groove > 0 else 0

def egg_derivative(y, A):
    """ 
    Calculates derivative of egg function at position y
    
    :param y: y position for which to find egg prime
    :param A: A width coefficient for egg function
    :return: Returns derivative of egg function at position y
    """
    return A**2 * A0 / (2*B * egg_function(y, A)) * (3/2 * y**2 / B**2 - 3 * y / B + 1)
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
    print("%f <= %f <= %f ?" % (max_y, y, 1))
    max_y = 0 #TODO remove this
    assert max_y <= y <= 1, "y out of range" #checks range (does not include perimeter below max point)
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

def accel(y, vel, Cd, density, height, width, groove_angle, n, sur_ten, capillary_angle, mass, depth):
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
    drag_force = 0.5 * Cd * frontal_area(height, width) * density * vel**2
    TODO = 0 #should be slope
    capillary_force = perimeter(y, width, height, groove_angle, n, depth) * math.cos(capillary_angle + math.atan(TODO)) * sur_ten
    force = -drag_force + capillary_force #drag pushes up and capillary force pulls down #TODO check signs
    return force / mass

def v_next(y_next, v_i, a_i, Cd, density, height, width, groove_angle, n, sur_ten, capillary_angle, mass, depth, dt):
    u_part = v_i + a_i / 2 * dt
    TODO = 0 #slope of zero implying it never changes - should be derivative or somethign similar
    z_part = perimeter(y_next, width, height, groove_angle, n, depth) * math.cos(capillary_angle + math.atan(TODO)) * sur_ten
    w = egg_function(y_next, width, height) * 2 #egg width
    A = -dt * Cd * tau * w**2 * density / (8 * mass) # A in quadratic formula Ax^2 + Bx + C
    B = -1
    C = u_part + dt / (2 * mass) * z_part
    v_f = (-B + math.sqrt(B**2 - 4 * A * C)) / (2 * A) #quadratic formula #TODO is it plus or minus?
    return v_f

def integrate(Cd, density, height, width, groove_angle, n, sur_ten, capillary_angle, mass, depth, dt=0.01, time=2.2):
    """ 
    Integrates starting with y0 at maximum of egg and v0 = 0
    :TODO params
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
        a_i = accel(y_i, v_i, Cd, density, height, width, groove_angle, n, sur_ten, capillary_angle, mass, depth)

        #leapfrog algorithm (modified to work better)
        data[0,i] = t_i + dt
        data[1,i] = y_f = y_i + v_i * dt + 0.5 * a_i * dt**2
        data[2,i] = v_next(y_f, v_i, a_i, Cd, density, height, width, groove_angle, n, sur_ten, capillary_angle, mass, depth, dt)
    return data

if __name__ == "__main__":
    print(tau)
    data = integrate(Cd=0.5, density=0.3, height=1, width=1, groove_angle=tau/8, n=8, sur_ten=0.5, capillary_angle=tau/20, mass=0.1, depth=0.05, dt=0.00001)
    print("depth: %f" % max(data[1]))
    plt.plot(data[0], data[1], 'r-')
    plt.plot(data[0], data[2], 'b-')
    plt.show()
