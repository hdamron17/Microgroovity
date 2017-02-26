#! /usr/bin/python
""" 
Generic numeric integrations for 2017 microgravity project
"""
import math


max_y_base = 1 - 1/math.sqrt(3) #max y (when multiplied by B in egg_function)

def egg_function(y, width, height):
    """ 
    Creates an egg with the desired parameters
    
    :param y: the desired Y value of the egg function
    :param width: coefficient which increases width of egg (formerly 'A')
    :param height: height of egg (egg is normally y in [0,1] but this changes 1)
    """
    #TODO figure out what other parameters are required for egg function
    max_y = max_y_base * height
    if y > max_y:
        pass #TODO James's ellipse formula
    else:
        A0 = 1 / (4 * max_y_base * (1-max_y_base) * (1-max_y_base/2))
        return width * math.sqrt(A0 / height * (1 - y/height) * (1 - y/(2*height)))

def groove_function(y, depth):
    """ 
    Creates a Groove for an egg with the desired parameters
    
    :param y: the Y value of the groove
    :param depth: depth coefficient of groove (#TODO the meaning of this is up to James)
    :return: Returns groove function x value at position y
    """
    pass #TODO James

def egg_derivative(y, A):
    """ 
    Calculates derivative of egg function at position y
    
    :param y: y position for which to find egg prime
    :param A: A width coefficient for egg function
    :return: Returns derivative of egg function at position y
    """
    return A / (2 * egg_function(y, A)) * (3/2 * y**2 - 3 * y + 1)

def perimeter(y, groove_angle, n, A):
    """ 
    Calculates the perimeter of an egg on a cross-section
    
    :param y: the Y value of the desired cross-section
    :param groove_angle: surface angle to center of egg of each groove
    :param n: number of grooves
    :return: Returns perimeter of egg intersection with plane at position y
    """
    assert max_y <= y <= 1, "y out of range" #checks range (does not include perimeter below max point)
    w = egg_function(y, A) * 2 #egg width
    d = w / 2 - groove_function(y) #depth = w/2 - groove function
    #separated into parts for typing simplicity
    #derivation by Gracen
    part1 = n * w / 2 * (math.tau / n - groove_angle)
    part2 = 2*n * math.sqrt(d**2 + 1/8 * w**2 * (1-math.cos(groove_angle)))
    return part1 + part2

def area(y):
    """ 
    Calculates the frontal area of an egg
    
    :param y: the Y value of a desired egg for calculation
    """
    

def accel(y, vel, Cd):
    """ 
    Calculates acceleration on a specific egg
    
    :param y: y position of egg
    :param vel: y velocity of egg
    :param Cd: drag coefficient
    :param density: density of outer fluid
    :TODO:
    """
    water_density = 1 #TODO density of water in correct units
    drag_force = 0.5 * Cd * area * water_density * vel**2

def depth():
    y0 = max_y

if __name__ == "__main__":
    print("TODO")
