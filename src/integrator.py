#! /usr/bin/python
""" 
Generic numeric integrations for 2017 microgravity project
"""


def egg_function(y, A):
    """
    Creates an egg with the desired parameters
    
    :param y: the desired Y value of the egg function
    :param A: coefficient which increases width of egg
    """
    #TODO figure out what other parameters are required for egg function
    pass #TODO

def groove_function(y):
    """
    Creates a Groove for an egg with the desired parameters
    
    :param y: the Y value of the groove
    :param TODO: other parameters governing grove depth
    :return: Returns groove function x value at position y
    """
    #TODO figure out what other parameters
    pass #TODO

def perimeter(y, groove_angle, n):
    """
    Calculates the perimeter of an egg on a cross-section
    
    :param y: the Y value of the desired cross-section
    :param groove_angle: surface angle to center of egg of each groove
    :param n: number of grooves
    :return: Returns perimeter of egg intersection with plane at position y
    """
    pass #TODO

def area(y):
    """
    Calculates the frontal area of an egg
    
    :param y: the Y value of a desired egg for calculation
    """
    pass #TODO

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
