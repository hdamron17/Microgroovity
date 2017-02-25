#! /usr/bin/python
""" 
Generic numeric integrations for 2017 microgravity project
"""


def egg_function(y):
    """
    Creates an egg with the desired parameters
    
    :param y: the desired Y value of the egg function
    """
    #TODO figure out what other parameters are required for egg function
    pass #TODO

def groove_function(y):
    """
    Creates a Groove for an egg with the desired parameters
    
    :param y: the Y value of the groove
    """
    #TODO figure out what other parameters
    pass #TODO

def perimeter(y):
    """
    Calculates the perimeter of an egg on a cross-section
    
    :param y: the Y value of the desired cross-section
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
