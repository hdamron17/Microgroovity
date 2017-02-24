#! /usr/bin/python
""" 
Generic numeric integrations for 2017 microgravity project
"""


def egg_function(y):
    #TODO figure out what other parameters are required for egg function
    pass #TODO

def groove_function(y):
    #TODO figure out what other parameters
    pass #TODO

def perimeter(y):
    pass #TODO

def area(y):
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
    
    drag_force = 0.5 * Cd * area * density * vel**2
