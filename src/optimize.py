#! /usr/bin/env python
import numpy as np
from scipy.optimize import minimize, brute, basinhopping
import integrator #numerical integration of egg diving
import sys
from skopt import gp_minimize, forest_minimize


def printer(p):
    print("%s -> %s" % (p, integrator.depth_wrapper(p, *integrator.example_params[7:])))

if __name__ == "__main__":
    #initial_guess = integrator.example_params[:6]
    args = integrator.example_params[6:]

    #print(initial_guess)
    #print(args)

    #names  = ['Cd',             'height',        'width',           'groove_angle',         'n',               'mass',           'depth']
    domains = [(0, 1),           (0, 0.1),        (0, 0.1),          (0, integrator.tau/20), (0, 20),           (0, 1),           (0, 0.05)]
    #values = [4.99980651e-01,   5.02834681e-02,  -2.53530591e+03,  7.85397870e-01,          7.99999999e+00,    1.00068048e-01,   1.01541872e-03] #from one optimization
    #values = [-6.99519307e-03,  7.09793499e-02,  7.09279929e-02,   -8.90682375e+07,         1.20004628e+01,    8.76918719e+01,   5.00218628e+00] #from Powell solver
    #values = [-3.98149522e-02,  9.04317205e-02,  3.85213854e-02,   6.84185878e-01,          1.06602889e+01,    1.38436173e-01,   8.96350513e-04] #from Nelder-Mead

    ### Then found out the return was inside the for loop and messed up checking domains

    #print(integrator.depth_wrapper(values, *args))

    initial_guesses = np.divide(np.genfromtxt("initials.csv", delimiter=",", comments="#"), [0.6,0.6,integrator.tau/2,20,1000,0.3])
    print(initial_guesses)
    output = ""

    type = "Nelder-Mead"
    for initial_guess in initial_guesses:
        print(initial_guess)
        done = minimize(integrator.depth_wrapper, initial_guess, args=args, method=type,  options={'maxiter':100000, 'disp':False})
        #done = basinhopping(integrator.depth_wrapper, initial_guess, minimizer_kwargs={'args':args, 'method':type}, niter=20000, disp=False)
        #done = forest_minimize(integrator.depth_outer_wrapper,
        #    [(0,0.06), (0,0.06), (0,3.14), (0,20), (500,1000), (0,0.03)],
        #x0=initial_guess, verbose=False, n_calls=200, n_random_starts=20, kappa=10)
        #print("%s\n-----\n%s\n=====\n" % (type, done))
        dive_depth = integrator.depth_wrapper(done['x'], *args)
        print("depth(%s) = %s" % ( np.multiply(done['x'],[0.6,0.6,integrator.tau/2,20,1000,0.3]), dive_depth))
        output += str(dive_depth) + "\n"
        if len(sys.argv) > 1 and sys.argv[1] == "plot":
            try: integrator.plot_dive(*(tuple(done['x']) + args))
            except: pass
    with open("outs.csv", "w+") as fp:
        fp.write(output)
