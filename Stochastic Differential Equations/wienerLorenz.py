import numpy as np
import random as r
import matplotlib.pyplot as plt
import math

# Compute the lorenz dynamical system to desired endtime

def calculateLorenz(seed, sigma, rho, beta, init_x, init_y, init_z, magnitude, endtime, delta_t = .001, compare = False,):
    
    r.seed(seed)

    sqrt_t = np.sqrt(delta_t)

    # Initialize x, y, z
    x = [init_x]
    y = [init_y]
    z = [init_z]

    # Compute
    i = 0
    computes = endtime/delta_t
    while i < computes:
        x.append((sigma * (y[i] - x[i]) + x[i] * magnitude * sqrt_t * r.normalvariate(0, delta_t)) * delta_t + x[i])
        y.append((x[i] * (rho - z[i]) - y[i]*(1 + (magnitude * sqrt_t * r.normalvariate(0, delta_t)))) * delta_t + y[i])
        z.append(((x[i] * y[i]) - (beta * z[i]) + z[i]*(magnitude * sqrt_t * r.normalvariate(0, delta_t))) * delta_t + z[i])
        i += 1

    return x, y, z

# Compare differenec with specific changed conditions

def lorenzPlot(sigma, rho, beta, init_x, init_y, init_z, magnitude, endtime, seed = None,  delta_t = .001, portrait = False):

    # Set up plot
    ax = plt.figure(figsize =(8,5)).add_subplot(projection='3d')

    # Simulate and plot
    x, y, z = calculateLorenz(seed, sigma, rho, beta, init_x, init_y, init_z, magnitude, endtime, delta_t)
    ax.plot(x, y, z, lw = .5)

    # Decorate the pLot
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title("Lorenz ODE")

    # Alternatively, just view the system
    if portrait == True:    
        ax.set_axis_off()
        plt.subplots_adjust(0, 0, 1, 1, 0, 0)
    plt.show()

    return

lorenzPlot(10, 28, 8/3, 1, 0, 0, 100000, 50, 42, .001, 
          portrait = True)

