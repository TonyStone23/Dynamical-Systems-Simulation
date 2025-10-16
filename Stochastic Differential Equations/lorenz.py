import numpy as np
import random as r
import matplotlib.pyplot as plt
import math

# Compute the lorenz dynamical system to desired endtime

def calculateLorenz(sigma, rho, beta, init_x, init_y, init_z, endtime, delta_t = .001, compare = False):

    # Initialize x, y, z
    x = [init_x]
    y = [init_y]
    z = [init_z]

    # Compute
    i = 0
    computes = endtime/delta_t
    while i < computes:
        x.append(sigma * (y[i] - x[i]) * delta_t + x[i])
        y.append((x[i] * (rho - z[i]) - y[i]) * delta_t + y[i])
        z.append(((x[i] * y[i]) - (beta * z[i]))* delta_t + z[i])
        i += 1

    return x, y, z

# Compare differenec with specific changed conditions

def lorenzPlot(sigma, rho, beta, init_x, init_y, init_z, endtime, delta_t = .001, portrait = False):

    # Set up plot
    ax = plt.figure(figsize =(8,5)).add_subplot(projection='3d')

    # Simulate and plot
    x, y, z = calculateLorenz(sigma, rho, beta, init_x, init_y, init_z, endtime, delta_t)
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

def lorenzMultiplot(sigma, rho, beta, init_x, init_y, init_z, endtime, delta_t = .001, change = 0, changes = 0, portrait = False):

    colors = ['plum', 'lemonchiffon', 'palegreen', 'palevioletred', 'peachpuff','lightcoral']

    # Set up plot
    ax = plt.figure(figsize =(8,5)).add_subplot(projection='3d')

    # Run the simulation and plot for each iteration of x, y, and z initial values
    for i in range(0, changes + 1):
        x, y, z = calculateLorenz(sigma, rho, beta, init_x, init_y, init_z, endtime, delta_t)
        ax.plot(x, y, z, lw = .5, color = colors[i])

        init_x += change
        init_y += change
        init_z += change

    # Decorate the plot
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

# Compare the lorenz plot with nudged constants.
def lorenzPlotConstants(sigma, rho, beta, init_x, init_y, init_z, endtime, delta_t = .001, nudge_sigma = None, nudge_rho = None, nudge_beta = None, portrait = False):
   
    plotTwo = False
    # Set up plot
    ax = plt.figure(figsize =(8,5)).add_subplot(projection='3d')

    # Calculate
    x, y, z = calculateLorenz(sigma, rho, beta, init_x, init_y, init_z, endtime, delta_t)
    ax.plot(x, y, z, lw = .5)

    # Nudge a variable, if at least one is nudged a new simulation will be plotted
    if nudge_sigma is not None:
        sigma = nudge_sigma
        plotTwo = True
    if nudge_rho is not None:
        rho = nudge_rho
        plotTwo = True
    if nudge_beta is not None:
        beta = nudge_beta
        plotTwo = True
    
    if plotTwo == True:
        x, y, z = calculateLorenz(sigma, rho, beta, init_x, init_y, init_z, endtime, delta_t)
        ax.plot(x, y, z, lw = .5, color = 'red')

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

lorenzPlot(10, 28, 8/3, 1, 0, 0, 50, .001, 
          portrait = True)

lorenzMultiplot(10, 28, 8/3, 1, 0, 0, 25, .001, 10, 5,
          portrait = True)


lorenzPlotConstants(10, 28, 8/3, 1, 0, 0, 50, .001,
          nudge_sigma = 10, nudge_rho= 20, nudge_beta = 8/3,
          portrait = True)

