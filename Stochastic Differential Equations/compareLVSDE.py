import numpy as np
import random as r
import matplotlib.pyplot as plt

def calculateSLV(x, y, a, b, d, c, sigma_x, sigma_y, timesteps, delta_t, seed = None):
    """
        x: inital prey population
        y: initial predator population
        a: prey growth rate
        b: prey death rate
        d: predator growth rate
        c: predator death rate
        sigma_x: intensity of noise for prey 
        sigma_y: intensity of noise for predator 
        timesteps: time passed
        delta_t: spcace between each point
    """
    # Set seed for reproducabiltiy
    if seed is not None:
        r.seed(seed)
    else:
        seed = r.randint(0, 10000000)
        r.seed(seed)

    # Predator and Prey Populations
    Xs = [x]
    Ys = [y]
    t_axis = []

    points = timesteps/delta_t
    i = 0
    root_delta_t = np.sqrt(delta_t)

    while i < points:
        Xs.append(Xs[i]*(a*delta_t - b*Ys[i]*delta_t + 1) + (sigma_x * Xs[i] * root_delta_t * r.normalvariate(0, 1) * delta_t))
        Ys.append(Ys[i]*(d*Xs[i]*delta_t - c*delta_t + 1) + (sigma_y * Ys[i] * root_delta_t * r.normalvariate(0, 1) * delta_t))
        i += 1
        
    # For Plotting against time
    t_axis = [j*delta_t for j in range(len(Xs))]

    return Xs, Ys, t_axis, seed

def SLV(x, y, alpha, beta, delta, gamma, sigma_xValues, sigma_yValues, timesteps = 100, delta_t =.001, seed = None):

    colors = ['blue', 'green', 'purple', 'red', 'Orange', 'Cyan', 'grey']

    allPrey = []
    allPredators = []
    allSeeds = []

    # Run simulation at differnet values of noise variables
    for s in range(0, len(sigma_xValues)):
        prey, predator, time, seed = calculateSLV(x, y, alpha, beta, delta, gamma, sigma_xValues[s], sigma_yValues[s], timesteps, delta_t, seed)
        allPrey.append(prey)
        allPredators.append(predator)
        allSeeds.append(seed)

    # Plot the Simulations against time

    fig, ax = plt.subplots(figsize = (12, 6))

    for i in range(len(allSeeds)):

        # Populate plot
        ax.plot(time, allPrey[i], color = colors[i], linestyle = 'solid', linewidth = .5)
        ax.plot(time, allPredators[i], color = colors[i], linestyle = 'solid', linewidth = .5)

    # Decorate plot
    fig.subplots_adjust(left = .175)
    for s in range(len(allSeeds)):
        fig.text(.01, .8 - (s * .1), f"{colors[s]}: sigma = {sigma_xValues[s]}", fontsize = 8)
    plt.title("Comparing Lotka-Volterra SDE Noise Constants")
    fig.text(.8, .95, f"{int(timesteps/delta_t)} outputs")
    fig.text(.8, .9, f"\u0394t = {delta_t}")
    plt.xlabel('Timesteps')
    plt.ylabel('Population')

    plt.show()

    return

initialPrey = 10
initialPredator = 2
alpha = 1.1
beta = .4
delta = .1
gamma = .4
sigma_xValues = [0, 50, 150, 200, 250, 300, 350]
sigma_yValues = [0, 50, 150, 200, 250, 300, 350]
timesteps = 50     

SLV(initialPrey, initialPredator,
    alpha,
    beta,
    delta,
    gamma,
    sigma_xValues,
    sigma_yValues,
    timesteps,
    delta_t = .0001,
    seed = 42)