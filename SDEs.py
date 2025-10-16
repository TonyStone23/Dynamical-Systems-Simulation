import numpy as np
import random as r
import matplotlib.pyplot as plt

def calculateSLV(x, y, a, b, d, c, sigma_x, sigma_y, timesteps, delta_t, smoothen, seed = None):
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
        Xs.append(Xs[i]*(a*delta_t - b*Ys[i]*delta_t + 1) + (sigma_x * Xs[i] * root_delta_t * r.normalvariate(0, delta_t) * delta_t))
        Ys.append(Ys[i]*(d*Xs[i]*delta_t - c*delta_t + 1) + (sigma_y * Ys[i] * root_delta_t * r.normalvariate(0, delta_t) * delta_t))
        i += 1
        
    # For Plotting against time
    t_axis = [j*smoothen for j in range(len(Xs))]

    return Xs, Ys, t_axis, seed

def SLV(x, y, alpha, beta, delta, gamma, sigma_x, sigma_y, timesteps = 100, delta_t =.001, seed = None, plot = False, line = True, phasePortrait = False):

    # Run the simulation
    prey, predator, time, seed = calculateSLV(x, y, alpha, beta, delta, gamma, sigma_x, sigma_y, timesteps, delta_t, delta_t, seed)

    # Plot the populations against time
    if plot == True:
        fig, ax = plt.subplots(figsize = (10, 5))

        # Populate plot
        ax.plot(time, prey, color = "cyan", linestyle = 'solid')
        ax.plot(time, predator, color = "Navy", linestyle = 'solid')

        if line == True:
            preyMax = max(prey)
            predatorMax = max(predator)
            preyLine = [preyMax for i in time]
            predatorLine = [predatorMax for i in time]
            ax.plot(time, preyLine, color = 'black', linewidth = .5)
            ax.plot(time, predatorLine, color = 'black', linewidth = .5)
        # Decorate plot
        plt.title("Lotka-Volterra SDEs")
        plt.legend(['Prey (x\u209C)', 'Predators (y\u209C)'], loc = (0, 1))
        fig.text(.7, .95, f"{int(timesteps/delta_t)} outputs")
        fig.text(.7, .9, f"\u0394t = {delta_t}  Seed: {seed} ")
        plt.xlabel('Timesteps')
        plt.ylabel('Population')
        
    # Still a work in progress
    if phasePortrait == True:

        fig, ax = plt.subplots()

        # Populate plot
        ax.plot(prey, predator, color = 'purple')

        # Decorate PLot
        plt.xlabel("Prey Population")
        plt.ylabel("Predator Population")
        plt.title("Phase Portrait of\nLotka Volterra SDEs")
        fig.text(.7, .95, f"{int(timesteps/delta_t)} outputs")
        fig.text(.7, .9, f"\u0394t = {delta_t}  Seed: {seed} ")

    plt.show()

    return

initialPrey = 10
initialPredator = 2
alpha = 1.1
beta = .4
delta = .1
gamma = .4
sigma_x = 150
sigma_y = 150
timesteps = 200     

SLV(initialPrey, initialPredator,
    alpha,
    beta,
    delta,
    gamma,
    sigma_x,
    sigma_y,
    timesteps,
    delta_t = .0001,
    seed = 42, 
    plot = True,
    phasePortrait=True,
    line = False)