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

def SLV(x, y, alpha, beta, delta, gamma, numSigma, sigmaIncrement, timesteps = 100, delta_t =.001, seed = None):

    xRange = 15
    yRange = 8

    fig, axs = plt.subplots(numSigma, numSigma, sharex = True, sharey = True, figsize = (12, 6))

    sigma_xValues = [s * sigmaIncrement for s in range (0, numSigma)]
    sigma_yValues = sigma_xValues

    # Run simulation at differnet values of noise variables
    for i in range(len(sigma_xValues)):
        for ii in range(len(sigma_yValues)):
            prey, predator, time, seed = calculateSLV(x, y, alpha, beta, delta, gamma, sigma_xValues[i], sigma_yValues[ii], timesteps, delta_t, seed)
            axs[i][ii].plot(prey, predator, linestyle = 'solid', linewidth = .3)
            if ii == (0):
                axs[i][ii].set_ylabel(f"Prey Sigma = {sigma_xValues[i]}", fontsize = 6)
            if i == (len(sigma_yValues)-1):
                axs[i][ii].set_xlabel(f"Predator Sigma = {sigma_yValues[ii]}", fontsize = 6)
    
    # Decorate plot
    fig.text(.8, .95, f"{int(timesteps/delta_t)} outputs | \u0394t = {delta_t}")
    plt.suptitle("Comparing Noise Constants of the LV SDEs")
    fig.supxlabel("Prey Population")
    fig.supylabel("Predator Population")
    plt.tight_layout()
    plt.show()

    return

initialPrey = 10
initialPredator = 2
alpha = 1.1
beta = .4
delta = .1
gamma = .4
numSigma = 4
sigmaIncrement = 100
timesteps = 50     

SLV(initialPrey, initialPredator,
    alpha,
    beta,
    delta,
    gamma,
    numSigma,
    sigmaIncrement,
    timesteps,
    delta_t = .0001,
    seed = 42)