import numpy as np
import matplotlib.pyplot as plt

# Calculate a given function X(t)
def calculate(X):
    """
    X = [S0, I0, R0, V0, beta, gamma, nu, iota, l, timesteps]
    """
    # Initial Values
    Ss = [X[0]]
    Is = [X[1]]
    Rs = [X[2]]
    Vs = [X[3]]
    # Constants
    beta = X[4]
    gamma = X[5]
    nu = X[6]
    iota = X[7]
    l = X[8]

    timesteps = X[9]

    # Calculations
    i = 0
    for i in range(0, timesteps - 1):
        Ss.append(Ss[i] * ((-1 * beta) * Is[i] - nu + 1) + l * Vs[i])
        Is.append(Is[i] * (beta * Ss[i] - gamma + 1) + iota * Vs[i])
        Rs.append(Rs[i] + (gamma * Is[i]))
        Vs.append(Vs[i] * ((-1 * iota) - l + 1) + nu * Ss[i])

        # List of timesteps for plotting
        ts = [j for j in range(0, timesteps)]

    return Ss, Is, Rs, Vs, ts

# Plot all proportions of X(t)
def plotAll(X):
    """
    Plots all of the functions against time
    """
    Ss, Is, Rs, Vs, ts = calculate(X)

    fig, ax = plt.subplots()
    ax.plot(ts, Ss, color = 'green', linestyle = 'solid')
    ax.plot(ts, Is, color = 'orange', linestyle = 'solid')
    ax.plot(ts, Rs, color = 'red', linestyle = 'solid')
    ax.plot(ts, Vs, color = 'blue', linestyle = 'solid')
    plt.legend(['S\u209C', 'I\u209C', 'R\u209C', 'V\u209C'])
    plt.title("SIRV Plot")
    plt.xlabel("Timesteps")
    plt.ylabel("Proportion")
    fig.subplots_adjust(bottom =.2)
    fig.text(.125, .03, f"Final Proportions: \nS = {Ss[-1]:.5f}, I = {Is[-1]:.5f}, R = {Rs[-1]:.5f}, V = {Vs[-1]:.5f}")
    plt.show()

    return

# Plot any two proportions of X(t) together
def plotTwo(X, first = None, second = None):
    """
    Plots two different functions of the SIRV Model

    For plots:
        0 = S
        1 = I
        2 = R
        3 = V
    """
    Ss, Is, Rs, Vs, ts = calculate(X)
    # List of values to be indexed for plotting
    SIRV = [Ss, Is, Rs, Vs, ts]
    labels = ['S\u209C', 'I\u209C', 'R\u209C', 'V\u209C']
    initialVal = ['Inital S = ', 'Inital I = ', 'Inital R = ', 'Inital V = ']
    finalVal = ['Final S = ', 'Final I = ', 'Final R = ', 'Final V = ']

    # Set up plot
    fig, ax = plt.subplots()
    ax.plot(SIRV[first], SIRV[second], color = 'purple')
    ax.set_xlabel(labels[first], fontsize = 12)
    ax.set_ylabel(labels[second], fontsize = 12)
    fig.subplots_adjust(bottom = .2)
    fig.text(.125, .1, f"{initialVal[first]}{SIRV[first][0]:.05f}")
    fig.text(.7, .1, f"{finalVal[first]}{SIRV[first][-1]:.05f}")
    fig.text(.125, .05, f"{initialVal[second]}{SIRV[second][0]:.05f}")
    fig.text(.7, .05, f"{finalVal[second]}{SIRV[second][-1]:.05f}")
    plt.title(f"Comparing {labels[first]} & {labels[second]}", fontsize = 14)
    fig.text(.125, .9, f"Timesteps = {len(SIRV[first])}")
    plt.show()

# Plot different values for vaccination rates and immunity loss to visualise 'improved vaccinations'
def analyzeVaccine(X, newNu, newLoss):
    # Set up plot
    fig, axs = plt.subplots(2, 2, figsize = (9, 9))
    plotColors = ['red', 'orange', 'green', 'blue', 'black', 'cyan']
    
    # New values for vaccination rate
    newNuX = [X] * len(newNu)
    for i in range(0, len(newNu)):
        newNuX[i][6] = newNu[i]

        Ss, Is, Rs, Vs, ts, = calculate(newNuX[i])
        axs[0][0].plot(Vs, Is, color = plotColors[i]) #top left
        axs[0][1].plot(Rs, Is, color = plotColors[i]) #top right

    #new values for vaccine immunities
    newLossX = [X] * len(newLoss)
    for i in range(0, len(newLoss)):
        newLossX[i][8] = newLoss[i]

        Ss, Is, Rs, Vs, ts, = calculate(newLossX[i])
        axs[1][0].plot(Vs, Is, color = plotColors[i]) #bottom left
        axs[1][1].plot(Rs, Is, color = plotColors[i]) #bottom right

    #top left
    axs[0][0].set_ylabel('Infected', fontsize = 10)
    axs[0][0].set_title(f"Different Vaccination Rates\n(Immunity Loss = {newNuX[0][8]})", fontsize = 8)
    axs[0][0].legend(newNu)
    #top right
    axs[0][1].set_title(f"Different Vaccination Rates\n(Immunity Loss = {newNuX[0][8]})", fontsize = 8)
    axs[0][1].legend(newNu)
    #bottom left
    axs[1][0].set_xlabel('Vaccinated', fontsize = 10)
    axs[1][0].set_ylabel('Infected', fontsize = 10)
    axs[1][0].set_title(f"Different Vaccine-Imunity Rates of loss\n(Vaccination Rate = {newLossX[0][6]})", fontsize = 8)
    axs[1][0].legend(newLoss)
    #bottom right
    axs[1][1].set_xlabel('Removed', fontsize = 10)
    axs[1][1].set_title(f"Different Vaccine-Imunity Rates of loss\n(Vaccination Rate = {newLossX[0][6]})", fontsize = 8)
    axs[1][1].legend(newLoss)
    plt.subplots_adjust(hspace = .25)
    plt.suptitle("Improvements in Vaccination")
    plt.show()

    return

# Main function
def SIRV(S0, I0, R0, V0, beta, gamma, nu, iota, l, timesteps, plotChoice1 = None, plotChoice2 = None, newNu = None, newLoss = None, mainPlot = True):
    """
    For choosing two proportions to plot:
        0 = S
        1 = I
        2 = R
        3 = V
    """
    # Format values into list for helper functions
    F = [S0, I0, R0, V0, beta, gamma, nu, iota, l, timesteps]

    # Plots appear when conditions are met
    if mainPlot is True:
        plotAll(F)
    if plotChoice1 is not None and plotChoice2 is not None:
        plotTwo(F, plotChoice1, plotChoice2)
    if newNu is not None and newLoss is not None:
        analyzeVaccine(F, newNu, newLoss)

    return
