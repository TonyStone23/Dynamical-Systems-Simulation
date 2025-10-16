import matplotlib.pyplot as plt
import numpy as np

def vectorFieldLV(a, b, d, c, timesteps, delta_t):

    # Set up chart
    fig, ax = plt.subplots()
    # Determine the amount of initial conditions, and their values
    x = [i for i in range(10, 41, 1)]
    y = [i for i in range(1, 21, 1)]
    X, Y = np.meshgrid(x, y)
    U = []
    V = []

    points = timesteps/delta_t
    
    # Run the simulation at each (x, y)
    i = 0
    for iii in range(len(y)):
        for ii in range(len(x)):
            # Initialize with a particular (x, y)
            Xs = [x[ii]]
            Ys = [y[iii]]
            while i < points:
                Xs.append(Xs[i]*(a*delta_t - b*Ys[i]*delta_t + 1))
                Ys.append(Ys[i]*(d*Xs[i]*delta_t - c*delta_t + 1))
                i += 1

            # Capture final values to create vectors
            U.append(Xs[-1])
            V.append(Ys[-1])

    # Populate plot
    magnitude = np.sqrt(np.array(U) ** 2 + np.array(V) ** 2)
    plt.quiver(X, Y, U, V, magnitude)

    # Decorate plot
    plt.title("Vector Field of\nPrey and Predator Populations")
    plt.xlabel("Prey")
    plt.ylabel("Predator")
    fig.text(.1, .9, f"{timesteps} timesteps")
    plt.show()
    return

def vectorFieldSIRV(S0, I0, R0, V0, beta, gamma, nu, iota, l, timesteps, plotNum):
    """
    plotNum = {1, 2, 3}
        1 = Susceptible & Infected
        2 = Susceptible & Vaccinated
        3 = Infected & Vaccinated
        4 = Susceptible & removed
    """
    
    # Proportions
    Ss = [S0]
    Is = [I0]
    Rs = [R0]
    Vs = [V0]

    # Set up chart
    fig, ax = plt.subplots()
    # Determine the amount of arrows
    x = [i/15 for i in range(0, 11, 1)]
    y = x
    X, Y = np.meshgrid(x, y)
    U = []
    V = []

    # Susceptible & Vaccinated
    if plotNum == 1:
        # Run the simulation for each (x, y)
        for iii in range(len(y)):
            for ii in range(len(x)):

                # Initialize with particular (x, y)
                Ss = [x[ii]]
                Vs = [y[iii]]
                for i in range(0, timesteps - 1):
                    Ss.append(Ss[i] * ((-1 * beta) * Is[i] - nu + 1) + l * Vs[i])
                    Is.append(Is[i] * (beta * Ss[i] - gamma + 1) + iota * Vs[i])
                    Rs.append(Rs[i] + (gamma * Is[i]))
                    Vs.append(Vs[i] * ((-1 * iota) - l + 1) + nu * Ss[i])
                
                # Capture the final values to create vectors
                U.append(Ss[-1])
                V.append(Vs[-1])

        # Populate plot
        magnitude = np.sqrt(np.array(U) ** 2 + np.array(V) ** 2)
        plt.quiver(X, Y, U, V, magnitude)

        # Decorate plot
        plt.xlabel("Susceptible")
        plt.ylabel("Vaccinated")
        plt.title("Susceptible & Vaccinated")

    # Infected & Vaccinated (See plotNum == 1 for comments)
    if plotNum == 2:
        for iii in range(len(y)):
            for ii in range(len(x)):
                Is = [x[ii]]
                Vs = [y[iii]]
                for i in range(0, timesteps - 1):
                    Ss.append(Ss[i] * ((-1 * beta) * Is[i] - nu + 1) + l * Vs[i])
                    Is.append(Is[i] * (beta * Ss[i] - gamma + 1) + iota * Vs[i])
                    Rs.append(Rs[i] + (gamma * Is[i]))
                    Vs.append(Vs[i] * ((-1 * iota) - l + 1) + nu * Ss[i])
                U.append(Is[-1])
                V.append(Vs[-1])
        magnitude = np.sqrt(np.array(U) ** 2 + np.array(V) ** 2)
        plt.quiver(X, Y, U, V, magnitude)
        plt.xlabel("Infected")
        plt.ylabel("Vaccinated")
        plt.title("Infected & Vaccinated")

    # Removed & Vaccinated (See plotNum == 1 for comments)
    if plotNum == 3:
        for iii in range(len(y)):
            for ii in range(len(x)):
                Rs = [x[ii]]
                Vs = [y[iii]]
                for i in range(0, timesteps - 1):
                    Ss.append(Ss[i] * ((-1 * beta) * Is[i] - nu + 1) + l * Vs[i])
                    Is.append(Is[i] * (beta * Ss[i] - gamma + 1) + iota * Vs[i])
                    Rs.append(Rs[i] + (gamma * Is[i]))
                    Vs.append(Vs[i] * ((-1 * iota) - l + 1) + nu * Ss[i])
                U.append(Rs[-1])
                V.append(Vs[-1])
        magnitude = np.sqrt(np.array(U) ** 2 + np.array(V) ** 2)
        plt.quiver(X, Y, U, V, magnitude)
        plt.xlabel("Removed")
        plt.ylabel("Vaccinated")
        plt.title("Removed & Vaccinated")
    
    fig.text(.1, .9, f"{timesteps} timesteps")
    plt.show()
    return

# Vector Field of Lotka Volterra Model
vectorFieldLV(1.1, .4, .1, .4, 30, .0001)

# Vector Field of SIRV Model
vectorFieldSIRV(0.95, 0.05, 0, 0, 0.2, .1, .05, .0001, .005, 30, 1)
"""
vectorFieldSIRV(0.95, 0.05, 0, 0, 0.2, .1, .05, .0001, .005, 30, 2)
vectorFieldSIRV(0.95, 0.05, 0, 0, 0.2, .1, .05, .0001, .005, 30, 3)
"""