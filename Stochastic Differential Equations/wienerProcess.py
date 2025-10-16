import random as r
import numpy as np
import matplotlib.pyplot as plt

#wienerProcess(timesteps=20, plots = 3, seed=42)

def wienerProcess(timesteps, delta_t, plots = 5, seed = 42):

    r.seed(42)
    computes = timesteps/delta_t

    fig, ax = plt.subplots()
    for plot in range(0, plots):
        W = [0]
        i = 0
        while i < computes:
            W.append(W[i]+ np.sqrt(delta_t) * r.normalvariate(0, delta_t))
            i += 1
            
            # Try list comprehension

        t_axis = [t*delta_t for t in range(0, len(W))]
        ax.plot(t_axis, W)

    plt.title("Wiener Process")
    fig.text(.7, .9, f"Seed: {seed} ")
    fig.text(.2, .9, f"{int(timesteps/delta_t)} computes")
    plt.xlabel('t')
    plt.ylabel('W(t)')
    plt.show()

    return

wienerProcess(20, .001)