import matplotlib.pyplot as plt

import numpy as np

def EulerMethod(r, V, x, y, endtime, delta_t = 1):
    Xs = [x]
    Ys = [y]
    i = 0
    points = endtime/delta_t

    #Indexing is now out determined by the amount of outputs, rather than the number of timesteps
    while i < points:

        nextX = Xs[i] + ((r/V) * (Xs[i] - Ys[i]))
        nextY = Ys[i] + ((r/V) * (Ys[i] - Xs[i]))
        Xs.append(nextX)
        Ys.append(nextY)
        i += 1

    #Correction of x-axis points for steps taken: the x-coordinate = Output Index * delta_t
    #EX: if delta_t = 0.5, Output #20 happens at timestep #10
    t_axis = [j*delta_t for j in range(len(Xs))]

    fig, ax = plt.subplots()
    ax.plot(t_axis, Xs, color = 'cyan')
    ax.plot(t_axis, Ys, color = 'navy')
    plt.legend(['X\u209C', 'Y\u209C'])
    plt.title(f"t = 0, 1, 2, ..., {int(points)}")
    plt.xlabel("t")
    plt.ylabel("")
    fig.text(.75, .95, f"{int(points)} outputs")
    fig.text(.75, .90, f"\u0394t = {delta_t}")
    plt.show()
    return

def LotkaVolterra(x, y, a, b, d, c, endtime, delta_t = 1):
    """
        x: inital prey population
        y: initial predator population
        a: prey growth rate
        b: prey death rate
        d: predator growth rate
        c: predator death rate
        timesteps: time passed
        delta_t: spcace between each point
    """
    Xs = [x]
    Ys = [y]
    points = endtime/delta_t
    i = 0

    while i < points:
        Xs.append(Xs[i]*(a*delta_t - b*Ys[i]*delta_t + 1))
        Ys.append(Ys[i]*(d*Xs[i]*delta_t - c*delta_t + 1))
        i += 1

    t_axis = [ii*delta_t for ii in range(0, len(Xs))]
    
    def LVPlotOne():
        fig, ax = plt.subplots(figsize = (10, 5))
        ax.plot(t_axis, Xs, color = "cyan", linestyle = 'solid')
        ax.plot(t_axis, Ys, color = "Navy", linestyle = 'solid')
        plt.legend(['Prey (x\u209C)', 'Predators (y\u209C)'], loc = (0, 1))
        plt.title("Lotka-Volterra ODEs")
        fig.text(.75, .95, f"{int(points)} Computes")
        fig.text(.75, .9, f"\u0394t = {delta_t}")
        plt.xlabel('Timesteps')
        plt.ylabel('Population')
        plt.show()
        return
    
    def LVPlotTwo():
        fig, ax = plt.subplots()
        ax.plot(Xs, Ys, color = "black", linestyle = 'solid')
        plt.title("Lotka-Volterra\nPhase Portrait")
        fig.text(.75, .95, f"{int(points)} outputs")
        fig.text(.75, .9, f"\u0394t = {delta_t}")
        plt.ylabel('Predator Population (y\u209C)')
        plt.xlabel('Prey Population (x\u209C)')
        plt.show()
        return
    
    LVPlotOne()
    LVPlotTwo()

    return

def SIR(s, i, b, c, endtime, delta_t = 1):
    n = s + i
    Ss = [s]
    Is = [i]
    Rs = [i]

    points = endtime/delta_t
    j = 0

    while j < points:
        nextS = Ss[j]*((-b/n)*Is[j]*delta_t + 1)
        nextI = Is[j]*((b/n)*Ss[j]*delta_t - c*delta_t + 1)
        nextR = n - (nextS + nextI)
        Ss.append(nextS)
        Is.append(nextI)
        Rs.append(nextR)
        j += 1

    t_axis = [j*delta_t for j in range(len(Rs))]

    def SIRPlotOne():
        fig, ax = plt.subplots()
        ax.plot(t_axis, Ss, color = 'green', linestyle = 'solid')
        ax.plot(t_axis, Is, color = 'orange', linestyle = 'solid')
        ax.plot(t_axis, Rs, color = 'red', linestyle = 'solid')
        plt.legend(['S\u209C', 'I\u209C', 'R\u209C'])
        plt.title("SIR Plot I")
        fig.text(.75, .95, f"{int(points)} outputs")
        fig.text(.75, .9, f"\u0394t = {delta_t}")
        plt.xlabel('t')
        plt.ylabel('Populations')
        plt.show()
        return
    
    def SIRPlotTwo():
        fig, axs = plt.subplots(1, 3, figsize = (9, 4))
        plt.suptitle("SIR Plot II", y = .9)
        fig.tight_layout(pad = 2)
        fig.text(.8, .90, f"{int(points)} outputs")
        fig.text(.8, .85, f"\u0394t = {delta_t}")
        fig.text(.1, .90, "Graphs of S, I, and R")
        fig.text(.1, .85, "Populations")
        axs[0].plot(Is, Ss, color = 'black')
        axs[0].set_xlabel('I\u209C')
        axs[0].set_ylabel('S\u209C')
        axs[1].plot(Ss, Rs, color = 'black')
        axs[1].set_xlabel('S\u209C')
        axs[1].set_ylabel('R\u209C')
        axs[2].plot(Is, Rs, color = 'black')
        axs[2].set_xlabel('I\u209C')
        axs[2].set_ylabel('R\u209C')
        plt.show()

    SIRPlotOne()
    SIRPlotTwo()
    return
    
def logisticMap(rLB, rUB, steps, x, computes, keep = 100):
    
    Rs = [r/steps for r in range(int(rLB*steps), int(rUB*steps))]
    fig, ax = plt.subplots(figsize = (12, 6))

    for r in Rs:
        Xs = [x]   
        i = 0
        while i < computes:
            Xs.append(r * Xs[i] * (1 - Xs[i]))
            i += 1
        rs = [r]*keep
        ax.scatter(rs, Xs[-keep:], color = 'grey', marker = '.', s =.02)
    
    plt.title("Bifurcation Diagram", fontsize = 20)
    fig.text(.1, .925, f"Steps = {steps}")
    fig.text(.1, .9, f"Computes = {computes}")
    plt.ylabel("x", fontsize = 15)
    plt.xlabel(f"{rLB} < r < {rUB}", fontsize = 15)
    plt.show()


def main(odes = True, LM = True, allLM = True):

    if odes is True:
        EulerMethod(2, 40, 10, 20, 20, .1)
        LotkaVolterra(10, 1, 1.1, .4, .1, .4, 30, .0001)
        SIR(997, 3, .4, .04, 100, .01)
    if LM is True:
        logisticMap(0, 1, 1000, .5, 200)
        logisticMap(1, 2, 1000, .5, 200)
        logisticMap(2, 3, 1000, .5, 200)
        logisticMap(3, 3.44949, 1000, .5, 200)
        logisticMap(3.44949, 3.56995, 1000, .5, 200)
    if LM is True and allLM is True:
        logisticMap(1, 4, 1000, .5, 200)
    






