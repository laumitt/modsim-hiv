import numpy as np
import matplotlib.pyplot as plt
import pandas
'''
R is the number of activated uninfected cells
L is the number of latently infected cells
E is the number of activated infected cells
V is the number of free virions
gamma is the rate at which uninfected cells arise
tau is the proportion of cells activated
mu is the HIV-independent death rate of uninfected cells
beta is the rate of infection of cells per virion
rho is the proportion of latently infected cells upon infection
alpha is the activation rate of latently infected cells
delta is the death rate of actively infected cells
pi is the rate of porduction of virions by an actively infected cells
sigma is the removal rate of cell-free virus
'''

sim_length = 120 # days
# parameters that can change - initial values given by paper
total_cells = 1000
gamma = 1.36
mu = 1.36 * 10**-3
tau = 0.2
beta = 0.00027
rho = 0.1
alpha = 3.6 * 10**-2
sigma = 2
delta = 0.33
pi = 100
deltat = 0.75 # days

class Plasma:
    def __init__(self):
        self.R = total_cells * tau
        self.L = 0
        self.E = 0
        self.V = 4 * 10**-7
        #self.unactivated_cells = total_cells - self.R - self.L - self.E
        print('initialized')

    def update(self, deltat):
        print('updating')
        global total_cells
        dRdt = gamma * tau - mu * self.R - beta * self.R * self.V
        dLdt = rho * beta * self.R * self.V - mu * self.L - alpha * self.L
        dEdt = (1 - rho) * beta * self.R * self.V + alpha * self.L - delta * self.E
        dVdt = pi * self.E - sigma * self.V
        print('dRdt {} dLdt {} dEdt {} dVdt {}'.format(dRdt, dLdt, dEdt, dVdt))

        self.R += dRdt * deltat
        self.L += dLdt * deltat
        self.E += dEdt * deltat
        self.V += dVdt * deltat
        print("values: R {}, L {}, E {}, V {}".format(self.R, self.L, self.E, self.V))
        return self.R, self.L, self.E, self.V

def main():
    num_steps = sim_length / deltat
    print(num_steps)
    plasma = Plasma()
    graph_array = [[-1, plasma.R, plasma.L, plasma.E, plasma.V]]
    time_list = []
    R_list = []
    L_list = []
    E_list = []
    V_list = []
    for i in range(int(num_steps)):
        R, L, E, V = plasma.update(deltat)
        time_list.append(i)
        R_list.append(R)
        L_list.append(L)
        E_list.append(E)
        V_list.append(V)
        graph_array.append([i, R, L, E, V])
    for x in time_list:
        time_list[x] = time_list[x] * deltat
    array_collection = [time_list, R_list, L_list, E_list, V_list]
    # code below adapted from : https://matplotlib.org/gallery/api/two_scales.html
    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('R', color=color)
    ax1.plot(array_collection[0], array_collection[1], color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:blue'
    ax2.set_ylabel('L and E', color=color)  # we already handled the x-label with ax1
    ax2.semilogy(array_collection[0], array_collection[2], color='blue', label='L')
    ax2.semilogy(array_collection[0], array_collection[3], color='green', label='E')
    ax2.tick_params(axis='y', labelcolor=color)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()
