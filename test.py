import math
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
import statistics as stat

video = True
avg = []


# define a dot product function used for the rotate operation
def v_dot(a): return lambda b: np.dot(a, b)


def endtoend():
    stappen = np.arange(10, 50, 5)
    square = np.sqrt(stappen)
    avgR = np.sqrt(avg)
    plt.plot(stappen, avgR, 'o-')
    plt.plot(stappen, square, '-')
    plt.legend(('Simulatie', 'sqrt(N)'))
    plt.xlabel("N")
    plt.ylabel("End to end")
    plt.grid(b=True, axis='both')
    plt.show()


class lattice_SAW:
    def __init__(self, N, l0, amino):
        self.N = N
        self.type = amino
        self.l0 = l0
        # initial configuration. Usually we just use a straight chain as inital configuration
        self.init_state = np.dstack((np.arange(N), np.zeros(N)))[0]
        self.state = self.init_state.copy()
        # 3 possible rotations: rotate angles(90,180,270)
        self.rotate_matrix = np.array([[[0, -1], [1, 0]],
                                       [[-1, 0], [0, -1]],
                                       [[0, 1], [1, 0]]])

    # define pivot algorithm process where t is the number of successful steps
    def walk(self, t, c):
        print(self.init_state)
        afstand = []
        for repeat in range(20):
            acpt = 0
            x = []
            y = []
            while acpt < t:
                pick_pivot = np.random.randint(1, self.N - 1)  # pick a pivot site

                old_chain = self.state[pick_pivot:]
                temp_chain = self.state[0:pick_pivot]

                # pick a symmetry operator
                symtry_oprtr = self.rotate_matrix[np.random.randint(len(self.rotate_matrix))]
                # new chain after symmetry operator
                new_chain = np.apply_along_axis(v_dot(symtry_oprtr), 1, temp_chain - self.state[pick_pivot]) \
                            + self.state[pick_pivot]

                # use cdist function of scipy package to calculate the pair-pair distance between old_chain and new_chain
                overlap = cdist(new_chain, old_chain)
                overlap = overlap.flatten()

                # determinte whether the new state is accepted or rejected
                if len(np.nonzero(overlap)[0]) != len(overlap):
                    continue
                else:
                    self.state = np.concatenate((new_chain, old_chain), axis=0)
                    afstand.append(((self.state[N - 1][0] - self.state[0][0]) ** 2)
                                   + ((self.state[N - 1][1] - self.state[0][1]) ** 2))

                    if video is True:
                        x = []
                        y = []
                        for i in range(N):
                            x.append(self.state[i][0])
                            y.append(self.state[i][1])
                        plt.plot(x, y, "-", markersize=3)

                        plt.ylabel('y-as')
                        plt.xlabel('x-as')

                        plt.xlim(-N, N)
                        plt.ylim(0.5 * -N, 0.5 * N)
                        plt.grid(b=True, axis='both')

                        plt.draw()
                        plt.pause(0.2)
                        plt.clf()

                    acpt += 1

            if video is False:
                for i in range(N):
                    x.append(self.state[i][0])
                    y.append(self.state[i][1])
                plt.plot(x, y, '-')
                plt.xlabel("x-as")
                plt.ylabel("y-as")
                plt.grid(b=True, axis='both')
                plt.show()

            for d in range(c * N):
                del afstand[0]

        avg.append(stat.mean(afstand))

        # place the center of mass of the chain on the origin
        # self.state = self.l0 * (self.state - np.int_(np.mean(self.state, axis=0)))


amino = "HHPHHHPH"
sample_size = 10
c = 1
N = 5  # number of monomers(number of steps)
l0 = 1  # bond length(step length)
t = 100  # number of pivot steps

chain = lattice_SAW(N, l0, amino)
chain.walk(t, c)

# for N in range(10, 50, 5):
#    chain = lattice_SAW(N, l0)
#    chain.walk(t, c)

# endtoend()
