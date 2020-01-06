import math
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
import statistics as stat

video = False
amino = "HHPHHHPHHHHHHHHHHHHPPPPHPHPPHPPHPHPPHPPHPHPPHPHPPPPHPHPHPPPHH"

class lattice_SAW:
    def __init__(self, amino):
        self.overlap_counter = 0
        self.amino = amino
        self.color = ['red' if x == "H" else 'blue' for x in self.amino]
        self.N = len(amino)
        self.state = np.array([[0, 0]])
        self.moves = [[1, 0], [0, 1], [-1, 0], [0, -1]]
        self.total_nodes = 1

    def walk(self):
        while self.total_nodes < self.N:
            new_node = self.state[-1] + self.moves[np.random.randint(len(self.moves))]
            overlap = False

            # check if this node overlaps the chain
            for node in self.state:
                if node[0] == new_node[0] and node[1] == new_node[1]:
                    overlap = True
                    self.overlap_counter += 1

            if self.overlap_counter > 100:
                print("Chain got stuck!")
                break

            # if new node doesnt overlap
            if overlap == False:
                self.state = np.vstack((self.state, new_node))
                self.total_nodes += 1

                if video is True:
                    x = []
                    y = []

                    for i in range(self.total_nodes):
                        x.append(self.state[i][0])
                        y.append(self.state[i][1])

                    plt.plot(x, y, "-", markersize=1, color='black', zorder=1)
                    plt.scatter(x, y, color=self.color[0:len(x)], zorder=2)
                    plt.ylabel('y-as')
                    plt.xlabel('x-as')
                    plt.grid(b=True, axis='both')
                    plt.draw()
                    plt.pause(0.001)
                    plt.clf()

        x = []
        y = []

        for i in range(self.total_nodes):
            x.append(self.state[i][0])
            y.append(self.state[i][1])

        plt.plot(x, y, "-", markersize=1, color='black', zorder=1)
        plt.scatter(x, y, color=self.color[0:len(x)], zorder=2)
        plt.ylabel('y-as')
        plt.xlabel('x-as')
        plt.grid(b=True, axis='both')
        plt.show()

chain = lattice_SAW(amino)
chain.walk()
