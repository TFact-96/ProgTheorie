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
        self.state = np.array([[0, 0]])
        self.moves = [[1, 0], [0, 1], [-1, 0], [0, -1]]
        self.move_code = [1, 2, -1, -2]
        self.output_moves = []
        self.total_nodes = 1

    def walk(self):
        while self.total_nodes < len(self.amino):
            # get random move of node
            random_move_index = np.random.randint(len(self.moves))
            new_node = self.state[-1] + self.moves[random_move_index]
            move_code = self.move_code[random_move_index]

            overlap = False
            # check if this node overlaps the chain
            for node in self.state:
                if node[0] == new_node[0] and node[1] == new_node[1]:
                    overlap = True
                    self.overlap_counter += 1

            # check if chain got stuck
            if self.overlap_counter > 100:
                print("Chain got stuck!")
                break

            # if new node doesnt overlap, add to state
            if overlap == False:
                self.overlap_counter = 0
                self.output_moves.append([self.amino[self.total_nodes - 1], move_code])
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

        # Print desired output moves
        for move in self.output_moves:
            print(move[0], move[1])

        stability = 0

        # check each H node
        for node in range(self.total_nodes):
            if self.amino[node] == "H":
                # compare with each other H node
                for compare_node in range(self.total_nodes):
                    # if index length greater than 1 (not neighbor nodes)
                    if self.amino[compare_node] == "H" and (abs(compare_node - node) > 1):
                        # calculate distance
                        dist = math.sqrt((self.state[node][0] - self.state[compare_node][0])**2 + (self.state[node][1] - self.state[compare_node][1])**2)

                        # H-bond if distance is 1
                        if dist == 1:
                            stability -= 0.5
                            plt.plot([self.state[node][0], self.state[compare_node][0]], [self.state[node][1], self.state[compare_node][1]], "--", markersize=1, color='orange', zorder=1)

        # Plot final result
        x = []
        y = []

        for i in range(self.total_nodes):
            x.append(self.state[i][0])
            y.append(self.state[i][1])

        plt.title("Amino acid chain")
        plt.plot(x, y, "-", linewidth=3, color='black', zorder=1, label=f"Amino chain (Stability: {stability})")
        plt.scatter(x, y, color=self.color[0:len(x)], zorder=2)
        plt.ylabel('y-as')
        plt.xlabel('x-as')
        plt.legend()
        plt.show()

chain = lattice_SAW(amino)
chain.walk()
