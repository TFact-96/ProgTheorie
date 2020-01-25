import math
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
import statistics as stat

# from nn import DQNAgent

video = False
amino = "HHPHHHPHPHPH"
optimalization_tries = 10


class lattice_SAW:
    def __init__(self, amino, optimalization_tries):
        self.overlap_counter = 0
        self.amino = amino
        self.color = ['red' if x == "H" else 'blue' for x in self.amino]
        self.state = np.array([[0, 0]])
        self.moves = [[1, 0], [0, 1], [-1, 0], [0, -1]]  # right, up, left, down
        self.move_code = [1, 2, -1, -2]
        self.output_moves = []
        self.node_number = 0
        self.optimalization_tries = optimalization_tries
        self.stability = 0

    def walk(self):
        while self.node_number < len(self.amino):
            new_node, move_code = self.bond(self.optimalization_tries)

            if not self.overlap(new_node):
                self.overlap_counter = 0
                self.output_moves.append([self.amino[self.node_number], move_code])
                self.state = np.vstack((self.state, new_node))
                self.node_number += 1

            # Plot final result
            x = []
            y = []

            for i in range(self.node_number):
                x.append(self.state[i][0])
                y.append(self.state[i][1])

            plt.title("Amino acid chain")
            plt.plot(x, y, "-", linewidth=3, color='black', zorder=1,
                     label=f"Amino chain (Stability: {self.stability})")
            plt.scatter(x, y, color=self.color[0:len(x)], zorder=2)
            plt.ylabel('y-as')
            plt.xlabel('x-as')
            plt.legend()
            plt.show()

    def move(self):
        # move the node
        random_move = str(input())
        if random_move == 'right':
            random_move_index = 0
        elif random_move == 'up':
            random_move_index = 1
        elif random_move == 'left':
            random_move_index = 2
        elif random_move == 'down':
            random_move_index = 3
        else:
            random_move_index = 7

        new_node = self.state[-1] + self.moves[random_move_index]
        move_code = self.move_code[random_move_index]

        return new_node, move_code

    def overlap(self, new_node):
        # check if chain got stuck
        if self.overlap_counter > 100:
            print("Chain got stuck!")
            return True

        # check if this node overlaps the chain
        for node in self.state:
            if node[0] == new_node[0] and node[1] == new_node[1]:
                self.overlap_counter += 1
                return True

        return False

    def bond(self):
        new_node, move_code = self.move()

        # Optimalize H-bonds
        if not self.overlap(new_node) and self.amino[self.node_number] == "H" and self.node_number > 2:
            # compare with each other H node
            for compare_node in range(self.node_number + 1):
                # if index length greater than 1 (not neighbor nodes) and other node also an H
                if self.amino[compare_node] == "H" and (abs(compare_node - self.node_number) > 1):
                    # calculate distance
                    dist = math.sqrt((new_node[0] - self.state[compare_node][0]) ** 2 + (
                            new_node[1] - self.state[compare_node][1]) ** 2)

                    # H-bond if distance is 1
                    if dist == 1:
                        self.stability -= 0.5
                        plt.plot([new_node[0], self.state[compare_node][0]],
                                 [new_node[1], self.state[compare_node][1]], "--", markersize=1, color='orange',
                                 zorder=1)

        return new_node, move_code


chain = lattice_SAW(amino, optimalization_tries)
chain.walk()