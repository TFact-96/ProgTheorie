import math
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
import statistics as stat

video = False
amino = "HHPHHHPHPHHHPHHHHHHHHHHHHHHHHHHHHHHHHHHHHHPHPHPHPPPHPPPHHHHHHHHHHHHHH"
optimalization_tries = 100

class lattice_SAW:
    def __init__(self, amino, optimalization_tries):
        self.overlap_counter = 0
        self.amino = amino
        self.color = ['red' if x == "H" else 'blue' for x in self.amino]
        self.state = np.array([[0, 0]])
        self.moves = [[1, 0], [0, 1], [-1, 0], [0, -1]]
        self.move_code = [1, 2, -1, -2]
        self.output_moves = []
        self.node_number = 0
        self.optimalization_tries = optimalization_tries
        self.stability = 0

    def walk(self):

        # while amount of nodes is smaller than the whole amino string length
        while len(self.state) < len(self.amino):

            # handle new node coords
            new_node, move_code = self.bond_optimalization(self.optimalization_tries)

            # if the node doesn't overlap its own chain
            if not self.check_overlap(new_node):

                # append node to data output with its fold code
                self.output_moves.append([self.amino[len(self.state)], move_code])

                # append the new node to the existing chain
                self.state = np.vstack((self.state, new_node))

        # Print desired output data
        for move in self.output_moves:
            print(move[0], move[1])

        x = []
        y = []

        # append all nodes in coord list
        for i in range(len(self.state)):
            x.append(self.state[i][0])
            y.append(self.state[i][1])

        # plot final result
        plt.title("Amino acid chain")
        #plt.plot([self.state[node][0], self.state[compare_node][0]], [self.state[node][1], self.state[compare_node][1]], "--", markersize=1, color='orange', zorder=1)
        plt.plot(x, y, "-", linewidth=3, color='black', zorder=1, label=f"Amino chain (Stability: {self.stability})")
        plt.scatter(x, y, color=self.color[0:len(x)], zorder=2)
        plt.ylabel('y-as')
        plt.xlabel('x-as')
        plt.legend()
        plt.show()

    def get_random_move(self):
        # get random move of new node added to chain
        random_move_index = np.random.randint(len(self.moves))
        # last old node coords + random move
        new_node = self.state[-1] + self.moves[random_move_index]
        # for data output
        move_code = self.move_code[random_move_index]

        return new_node, move_code

    def check_overlap(self, new_node):
        # check if chain got stuck, for if the program already tried random moves 100 times and none accepted.
        if self.overlap_counter > 100:
            print("Chain got stuck!")
            return True

        # check if this node overlaps the chain (new node overlaps other node)
        for node in self.state:
            if node[0] == new_node[0] and node[1] == new_node[1]:
                self.overlap_counter += 1
                return True

        # reset overlap_counter if new node is added
        self.overlap_counter = 0
        return False

    def bond_optimalization(self, tries):
        dist = 2
        i = 0

        # try to repeat 'i' random moves, until another H-node is found that is 1 distance away from own H-node
        # if not found after i times, just add the node with a random move.
        while dist > 1 and i < tries:

            # get new node coords
            new_node, move_code = self.get_random_move()

            # Check for H-bonds (new node doesn't overlap own chain, new node is H, more than 2 nodes already exist)
            if self.amino[len(self.state)] == "H" and not self.check_overlap(new_node):

                # compare own H-node with all other nodes already put into the chain
                for compare_node in range(len(self.state)):

                    # if index length greater than 1 (they're not neighbor nodes) and other node also an H
                    if self.amino[compare_node] == "H" and (abs(compare_node - len(self.state)) > 1):

                        # calculate distance between these nodes
                        dist = math.sqrt((new_node[0] - self.state[compare_node][0])**2 + (new_node[1] - self.state[compare_node][1])**2)

                        # H-bond if distance is 1
                        if dist <= 1:
                            print(f"node1: {len(self.state)}, node2: {compare_node}")

                            # decrease stability
                            self.stability -= 1

                            # plot this H-bond between the two nodes
                            plt.plot([new_node[0], self.state[compare_node][0]], [new_node[1], self.state[compare_node][1]], "--", markersize=1, color='orange', zorder=1)
                            break
            i += 1

        return new_node, move_code

chain = lattice_SAW(amino, optimalization_tries)
chain.walk()
