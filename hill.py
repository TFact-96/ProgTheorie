import numpy as np
import math
import copy
import random
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from matplotlib import style
import itertools


class Node:
    # object for an node in the amino chain
    def __init__(self, x, y):

        # coords
        self.x = x
        self.y = y
        self.n = 0

        # atom type (default P)
        self.type = "P"

        # the fold that this atom makes to the next atom (default 0)
        self.fold_code = 0

        # All neighbours of a single node
        self.neighbours = []

    # print the type if printing the object
    def __repr__(self):
        return f"{self.type}"


class NodeChain:
    def __init__(self, amino):

        self.stability = 0
        self.old_stability = 0
        self.amino = amino
        self.first_node = Node(0, 0)
        self.first_node.type = amino[0]
        self.chain = {0: self.first_node}

        self.cc_bonds = []
        self.ch_bonds = []

        # available moves
        self.moves = [[1, 0], [0, 1], [-1, 0], [0, -1]]

        self.diagonal_moves = [[1, 1], [1, -1], [-1, 1], [-1, -1]]

    def link_neighbours(self):
        temp_stability = 0

        for node_key, node in self.chain.items():

            for node_neighbour in node.neighbours:

                if node.type == "C" and node_neighbour.type == "C":
                    temp_stability -= 5
                    self.cc_bonds.append(
                        [[node.x, node_neighbour.x], [node.y, node_neighbour.y]]
                    )
                elif node.type == "C":
                    temp_stability -= 1
                    self.ch_bonds.append(
                        [[node.x, node_neighbour.x], [node.y, node_neighbour.y]]
                    )
                else:
                    temp_stability -= 1
                    self.ch_bonds.append(
                        [[node.x, node_neighbour.x], [node.y, node_neighbour.y]]
                    )

        self.stability = temp_stability

    def update_neighbours(self):

        for a, b in itertools.combinations(list(self.chain.items()), 2):
            node_1_index = a[0]
            node_2_index = b[0]
            node_1 = a[1]
            node_2 = b[1]

            # Deduct coordinates to create vector and check vector length
            if (
                np.linalg.norm(
                    np.array([node_1.x, node_1.y]) - np.array([node_2.x, node_2.y])
                )
                == 1.0
            ) and (abs(node_1_index - node_2_index) != 1):

                node_1.neighbours.append(node_2)

        self.link_neighbours()

    def reset_neighbours(self):

        for node_key, node in self.chain.items():
            node.neighbours = []

    def create_chain(self):
        index = 0

        # Create a random chain
        while len(self.chain) < len(self.amino):
            current_node = self.chain[index]
            random_move = np.random.randint(len(self.moves))

            new_x = current_node.x + self.moves[random_move][0]
            new_y = current_node.y + self.moves[random_move][1]

            if self.no_overlap(new_x, new_y):
                index += 1
                self.chain[index] = Node(new_x, new_y)
                self.chain[index].n = index
                self.chain[index].type = self.amino[index]

    def no_overlap(self, x, y):

        for node_key, node in self.chain.items():

            if [node.x, node.y] == [x, y]:
                return False

        return True

    def plot_chain(self):
        plot_x = []
        plot_y = []
        color = []

        for node_key, node in self.chain.items():

            plot_x.append(node.x)
            plot_y.append(node.y)
            node_type = node.type

            if node_type == "C":
                color.append("yellow")
            elif node_type == "H":
                color.append("red")
            else:
                color.append("blue")

        # for cc_bond in self.cc_bonds:
        #     plt.plot(cc_bond[0], cc_bond[1], "y--")

        # for ch_bond in self.ch_bonds:
        #     plt.plot(ch_bond[0], ch_bond[1], "r--")

        # plot atomtype name at its node coord
        for i in range(len(plot_x)):
            plt.text(plot_x[i] + 0.05, plot_y[i] + 0.05, self.chain[i].type)

        plt.plot(plot_x, plot_y)
        plt.scatter(plot_x, plot_y, color=color)
        plt.show()

    # checks if a point has a node or not
    def check_point(self, array):
        for index, nodes in self.chain.items():
            if str(array) == str(np.array([nodes.x, nodes.y])):
                return False
        return True

    def pull_move(self, node):

        node_i_coords = np.array([node.x, node.y])
        node_i1 = self.chain[int(node.n) + 1]
        node_i1_coords = np.array([node_i1.x, node_i1.y])
        vector1 = node_i1_coords - node_i_coords

        chain_copy = copy.deepcopy(self.chain)
        self.old_stability = self.stability

        checker = {
            "[1, 1]": [True, [1, 1]],
            "[1, -1]": [True, [1, -1]],
            "[-1, 1]": [True, [-1, 1]],
            "[-1, -1]": [True, [-1, -1]],
        }

        for d_move in self.diagonal_moves:
            if not self.check_point(node_i_coords + np.array(d_move)):
                checker[str(d_move)][0] = False

        for d_check, d_bool in checker.items():
            L = node_i_coords + np.array(d_bool[1])
            C = L - vector1

            # checks pull move requirements
            if (
                (d_bool[0] == True)
                and (np.linalg.norm(L - node_i1_coords) == 1.0)
                and (self.check_point(C))
            ):
                # residue of chain follows in footsteps
                for index in range(int(node.n - 1)):

                    residue_node = self.chain[index]
                    residue_node_next = self.chain[index + 2]

                    residue_node.x = residue_node_next.x
                    residue_node.y = residue_node_next.y

                # node moves to L
                node.x = L[0]
                node.y = L[1]

                # Previous node moves to C
                previous_node = self.chain[int(node.n) - 1]
                previous_node.x = C[0]
                previous_node.y = C[1]
                break

        self.reset_neighbours()
        self.update_neighbours()

        if self.stability > self.old_stability:
            self.chain = chain_copy

            self.reset_neighbours()
            self.update_neighbours()

            print(self.stability)

    def random_pull(self):

        for x in range(100):
            d = random.randint(1, 10)
            self.pull_move(self.chain[d])


k = NodeChain("CHHCHHCHHCHC")
k.create_chain()
k.random_pull()
k.plot_chain()
