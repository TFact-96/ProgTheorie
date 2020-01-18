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
        self.chain = self.create_chain()

        self.best_chain = []
        self.best_hh_bonds = []

        # available moves
        self.moves = [[1, 0], [0, 1], [-1, 0], [0, -1]]

        self.diagonal_moves = [[1, 1], [1, -1], [-1, 1], [-1, -1]]

    def link_neighbours(self, chain):
        temp_stability = 0
        hh_bonds = []

        for node_key, node in chain.items():

            for node_neighbour in node.neighbours:

                if node.type == "H" and node_neighbour.type == "H":
                    temp_stability -= 1
                    hh_bonds.append(
                        [[node.x, node_neighbour.x], [node.y, node_neighbour.y]]
                    )

        return temp_stability, hh_bonds

    def update_neighbours(self, chain):

        for a, b in itertools.combinations(list(chain.items()), 2):
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

        stability, hh_bonds = self.link_neighbours(chain)
        return stability, hh_bonds

    def reset_neighbours(self, chain):

        self.hh_bonds = []
        for node_key, node in chain.items():
            node.neighbours = []

    def chain_stuck(self, x, y, chain):
        if (
            (not self.no_overlap(x + 1, y, chain))
            and (not self.no_overlap(x, y + 1, chain))
            and (not self.no_overlap(x - 1, y, chain))
            and (not self.no_overlap(x, y - 1, chain))
        ):
            return True

        return False

    def create_chain(self):
        index = 0
        self.moves = [[1, 0], [0, 1], [-1, 0], [0, -1]]
        first_node = Node(0, 0)
        first_node.type = self.amino[0]
        chain = {0: first_node}

        # Create a random chain
        while len(chain) < len(self.amino):
            current_node = chain[index]
            random_move = np.random.randint(len(self.moves))

            new_x = current_node.x + self.moves[random_move][0]
            new_y = current_node.y + self.moves[random_move][1]

            if self.no_overlap(new_x, new_y, chain):
                index += 1
                chain[index] = Node(new_x, new_y)
                chain[index].n = index
                chain[index].type = self.amino[index]

            if self.chain_stuck(new_x, new_y, chain):
                chain = {0: first_node}
                index = 0

        return chain

    def no_overlap(self, x, y, chain):

        for node_key, node in chain.items():

            if [node.x, node.y] == [x, y]:
                return False

        return True

    def plot_chain(self):
        plot_x = []
        plot_y = []
        color = []

        for node_key, node in self.best_chain.items():

            plot_x.append(node.x)
            plot_y.append(node.y)
            node_type = node.type

            if node_type == "C":
                color.append("yellow")
            elif node_type == "H":
                color.append("red")
            else:
                color.append("blue")

        for hh_bond in self.best_hh_bonds:
            plt.plot(hh_bond[0], hh_bond[1], "y--")

        # for ch_bond in self.ch_bonds:
        #     plt.plot(ch_bond[0], ch_bond[1], "r--")

        # plot atomtype name at its node coord
        for i in range(len(plot_x)):
            plt.text(plot_x[i] + 0.05, plot_y[i] + 0.05, self.chain[i].type)

        plt.plot(plot_x, plot_y)
        plt.scatter(plot_x, plot_y, color=color)
        plt.show()

    # checks if a point has a node or not
    def check_point(self, array, chain):
        for index, nodes in chain.items():
            if str(array) == str(np.array([nodes.x, nodes.y])):
                return False
        return True

    def create_vectors(self, node, chain):
        node_i_coords = np.array([node.x, node.y])
        node_i1 = chain[int(node.n) + 1]
        node_i1_coords = np.array([node_i1.x, node_i1.y])
        vector1 = node_i1_coords - node_i_coords

        return node_i_coords, node_i1, node_i1_coords, vector1

    def pull_move(self, node, chain):

        node_i_coords, node_i1, node_i1_coords, vector1 = self.create_vectors(
            node, chain
        )

        checker = {
            "[1, 1]": [True, [1, 1]],
            "[1, -1]": [True, [1, -1]],
            "[-1, 1]": [True, [-1, 1]],
            "[-1, -1]": [True, [-1, -1]],
        }

        viable_moves = []

        for d_move in self.diagonal_moves:
            if not self.check_point(node_i_coords + np.array(d_move), chain):
                checker[str(d_move)][0] = False

        for checker_key, d_bool in checker.items():
            if d_bool[0] == True:
                viable_moves.append(d_bool[1])

        if viable_moves != []:
            random_move = random.choice(viable_moves)

            L = node_i_coords + np.array(random_move)
            C = L - vector1

            # checks pull move requirements
            if (np.linalg.norm(L - node_i1_coords) == 1.0) and (
                self.check_point(C, chain)
            ):
                # residue of chain follows in footsteps
                for index in range(int(node.n - 1)):

                    residue_node = chain[index]
                    residue_node_next = chain[index + 2]

                    residue_node.x = residue_node_next.x
                    residue_node.y = residue_node_next.y

                # node moves to L
                node.x = L[0]
                node.y = L[1]

                # Previous node moves to C
                previous_node = chain[int(node.n) - 1]
                previous_node.x = C[0]
                previous_node.y = C[1]

    def random_pull(self, max_iteration):

        # Current protein chain
        self.current_hilltop = self.create_chain()
        self.current_stability, self.current_hh = self.update_neighbours(
            self.current_hilltop
        )

        # Save as best chain
        self.best_c = self.current_hilltop
        self.best_stab_c = self.current_stability
        self.hh_bonds = self.current_hh
        self.best_chain.append([self.best_c, self.best_stab_c, self.hh_bonds])

        for iteration in range(max_iteration):
            best_c_found = False
            print(iteration)

            for it in range(100):
                print(it)
                for index in range(1, len(self.current_hilltop) - 1):
                    self.pull_move(self.current_hilltop[index], self.current_hilltop)
                    self.reset_neighbours(self.current_hilltop)
                    temp_stability, temp_hh = self.update_neighbours(
                        self.current_hilltop
                    )

                    if temp_stability < self.best_stab_c:
                        self.best_c = self.current_hilltop
                        self.best_stab_c = temp_stability
                        self.hh_bonds = temp_hh
                        best_c_found = True

            if best_c_found:
                self.current_hilltop = self.best_c
                self.current_stability = self.best_stab_c
                self.current_hh = self.hh_bonds
                best_c_found = False

            else:
                self.best_chain.append([self.best_c, self.best_stab_c, self.hh_bonds])
                self.current_hilltop = self.create_chain()
                self.current_stability, self.current_hh = self.update_neighbours(
                    self.current_hilltop
                )

                self.best_c = self.current_hilltop
                self.best_stab_c = self.current_stability
                self.hh_bonds = self.current_hh

    def find_best_c(self):
        for list_ in self.best_chain:
            print(list_[1])


k = NodeChain("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
k.random_pull(100)
k.find_best_c()

