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


class Grid_point:
    def __init__(self, filled, coords):
        self.filled = filled
        self.coords = coords
        self.nodes = []

    def add_node(self, node):
        self.nodes.append(node)
        self.filled = True

    def remove_node(self, node):
        self.nodes.remove(node)
        if self.nodes == []:
            self.filled = False

    def __repr__(self):
        return f"{self.nodes}"


class Grid:
    def __init__(self, amino):
        self.grid = {}
        self.grid_chain = {}
        self.hh_bonds = []

        self.amino = amino
        self.best_chain = {}
        self.diagonal_moves = [[1, 1], [1, -1], [-1, 1], [-1, -1]]
        self.moves = [[1, 0], [0, 1], [-1, 0], [0, -1]]

        self.stability = 0

    # Create n x n grid
    def create_grid(self, n):
        for y in range(-n, n + 1):
            for x in range(-n, n + 1):
                self.grid[f"{x, y}"] = Grid_point(False, [x, y])

    def clear_grid(self):
        for key, value in self.grid.items():
            self.grid[key].filled = False
            self.grid[key].nodes = []
            self.stability = 0
            self.hh_bonds = []

    def add_neighbours(self, node):
        x = node.x
        y = node.y
        n = node.n

        if node.type == "P":
            return 0

        if (
            self.grid[f"{x + 1, y}"].filled
            and (abs(self.grid[f"{x + 1, y}"].nodes[0].n - n) > 1)
            and (self.grid[f"{x + 1, y}"].nodes[0].type == "H")
        ):
            if [[x, x + 1], [y, y]] and [[x + 1, x], [y, y]] not in self.hh_bonds:
                self.hh_bonds.append([[x, x + 1], [y, y]])

        if (
            self.grid[f"{x, y + 1}"].filled
            and (abs(self.grid[f"{x, y + 1}"].nodes[0].n - n) > 1)
            and (self.grid[f"{x, y + 1}"].nodes[0].type == "H")
        ):
            if [[x, x], [y, y + 1]] and [[x, x], [y + 1, y]] not in self.hh_bonds:
                self.hh_bonds.append([[x, x], [y, y + 1]])

        if (
            self.grid[f"{x - 1, y}"].filled
            and (abs(self.grid[f"{x - 1, y}"].nodes[0].n - n) > 1)
            and (self.grid[f"{x - 1, y}"].nodes[0].type == "H")
        ):
            if [[x, x - 1], [y, y]] and [[x - 1, x], [y, y]] not in self.hh_bonds:
                self.hh_bonds.append([[x, x - 1], [y, y]])

        if (
            self.grid[f"{x, y - 1}"].filled
            and (abs(self.grid[f"{x, y - 1}"].nodes[0].n - n) > 1)
            and (self.grid[f"{x, y - 1}"].nodes[0].type == "H")
        ):
            if [[x, x], [y, y - 1]] and [[x, x], [y - 1, y]] not in self.hh_bonds:
                self.hh_bonds.append([[x, x], [y, y - 1]])

    def update_neighbours(self):
        self.hh_bonds = []
        for key, value in self.grid_chain.items():
            self.add_neighbours(self.grid[value[0]].nodes[0])

        self.stability = -len(self.hh_bonds)

    def add_point(self, node, n):
        x = node.x
        y = node.y
        self.grid[f"{x, y}"].add_node(node)
        self.grid_chain[n] = [f"{x, y}", [x, y]]

    def clear_point(self, node, n):
        x = node.x
        y = node.y
        self.grid[f"{x, y}"].remove_node(node)

    def transfer_point(self, node1, x2, y2):
        n = node1.n
        self.clear_point(node1, n)
        node1.x = x2
        node1.y = y2
        self.add_point(node1, n)

    def overlap(self, x, y):
        if self.grid[f"{x, y}"].filled:
            return True

        return False

    def chain_stuck(self, x, y):

        if (
            self.overlap(x + 1, y)
            and self.overlap(x - 1, y)
            and self.overlap(x, y + 1)
            and self.overlap(x, y - 1)
        ):
            return True

        return False

    def check_diagonals(self, x, y):
        available_moves = []

        if not self.grid[f"{x + 1, y + 1}"].filled:
            available_moves.append([1, 1])

        if not self.grid[f"{x + 1, y - 1}"].filled:
            available_moves.append([1, -1])

        if not self.grid[f"{x - 1, y + 1}"].filled:
            available_moves.append([-1, 1])

        if not self.grid[f"{x - 1, y - 1}"].filled:
            available_moves.append([-1, -1])

        return available_moves

    def create_chain(self):
        index = 0

        # Create grid
        self.create_grid(36)
        first_node = Node(0, 0)
        first_node.type = self.amino[0]
        self.add_point(first_node, index)

        while index < len(self.amino) - 1:
            current_node_key = self.grid_chain[index][0]
            current_node = self.grid[current_node_key].nodes[0]

            random_move = np.random.randint(len(self.moves))

            new_x = current_node.x + self.moves[random_move][0]
            new_y = current_node.y + self.moves[random_move][1]

            if not self.overlap(new_x, new_y):
                index += 1
                new_node = Node(new_x, new_y)
                new_node.n = index
                new_node.type = self.amino[index]
                self.add_point(new_node, index)

            if self.chain_stuck(new_x, new_y):
                index = 0
                self.clear_grid()
                self.grid_chain = {}
                self.add_point(first_node, index)

    def create_vectors(self, node):
        node_i_coords = [node.x, node.y]
        node_i1 = self.grid_chain[int(node.n) + 1]
        node_i1_coords = np.array([node_i1[1][0], node_i1[1][1]])
        vector1 = node_i1_coords - np.array(node_i_coords)

        return node_i_coords, node_i1_coords, vector1

    def check_requirements(
        self, available_moves, vector1, node_i_coords, node_i1_coords
    ):

        for move in available_moves:
            L = node_i_coords + np.array(move)
            C = L - vector1

            if (
                (not self.overlap(L[0], L[1]))
                and (not self.overlap(C[0], C[1]))
                and (np.linalg.norm(L - node_i1_coords) == 1.0)
            ):
                return L, C, True

        return 0, 0, False

    def pull_move(self, node):

        node_i_coords, node_i1_coords, vector1 = self.create_vectors(node)

        available_moves = self.check_diagonals(node_i_coords[0], node_i_coords[1])

        L, C, check = self.check_requirements(
            available_moves, vector1, node_i_coords, node_i1_coords
        )

        if check:

            # residue of chain follows in footsteps
            for index in range(int(node.n - 1)):

                residue_node_key = self.grid_chain[index][0]
                residue_node = self.grid[residue_node_key].nodes[0]

                residue_node_next_key = self.grid_chain[index + 2][0]
                residue_node_next = self.grid[residue_node_next_key].nodes[0]

                self.transfer_point(
                    residue_node, residue_node_next.x, residue_node_next.y
                )

            # node moves to L
            self.transfer_point(node, L[0], L[1])

            # Previous node moves to C
            previous_node_key = self.grid_chain[int(node.n) - 1][0]
            previous_node = self.grid[previous_node_key].nodes[0]

            self.transfer_point(previous_node, C[0], C[1])
        self.update_neighbours()

    def hill_climber(self, max_iteration):
        # Current protein chain
        self.create_chain()
        self.current_hilltop = self.grid_chain
        self.current_grid = self.grid

        self.update_neighbours()
        self.current_stability = self.stability
        self.current_hh = self.hh_bonds

        # Save as best chain
        self.best_c = self.current_hilltop
        self.best_stab_c = self.current_stability
        self.best_grid = self.current_grid
        self.best_hh_bonds = self.current_hh
        self.best_chain[self.best_stab_c] = [self.best_c, self.best_hh_bonds]

        for iteration in range(max_iteration):
            best_c_found = False
            print(iteration)

            for it in range(100):
                for index in range(1, len(self.current_hilltop) - 1):

                    self.pull_move(
                        self.current_grid[self.current_hilltop[index][0]].nodes[0]
                    )
                    self.update_neighbours()
                    temp_stability, temp_hh, temp_grid = (
                        self.stability,
                        self.hh_bonds,
                        self.current_grid,
                    )

                    if temp_stability < self.best_stab_c:
                        self.best_c = self.grid_chain
                        self.best_stab_c = temp_stability
                        self.best_hh_bonds = temp_hh
                        self.best_grid = temp_grid
                        best_c_found = True

            if best_c_found:
                self.current_hilltop = self.best_c
                self.current_stability = self.best_stab_c
                self.current_hh = self.best_hh_bonds
                self.current_grid = self.best_grid
                best_c_found = False

            else:
                self.best_chain[self.best_stab_c] = [self.best_c, self.best_hh_bonds]
                self.create_chain()
                self.current_hilltop = self.grid_chain
                self.current_grid = self.grid
                self.current_hh = self.hh_bonds
                self.update_neighbours()

                self.current_stability = self.stability
                self.current_hh = self.hh_bonds

                self.best_c = self.current_hilltop
                self.best_stab_c = self.current_stability
                self.hh_bonds = self.current_hh

    def plot_chain(self):
        x = []
        y = []
        index = 0
        for key, l in self.best_c.items():
            x.append(l[1][0])
            y.append(l[1][1])

        for hh_bond in self.hh_bonds:
            plt.plot(hh_bond[0], hh_bond[1], "y--")

        plt.plot(x, y, "ro-")
        for x_p, y_p in zip(x, y):

            if self.amino[index] == "H":
                plt.plot(x_p, y_p, "bo-")
            else:
                plt.plot(x_p, y_p, "ro-")
            index += 1

        plt.show()

    def find_best_c(self):
        best_chain_key = min(self.best_chain.keys())
        best_chain_double = self.best_chain[best_chain_key]

        best_chain = best_chain_double[0]
        best_hh_bonds = best_chain_double[1]

        x = []
        y = []

        for key, value in best_chain.items():
            x.append(value[1][0])
            y.append(value[1][1])

        for hh_bond in best_hh_bonds:
            plt.plot(hh_bond[0], hh_bond[1], "y--")

        plt.plot(x, y, "bo")
        plt.show()


k = Grid("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
k.hill_climber(2)
k.find_best_c()

