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


class Grid:
    def __init__(self):
        self.grid = {}
        self.grid_chain = []

    # Create n x n grid
    def create_grid(self, n):
        for y in range(-n, n + 1):
            for x in range(-n, n + 1):
                self.grid[f"{x, y}"] = [False, [x, y], 0]

    def clear_grid(self):
        for key, value in self.grid.items():
            if key == f"{0, 0}":
                continue
            self.grid[key][0] = False

    def add_point(self, x, y):

        if self.grid[f"{x, y}"][0]:
            self.grid[f"{x, y}"][2] += 1

        else:
            self.grid[f"{x, y}"][0] = True
            self.grid_chain.append(f"{x, y}")

    def clear_point(self, x, y):

        if self.grid[f"{x, y}"][2] > 0:
            self.grid[f"{x, y}"][2] -= 1

        else:
            self.grid[f"{x, y}"][0] = False
            self.grid_chain.remove(f"{x, y}")

    def overlap(self, x, y):
        if self.grid[f"{x, y}"][0]:
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

        if not self.grid[f"{x + 1, y + 1}"][0]:
            available_moves.append([1, 1])

        if not self.grid[f"{x + 1, y - 1}"][0]:
            available_moves.append([1, -1])

        if not self.grid[f"{x - 1, y + 1}"][0]:
            available_moves.append([-1, 1])

        if not self.grid[f"{x - 1, y - 1}"][0]:
            available_moves.append([-1, -1])

        return available_moves

    def __repr__(self):
        return f"{self.grid}"


class NodeChain:
    def __init__(self, amino):

        self.amino = amino
        self.best_chain = []
        self.diagonal_moves = [[1, 1], [1, -1], [-1, 1], [-1, -1]]

    def create_chain(self):

        index = 0
        self.moves = [[1, 0], [0, 1], [-1, 0], [0, -1]]

        # Create grid
        grid = Grid()
        grid.create_grid(36)

        first_node = Node(0, 0)
        first_node.type = self.amino[0]
        chain = {0: first_node}
        grid.add_point(0, 0)

        # Create a random chain
        while len(chain) < len(self.amino):
            current_node = chain[index]
            random_move = np.random.randint(len(self.moves))

            new_x = current_node.x + self.moves[random_move][0]
            new_y = current_node.y + self.moves[random_move][1]

            if not grid.overlap(new_x, new_y):
                index += 1
                chain[index] = Node(new_x, new_y)
                chain[index].n = index
                chain[index].type = self.amino[index]
                grid.add_point(new_x, new_y)

            if grid.chain_stuck(new_x, new_y):
                chain = {0: first_node}
                index = 0
                grid.clear_grid()

        return chain, grid

    def create_vectors(self, node, chain):
        node_i_coords = [node.x, node.y]
        node_i1 = chain[int(node.n) + 1]
        node_i1_coords = np.array([node_i1.x, node_i1.y])
        vector1 = node_i1_coords - np.array(node_i_coords)

        return node_i_coords, node_i1, node_i1_coords, vector1

    def check_requirements(
        self, available_moves, vector1, node_i_coords, node_i1_coords, grid
    ):

        for move in available_moves:
            L = node_i_coords + np.array(move)
            C = L - vector1

            if (
                (not grid.overlap(L[0], L[1]))
                and (not grid.overlap(C[0], C[1]))
                and (np.linalg.norm(L - node_i1_coords) == 1.0)
            ):
                return L, C, True

        return 0, 0, False

    def pull_move(self, node, chain, grid):

        self.grid = grid

        node_i_coords, node_i1, node_i1_coords, vector1 = self.create_vectors(
            node, chain
        )

        available_moves = self.grid.check_diagonals(node_i_coords[0], node_i_coords[1])

        L, C, check = self.check_requirements(
            available_moves, vector1, node_i_coords, node_i1_coords, grid
        )

        if check:

            # residue of chain follows in footsteps
            for index in range(int(node.n - 1)):

                residue_node = chain[index]
                residue_node_next = chain[index + 2]

                self.grid.clear_point(residue_node.x, residue_node.y)
                self.grid.add_point(residue_node_next.x, residue_node_next.y)

                residue_node.x = residue_node_next.x
                residue_node.y = residue_node_next.y

            # node moves to L
            self.grid.clear_point(node.x, node.y)
            self.grid.add_point(L[0], L[1])

            node.x = L[0]
            node.y = L[1]

            # Previous node moves to C
            previous_node = chain[int(node.n) - 1]

            self.grid.clear_point(previous_node.x, previous_node.y)
            self.grid.add_point(C[0], C[1])

            previous_node.x = C[0]
            previous_node.y = C[1]

        return chain, grid

    def update_neighbours(self, chain):
        temp_stability = 0
        hh_bonds = []

        i = 0
        for a, b in itertools.combinations(list(chain.items()), 2):
            i += 1
            print(i)
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
                temp_stability -= 1

        stability = temp_stability
        return stability, hh_bonds

    def reset_neighbours(self, chain):

        self.hh_bonds = []
        for node_key, node in chain.items():
            node.neighbours = []

    def hill_climber(self, max_iteration):
        # Current protein chain
        self.current_hilltop, grid = self.create_chain()
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

            for it in range(10):
                for index in range(1, len(self.current_hilltop) - 1):
                    self.pull_move(
                        self.current_hilltop[index], self.current_hilltop, grid
                    )
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
                self.current_hilltop, grid = self.create_chain()
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
k.hill_climber(100)
k.find_best_c()

