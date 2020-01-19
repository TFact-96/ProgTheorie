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

    def clear_point(self, x, y):

        if self.grid[f"{x, y}"][2] > 0:
            self.grid[f"{x, y}"][2] -= 1

        else:
            self.grid[f"{x, y}"][0] = False

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
        grid.create_grid(15)

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

        return L, C, False

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

    def test_move(self):

        beforepullx = []
        beforepully = []
        afterpullx = []
        afterpully = []

        beforepullxg = []
        beforepullyg = []
        afterpullxg = []
        afterpullyg = []

        self.chain, grid = self.create_chain()

        for a, b in self.chain.items():
            beforepullx.append(b.x)
            beforepully.append(b.y)

        for a, b in grid.grid.items():

            if b[0]:
                beforepullxg.append(b[1][0])
                beforepullyg.append(b[1][1])

        for index in range(1, 14):
            self.chain, grid = self.pull_move(self.chain[index], self.chain, grid)

        for a, b in self.chain.items():
            afterpullx.append(b.x)
            afterpully.append(b.y)

        for a, b in grid.grid.items():
            if b[0]:
                afterpullxg.append(b[1][0])
                afterpullyg.append(b[1][1])

        print(grid)
        plt.subplot(2, 2, 1)
        plt.plot(beforepullx, beforepully, "ro-")

        plt.subplot(2, 2, 2)
        plt.plot(afterpullx, afterpully, "bo-")

        plt.subplot(2, 2, 3)
        plt.plot(beforepullxg, beforepullyg, "ro")

        plt.subplot(2, 2, 4)
        plt.plot(afterpullxg, afterpullyg, "bo")

        plt.show()


k = NodeChain("HPHPPHPHPHPHPHP")
k.test_move()

