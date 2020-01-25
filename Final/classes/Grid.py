import numpy as np
import random
from algorithms.CalcUpperbound import calc_upperbound
from classes.GridPoint import GridPoint
from classes.Node import Node

class Grid:
    def __init__(self, amino):

        self.amino = amino
        self.best_chain = {}
        
        # available moves
        self.moves = [[1, 0, 0], [0, 1, 0], [0, 0, 1],
                    [-1, 0, 0], [0, -1, 0], [0, 0, -1]]
                    
        # diagonal moves for chain pulling
        self.diagonal_moves = [[-1, 1, 0], [-1, -1, 0], [-1, 0, 1], [-1, 0, -1],
                                [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1],
                                [1, -1, 0], [1, 1, 0], [1, 0, 1], [1, 0, -1]]

        # estimate of lowest possible stability
        self.pivot_upperbound = int(calc_upperbound(amino))
        
    # Create n x n grid
    def create_grid(self, n):
        grid = {}
        for z in range(-n, n + 1):
            for y in range(-n, n + 1):
                for x in range(-n, n + 1):
                    grid[f"{x, y, z}"] = GridPoint(False, [x, y, z])

        return grid

    def clear_grid(self, grid):
        for key, value in grid.items():
            grid[key].filled = False
            grid[key].nodes = []

    def add_neighbours(self, node, grid, hh_bonds):
        x = node.x
        y = node.y
        z = node.z
        n = node.n

        if node.type == "P":
            return False

        if (
            grid[f"{x + 1, y, z}"].filled
            and (abs(grid[f"{x + 1, y, z}"].nodes[0].n - n) > 1)
            and (grid[f"{x + 1, y, z}"].nodes[0].type == "H")
        ):
            if [[x, x + 1], [y, y], [z, z]] and [[x + 1, x], [y, y], [z, z]] not in hh_bonds:
                hh_bonds.append([[x, x + 1], [y, y], [z, z]])

        if (
            grid[f"{x, y + 1, z}"].filled 
            and (abs(grid[f"{x, y + 1, z}"].nodes[0].n - n) > 1) 
            and (grid[f"{x, y + 1, z}"].nodes[0].type == "H")
        ):
            if [[x, x], [y, y + 1], [z, z]] and [[x, x], [y + 1, y], [z, z]] not in hh_bonds:
                hh_bonds.append([[x, x], [y, y + 1], [z, z]])

        if (
            grid[f"{x, y, z + 1}"].filled 
            and (abs(grid[f"{x, y, z + 1}"].nodes[0].n - n) > 1) 
            and (grid[f"{x, y, z + 1}"].nodes[0].type == "H")
        ):
            if [[x, x], [y, y], [z + 1, z]] and [[x, x], [y, y], [z, z + 1]] not in hh_bonds:
                hh_bonds.append([[x, x], [y, y], [z, z + 1]])

        if (
            grid[f"{x - 1, y, z}"].filled
            and (abs(grid[f"{x - 1, y, z}"].nodes[0].n - n) > 1)
            and (grid[f"{x - 1, y, z}"].nodes[0].type == "H")
        ):
            if [[x, x - 1], [y, y], [z, z]] and [[x - 1, x], [y, y], [z, z]] not in hh_bonds:
                hh_bonds.append([[x, x - 1], [y, y], [z, z]])

        if (
            grid[f"{x, y - 1, z}"].filled
            and (abs(grid[f"{x, y - 1, z}"].nodes[0].n - n) > 1)
            and (grid[f"{x, y - 1, z}"].nodes[0].type == "H")
        ):
            if [[x, x], [y, y - 1], [z, z]] and [[x, x], [y - 1, y], [z, z]] not in hh_bonds:
                hh_bonds.append([[x, x], [y, y - 1], [z, z]])

        if (
            grid[f"{x, y, z - 1}"].filled
            and (abs(grid[f"{x, y, z - 1}"].nodes[0].n - n) > 1) 
            and (grid[f"{x, y, z - 1}"].nodes[0].type == "H")
        ):
            if [[x, x], [y, y], [z - 1, z]] and [[x, x], [y, y], [z, z - 1]] not in hh_bonds:
                hh_bonds.append([[x, x], [y, y], [z, z - 1]])

    def update_neighbours(self, grid, grid_chain):
        hh_bonds = []
        for key, value in grid_chain.items():
            self.add_neighbours(grid[value[0]].nodes[0], grid, hh_bonds)

        return -len(hh_bonds), hh_bonds

    def add_point(self, node, n, grid, grid_chain):
        x = node.x
        y = node.y
        z = node.z
        grid[f"{x, y, z}"].add_node(node)
        grid_chain[n] = [f"{x, y, z}", [x, y, z]]

    def clear_point(self, node, n, grid):
        x = node.x
        y = node.y
        z = node.z
        grid[f"{x, y, z}"].remove_node(node)

    def transfer_point(self, node1, x2, y2, z2, grid, grid_chain):
        n = node1.n
        self.clear_point(node1, n, grid)
        node1.x = x2
        node1.y = y2
        node1.z = z2
        self.add_point(node1, n, grid, grid_chain)

    def overlap(self, x, y, z, grid):
        if grid[f"{x, y, z}"].filled:
            return True

        return False

    def chain_stuck(self, x, y, z, grid):

        if (
            self.overlap(x + 1, y, z, grid)
            and self.overlap(x - 1, y, z, grid)
            and self.overlap(x, y + 1, z, grid)
            and self.overlap(x, y - 1, z, grid)
            and self.overlap(x, y, z + 1, grid)
            and self.overlap(x, y, z - 1, grid)
        ):
            return True

        return False

    def check_diagonals(self, x, y, z, grid):
        # diagonal moves for chain pulling
        #self.diagonal_moves = [[-1, 1, 0], [-1, -1, 0], [-1, 0, 1], [-1, 0, -1],
        #                        [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1],
        #                        [1, -1, 0], [1, 1, 0], [1, 0, 1], [1, 0, -1]]

        available_moves = []

        if not grid[f"{x + 1, y + 1, z}"].filled:
            available_moves.append([1, 1, 0])

        if not grid[f"{x + 1, y - 1, z}"].filled:
            available_moves.append([1, -1, 0])

        if not grid[f"{x - 1, y + 1 , z}"].filled:
            available_moves.append([-1, 1, 0])

        if not grid[f"{x - 1, y - 1, z}"].filled:
            available_moves.append([-1, -1, 0])

        if not grid[f"{x, y + 1, z + 1}"].filled:
            available_moves.append([0, 1, 1])

        if not grid[f"{x, y + 1, z - 1}"].filled:
            available_moves.append([0, 1, -1])

        if not grid[f"{x, y - 1, z + 1}"].filled:
            available_moves.append([0, -1, 1])

        if not grid[f"{x, y - 1, z - 1}"].filled:
            available_moves.append([0, -1, -1])

        if not grid[f"{x + 1, y, z + 1}"].filled:
            available_moves.append([1, 0, 1])

        if not grid[f"{x + 1, y, z - 1}"].filled:
            available_moves.append([1, 0, -1])

        if not grid[f"{x - 1, y, z + 1}"].filled:
            available_moves.append([-1, 0, 1])

        if not grid[f"{x - 1, y, z - 1}"].filled:
            available_moves.append([-1, 0, -1])

        return available_moves

    def create_chain(self):
        index = 0
        grid_chain = {}

        # Create grid
        # too big?
        grid = self.create_grid(int(len(self.amino) / 2))
        first_node = Node(0, 0, 0)
        first_node.type = self.amino[0]
        self.add_point(first_node, index, grid, grid_chain)

        while index < len(self.amino) - 1:
            current_node_key = grid_chain[index][0]
            current_node = grid[current_node_key].nodes[0]

            random_move = np.random.randint(len(self.moves))

            new_x = current_node.x + self.moves[random_move][0]
            new_y = current_node.y + self.moves[random_move][1]
            new_z = current_node.z + self.moves[random_move][2]

            if not self.overlap(new_x, new_y, new_z, grid):
                index += 1
                new_node = Node(new_x, new_y, new_z)
                new_node.n = index
                new_node.type = self.amino[index]
                self.add_point(new_node, index, grid, grid_chain)

            if self.chain_stuck(new_x, new_y, new_z, grid):
                index = 0
                self.clear_grid(grid)
                grid_chain = {}
                self.add_point(first_node, index, grid, grid_chain)

        return grid_chain, grid

    def create_vectors(self, node, grid_chain):
        node_i_coords = [node.x, node.y, node.z]
        node_i1 = grid_chain[int(node.n) + 1]
        node_i1_coords = np.array([node_i1[1][0], node_i1[1][1], node_i1[1][2]])
        vector1 = node_i1_coords - np.array(node_i_coords)

        return node_i_coords, node_i1_coords, vector1

    def check_requirements(self, available_moves, vector1, node_i_coords, node_i1_coords, grid):
        viable_moves = []
        found = False
        
        for move in available_moves:
            L = node_i_coords + np.array(move)
            C = L - vector1

            if (
                (not self.overlap(L[0], L[1], L[2], grid))
                and (not self.overlap(C[0], C[1], C[2], grid))
                and (np.linalg.norm(L - node_i1_coords) == 1.0)
            ):
                viable_moves.append([L, C])
                found = True

        if not found:
            return 0, 0, False

        random_choice = random.choice(viable_moves)
        return random_choice[0], random_choice[1], True

    def pull_move(self, node, grid, grid_chain):

        node_i_coords, node_i1_coords, vector1 = self.create_vectors(node, grid_chain)

        available_moves = self.check_diagonals(node_i_coords[0], node_i_coords[1], node_i_coords[2], grid)

        L, C, check = self.check_requirements(
            available_moves, vector1, node_i_coords, node_i1_coords, grid
        )

        if check:

            # residue of chain follows in footsteps
            for index in range(int(node.n - 1)):

                residue_node_key = grid_chain[index][0]
                residue_node = grid[residue_node_key].nodes[0]

                residue_node_next_key = grid_chain[index + 2][0]
                residue_node_next = grid[residue_node_next_key].nodes[0]

                self.transfer_point(
                    residue_node,
                    residue_node_next.x,
                    residue_node_next.y,
                    residue_node_next.z,
                    grid,
                    grid_chain,
                )

            # node moves to L
            self.transfer_point(node, L[0], L[1], L[2], grid, grid_chain)

            # Previous node moves to C
            previous_node_key = grid_chain[int(node.n) - 1][0]
            previous_node = grid[previous_node_key].nodes[0]

            self.transfer_point(previous_node, C[0], C[1], C[2], grid, grid_chain)

        return grid, grid_chain