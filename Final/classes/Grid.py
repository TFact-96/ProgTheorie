import numpy as np
import random
from algorithms.CalcUpperbound import calc_upperbound
from classes.GridPoint import GridPoint
from classes.Node import Node

class Grid:
    def __init__(self, amino):

        self.amino = amino
        
        # available moves
        self.moves = [[1, 0, 0], [0, 1, 0], [0, 0, 1],
                    [-1, 0, 0], [0, -1, 0], [0, 0, -1]]
                    
        # diagonal moves for chain pulling
        self.diagonal_moves = [[-1, 1, 0], [-1, -1, 0], [-1, 0, 1], [-1, 0, -1],
                                [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1],
                                [1, -1, 0], [1, 1, 0], [1, 0, 1], [1, 0, -1]]
                                        
        # gridpoints where chain can be in
        self.grid = {}
        
        # chain itself
        self.grid_chain = {}
        
        # bonds
        self.hh_bonds = []
        self.ch_bonds = []
        self.cc_bonds = []
        
        # stability of the chain
        self.stability = 0
                
    # Create n x n grid
    def create_grid(self, n):
        for z in range(-n, n + 1):
            for y in range(-n, n + 1):
                for x in range(-n, n + 1):
                    self.grid[f"{x, y, z}"] = GridPoint(False, [x, y, z])

    def clear_grid(self):
        for key, value in self.grid.items():
            self.grid[key].filled = False
            self.grid[key].nodes = []

    def add_neighbours(self, node):
        x = node.x
        y = node.y
        z = node.z
        n = node.n

        if node.type == "P":
            return False
        
        for move in self.moves:
            # coords of neighbor
            x2, y2, z2 = x + move[0], y + move[1], z + move[2]
            
            # neighbor node in grid
            neighbor_grid = self.grid[f"{x2, y2, z2}"]
            
            # if neighbor exists and not next to eachother in chain
            if neighbor_grid.filled and (abs(neighbor_grid.nodes[0].n - n) > 1):
                
                # get the node from the grid
                neighbor_node = neighbor_grid.nodes[0]
                
                # create bond lines between the nodes
                bond_line = [[x, x2], [y, y2], [z, z2]]
                inverse_bond_line = [[x2, x], [y2, y], [z2, z]]
                
                # H-H bonds
                if node.type == "H" and (neighbor_node.type == "H"):
                    # prevent double bond counting
                    if bond_line and inverse_bond_line not in self.hh_bonds:
                        self.hh_bonds.append(bond_line)
                
                # H-C bonds
                if (
                    (node.type == "H" and neighbor_node.type == "C")
                    or (node.type == "C" and neighbor_node.type == "H")
                ):
                    # prevent double bond counting
                    if bond_line and inverse_bond_line not in self.ch_bonds:
                        self.ch_bonds.append(bond_line)
                        
                # C-C bonds
                if node.type == "C" and (neighbor_node.type == "C"):
                    # prevent double bond counting
                    if bond_line and inverse_bond_line not in self.cc_bonds:
                        self.cc_bonds.append(bond_line)

    def update_neighbours(self):
        # reset bonds
        self.hh_bonds = []
        self.ch_bonds = []
        self.cc_bonds = []

        for key, value in self.grid_chain.items():
            self.add_neighbours(self.grid[value[0]].nodes[0])
        
        self.stability = - (len(self.hh_bonds) + len(self.ch_bonds) + (5 * len(self.cc_bonds)))
        
    def add_point(self, node, n):
        x = node.x
        y = node.y
        z = node.z
        self.grid[f"{x, y, z}"].add_node(node)
        self.grid_chain[n] = [f"{x, y, z}", [x, y, z]]

    def clear_point(self, node, n):
        x = node.x
        y = node.y
        z = node.z
        self.grid[f"{x, y, z}"].remove_node(node)

    def transfer_point(self, node1, x2, y2, z2):
        n = node1.n
        self.clear_point(node1, n)
        node1.x = x2
        node1.y = y2
        node1.z = z2
        self.add_point(node1, n)

    def overlap(self, x, y, z):
        if self.grid[f"{x, y, z}"].filled:
            return True

        return False

    def chain_stuck(self, x, y, z):
        overlap_counter = 0
        
        for neighbor_coords in self.moves:
            if (
                self.overlap(x + neighbor_coords[0], y + neighbor_coords[1], z + neighbor_coords[2])
            ):
                overlap_counter += 1

        if overlap_counter == len(self.moves):
            return True
            
        return False

    def check_diagonals(self, x, y, z):
        # diagonal moves for chain pulling
        available_moves = []
        
        for move in self.diagonal_moves:
            if not self.grid[f"{x + move[0], y + move[1], z + move[2]}"].filled:
                available_moves.append(move)

        return available_moves

    def create_chain(self):
        index = 0

        # Create grid
        # too big?
        self.create_grid(int(len(self.amino) / 2))
        first_node = Node(0, 0, 0)
        first_node.type = self.amino[0]
        self.add_point(first_node, index)

        while index < len(self.amino) - 1:
            current_node_key = self.grid_chain[index][0]
            current_node = self.grid[current_node_key].nodes[0]

            random_move = np.random.randint(len(self.moves))

            new_x = current_node.x + self.moves[random_move][0]
            new_y = current_node.y + self.moves[random_move][1]
            new_z = current_node.z + self.moves[random_move][2]

            if not self.overlap(new_x, new_y, new_z):
                index += 1
                new_node = Node(new_x, new_y, new_z)
                new_node.n = index
                new_node.type = self.amino[index]
                self.add_point(new_node, index)

            if self.chain_stuck(new_x, new_y, new_z):
                index = 0
                self.clear_grid()
                self.grid_chain = {}
                self.add_point(first_node, index)

    def create_vectors(self, node):
        node_i_coords = [node.x, node.y, node.z]
        node_i1 = self.grid_chain[int(node.n) + 1]
        node_i1_coords = np.array([node_i1[1][0], node_i1[1][1], node_i1[1][2]])
        vector1 = node_i1_coords - np.array(node_i_coords)

        return node_i_coords, node_i1_coords, vector1

    def check_requirements(self, available_moves, vector1, node_i_coords, node_i1_coords):
        viable_moves = []
        found = False
        
        for move in available_moves:
            L = node_i_coords + np.array(move)
            C = L - vector1

            if (
                (not self.overlap(L[0], L[1], L[2]))
                and (not self.overlap(C[0], C[1], C[2]))
                and (np.linalg.norm(L - node_i1_coords) == 1.0)
            ):
                viable_moves.append([L, C])
                found = True

        if not found:
            return 0, 0, False

        random_choice = random.choice(viable_moves)
        return random_choice[0], random_choice[1], True

    def pull_move(self, node):

        node_i_coords, node_i1_coords, vector1 = self.create_vectors(node)

        available_moves = self.check_diagonals(node_i_coords[0], node_i_coords[1], node_i_coords[2])

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
                    residue_node,
                    residue_node_next.x,
                    residue_node_next.y,
                    residue_node_next.z,
                )

            # node moves to L
            self.transfer_point(node, L[0], L[1], L[2])

            # Previous node moves to C
            previous_node_key = self.grid_chain[int(node.n) - 1][0]
            previous_node = self.grid[previous_node_key].nodes[0]

            self.transfer_point(previous_node, C[0], C[1], C[2])