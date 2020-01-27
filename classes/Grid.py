import random
from classes.GridPoint import GridPoint
from classes.Node import Node

# Creates a Grid the size of NxN where N is the stringlength of the protein.
# Can add Nodes onto it, clear Nodes, transfer Nodes to new coords,
# clear whole grid, check for overlap of coordinates,
# check if all points around a coord are filled (stuck chain),
# update a node's possible bonds, or update all nodes' bonds.
class Grid:
    def __init__(self, protein):

        # string of all amino nodes
        self.protein = protein

        # available moves from a Node to another Node
        self.moves = [[1, 0, 0], [0, 1, 0], [0, 0, 1],
                    [-1, 0, 0], [0, -1, 0], [0, 0, -1]]

        # Diagonal moves for chain pulling
        self.diagonal_moves = [[-1, 1, 0], [-1, -1, 0], [-1, 0, 1], [-1, 0, -1],
                                [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1],
                                [1, -1, 0], [1, 1, 0], [1, 0, 1], [1, 0, -1]]

        # gridpoints where chain can be in
        self.grid = {}

        # only filled gridpoints for quicker deepcopying
        self.filled_gridpoints = {}

        # chain itself
        self.grid_chain = {}

        # Stores the bond-line coordinates between nodes that are bonded.
        self.hh_bonds = []
        self.ch_bonds = []
        self.cc_bonds = []

        # stability weight of each bond combination
        self.hh_weight = -1
        self.ch_weight = -1
        self.cc_weight = -5

        # stability of the chain (length of the bond-lists combined)
        self.stability = 0

    # Create n x n grid
    def create_grid(self, n):
        for z in range(-n, n + 1):
            for y in range(-n, n + 1):
                for x in range(-n, n + 1):
                    self.grid[f"{x, y, z}"] = GridPoint(False, [x, y, z])

    # Clear all filled points from the grid
    def clear_grid(self):
        for key, value in self.grid.items():
            self.grid[key].filled = False
            self.grid[key].nodes = []

    # reducing deepcopy time by making dict of only the filled gridpoints
    def set_filled_gridpoints_from_grid(self):
        self.filled_gridpoints = {}

        for key, value in self.grid.items():
            if value.filled:
                self.filled_gridpoints[key] = value

    # merge the filled gridpoints onto the grid (deepcopy reduce time)
    def set_grid_from_filled_gridpoints(self):
        # clear grid first
        self.clear_grid()

        # put filled gridpoints into it
        for key, value in self.filled_gridpoints.items():
            self.grid[key] = self.filled_gridpoints[key]

    # checking the neighbors of a node and adding bonds
    def update_bonds_of_node(self, node):
        x = node.x
        y = node.y
        z = node.z
        n = node.n

        if node.type == "P":
            return False

        # neighbor of node = 1 move away
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

    # updating bonds for all nodes and calculating stability with it
    def update_all_bonds(self):
        # reset bonds
        self.hh_bonds = []
        self.ch_bonds = []
        self.cc_bonds = []

        # finding bonds for each node
        for key, value in self.grid_chain.items():
            self.update_bonds_of_node(self.grid[value[0]].nodes[0])

        # stability is the amount of bonds times its individual weight
        self.stability = (len(self.hh_bonds) * self.hh_weight) + (len(self.ch_bonds) * self.ch_weight) + (len(self.cc_bonds) * self.cc_weight)

    # adding a node onto the grid
    def add_point(self, node, n):
        x = node.x
        y = node.y
        z = node.z
        self.grid[f"{x, y, z}"].add_node(node)
        self.grid_chain[n] = [f"{x, y, z}", [x, y, z]]

    # clearing a node from the grid
    def clear_point(self, node, n):
        x = node.x
        y = node.y
        z = node.z
        self.grid[f"{x, y, z}"].remove_node(node)

    # transferring a node to new coordinates
    def transfer_point(self, node1, x2, y2, z2):
        n = node1.n
        self.clear_point(node1, n)
        node1.x = x2
        node1.y = y2
        node1.z = z2
        self.add_point(node1, n)

    # True if given coordinates are already filled
    # --> overlap of nodes
    def overlap(self, x, y, z):
        if self.grid[f"{x, y, z}"].filled:
            return True

        return False

    # True if all coordinates around a point are already filled
    # --> chain is stuck
    def chain_stuck(self, x, y, z):
        overlap_counter = 0

        for neighbor_coords in self.moves:
            if (
                self.overlap(x + neighbor_coords[0], y + neighbor_coords[1], z + neighbor_coords[2])
            ):
                overlap_counter += 1

        # all possible moves are overlapping --> chain stuck
        if overlap_counter == len(self.moves):
            return True

        return False
