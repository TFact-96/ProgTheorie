import random
import math
from classes.Atom import Atom3D, Atom2D
import numpy as np
import itertools
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from matplotlib import style
import copy
import matplotlib.pyplot as plt

class AminoLattice:
    def __init__(self, amino, ThreeD):
        # 2D or 3D
        self.ThreeD = ThreeD
        
        # for stuck nodes
        self.overlap_counter = 0
        self.chain_stuck = False

        # amino string
        self.amino = amino

        # create atom object for initial node at (0,0) with first string char as type
        if self.ThreeD:
            self.first_node = Atom3D(0, 0, 0)
        else:
            self.first_node = Atom2D(0, 0)
            
        self.first_node.type = self.amino[0]

        # whole chain array of atom nodes with initial node as start        
        self.chain = {0: self.first_node}
        

        # available moves
        if self.ThreeD:
            self.moves = [[1, 0, 0], [0, 1, 0], [0, 0, 1],
                        [-1, 0, 0], [0, -1, 0], [0, 0, -1]]
        else:
            self.moves = [[1, 0], [0, 1], [-1, 0], [0, -1]]
            self.diagonal_moves = [[1, 1], [1, -1], [-1, 1], [-1, -1]]

        if self.ThreeD:
            # folding code corresponding to move index
            self.fold_code_to_index = {"1":0, "2":1, "3":2, "-1":3, "-2":4, "-3":5}
            self.index_to_fold_code = {0:"1", 1:"2", 2:"3", 3:"-1", 4:"-2", 5:"-3"}
        else:
            self.index_to_fold_code = {0:"1", 1:"-1", 2:"2", 3:"-2"}
            self.fold_code_to_index = {"1":0, "-1":1, "2":3, "-2":4,}


        # different types of bonds between nodes (consists of coordinates between nodes)
        self.hh_bonds = []
        self.cc_bonds = []
        self.ch_bonds = []

        # general chain stability
        self.stability = 0
        self.old_stability = 0

    ###################################### Returns next Atom object with non-self-overlapping coords. Returns none if stuck.
    def generate_atom_random_move(self):
        # if more than 100 random moves are tried (see check_overlap()); almost 100% chance that chain is stuck; Exit program
        if self.overlap_counter > 100:
            self.chain_stuck = True
            return

        # get random move of new node added to chain
        random_move_index = np.random.randint(len(self.moves))

        # get fold code that creates this move
        fold_code = self.index_to_fold_code[random_move_index]

        # last atom needed to generate next one
        last_atom = self.chain[len(self.chain) - 1]

        # last atom coords + random move
        if self.ThreeD:
            new_x = last_atom.x + self.moves[random_move_index][0]
            new_y = last_atom.y + self.moves[random_move_index][1]
            new_z = last_atom.z + self.moves[random_move_index][2]
            new_coords = [new_x, new_y, new_z]
        else:
            new_x = last_atom.x + self.moves[random_move_index][0]
            new_y = last_atom.y + self.moves[random_move_index][1]
            new_coords = [new_x, new_y]

        # if the new node overlaps its own chain; try new move.
        if self.check_node_overlap(new_coords):
            return self.generate_atom_random_move()

        # no overlap: accepted, create an atom object at this coord
        return self.create_atom_object(new_coords, last_atom, fold_code)

    ###################################### Checking if new coords overlap the already generated chain
    def check_node_overlap(self, new_coords):
        # check if this node overlaps the chain (new node overlaps any other node)
        for node_key, node in self.chain.items():
            if self.ThreeD:
                if [node.x, node.y, node.z] == [new_coords[0], new_coords[1], new_coords[2]]:
                    self.overlap_counter += 1
                    return True
            else:
                if [node.x, node.y] == [new_coords[0], new_coords[1]]:
                    self.overlap_counter += 1
                    return True
                    
        # reset overlap_counter if new node is added
        self.overlap_counter = 0
        return False

    ####################################### creating a next atom object based on the last atom and fold code to get to new coords
    def create_atom_object(self, new_coords, last_atom, fold_code):
        # make Atom object
        if self.ThreeD:
            new_node = Atom3D(new_coords[0], new_coords[1], new_coords[2])
        else:
            new_node = Atom2D(new_coords[0], new_coords[1])

        # set atomnumber
        new_node.n = last_atom.n + 1

        # Get atom type from new char index of amino-string
        new_node.type = self.amino[new_node.n]

        # set the folding move from last atom that creates this one
        last_atom.fold_code = fold_code

        return new_node

    ###################################### Calculates current HH/CH/CC bonds (and coords), and current stability based on current chain
    def get_stability_and_bonds(self, only_stability):
        # create list of all bondable atoms in the chain
        bondable_atoms = [
            atom for atom_key, atom in self.chain.items() if atom.type == "C" or atom.type == "H"
        ]

        stability = 0
        hh_bonds = []
        ch_bonds = []
        cc_bonds = []

        # compare atom with each other atom in chain
        for atom, compare_atom in itertools.combinations(bondable_atoms, 2):
            if self.ThreeD:
                # if theyre not neighbor nodes
                if (abs(compare_atom.n - atom.n) > 1):
                    dist = math.sqrt(
                        (atom.x - compare_atom.x)**2
                        + (atom.y - compare_atom.y)**2
                        + (atom.z - compare_atom.z)**2
                    )

                    # if distance is 1 nontheless => bond depending on atom types
                    if dist <= 1:
                        if atom.type == "H" and compare_atom.type == "H":
                            stability -= 1

                            if not only_stability:
                                hh_bonds.append(
                                    [[atom.x, compare_atom.x],
                                    [atom.y, compare_atom.y],
                                    [atom.z, compare_atom.z]]
                                )

                        if (atom.type == "H" and compare_atom.type == "C") or (atom.type == "C" and compare_atom.type == "H"):
                            stability -= 1

                            if not only_stability:
                                ch_bonds.append(
                                    [[atom.x, compare_atom.x],
                                    [atom.y, compare_atom.y],
                                    [atom.z, compare_atom.z]]
                                )

                        if atom.type == "C" and compare_atom.type == "C":
                            stability -= 5

                            if not only_stability:
                                cc_bonds.append(
                                    [[atom.x, compare_atom.x],
                                    [atom.y, compare_atom.y],
                                    [atom.z, compare_atom.z]]
                                )
            else:
                # Otherwise 2D
                if (abs(compare_atom.n - atom.n) > 1):
                    dist = math.sqrt(
                        (atom.x - compare_atom.x)**2
                        + (atom.y - compare_atom.y)**2
                    )

                    # if distance is 1 nontheless => bond depending on atom types
                    if dist <= 1:
                        if atom.type == "H" and compare_atom.type == "H":
                            stability -= 1

                            if not only_stability:
                                hh_bonds.append(
                                    [[atom.x, compare_atom.x],
                                    [atom.y, compare_atom.y]])

                        if (atom.type == "H" and compare_atom.type == "C") or (atom.type == "C" and compare_atom.type == "H"):
                            stability -= 1

                            if not only_stability:
                                ch_bonds.append(
                                    [[atom.x, compare_atom.x],
                                    [atom.y, compare_atom.y]]
                                    )

                        if atom.type == "C" and compare_atom.type == "C":
                            stability -= 5

                            if not only_stability:
                                cc_bonds.append(
                                    [[atom.x, compare_atom.x],
                                    [atom.y, compare_atom.y]]
                                )

        if only_stability == True:
            return stability

        return stability, hh_bonds, ch_bonds, cc_bonds

    ######################################### Sets the current stability and bonds attributes based on the current chain
    def set_stability_and_bonds(self):
        stability, hh_bonds, ch_bonds, cc_bonds = self.get_stability_and_bonds(False)
        self.stability = stability
        self.hh_bonds = hh_bonds
        self.ch_bonds = ch_bonds
        self.cc_bonds = cc_bonds

        return
        
        
    ############################################## MEHMET'S HILL CLIMB, only 2D
    def link_neighbours(self):
        temp_stability = 0
        
        for node_key, node in self.chain.items():

            for node_neighbour in node.neighbours:

                if node.type == "C" and node_neighbour.type == "C":
                    temp_stability -= 5
                    self.cc_bonds.append(
                        [[node.x, node_neighbour.x], [node.y, node_neighbour.y]]
                    )
                if (node.type == "C" and node_neighbour.type == "H") or (node.type == "H" and node_neighbour.type == "C"):
                    temp_stability -= 1
                    self.ch_bonds.append(
                        [[node.x, node_neighbour.x], [node.y, node_neighbour.y]]
                    )
                if node.type == "H" and node_neighbour.type == "H":
                    temp_stability -= 1
                    self.hh_bonds.append(
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
                == 1.0 ) and (abs(node_1_index - node_2_index) != 1):

                node_1.neighbours.append(node_2)

        self.link_neighbours()

    def reset_neighbours(self):
        for node_key, node in self.chain.items():
            node.neighbours = []

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

    def random_pull(self, pull_times_per_chain):                
        for x in range(pull_times_per_chain):
            d = random.randint(1, len(self.chain) - 2)
            self.pull_move(self.chain[d])
        
