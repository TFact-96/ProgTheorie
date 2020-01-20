import random
import math
from classes.Amino import Amino
import numpy as np
import itertools
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from matplotlib import style
import copy
import matplotlib.pyplot as plt

class ChainLattice:
    def __init__(self, protein):
        # for stuck aminos
        self.overlap_counter = 0
        self.state_stuck = False

        # protein string
        self.protein = protein

        # create amino object for initial amino at (0,0) with first string char as type
        self.first_amino = Amino(0, 0, 0)
        self.first_amino.type = self.protein[0]

        # whole chain array of amino aminos with initial amino as start
        self.state = {0: self.first_amino}

        # available moves
        self.moves = [[1, 0, 0], [0, 1, 0], [0, 0, 1],
                    [-1, 0, 0], [0, -1, 0], [0, 0, -1]]
                    
        # diagonal moves for chain pulling
        self.diagonal_moves = [[-1, 1, 0], [-1, -1, 0], [-1, 0, 1], [-1, 0, -1],
                                [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1],
                                [1, -1, 0], [1, 1, 0], [1, 0, 1], [1, 0, -1]]

        # folding code corresponding to move index
        self.fold_code_to_index = {"1":0, "2":1, "3":2, "-1":3, "-2":4, "-3":5}
        self.index_to_fold_code = {0:"1", 1:"2", 2:"3", 3:"-1", 4:"-2", 5:"-3"}


        # different types of bonds between aminos (consists of coordinates between aminos)
        self.hh_bonds = []
        self.cc_bonds = []
        self.ch_bonds = []

        # general state stability
        self.stability = 0
        self.old_stability = 0

        # for Mehmet's hillclimb
        self.best_chain = []
        self.best_hh_bonds = []
        self.best_ch_bonds = []
        self.best_cc_bonds = []


    ###################################### Returns next amino object with non-self-overlapping coords. Returns none if stuck.
    def generate_amino_random_move(self):
        # get random move of new amino added to state
        random_move_index = np.random.randint(len(self.moves))

        # get fold code that creates this move
        fold_code = self.index_to_fold_code[random_move_index]

        # last amino needed to generate next one
        last_amino = self.state[len(self.state) - 1]

        # last amino coords + random move
        new_x = last_amino.x + self.moves[random_move_index][0]
        new_y = last_amino.y + self.moves[random_move_index][1]
        new_z = last_amino.z + self.moves[random_move_index][2]
        new_coords = [new_x, new_y, new_z]            

        # if the new amino overlaps its own state; try new move.
        if self.check_amino_overlap(new_coords):
            return self.generate_amino_random_move()

        # no overlap: accepted, create an amino object at this coord
        return self.create_amino_object(new_coords, last_amino, fold_code)

    ###################################### Checking if new coords overlap the already generated state
    def check_amino_overlap(self, new_coords):
        # check if this amino overlaps the state (new amino overlaps any other amino)
        for amino_key, amino in self.state.items():
            if [amino.x, amino.y, amino.z] == [new_coords[0], new_coords[1], new_coords[2]]:
                return True

        return False

    ####################################### Check if chain is stuck
    def chain_stuck(self, x, y, z):
        if (
            (not self.check_amino_overlap([x + 1, y, z]))
            and (not self.check_amino_overlap([x - 1, y, z]))
            and (not self.check_amino_overlap([x, y + 1, z]))
            and (not self.check_amino_overlap([x, y - 1, z]))
            and (not self.check_amino_overlap([x, y, z + 1]))
            and (not self.check_amino_overlap([x, y, z - 1]))
        ):
            return True

        return False
    ####################################### creating a next amino object based on the last amino and fold code to get to new coords
    def create_amino_object(self, new_coords, last_amino, fold_code):
        # make amino object
        new_amino = Amino(new_coords[0], new_coords[1], new_coords[2])

        # set aminonumber
        new_amino.n = last_amino.n + 1

        # Get amino type from new char index of protein-string
        new_amino.type = self.protein[new_amino.n]

        # set the folding move from last amino that creates this one
        last_amino.fold_code = fold_code

        return new_amino

    ###################################### Calculates current HH/CH/CC bonds (and coords), and current stability based on current state
    def get_stability_and_bonds(self, only_stability):
        # create list of all bondable aminos in the state
        bondable_aminos = [
            amino for amino_key, amino in self.state.items() if amino.type == "C" or amino.type == "H"
        ]

        stability = 0
        hh_bonds = []
        ch_bonds = []
        cc_bonds = []

        # compare amino with each other amino in state
        for amino, compare_amino in itertools.combinations(bondable_aminos, 2):
            # if theyre not neighbor aminos
            if (abs(compare_amino.n - amino.n) > 1):
                dist = math.sqrt(
                    (amino.x - compare_amino.x)**2
                    + (amino.y - compare_amino.y)**2
                    + (amino.z - compare_amino.z)**2
                )

                # if distance is 1 nontheless => bond depending on amino types
                if dist <= 1:
                    if amino.type == "H" and compare_amino.type == "H":
                        stability -= 1

                        if not only_stability:
                            hh_bonds.append(
                                [[amino.x, compare_amino.x],
                                [amino.y, compare_amino.y],
                                [amino.z, compare_amino.z]]
                            )

                    if (amino.type == "H" and compare_amino.type == "C") or (amino.type == "C" and compare_amino.type == "H"):
                        stability -= 1

                        if not only_stability:
                            ch_bonds.append(
                                [[amino.x, compare_amino.x],
                                [amino.y, compare_amino.y],
                                [amino.z, compare_amino.z]]
                            )

                    if amino.type == "C" and compare_amino.type == "C":
                        stability -= 5

                        if not only_stability:
                            cc_bonds.append(
                                [[amino.x, compare_amino.x],
                                [amino.y, compare_amino.y],
                                [amino.z, compare_amino.z]]
                            )

        if only_stability == True:
            return stability

        return stability, hh_bonds, ch_bonds, cc_bonds

    ######################################### Sets the current stability and bonds attributes based on the current state
    def set_stability_and_bonds(self):
        stability, hh_bonds, ch_bonds, cc_bonds = self.get_stability_and_bonds(False)
        self.stability = stability
        self.hh_bonds = hh_bonds
        self.ch_bonds = ch_bonds
        self.cc_bonds = cc_bonds
        return


    ############################################## MEHMET'S PULL-NODE ALGORITHM
    # checks if a point has a amino or not
    def check_point(self, array, chain):
        for index, nodes in chain.items():
            if str(array) == str(np.array([nodes.x, nodes.y, nodes.z])):
                return False
        return True
            
    def create_vectors(self, node, chain):
        node_i_coords = np.array([node.x, node.y, node.z])
        node_i1 = chain[int(node.n) + 1]
        node_i1_coords = np.array([node_i1.x, node_i1.y, node_i1.z])
        vector1 = node_i1_coords - node_i_coords

        return node_i_coords, node_i1, node_i1_coords, vector1

    def pull_move(self, node, chain):

        node_i_coords, node_i1, node_i1_coords, vector1 = self.create_vectors(
            node, chain
        )

        checker = {
            "[-1, 1, 0]": [True, [-1, 1, 0]],
            "[-1, -1, 0]": [True, [-1, -1, 0]],
            "[-1, 0, 1]": [True, [-1, 0, 1]],
            "[-1, 0, -1]": [True, [-1, 0, -1]],
            "[0, 1, 1]": [True, [0, 1, 1]],
            "[0, 1, -1]": [True, [0, 1, -1]],
            "[0, -1, 1]": [True, [0, -1, 1]],
            "[0, -1, -1]": [True, [0, -1, -1]],
            "[1, -1, 0]": [True, [1, -1, 0]],
            "[1, 1, 0]": [True, [1, 1, 0]],
            "[1, 0, 1]": [True, [1, 0, 1]],
            "[1, 0, -1]": [True, [1, 0, -1]],
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
                    residue_node.z = residue_node_next.z

                # node moves to L
                node.x = L[0]
                node.y = L[1]
                node.z = L[2]

                # Previous node moves to C
                previous_node = chain[int(node.n) - 1]
                previous_node.x = C[0]
                previous_node.y = C[1]
                previous_node.z = C[2]