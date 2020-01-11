import random
import math
from classes.Atom import Atom
import numpy as np
import itertools

class AminoLattice:
    def __init__(self, amino):
        # for stuck nodes
        self.overlap_counter = 0
        self.chain_stuck = False

        # amino string
        self.amino = amino

        # create atom object for initial node at (0,0) with first string char as type
        self.first_node = Atom(0, 0, 0)
        self.first_node.type = self.amino[0]

        # whole chain array of atom nodes with initial node as start
        self.chain = [self.first_node]

        # available moves
        self.moves = [[1, 0, 0], [0, 1, 0], [0, 0, 1],
                      [-1, 0, 0], [0, -1, 0], [0, 0, -1]]
        # folding code corresponding to move index
        self.move_code = [1, 2, 3, -1, -2, -3]

        # already used moves
        self.leftover_moves = self.moves

        # different types of bonds between nodes (consists of coordinates between nodes)
        self.hh_bonds = []
        self.cc_bonds = []
        self.ch_bonds = []

        # general chain stability
        self.stability = 0

    ###################################### Returns next Atom object with non-self-overlapping coords. Returns none if stuck.
    def generate_random_valid_node(self):
        # if more than 100 random moves are tried (see check_overlap()); almost 100% chance that chain is stuck; Exit program
        if self.overlap_counter > 100:
            self.chain_stuck = True
            return

        # get random move of new node added to chain
        random_move_index = np.random.randint(len(self.leftover_moves))

        # last atom needed to generate next one
        last_atom = self.chain[-1]

        # last atom coords + random move
        new_x = last_atom.x + self.leftover_moves[random_move_index][0]
        new_y = last_atom.y + self.leftover_moves[random_move_index][1]
        new_z = last_atom.z + self.leftover_moves[random_move_index][2]
        new_coords = [new_x, new_y, new_z]

        # get fold code that creates this move
        fold_code = self.move_code[random_move_index]

        # if the new node overlaps its own chain; try new move.
        if self.check_node_overlap(new_coords):
            self.leftover_moves.pop(random_move_index)
            return self.generate_random_valid_node()

        # if doesnt overlap; accepted!
        # reset leftover_moves to all moves (for next new node generation)
        self.leftover_moves = self.moves

        # make Atom object
        new_node = Atom(new_coords[0], new_coords[1], new_coords[2])

        # set atomnumber
        new_node.n = last_atom.n + 1

        # Get atom type from new char index of amino-string
        new_node.type = self.amino[new_node.n]

        # set the folding move from last atom that creates this one
        last_atom.fold_code = fold_code

        return new_node

    ###################################### Checking if new coords overlap the already generated chain
    def check_node_overlap(self, new_coords):
        # check if this node overlaps the chain (new node overlaps any other node)
        for node in self.chain:
            if node.x == new_coords[0] and node.y == new_coords[1] and node.z == new_coords[2]:
                self.overlap_counter += 1
                return True

        # reset overlap_counter if new node is added
        self.overlap_counter = 0
        return False

    ###################################### Calculates current HH/CH/CC bonds (and coords), and current stability based on current chain
    def get_stability_and_bonds(self, only_stability):
        # create list of all bondable atoms in the chain
        bondable_atoms = [
            atom for atom in self.chain if atom.type == "C" or atom.type == "H"
        ]

        atom_nr = 0
        stability = 0
        hh_bonds = []
        ch_bonds = []
        cc_bonds = []

        # compare atom with each other atom in chain
        for atom, compare_atom in itertools.combinations(bondable_atoms, 2):

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

                    if (atom.type == "H" and compare_atom.type == "C") or (
                        atom.type == "C" and compare_atom.type == "H"
                    ):
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

        # delete this index for double count prevention
        atom_nr += 1

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
