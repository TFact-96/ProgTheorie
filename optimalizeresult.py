import math
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
import statistics as stat

amino = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
optimalization_tries = 100

class Atom:
    # object for an atom in the amino chain
    def __init__(self, x, y):
        # nth atom
        self.n = 0

        # coords
        self.x = x
        self.y = y

        # atom type (C, H, P)
        self.type = type

        # the fold that this atom makes to the next atom (default 0)
        self.fold_code = 0

        # default color (P)
        self.color = 'blue'

    # print the type if printing the object
    def __str__(self):
        return f"{self.type}"

class AminoLattice:
    def __init__(self, amino):
        self.overlap_counter = 0
        self.optimalization_tries = optimalization_tries
        self.amino = amino

        # create atom object for initial node at (0,0) with first string char as type
        self.first_node = Atom(0, 0)
        self.first_node.type = self.amino[0]

        # put this in a chain array
        self.chain = [self.first_node]

        # available moves
        self.moves = [[1, 0], [0, 1], [-1, 0], [0, -1]]
        # folding code corresponding to move index
        self.move_code = [1, 2, -1, -2]

        # different types of bonds between nodes
        self.hh_bonds = []
        self.cc_bonds = []
        self.ch_bonds = []

        # general chain stability
        self.stability = 0

    def generate_nodes(self):
        # while amount of nodes is smaller than the whole amino string length
        while len(self.chain) < len(self.amino):
            # generate a random move
            #new_node = self.get_next_node()

            # generate a optimalized move (new node always makes bond-making moves if possible)
            new_node = self.get_optimal_node()

            # append the new node to the existing chain
            self.chain.append(new_node)

        # prepare for plotting
        x = []
        y = []
        color = []

        # append all nodes in coord list
        for node in self.chain:
            x.append(node.x)
            y.append(node.y)

            # set node colors
            if node.type == "H":
                node.color = 'red'
            elif node.type == "C":
                node.color = 'orange'

            color.append(node.color)

            # Print the type and fold code in terminal while we're at it.
            print(node.type, node.fold_code)

        # calculate the bonds in this chain
        self.stability, self.hh_bonds, self.ch_bonds, self.cc_bonds = self.calculate_bonds(False)

        # plot the bond lines
        for hh_bond in self.hh_bonds:
            plt.plot(hh_bond[0], hh_bond[1], "--", markersize=1, color='yellow', zorder=1)

        for ch_bond in self.ch_bonds:
            plt.plot(ch_bond[0], ch_bond[1], "--", markersize=1, color='green', zorder=1)

        for cc_bond in self.cc_bonds:
            plt.plot(cc_bond[0], cc_bond[1], "--", markersize=1, color='pink', zorder=1)

        # plot the chain itself
        plt.plot(x, y, "-", linewidth=3, color='black', zorder=1, label=f"Amino chain (Stability: {self.stability})")
        plt.scatter(x, y, color=color, zorder=2)

        # rest info and show plot
        plt.title("Amino acid chain")
        plt.ylabel('y-as')
        plt.xlabel('x-as')
        plt.legend()
        plt.show()

    def get_next_node(self):
        # get random move of new node added to chain
        random_move_index = np.random.randint(len(self.moves))
        # last atom needed to generate next one
        last_atom = self.chain[-1]

        # last atom coords + random move
        new_x = last_atom.x + self.moves[random_move_index][0]
        new_y = last_atom.y + self.moves[random_move_index][1]
        new_coords = [new_x, new_y]

        # get fold code that creates this move
        fold_code = self.move_code[random_move_index]

        # if more than 100 random moves are tried (see check_overlap()); almost 100% chance that chain is stuck; Exit program
        if self.overlap_counter > 100:
            print("Chain got stuck!")
            exit(0)

        # if the new node overlaps its own chain; try new move.
        if self.check_overlap(new_coords):
            return self.get_next_node()

        # if doesnt overlap; accepted!
        # make Atom object with this
        new_node = Atom(new_coords[0], new_coords[1])

        # set atomnumber
        new_node.n = last_atom.n + 1

        # Get atom type from new char index of amino-string
        new_node.type = self.amino[new_node.n]

        # set the folding move from last atom that creates this one
        last_atom.fold_code = fold_code

        return new_node

    def check_overlap(self, new_coords):
        # check if this node overlaps the chain (new node overlaps any other node)
        for node in self.chain:
            if node.x == new_coords[0] and node.y == new_coords[1]:
                self.overlap_counter += 1
                return True

        # reset overlap_counter if new node is added
        self.overlap_counter = 0
        return False

    def calculate_bonds(self, only_stability):
        # create list of all bondable atoms in the chain
        bondable_atoms = [atom for atom in self.chain if atom.type == "C" or atom.type == "H"]
        atom_nr = 0
        stability = 0
        hh_bonds = []
        ch_bonds = []
        cc_bonds = []

        # check each atom
        for atom in bondable_atoms:
            # compare with each other atom in chain
            for compare_atom in bondable_atoms:
                # if theyre not neighbor nodes
                if (abs(compare_atom.n - atom.n) > 1):
                    dist = math.sqrt((atom.x - compare_atom.x)**2 + (atom.y - compare_atom.y)**2)

                    # if distance is 1 nontheless => bond depending on atom types
                    if dist <= 1:
                        if atom.type == "H" and compare_atom.type == "H":
                            stability -= 1
                            hh_bonds.append([[atom.x, compare_atom.x], [atom.y, compare_atom.y]])

                        if atom.type == "H" and compare_atom.type == "C":
                            stability -= 1
                            ch_bonds.append([[atom.x, compare_atom.x], [atom.y, compare_atom.y]])

                        if atom.type == "C" and compare_atom.type == "C":
                            stability -= 5
                            cc_bonds.append([[atom.x, compare_atom.x], [atom.y, compare_atom.y]])

            # delete this index for double count prevention
            bondable_atoms.pop(atom_nr)
            atom_nr += 1

        if only_stability == False:
            return stability, hh_bonds, ch_bonds, cc_bonds

        return stability

    def get_optimal_node(self):
        stability = self.calculate_bonds(True)
        new_stability = stability
        i = 0

        # if a move makes stability change; keep that move
        while i < self.optimalization_tries and stability == new_stability:
            # get a random new atom
            new_atom = self.get_next_node()

            # append to chain for testing
            self.chain.append(new_atom)

            # calculate new stability with this move
            new_stability = self.calculate_bonds(True)

            # pop atom (otherwise it keeps adding)
            self.chain.pop(-1)

            # break if stability changed
            if new_stability < stability:
                break

            i += 1

        # return the new atom
        return new_atom

chain = AminoLattice(amino)
chain.generate_nodes()
