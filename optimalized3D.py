import math
import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Atom:
    # object for an atom in the amino chain
    def __init__(self, x, y, z):
        # nth atom in the chain
        self.n = 0

        # coords
        self.x = x
        self.y = y
        self.z = z

        # atom type (default P)
        self.type = "P"

        # the fold that this atom makes to the next atom (default 0)
        self.fold_code = 0

        # default color (P)
        self.color = 'blue'

    # print the type if printing the object
    def __str__(self):
        return f"{self.type}"


class AminoLattice:
    def __init__(self, amino, random_move, optimal_move):
        # either use random moves per new atom added or an optimal move
        self.random_move = random_move
        self.optimal_move = optimal_move

        # amount of tries to make atom find an atom bond with its random move
        self.optimalization_tries = optimalization_tries

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

        # different types of bonds between nodes (consists of coordinates between nodes)
        self.hh_bonds = []
        self.cc_bonds = []
        self.ch_bonds = []

        # general chain stability
        self.stability = 0

    ###################################### Generating all atoms one by one onto the lattice
    def generate_chain(self):
        # while amount of nodes already generated is smaller than the whole amino string length
        while len(self.chain) < len(self.amino) and not self.chain_stuck:

            if self.random_move == True:
                # get a random move per new atom added to chain
                new_node = self.get_next_node()
            elif self.optimal_move == True:
                # generate a optimalized move (new node always makes bond-making moves if possible)
                new_node = self.get_optimal_node()

            # append the new node to the existing chain
            self.chain.append(new_node)

    ###################################### Returns next Atom object with non-self-overlapping coords. Returns none if stuck.
    def get_next_node(self):
        # get random move of new node added to chain
        random_move_index = np.random.randint(len(self.moves))
        # last atom needed to generate next one
        last_atom = self.chain[-1]

        # last atom coords + random move
        new_x = last_atom.x + self.moves[random_move_index][0]
        new_y = last_atom.y + self.moves[random_move_index][1]
        new_z = last_atom.z + self.moves[random_move_index][2]
        new_coords = [new_x, new_y, new_z]

        # get fold code that creates this move
        fold_code = self.move_code[random_move_index]

        # if more than 100 random moves are tried (see check_overlap()); almost 100% chance that chain is stuck; Exit program
        if self.overlap_counter > 100:
            self.chain_stuck = True
            return

        # if the new node overlaps its own chain; try new move.
        if self.check_overlap(new_coords):
            return self.get_next_node()

        # if doesnt overlap; accepted!
        # make Atom object with this
        new_node = Atom(new_coords[0], new_coords[1], new_coords[2])

        # set atomnumber
        new_node.n = last_atom.n + 1

        # Get atom type from new char index of amino-string
        new_node.type = self.amino[new_node.n]

        # set the folding move from last atom that creates this one
        last_atom.fold_code = fold_code

        return new_node

    ###################################### Checking if new coords overlap the already generated chain
    def check_overlap(self, new_coords):
        # check if this node overlaps the chain (new node overlaps any other node)
        for node in self.chain:
            if node.x == new_coords[0] and node.y == new_coords[1] and node.z == new_coords[2]:
                self.overlap_counter += 1
                return True

        # reset overlap_counter if new node is added
        self.overlap_counter = 0
        return False

    ###################################### Calculates current HH/CH/CC bonds (and coords), and current stability
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
                    dist = math.sqrt((atom.x - compare_atom.x)**2 + (atom.y - compare_atom.y)**2 + (atom.z - compare_atom.z)**2)

                    # if distance is 1 nontheless => bond depending on atom types
                    if dist <= 1:
                        if atom.type == "H" and compare_atom.type == "H":
                            stability -= 1

                            if not only_stability:
                                hh_bonds.append([[atom.x, compare_atom.x], [atom.y, compare_atom.y], [atom.z, compare_atom.z]])

                        if (atom.type == "H" and compare_atom.type == "C") or (atom.type == "C" and compare_atom.type == "H"):
                            stability -= 1

                            if not only_stability:
                                ch_bonds.append([[atom.x, compare_atom.x], [atom.y, compare_atom.y], [atom.z, compare_atom.z]])

                        if atom.type == "C" and compare_atom.type == "C":
                            stability -= 5

                            if not only_stability:
                                cc_bonds.append([[atom.x, compare_atom.x], [atom.y, compare_atom.y], [atom.z, compare_atom.z]])

            # delete this index for double count prevention
            bondable_atoms.pop(atom_nr)
            atom_nr += 1

        if only_stability == True:
            return stability

        return stability, hh_bonds, ch_bonds, cc_bonds

    ###################################### OPTIMALIZATION ALGORITHM (The more C's and H's, the better it works)
    def get_optimal_node(self):
        stability = self.calculate_bonds(True)
        new_stability = stability
        i = 0

        # if a move makes stability change; keep that move
        while i < self.optimalization_tries and stability == new_stability:
            # get a random new atom
            new_atom = self.get_next_node()

            if new_atom:
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
    ###################################### END OPTIMALIZATION ALGORITHM

    ###################################### Returns resulting stability of whole chain and a tuple list of fold codes per atom
    def get_chain_data(self, print_data):
        # cant get data if chain was stuck
        if self.chain_stuck:
            return

        stability = self.calculate_bonds(True)
        moves = [[node.type, node.fold_code] for node in self.chain]

        if print_data:
            print(f"\nThe stability of this amino-acid is {stability}.\n")
            print(f"Its moves are:")

            for node in moves:
                print(f"{node[0]} {node[1]}")

        return stability, moves

    ###################################### Preparing lists for plotting amino chain
    def get_plot_data(self):
        # prepare for plotting
        x = []
        y = []
        z = []
        color = []

        # append all nodes in coord list
        for node in self.chain:
            x.append(node.x)
            y.append(node.y)
            z.append(node.z)

            # set node colors
            if node.type == "H":
                node.color = 'red'
            elif node.type == "C":
                node.color = 'orange'

            color.append(node.color)

        # calculate the bonds in this chain
        self.stability, self.hh_bonds, self.ch_bonds, self.cc_bonds = self.calculate_bonds(False)

        return x, y, z, color

    ###################################### Plotting the chain
    def plot_chain(self):
        x, y, z, color = self.get_plot_data()

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # plot the bond lines
        for hh_bond in self.hh_bonds:
            ax.plot3D(hh_bond[0], hh_bond[1], hh_bond[2], "--", markersize=1, color='red', zorder=1)

        for ch_bond in self.ch_bonds:
            ax.plot3D(ch_bond[0], ch_bond[1], ch_bond[2], "--", markersize=1, color='orange', zorder=1)

        for cc_bond in self.cc_bonds:
            ax.plot3D(cc_bond[0], cc_bond[1], cc_bond[2], "--", markersize=1, color='yellow', zorder=1)

        # plot the chain itself
        ax.plot3D(x, y, z, "-", linewidth=3, color='black', zorder=1, label=f"Amino chain (Stability: {self.stability})")
        ax.scatter3D(x, y, z, color=color, zorder=2)

        # plot atomtype name at its node coord
        for i in range(len(x)):
            ax.text(x[i] + 0.05, y[i] + 0.05, z[i] + 0.05, self.chain[i].type)

        plt.show()

        # rest info and show plot

###################################### Just get one calculation of the chain
def calculate_one_chain(random_move, optimal_move):
    chain = AminoLattice(amino, random_move, optimal_move)
    chain.generate_chain()

    if chain.chain_stuck:
        print("Chain got stuck! Try again.")
    else:
        chain.plot_chain()

###################################### (bruteforce) ALGORITHM
def iterate_generations_of_chains(random_move, optimal_move):
    print("\nBrute force generating chains:")
    print(f"Generating {iterations} chains...")
    best_stability = 0
    best_chain = AminoLattice(amino, random_move, optimal_move)

    # try n iterations for best stability
    for i in range(iterations):
        chain = AminoLattice(amino, random_move, optimal_move)
        chain.generate_chain()

        # only count non-stuck chains
        if not chain.chain_stuck:
            stability, moves = chain.get_chain_data(False)

            # if this generation is a new record
            if stability <= best_stability:
                best_chain = chain
                best_stability = stability
                print(f"Generation {i}: Stability {best_stability}.")


    # print this chain
    best_chain.get_chain_data(True)

    # plot this chain
    best_chain.plot_chain()
###################################### END (bruteforce) ALGORITHM

###################################### Main structure
if __name__ == "__main__":
    amino = input("\nEnter desired amino-chain (C's, H's and P's): ")
    str = input("Random generation or optimal generation? (y = optimal / n = random): ")
    str_brute = input("Brute force generation to find optimal amino fold? (y/n): ")

    if (str != "y" and str != "n") or (str_brute != "y" and str_brute != "n"):
        print("Only answer with y or n please.")
        exit(0)

    if str == "n":
        if str_brute == "y":
            iterations = int(input("How many chain generations for brute forcing?: "))
            iterate_generations_of_chains(True, False)
        else:
            calculate_one_chain(True, False)

    if str == "y":
        if str_brute == "y":
            iterations = int(input("How many chain generations for brute forcing?: "))
            optimalization_tries = int(input("How many times should a node try for an optimal move?: "))
            iterate_generations_of_chains(False, True)
        else:
            optimalization_tries = int(input("How many times should a node try for an optimal move?: "))
            calculate_one_chain(False, True)
###################################### End main structure
