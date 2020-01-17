import random
import math
from classes.Amino import Amino2D, Amino3D
import numpy as np
import itertools
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from matplotlib import style
import copy
import matplotlib.pyplot as plt

class ChainLattice:
    def __init__(self, protein, ThreeD):
        # 2D or 3D
        self.ThreeD = ThreeD

        # for stuck aminos
        self.overlap_counter = 0
        self.state_stuck = False

        # protein string
        self.protein = protein

        # create amino object for initial amino at (0,0) with first string char as type
        if self.ThreeD:
            self.first_amino = Amino3D(0, 0, 0)
        else:
            self.first_amino = Amino2D(0, 0)

        self.first_amino.type = self.protein[0]

        # whole chain array of amino aminos with initial amino as start
        self.state = {0: self.first_amino}


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
            self.index_to_fold_code = {0:"1", 1:"2", 2:"-1", 3:"-2"}
            self.fold_code_to_index = {"1":0, "2":1, "-1":2, "-2":3}


        # different types of bonds between aminos (consists of coordinates between aminos)
        self.hh_bonds = []
        self.cc_bonds = []
        self.ch_bonds = []

        # general state stability
        self.stability = 0
        self.old_stability = 0

    ###################################### Returns next amino object with non-self-overlapping coords. Returns none if stuck.
    def generate_amino_random_move(self):
        # if more than 100 random moves are tried (see check_overlap()); almost 100% chance that state is stuck; Exit program
        if self.overlap_counter > 100:
            self.state_stuck = True
            return

        # get random move of new amino added to state
        random_move_index = np.random.randint(len(self.moves))

        # get fold code that creates this move
        fold_code = self.index_to_fold_code[random_move_index]

        # last amino needed to generate next one
        last_amino = self.state[len(self.state) - 1]

        # last amino coords + random move
        if self.ThreeD:
            new_x = last_amino.x + self.moves[random_move_index][0]
            new_y = last_amino.y + self.moves[random_move_index][1]
            new_z = last_amino.z + self.moves[random_move_index][2]
            new_coords = [new_x, new_y, new_z]
        else:
            new_x = last_amino.x + self.moves[random_move_index][0]
            new_y = last_amino.y + self.moves[random_move_index][1]
            new_coords = [new_x, new_y]

        # if the new amino overlaps its own state; try new move.
        if self.check_amino_overlap(new_coords):
            return self.generate_amino_random_move()

        # no overlap: accepted, create an amino object at this coord
        return self.create_amino_object(new_coords, last_amino, fold_code)

    ###################################### Checking if new coords overlap the already generated state
    def check_amino_overlap(self, new_coords):
        # check if this amino overlaps the state (new amino overlaps any other amino)
        for amino_key, amino in self.state.items():
            if self.ThreeD:
                if [amino.x, amino.y, amino.z] == [new_coords[0], new_coords[1], new_coords[2]]:
                    self.overlap_counter += 1
                    return True
            else:
                if [amino.x, amino.y] == [new_coords[0], new_coords[1]]:
                    self.overlap_counter += 1
                    return True

        # reset overlap_counter if new amino is added
        self.overlap_counter = 0
        return False

    ####################################### creating a next amino object based on the last amino and fold code to get to new coords
    def create_amino_object(self, new_coords, last_amino, fold_code):
        # make amino object
        if self.ThreeD:
            new_amino = Amino3D(new_coords[0], new_coords[1], new_coords[2])
        else:
            new_amino = Amino2D(new_coords[0], new_coords[1])

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
            if self.ThreeD:
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
            else:
                # Otherwise 2D
                if (abs(compare_amino.n - amino.n) > 1):
                    dist = math.sqrt(
                        (amino.x - compare_amino.x)**2
                        + (amino.y - compare_amino.y)**2
                    )

                    # if distance is 1 nontheless => bond depending on amino types
                    if dist <= 1:
                        if amino.type == "H" and compare_amino.type == "H":
                            stability -= 1

                            if not only_stability:
                                hh_bonds.append(
                                    [[amino.x, compare_amino.x],
                                    [amino.y, compare_amino.y]])

                        if (amino.type == "H" and compare_amino.type == "C") or (amino.type == "C" and compare_amino.type == "H"):
                            stability -= 1

                            if not only_stability:
                                ch_bonds.append(
                                    [[amino.x, compare_amino.x],
                                    [amino.y, compare_amino.y]]
                                    )

                        if amino.type == "C" and compare_amino.type == "C":
                            stability -= 5

                            if not only_stability:
                                cc_bonds.append(
                                    [[amino.x, compare_amino.x],
                                    [amino.y, compare_amino.y]]
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


    ############################################## MEHMET'S HILL CLIMB, only 2D
    # checks if a point has a amino or not
    def check_point(self, array):
        for index, aminos in self.state.items():
            if str(array) == str(np.array([aminos.x, aminos.y])):
                return False
        return True

    def pull_move(self, amino):

        amino_i_coords = np.array([amino.x, amino.y])
        amino_i1 = self.state[int(amino.n) + 1]
        amino_i1_coords = np.array([amino_i1.x, amino_i1.y])
        vector1 = amino_i1_coords - amino_i_coords

        state_copy = copy.deepcopy(self.state)
        self.old_stability = self.stability

        checker = {
            "[1, 1]": [True, [1, 1]],
            "[1, -1]": [True, [1, -1]],
            "[-1, 1]": [True, [-1, 1]],
            "[-1, -1]": [True, [-1, -1]],
        }

        for d_move in self.diagonal_moves:
            if not self.check_point(amino_i_coords + np.array(d_move)):
                checker[str(d_move)][0] = False

        for d_check, d_bool in checker.items():
            L = amino_i_coords + np.array(d_bool[1])
            C = L - vector1

            # checks pull move requirements
            if (
                (d_bool[0] == True)
                and (np.linalg.norm(L - amino_i1_coords) == 1.0)
                and (self.check_point(C))
            ):
                # residue of state follows in footsteps
                for index in range(int(amino.n - 1)):

                    residue_amino = self.state[index]
                    residue_amino_next = self.state[index + 2]

                    residue_amino.x = residue_amino_next.x
                    residue_amino.y = residue_amino_next.y

                # amino moves to L
                amino.x = L[0]
                amino.y = L[1]

                # Previous amino moves to C
                previous_amino = self.state[int(amino.n) - 1]
                previous_amino.x = C[0]
                previous_amino.y = C[1]
                break

        # calculate self.stability, calculate all bonds and put coords in cc/ch/hh_bonds list. All with current self.state
        self.set_stability_and_bonds()

        if self.stability > self.old_stability:
            self.state = state_copy
            # reset stability and bonds
            self.set_stability_and_bonds()

        elif self.stability < self.old_stability:
            print(f"Pull on amino {amino.n} made better stability: {self.stability}")

    def random_pull(self, pull_times_per_chain):
        for x in range(pull_times_per_chain):
            d = random.randint(1, len(self.state) - 2)
            self.pull_move(self.state[d])
