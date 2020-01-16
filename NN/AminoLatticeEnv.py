import random
from classes.Atom import Atom
import numpy as np
import itertools
import gym
import math
from gym import spaces, logger
from gym.utils import seeding

class AminoLatticeEnv(gym.Env):
    def __init__(self, amino):
        # for stuck nodes
        self.action_overlap = False

        # amino string
        self.amino = amino

        # create atom object for initial node at (0,0) with first string char as type
        self.first_node = Atom(0, 0, 0)
        self.first_node.type = self.amino[0]

        # whole chain array of atom nodes with initial node as start
        self.state = [self.first_node]

        # available moves
        self.moves = [[1, 0, 0], [0, 1, 0], [0, 0, 1],
                      [-1, 0, 0], [0, -1, 0], [0, 0, -1]]

        # folding code corresponding to move index
        self.fold_code_to_index = {"1":0, "2":1, "3":2, "-1":3, "-2":4, "-3":5}
        self.index_to_fold_code = {0:"1", 1:"2", 2:"3", 3:"-1", 4:"-2", 5:"-3"}

        self.actions = []

        # different types of bonds between nodes (consists of coordinates between nodes)
        self.hh_bonds = []
        self.cc_bonds = []
        self.ch_bonds = []

        # general chain stability
        self.stability = 0

        # Grid attributes
        self.grid_length = 2 * len(amino) + 1
        self.midpoint = (len(amino), len(amino))
        self.grid = np.zeros(shape=(self.grid_length, self.grid_length), dtype=int)

        # Define action-observation spaces
        self.action_space = spaces.Discrete(6)
        self.observation_space = spaces.Box(low=-2, high=1,
                                            shape=(self.grid_length, self.grid_length),
                                            dtype=int)
        self.last_action = None
        self.prev_actions = []
        
    # random seed
    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
        
    # reset the chain state
    def reset(self):
        self.first_node = Atom(0, 0, 0)
        self.first_node.type = self.amino[0]
        
        # whole chain array of atom nodes with initial node as start
        self.state = [self.first_node]

        self.overlap_counter = 0
        self.hh_bonds = []
        self.hc_bonds = []
        self.cc_bonds = []

        self.stability = self.get_stability_and_bonds(True)
        observation = np.array([len(self.state), len(self.amino), self.overlap_counter, self.stability])
        done = True if (len(self.state) == len(self.amino)) else False
        reward = self.calculate_reward()
        
        info = {
            'chain_length': len(self.state),
            'seq_length': len(self.amino),
            'actions': [self.moves[i] for i in self.actions],
            'state_chain': self.state
        }
        return observation, reward, done, info        

    ###################################### Returns next Atom object with non-self-overlapping coords. Returns none if stuck.
    def step(self, action):
        # get fold code that creates this move
        fold_code = self.index_to_fold_code[action]

        # last atom needed to generate next one
        last_atom = self.state[-1]

        # last atom coords + random move
        new_x = last_atom.x + self.moves[action][0]
        new_y = last_atom.y + self.moves[action][1]
        new_z = last_atom.z + self.moves[action][2]
        new_coords = [new_x, new_y, new_z]

        # if the new node overlaps its own chain; try new move.
        if self.check_node_overlap(new_coords):
            # constraint relaxation, punish if chain self overlaps
            self.stability += 100
            self.overlap_counter += 1

        # make Atom object
        new_node = Atom(new_coords[0], new_coords[1], new_coords[2])

        # set atomnumber
        new_node.n = last_atom.n + 1

        # Get atom type from new char index of amino-string
        new_node.type = self.amino[new_node.n]

        # set the folding move from last atom that creates this one
        last_atom.fold_code = fold_code
        
        # append to state
        self.state.append(new_node)
        
        self.stability = self.get_stability_and_bonds(True)
        observation = np.array([len(self.state), len(self.amino), self.overlap_counter, self.stability])
        done = True if (len(self.state) == len(self.amino)) else False
        reward = self.calculate_reward()
        
        info = {
            'chain_length': len(self.state),
            'seq_length': len(self.amino),
            'actions': [self.moves[i] for i in self.actions],
            'state_chain': self.state
        }

        return observation, reward, done, info

    ###################################### Checking if new coords overlap the already generated chain
    def check_node_overlap(self, new_coords):
        # check if this node overlaps the chain (new node overlaps any other node)
        for node in self.state:
            if node.x == new_coords[0] and node.y == new_coords[1] and node.z == new_coords[2]:
                self.overlap_counter += 1
                return True

        # reset overlap_counter if new node is added
        return False

    ###################################### Calculates current HH/CH/CC bonds (and coords), and current stability based on current chain
    def get_stability_and_bonds(self, only_stability):
        # create list of all bondable atoms in the chain
        bondable_atoms = [
            atom for atom in self.state if atom.type == "C" or atom.type == "H"
        ]

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

        if only_stability == True:
            return stability

        return stability, hh_bonds, ch_bonds, cc_bonds

    def calculate_reward(self):
        stability = self.get_stability_and_bonds(True)
        
        reward = - stability - (self.overlap_counter * 100) + (len(self.state))
        return reward
        
    ######################################### Sets the current stability and bonds attributes based on the current chain
    def set_stability_and_bonds(self):
        stability, hh_bonds, ch_bonds, cc_bonds = self.get_stability_and_bonds(False)
        self.stability = stability
        self.hh_bonds = hh_bonds
        self.ch_bonds = ch_bonds
        self.cc_bonds = cc_bonds

        return
