# -*- coding: utf-8 -*-

"""
Implements the 2D Lattice Environment
"""
# Import gym modules
import sys
from math import floor
from collections import OrderedDict
import matplotlib.pyplot as plt
from six import StringIO

import gym
from gym import (spaces, utils, logger)
import numpy as np

# Human-readable
ACTION_TO_STR = {
    0: 'L', 1: 'D',
    2: 'U', 3: 'R'}

POLY_TO_INT = {
    'H': 1, 'P': -1
}


class Lattice2DEnv(gym.Env):
    def __init__(self, seq="HPPHHHP", collision_penalty=-2, trap_penalty=0.5):

        self.seq = seq.upper()

        self.collision_penalty = collision_penalty
        self.trap_penalty = trap_penalty

        self.state = OrderedDict({(0, 0): self.seq[0]})
        self.actions = []
        self.collisions = 0
        self.trapped = 0

        # Grid attributes
        self.grid_length = 2 * len(seq) + 1
        self.midpoint = (len(seq), len(seq))
        self.grid = np.zeros(shape=(self.grid_length, self.grid_length), dtype=int)

        # Automatically assign first element into grid
        self.grid[self.midpoint] = POLY_TO_INT[self.seq[0]]

        # Define action-observation spaces
        self.action_space = spaces.Discrete(4)
        self.observation_space = spaces.Box(low=-2, high=1,
                                            shape=(self.grid_length, self.grid_length),
                                            dtype=int)
        self.last_action = None

    def step(self, action):
        if not self.action_space.contains(action):
            raise ValueError("%r (%s) invalid" % (action, type(action)))

        self.last_action = action
        is_trapped = False  # Trap signal
        collision = False  # Collision signal

        # Obtain coordinate of previous polymer
        x, y = next(reversed(self.state))
        # Get all adjacent coords and next move based on action
        adj_coords = self._get_adjacent_coords((x, y))
        next_move = adj_coords[action]
        # Detects for collision or traps in the given coordinate
        idx = len(self.state)
        if set(adj_coords.values()).issubset(self.state.keys()):
            self.trapped += 1
            is_trapped = True
        elif next_move in self.state:
            self.collisions += 1
            collision = True
        else:
            self.actions.append(action)
            try:
                self.state.update({next_move: self.seq[idx]})
            except IndexError:
                logger.error('All molecules have been placed! Nothing can be added to the protein chain.')
                raise

        # Set-up return values
        obs = np.array([len(self.state), len(self.seq), self.collisions, is_trapped])
        done = True if (len(self.state) == len(self.seq) or is_trapped) else False
        reward = self._compute_reward(is_trapped, collision, done)
        info = {
            'chain_length': len(self.state),
            'seq_length': len(self.seq),
            'collisions': self.collisions,
            'actions': [ACTION_TO_STR[i] for i in self.actions],
            'is_trapped': is_trapped,
            'state_chain': self.state
        }

        return obs, reward, done, info

    def reset(self):
        """Resets the environment"""
        self.state = OrderedDict({(0, 0): self.seq[0]})
        self.actions = []
        self.collisions = 0
        self.trapped = 0
        self.grid = np.zeros(shape=(self.grid_length, self.grid_length), dtype=int)
        # Automatically assign first element into grid
        self.grid[self.midpoint] = POLY_TO_INT[self.seq[0]]

        obs = np.array([len(self.state), len(self.seq), self.collisions, self.trapped])

        return obs

    def render(self, mode='human'):
        """Renders the environment"""

        outfile = StringIO() if mode == 'ansi' else sys.stdout
        desc = self.grid.astype(str)

        # Convert everything to human-readable symbols
        desc[desc == '0'] = '*'
        desc[desc == '1'] = 'H'
        desc[desc == '-1'] = 'P'

        # Obtain all x-y indices of elements
        x_free, y_free = np.where(desc == '*')
        x_h, y_h = np.where(desc == 'H')
        x_p, y_p = np.where(desc == 'P')

        # Decode if possible
        desc.tolist()
        try:
            desc = [[c.decode('utf-8') for c in line] for line in desc]
        except AttributeError:
            pass

        # All unfilled spaces are gray
        for unfilled_coords in zip(x_free, y_free):
            desc[unfilled_coords] = utils.colorize(desc[unfilled_coords], "gray")

        # All hydrophobic molecules are bold-green
        for hmol_coords in zip(x_h, y_h):
            desc[hmol_coords] = utils.colorize(desc[hmol_coords], "green", bold=True)

        # All polar molecules are cyan
        for pmol_coords in zip(x_p, y_p):
            desc[pmol_coords] = utils.colorize(desc[pmol_coords], "cyan")

        # Provide prompt for last action
        if self.last_action is not None:
            outfile.write("  ({})\n".format(["Left", "Down", "Up", "Right"][self.last_action]))
        else:
            outfile.write("\n")

        # Draw desc
        outfile.write("\n".join(''.join(line) for line in desc)+"\n")

        if mode != 'human':
            return outfile

    def _get_adjacent_coords(self, coords):
        x, y = coords
        adjacent_coords = {
            0: (x - 1, y),
            1: (x, y - 1),
            2: (x, y + 1),
            3: (x + 1, y),
        }

        return adjacent_coords

    def _draw_grid(self, chain):

        for coord, poly in chain.items():
            trans_x, trans_y = tuple(sum(x) for x in zip(self.midpoint, coord))
            # Recall that a numpy array works by indexing the rows first
            # before the columns, that's why we interchange.
            self.grid[(trans_y, trans_x)] = POLY_TO_INT[poly]

        return np.flipud(self.grid)

    def _compute_reward(self, is_trapped, collision, done):

        state_reward = self._compute_free_energy(self.state) if done else 0
        collision_penalty = self.collision_penalty if collision else 0
        actual_trap_penalty = -floor(len(self.seq) * self.trap_penalty) if is_trapped else 0

        # Compute reward at timestep, the state_reward is originally
        # negative (Gibbs), so we invert its sign.
        reward = - state_reward + collision_penalty + actual_trap_penalty

        return reward

    def _compute_free_energy(self, chain):

        h_polymers = [x for x in chain if chain[x] == 'H']
        h_pairs = [(x, y) for x in h_polymers for y in h_polymers]

        # Compute distance between all hydrophobic pairs
        h_adjacent = []
        for pair in h_pairs:
            dist = np.linalg.norm(np.subtract(pair[0], pair[1]))
            if dist == 1.0:  # adjacent pairs have a unit distance
                h_adjacent.append(pair)

        # Get the number of consecutive H-pairs in the string,
        # these are not included in computing the energy
        h_consecutive = 0
        for i in range(1, len(self.state)):
            if (self.seq[i] == 'H') and (self.seq[i] == self.seq[i - 1]):
                h_consecutive += 1

        # Remove duplicate pairs of pairs and subtract the
        # consecutive pairs
        nb_h_adjacent = len(h_adjacent) / 2
        gibbs_energy = nb_h_adjacent - h_consecutive
        reward = - gibbs_energy
        return int(reward)