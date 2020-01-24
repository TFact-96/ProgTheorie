import numpy as np
import math
import copy
import random
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D
from upperbound import calc_upperbound
import matplotlib.animation as animation
from matplotlib import style
import itertools


class Node:
    # object for an node in the amino chain
    def __init__(self, x, y, z):

        # coords
        self.x = x
        self.y = y
        self.z = z
        self.n = 0

        # atom type (default P)
        self.type = "P"

        # the fold that this atom makes to the next atom (default 0)
        self.fold_code = 0

        # All neighbours of a single node
        self.neighbours = []

    # print the type if printing the object
    def __repr__(self):
        return f"{self.type}"


class Grid_point:
    def __init__(self, filled, coords):
        self.filled = filled
        self.coords = coords
        self.nodes = []

    def add_node(self, node):
        self.nodes.append(node)
        self.filled = True

    def remove_node(self, node):
        self.nodes.remove(node)
        if self.nodes == []:
            self.filled = False

    def __repr__(self):
        return f"{self.nodes}"


class Grid:
    def __init__(self, amino):

        self.amino = amino
        self.best_chain = {}
        
        # available moves
        self.moves = [[1, 0, 0], [0, 1, 0], [0, 0, 1],
                    [-1, 0, 0], [0, -1, 0], [0, 0, -1]]
                    
        # diagonal moves for chain pulling
        self.diagonal_moves = [[-1, 1, 0], [-1, -1, 0], [-1, 0, 1], [-1, 0, -1],
                                [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1],
                                [1, -1, 0], [1, 1, 0], [1, 0, 1], [1, 0, -1]]

        self.pivot_upperbound = calc_upperbound(amino)
        
    # Create n x n grid
    def create_grid(self, n):
        grid = {}
        for z in range(-n, n + 1):
            for y in range(-n, n + 1):
                for x in range(-n, n + 1):
                    grid[f"{x, y, z}"] = Grid_point(False, [x, y, z])

        return grid

    def clear_grid(self, grid):
        for key, value in grid.items():
            grid[key].filled = False
            grid[key].nodes = []

    def add_neighbours(self, node, grid, hh_bonds):
        x = node.x
        y = node.y
        z = node.z
        n = node.n

        if node.type == "P":
            return False

        if (
            grid[f"{x + 1, y, z}"].filled
            and (abs(grid[f"{x + 1, y, z}"].nodes[0].n - n) > 1)
            and (grid[f"{x + 1, y, z}"].nodes[0].type == "H")
        ):
            if [[x, x + 1], [y, y], [z, z]] and [[x + 1, x], [y, y], [z, z]] not in hh_bonds:
                hh_bonds.append([[x, x + 1], [y, y], [z, z]])

        if (
            grid[f"{x, y + 1, z}"].filled
            and (abs(grid[f"{x, y + 1, z}"].nodes[0].n - n) > 1)
            and (grid[f"{x, y + 1, z}"].nodes[0].type == "H")
        ):
            if [[x, x], [y, y + 1], [z, z]] and [[x, x], [y + 1, y], [z, z]] not in hh_bonds:
                hh_bonds.append([[x, x], [y, y + 1], [z, z]])

        if (
            grid[f"{x, y, z + 1}"].filled
            and (abs(grid[f"{x, y, z + 1}"].nodes[0].n - n) > 1)
            and (grid[f"{x, y, z + 1}"].nodes[0].type == "H")
        ):
            if [[x, x], [y, y], [z + 1, z]] and [[x, x], [y, y], [z + 1, z]] not in hh_bonds:
                hh_bonds.append([[x, x], [y, y], [z + 1, z]])

        if (
            grid[f"{x - 1, y, z}"].filled
            and (abs(grid[f"{x - 1, y, z}"].nodes[0].n - n) > 1)
            and (grid[f"{x - 1, y, z}"].nodes[0].type == "H")
        ):
            if [[x, x - 1], [y, y], [z, z]] and [[x - 1, x], [y, y], [z, z]] not in hh_bonds:
                hh_bonds.append([[x, x - 1], [y, y], [z, z]])

        if (
            grid[f"{x, y - 1, z}"].filled
            and (abs(grid[f"{x, y - 1, z}"].nodes[0].n - n) > 1)
            and (grid[f"{x, y - 1, z}"].nodes[0].type == "H")
        ):
            if [[x, x], [y, y - 1], [z, z]] and [[x, x], [y - 1, y], [z, z]] not in hh_bonds:
                hh_bonds.append([[x, x], [y, y - 1], [z, z]])

        if (
            grid[f"{x, y, z - 1}"].filled
            and (abs(grid[f"{x, y, z - 1}"].nodes[0].n - n) > 1)
            and (grid[f"{x, y, z - 1}"].nodes[0].type == "H")
        ):
            if [[x, x], [y, y], [z - 1, z]] and [[x, x], [y, y], [z - 1, z]] not in hh_bonds:
                hh_bonds.append([[x, x], [y, y], [z - 1, z]])

    def update_neighbours(self, grid, grid_chain):
        hh_bonds = []
        for key, value in grid_chain.items():
            self.add_neighbours(grid[value[0]].nodes[0], grid, hh_bonds)

        return -len(hh_bonds), hh_bonds

    def add_point(self, node, n, grid, grid_chain):
        x = node.x
        y = node.y
        z = node.z
        grid[f"{x, y, z}"].add_node(node)
        grid_chain[n] = [f"{x, y, z}", [x, y, z]]

    def clear_point(self, node, n, grid):
        x = node.x
        y = node.y
        z = node.z
        grid[f"{x, y, z}"].remove_node(node)

    def transfer_point(self, node1, x2, y2, z2, grid, grid_chain):
        n = node1.n
        self.clear_point(node1, n, grid)
        node1.x = x2
        node1.y = y2
        node1.z = z2
        self.add_point(node1, n, grid, grid_chain)

    def overlap(self, x, y, z, grid):
        if grid[f"{x, y, z}"].filled:
            return True

        return False

    def chain_stuck(self, x, y, z, grid):

        if (
            self.overlap(x + 1, y, z, grid)
            and self.overlap(x - 1, y, z, grid)
            and self.overlap(x, y + 1, z, grid)
            and self.overlap(x, y - 1, z, grid)
            and self.overlap(x, y, z + 1, grid)
            and self.overlap(x, y, z - 1, grid)
        ):
            return True

        return False

    def check_diagonals(self, x, y, z, grid):
        # diagonal moves for chain pulling
        #self.diagonal_moves = [[-1, 1, 0], [-1, -1, 0], [-1, 0, 1], [-1, 0, -1],
        #                        [0, 1, 1], [0, 1, -1], [0, -1, 1], [0, -1, -1],
        #                        [1, -1, 0], [1, 1, 0], [1, 0, 1], [1, 0, -1]]

        available_moves = []

        if not grid[f"{x + 1, y + 1, z}"].filled:
            available_moves.append([1, 1, 0])

        if not grid[f"{x + 1, y - 1, z}"].filled:
            available_moves.append([1, -1, 0])

        if not grid[f"{x - 1, y + 1 , z}"].filled:
            available_moves.append([-1, 1, 0])

        if not grid[f"{x - 1, y - 1, z}"].filled:
            available_moves.append([-1, -1, 0])

        if not grid[f"{x, y + 1, z + 1}"].filled:
            available_moves.append([0, 1, 1])

        if not grid[f"{x, y + 1, z - 1}"].filled:
            available_moves.append([0, 1, -1])

        if not grid[f"{x, y - 1, z + 1}"].filled:
            available_moves.append([0, -1, 1])

        if not grid[f"{x, y - 1, z - 1}"].filled:
            available_moves.append([0, -1, -1])

        if not grid[f"{x + 1, y, z + 1}"].filled:
            available_moves.append([1, 0, 1])

        if not grid[f"{x + 1, y, z - 1}"].filled:
            available_moves.append([1, 0, -1])

        if not grid[f"{x - 1, y, z + 1}"].filled:
            available_moves.append([-1, 0, 1])

        if not grid[f"{x - 1, y, z - 1}"].filled:
            available_moves.append([-1, 0, -1])

        return available_moves

    def create_chain(self):
        index = 0
        grid_chain = {}

        # Create grid
        grid = self.create_grid(36)
        first_node = Node(0, 0, 0)
        first_node.type = self.amino[0]
        self.add_point(first_node, index, grid, grid_chain)

        while index < len(self.amino) - 1:
            current_node_key = grid_chain[index][0]
            current_node = grid[current_node_key].nodes[0]

            random_move = np.random.randint(len(self.moves))

            new_x = current_node.x + self.moves[random_move][0]
            new_y = current_node.y + self.moves[random_move][1]
            new_z = current_node.z + self.moves[random_move][2]

            if not self.overlap(new_x, new_y, new_z, grid):
                index += 1
                new_node = Node(new_x, new_y, new_z)
                new_node.n = index
                new_node.type = self.amino[index]
                self.add_point(new_node, index, grid, grid_chain)

            if self.chain_stuck(new_x, new_y, new_z, grid):
                index = 0
                self.clear_grid(grid)
                grid_chain = {}
                self.add_point(first_node, index, grid, grid_chain)

        return grid_chain, grid

    def create_vectors(self, node, grid_chain):
        node_i_coords = [node.x, node.y, node.z]
        node_i1 = grid_chain[int(node.n) + 1]
        node_i1_coords = np.array([node_i1[1][0], node_i1[1][1], node_i1[1][2]])
        vector1 = node_i1_coords - np.array(node_i_coords)

        return node_i_coords, node_i1_coords, vector1

    def check_requirements(self, available_moves, vector1, node_i_coords, node_i1_coords, grid):
        viable_moves = []
        found = False
        
        for move in available_moves:
            L = node_i_coords + np.array(move)
            C = L - vector1

            if (
                (not self.overlap(L[0], L[1], L[2], grid))
                and (not self.overlap(C[0], C[1], C[2], grid))
                and (np.linalg.norm(L - node_i1_coords) == 1.0)
            ):
                viable_moves.append([L, C])
                found = True

        if not found:
            return 0, 0, False

        random_choice = random.choice(viable_moves)
        return random_choice[0], random_choice[1], True

    def pull_move(self, node, grid, grid_chain):

        node_i_coords, node_i1_coords, vector1 = self.create_vectors(node, grid_chain)

        available_moves = self.check_diagonals(node_i_coords[0], node_i_coords[1], node_i_coords[2], grid)

        L, C, check = self.check_requirements(
            available_moves, vector1, node_i_coords, node_i1_coords, grid
        )

        if check:

            # residue of chain follows in footsteps
            for index in range(int(node.n - 1)):

                residue_node_key = grid_chain[index][0]
                residue_node = grid[residue_node_key].nodes[0]

                residue_node_next_key = grid_chain[index + 2][0]
                residue_node_next = grid[residue_node_next_key].nodes[0]

                self.transfer_point(
                    residue_node,
                    residue_node_next.x,
                    residue_node_next.y,
                    residue_node_next.z,
                    grid,
                    grid_chain,
                )

            # node moves to L
            self.transfer_point(node, L[0], L[1], L[2], grid, grid_chain)

            # Previous node moves to C
            previous_node_key = grid_chain[int(node.n) - 1][0]
            previous_node = grid[previous_node_key].nodes[0]

            self.transfer_point(previous_node, C[0], C[1], C[2], grid, grid_chain)

        return grid, grid_chain

    def hill_climber(self, max_iteration):
        
        # Current protein chain
        current_hilltop, grid = self.create_chain()        
        current_stability = self.update_neighbours(grid, current_hilltop)[0]
        
        # Save as best chain
        best_c = current_hilltop
        best_stab_c = current_stability
        best_grid = grid
        self.best_chain[best_stab_c] = [best_c, best_grid]

        for iteration in range(max_iteration):
            best_c_found = False
            print(iteration)
            print(self.pivot_upperbound)

            for it in range(100):
                print(f"second: {it}")
                for index in range(1, len(current_hilltop) - 1):
                    print(f"deepest: {index}")
                    temp_grid, temp_chain = self.pull_move(
                        grid[current_hilltop[index][0]].nodes[0], grid, current_hilltop,
                    )

                    temp_stability = self.update_neighbours(temp_grid, temp_chain)[0]

                    if temp_stability < best_stab_c and temp_stability < self.pivot_upperbound:
                        print("pullmove made better stab!")
                        best_c = copy.deepcopy(temp_chain)
                        best_grid = copy.deepcopy(temp_grid)
                        best_stab_c = temp_stability
                        best_c_found = True

            if best_c_found:
                current_hilltop = copy.deepcopy(best_c)
                current_stability = best_stab_c
                grid = copy.deepcopy(best_grid)
                best_c_found = False

            else:
                self.best_chain[best_stab_c] = [best_c, best_grid]

                current_hilltop, grid = self.create_chain()
                current_stability = self.update_neighbours(grid, current_hilltop)[0]

                best_c = current_hilltop
                best_stab_c = current_stability
                best_grid = grid

    def find_best_c(self):
        best_chain_key = min(self.best_chain.keys())
        best_chain_double = self.best_chain[best_chain_key]

        best_chain = best_chain_double[0]
        best_grid = best_chain_double[1]

        best_stability, best_hh = self.update_neighbours(best_grid, best_chain)
        print(best_stability)
        
        # prepare for plotting
        x = []
        y = []
        z = []
        color = []

        # get coords and node colors of chain
        for key, value in best_chain.items():
            node = best_grid[value[0]].nodes[0]

            x.append(node.x)
            y.append(node.y)
            z.append(node.z)

            # set amino colors
            if node.type == "H":
                color.append('red')
                
            elif node.type == "C":
                color.append('yellow')
            
            else:
                color.append('blue')


        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # plot the bond lines
        for hh_bond in best_hh:
            ax.plot3D(
                hh_bond[0], hh_bond[1], hh_bond[2],
                "--",
                markersize=1,
                color='red',
                zorder=-1
            )

        # plot text at start and ending atom
        ax.text(
             x[0], y[0], z[0] + 0.1,
             "Start"
        )

        ax.text(
             x[-1], y[-1], z[-1] + 0.1,
             "Finish"
        )

        # plot the chain itself
        ax.plot3D(
            x, y, z,
            "-",
            linewidth=3,
            color='black',
            zorder=0,
        )

        ax.scatter3D(
            x, y, z,
            color=color,
            edgecolor='black',
            s=100,
            depthshade=False
        )

        custom_legend = [
                            Line2D([0], [0], color='black', lw=2, label="Protein chain"),
                            Line2D([0], [0], marker='o', color='black', label='H-amino',
                                   markerfacecolor='red', markersize=8),
                            Line2D([0], [0], marker='o', color='black', label='P-amino',
                                   markerfacecolor='blue', markersize=8),
                            Line2D([0], [0], marker='o', color='black', label='C-amino',
                                   markerfacecolor='yellow', markersize=8),
                            Line2D([0], [0], linestyle="--", color='red', lw=1, label="H-H (Stability -1)"),
                            Line2D([0], [0], linestyle="--", color='orange', lw=1, label="H-C (Stability -1)"),
                            Line2D([0], [0], linestyle="--", color='yellow', lw=1, label="C-C (Stability -5)"),

                        ]

        # Hide grid lines
        ax.grid(False)

        # Hide axes ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        # plt axis and grid off
        plt.axis('off')
        plt.grid(b=None)

        ax.legend(handles=custom_legend)
        ax.set_title(f"Protein chain (Stability: {best_stability})")
        plt.show()                



k = Grid("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
k.hill_climber(1)
k.find_best_c()
