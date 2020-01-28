import copy
import random
from algorithms.CalcUpperbound import calc_upperbound
from algorithms.RandomChain import random_chain
from algorithms.PullMove import pull_move
from classes.Grid import Grid


class Hill_climber:
    """
    The hillclimber algorithm. Creates a list of best chains with arguments: 
    amount of resets, amount of stability checks, amount of pulls per stability change
    """

    def __init__(
        self,
        protein,
        amount_of_reset_checks,
        amt_stab_change_checks,
        amt_pulls_per_stab_change_check,
    ):
        self.local_minimum_chains = {}
        self.stability_over_time = []
        self.chain_nr = 1
        self.protein = protein
        self.amount_of_reset_checks = amount_of_reset_checks
        self.amt_stab_change_checks = amt_stab_change_checks
        self.amt_pulls_per_stab_change_check = amt_pulls_per_stab_change_check

        # set minimal stability the chain has to get
        naive_upperbound = calc_upperbound(protein)

    def create_grid_object(self, protein):
        # create initial grid object with random chain in it
        grid_object = Grid(protein)
        grid_object = random_chain(grid_object)

        return grid_object

    def get_current_grid_chain(self, grid_object):
        # only get grid and chain
        best_current_grid = copy.deepcopy(grid_object.grid)
        best_current_chain = copy.deepcopy(grid_object.grid_chain)

        return best_current_grid, best_current_chain

    def pull_move_node(self, grid_object, index):
        # get node object
        node_coords = grid_object.grid_chain[index][0]
        node = grid_object.grid[node_coords].nodes[0]

        # perform a pullmove on this node and update stability and bonds
        pull_move(grid_object, node)

        # update new bonds and stability of this chain
        grid_object.update_all_bonds()

    def copy_best(self, grid_object):
        best_current_grid = copy.deepcopy(grid_object.grid)
        best_current_chain = copy.deepcopy(grid_object.grid_chain)
        best_stability = copy.copy(grid_object.stability)
        better_stab_found = True

        return best_current_grid, best_current_chain, best_stability, better_stab_found

    def stability_check(self, grid_object, best_stability):
        # update best current chain if the stability is better
        if grid_object.stability < best_stability:

            # only get grid and chain for lower computing time
            (
                best_current_grid,
                best_current_chain,
                best_stability,
                better_stab_found,
            ) = self.copy_best(grid_object)

            return (
                best_current_grid,
                best_current_chain,
                best_stability,
                better_stab_found,
            )

    def better_stab_check(
        self,
        better_stab_found,
        best_stability,
        best_current_grid,
        best_current_chain,
        grid_object,
    ):
        # if it didnt find any upgrades after the amount of random pulls of the whole chain, could be local minimum
        if not better_stab_found:
            # terminal info
            print(
                f"No stability change after {self.amt_stab_change_checks} checks. Saving chain {self.chain_nr} as local minima...\n"
            )

            # append to local minima list
            self.local_minimum_chains[copy.copy(best_stability)] = [
                copy.deepcopy(best_current_grid),
                copy.deepcopy(best_current_chain),
            ]

            # reset to new chain
            grid_object = random_chain(grid_object)
            self.chain_nr += 1

            # set this as best
            best_stability = copy.deepcopy(grid_object.stability)
            best_current_grid = copy.deepcopy(grid_object.grid)
            best_current_chain = copy.deepcopy(grid_object.grid_chain)

            return best_stability, best_current_grid, best_current_chain

    def iterator(self, grid_object):

        # Save as initial local minimum
        best_stability = copy.copy(grid_object.stability)
        best_current_grid, best_current_chain = self.get_current_grid_chain(grid_object)

        # how many checks if the chain is a local minimum
        for iteration in range(self.amount_of_reset_checks):
            better_stab_found = False

            # amount of stab change checks before resetting if it didn't change
            for _ in range(self.amt_stab_change_checks):
                self.stability_over_time.append(best_stability)

                # amount of random pulls before checking if new chain has better stability
                for _ in range(self.amt_pulls_per_stab_change_check):
                    # pulling random node
                    index = random.randint(1, len(grid_object.protein) - 2)

                    self.pull_move_node(grid_object, index)

                self.stability_check(grid_object, best_stability)

            print(
                f"Random chain {self.chain_nr}: Nodes pulled: {self.amt_stab_change_checks * self.amt_pulls_per_stab_change_check}. Stability: {best_stability}"
            )

            (
                best_stability,
                best_current_grid,
                best_current_chain,
            ) = self.better_stab_check(
                better_stab_found,
                best_stability,
                best_current_grid,
                best_current_chain,
                grid_object,
            )

    def hill_start(self):

        grid_object = self.create_grid_object(self.protein)

        self.iterator(grid_object)

        return self.local_minimum_chains, self.stability_over_time
