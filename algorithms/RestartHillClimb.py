import copy
import random
from algorithms.CalcUpperbound import calc_upperbound
from algorithms.RandomChain import random_chain
from algorithms.PullMove import pull_move
from classes.Grid import Grid


class Hill_climber:
    """
    The hillclimber algorithm. Creates a list of best chains.
    :param amount_of_reset_checks: amount of resets
    :param amt_stab_change_checks: amount of stability checks
    :param amt_pulls_per_stab_change_check: amount of pulls per stability change
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
        """
        Create initial grid object with random chain in it.
        :param protein: protein chain
        :return: grid object
        """
        grid_object = Grid(protein)
        grid_object = random_chain(grid_object)

        return grid_object

    def get_current_grid_chain(self, grid_object):
        """
        Deepcopies grid and chain.
        :param grid_object: grid object
        :return: best current grid, best current chain
        """
        best_current_grid = copy.deepcopy(grid_object.grid)
        best_current_chain = copy.deepcopy(grid_object.grid_chain)

        return best_current_grid, best_current_chain

    def pull_move_node(self, grid_object, index):
        """
        Performs pullmove on node.
        :param grid_object: grid object
        :param index: index
        """

        node_coords = grid_object.grid_chain[index][0]
        node = grid_object.grid[node_coords].nodes[0]

        pull_move(grid_object, node)

        # update new bonds and stability of this chain
        grid_object.update_all_bonds()

    def copy_best(self, grid_object):
        """
        Deepcopies best chain and grid.
        :param grid_object: grid object
        """
        best_current_grid = copy.deepcopy(grid_object.grid)
        best_current_chain = copy.deepcopy(grid_object.grid_chain)
        best_stability = copy.copy(grid_object.stability)
        better_stab_found = True

        return best_current_grid, best_current_chain, best_stability, better_stab_found

    def print_message(self, message_num, best_stability):
        """
        Prints messages.
        :param message_num: message number
        :param best_stability: best stability
        """
        if message_num == 1:
            print(
                f"Random chain {self.chain_nr}: Nodes pulled: {self.amt_stab_change_checks * self.amt_pulls_per_stab_change_check}. Stability: {best_stability}"
            )
        elif message_num == 2:
            print(
                f"No stability change after {self.amt_stab_change_checks} checks. Saving chain {self.chain_nr} as local minima...\n"
            )

    def save_local_minima(
        self, best_stability, best_current_grid, best_current_chain, grid_object
    ):
        """
        Saves local minima to list and creates new grid object.
        :param best_stability: best stability
        :param best_current_grid: best current grid
        :param best_current_chain: best current chain
        :param grid_object: grid object
        :return: grid object
        """
        self.local_minimum_chains[copy.copy(best_stability)] = [
            copy.deepcopy(best_current_grid),
            copy.deepcopy(best_current_chain),
        ]

        grid_object = random_chain(grid_object)
        self.chain_nr += 1
        return grid_object

    def hill_start(self):
        """
        The main hillclimber logic which performs the restart hill climber
        """

        grid_object = self.create_grid_object(self.protein)

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

                # update best current chain if the stability is better
                if grid_object.stability < best_stability:

                    # only get grid and chain for lower computing time
                    (
                        best_current_grid,
                        best_current_chain,
                        best_stability,
                        better_stab_found,
                    ) = self.copy_best(grid_object)

            self.print_message(1, best_stability)

            # if it didnt find any upgrades after the amount of random pulls of the whole chain, could be local minimum
            if not better_stab_found:
                # terminal info
                self.print_message(2, best_stability)

                # append to local minima list
                grid_object = self.save_local_minima(
                    best_stability, best_current_grid, best_current_chain, grid_object
                )

                # set this as best
                best_stability = copy.deepcopy(grid_object.stability)
                best_current_grid = copy.deepcopy(grid_object.grid)
                best_current_chain = copy.deepcopy(grid_object.grid_chain)

        return self.local_minimum_chains, self.stability_over_time
