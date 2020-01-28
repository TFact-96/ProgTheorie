import numpy as np  # Protein folding calculations


class Amino:
    def __init__(self, amino, index):
        self.index = index
        self.type = amino
        self.colour = self.Colour()
        self.jointAngle = ''  # Staight (S), clock wise (R) and counter clockwise (L)
        self.pos = [index, 0]

    def Colour(self):
        """
        Sets colour for plotting.
        :return: colour of amino's
        """
        if self.type == 'H':
            return 'red'
        elif self.type == 'C':
            return 'yellow'
        elif self.type == 'P':
            return 'blue'
        else:
            return 'grey'


class Protein:  # The state of the protein
    def __init__(self, proteinSequence):
        self.state = []
        for index, aminoAcidType in enumerate(proteinSequence):
            self.state.append(Amino(aminoAcidType, index))
        self.stability = 0
        rows, cols = (len(self.state) + 1, len(self.state) + 1)
        self.stateField = [[0 for i in range(-cols, cols)] for j in range(-rows, rows)]
        self.stateField[0][0] = self.state[0].type
        if len(self.state) >= 3:
            self.actions = ['S', 'L', 'R']
        else:
            self.actions = []
        self.nextAminoFold = 0
        self.direction = 'E'
        self.rewardAction = {'S': 0, 'L': 0, 'R': 0}  # The reward for each future action
        self.deadEnds = 0

    def getStability(self):
        return self.stability

    def getCompactState(self):
        """
        Returns current sequence of protein.
        :return: current sequence of protein
        """
        compactState = ''
        for amino in self.state:
            compactState += amino.jointAngle
        return compactState

    def CalculateReward(self, p_x, p_y, direction):
        """
        Calculate reward of current state.
        :param p_x: last amino x position
        :param p_y: last amino y position
        :param direction: direction of movement
        :return: reward
        """
        reward = 0
        nextAminioType = self.state[self.nextAminoFold + 1].type
        # Checking around amino acid of action for new bonds
        if nextAminioType == 'H':
            if self.stateField[p_x][p_y - 1] == 'H' and not direction == 'North':
                reward += 1
            if self.stateField[p_x][p_y - 1] == 'C' and not direction == 'North':
                reward += 1
            if self.stateField[p_x][p_y + 1] == 'H' and not direction == 'South':
                reward += 1
            if self.stateField[p_x][p_y + 1] == 'C' and not direction == 'South':
                reward += 1
            if self.stateField[p_x - 1][p_y] == 'H' and not direction == 'East':
                reward += 1
            if self.stateField[p_x - 1][p_y] == 'C' and not direction == 'East':
                reward += 1
            if self.stateField[p_x + 1][p_y] == 'H' and not direction == 'West':
                reward += 1
            if self.stateField[p_x + 1][p_y] == 'C' and not direction == 'West':
                reward += 1
        elif nextAminioType == 'C':
            if self.stateField[p_x][p_y - 1] == 'C' and not direction == 'North':
                reward += 5
            if self.stateField[p_x][p_y - 1] == 'H' and not direction == 'North':
                reward += 1
            if self.stateField[p_x][p_y + 1] == 'C' and not direction == 'South':
                reward += 5
            if self.stateField[p_x][p_y + 1] == 'H' and not direction == 'South':
                reward += 1
            if self.stateField[p_x - 1][p_y] == 'C' and not direction == 'East':
                reward += 5
            if self.stateField[p_x - 1][p_y] == 'H' and not direction == 'East':
                reward += 1
            if self.stateField[p_x + 1][p_y] == 'C' and not direction == 'West':
                reward += 5
            if self.stateField[p_x + 1][p_y] == 'H' and not direction == 'West':
                reward += 1

        return reward

    def CheckFreeSpot(self, x, y):
        """
        Checks if position next to amino is free.
        :param x: x coordinate of amino
        :param y: y coordinate of amino
        :return:
        """
        return self.stateField[x][y] == 0

    def ProcessAction(self, action):
        """
        Updates current amino chain with new action given by the agent.
        :param action: action given by agent
        :return: reward in new state
        """
        # Set fold of next amino acid
        self.nextAminoFold += 1
        i = self.nextAminoFold
#        print(f"Processing action {action} at amino {self.nextAminoFold}")

        # Place amino acid according to current facing direction or last action
        if self.direction == 'E':  # East
            self.state[i].pos[0] = self.state[i - 1].pos[0] + 1
            self.state[i].pos[1] = self.state[i - 1].pos[1]
        elif self.direction == 'W':  # West
            self.state[i].pos[0] = self.state[i - 1].pos[0] - 1
            self.state[i].pos[1] = self.state[i - 1].pos[1]
        elif self.direction == 'N':  # North
            self.state[i].pos[0] = self.state[i - 1].pos[0]
            self.state[i].pos[1] = self.state[i - 1].pos[1] + 1
        elif self.direction == 'S':  # South
            self.state[i].pos[0] = self.state[i - 1].pos[0]
            self.state[i].pos[1] = self.state[i - 1].pos[1] - 1
        else: print(f"Direction {self.direction} doesn't exist")

        # Add amino acid type to fieldState
        p_x = self.state[i].pos[0]
        p_y = self.state[i].pos[1]
        if self.CheckFreeSpot(p_x, p_y):
            self.stateField[self.state[i].pos[0]][self.state[i].pos[1]] = self.state[i].type
        else:
            print(f'Overlap at {i}! Should not happen!')

        self.state[i].jointAngle = action
        # No more actions possible after we reached the one but last amino acid
        self.actions = []
        reward = 0
        if self.nextAminoFold >= len(self.state) - 1:
            print("No more aminos")
        else:
            # Direction and possible reward for next amino acid
            if action == 'L':
                if self.direction == 'N': self.direction = 'W'
                elif self.direction == 'E': self.direction = 'N'
                elif self.direction == 'S': self.direction = 'E'
                elif self.direction == 'W': self.direction = 'S'
            if action == 'R':
                if self.direction == 'N': self.direction = 'E'
                elif self.direction == 'E': self.direction = 'S'
                elif self.direction == 'S': self.direction = 'W'
                elif self.direction == 'W': self.direction = 'N'

            if self.direction == 'E': # Next east
                reward = self.CalculateReward(p_x + 1, p_y, 'East')
            elif self.direction == 'W':  # Next west
                reward = self.CalculateReward(p_x - 1, p_y, 'West')
            elif self.direction == 'N':  # Next north
                reward = self.CalculateReward(p_x, p_y + 1, 'North')
            elif self.direction == 'S': # Next south
                reward = self.CalculateReward(p_x, p_y - 1, 'South')

            justPlacedAminoAcid = self.nextAminoFold
            if justPlacedAminoAcid < len(self.state) - 2:
                for action in self.rewardAction:
                    action = 1  # If this action is an overlap, you'll see it being 1,
                    # all will be overwritten by CalculateReward, except overlaps

                if self.direction == 'E':  # Next east
                    p_x += 1
                    if self.CheckFreeSpot(p_x + 1, p_y):
                        self.actions.append('S')
                    if self.CheckFreeSpot(p_x, p_y + 1):
                        self.actions.append('L')
                    if self.CheckFreeSpot(p_x, p_y - 1):
                        self.actions.append('R')
                elif self.direction == 'W':  # Next west
                    p_x -= 1
                    if self.CheckFreeSpot(p_x - 1, p_y):
                        self.actions.append('S')
                    if self.CheckFreeSpot(p_x, p_y - 1):
                        self.actions.append('L')
                    if self.CheckFreeSpot(p_x, p_y + 1):
                        self.actions.append('R')
                elif self.direction == 'N':  # Next north
                    p_y += 1
                    if self.CheckFreeSpot(p_x, p_y + 1):
                        self.actions.append('S')
                    if self.CheckFreeSpot(p_x + 1, p_y):
                        self.actions.append('R')
                    if self.CheckFreeSpot(p_x - 1, p_y):
                        self.actions.append('L')
                elif self.direction == 'S': # Next south
                    p_y -= 1
                    if self.CheckFreeSpot(p_x, p_y - 1):
                        self.actions.append('S')
                    if self.CheckFreeSpot(p_x - 1, p_y):
                        self.actions.append('R')
                    if self.CheckFreeSpot(p_x + 1, p_y):
                        self.actions.append('L')
                if len(self.actions) < 1:
                    # Dead end
                    reward = -1
#                    print(f'Dead end state for amino {i}\t(state: { {self.getCompactState()}})')
            else:
                # Place last amino acid according to current facing direction or last action
                i += 1
                if self.direction == 'E':  # East
                    self.state[i].pos[0] = self.state[i - 1].pos[0] + 1
                    self.state[i].pos[1] = self.state[i - 1].pos[1]
                elif self.direction == 'W':  # West
                    self.state[i].pos[0] = self.state[i - 1].pos[0] - 1
                    self.state[i].pos[1] = self.state[i - 1].pos[1]
                elif self.direction == 'N':  # North
                    self.state[i].pos[0] = self.state[i - 1].pos[0]
                    self.state[i].pos[1] = self.state[i - 1].pos[1] + 1
                elif self.direction == 'S': # South
                    self.state[i].pos[0] = self.state[i - 1].pos[0]
                    self.state[i].pos[1] = self.state[i - 1].pos[1] - 1
        self.stability += reward
        return reward
