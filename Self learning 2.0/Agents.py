import random
import numpy as np
import csv


class PassiveTDAgent:
    def __init__(self, discountFactor=1):
        self.utilityTable = {}
        self.policyStates = {}
        self.frequencyTable = {}
        self.previousState = ''
        self.previousAction = ''
        self.previousReward = 0
        self.dF = discountFactor
        self.trials = 0

    def PerceiveAndAct(self, actions):
        unchanged = True
        while unchanged:
            for policy in self.policyStates:
                actionUtilities = {}
                try:
                    actionUtilities['S'] = (self.utilityTable[self.previousState + 'S'])
                except:
                    actionUtilities['S'] = 0
                try:
                    actionUtilities['L'] = (self.utilityTable[self.previousState + 'L'])
                except:
                    actionUtilities['L'] = 0
                try:
                    actionUtilities['R'] = (self.utilityTable[self.previousState + 'R'])
                except:
                    actionUtilities['R'] = 0

                bestUtility = min(actionUtilities)

                if actionUtilities[bestUtility] < self.utilityTable[self.previousState]:
                    unchanged = False
        return policy

    def ProcessReward(self, currentReward, stateAfterAction):
        # Check if we have a new state
        if stateAfterAction not in self.utilityTable.keys():
            self.utilityTable[stateAfterAction] = currentReward
            self.frequencyTable[stateAfterAction] = 0

        if self.previousState != '':
            self.frequencyTable[self.previousState] += 1
            self.utilityTable[self.previousState] += self.StepSize() * (
                    self.previousReward + self.dF * self.utilityTable[stateAfterAction] -
                    self.utilityTable[self.previousState])

        self.previousState = stateAfterAction
        self.previousReward = currentReward

    def StepSize(self):
        return 60 / (59 + self.frequencyTable[self.previousState])

    def ActionPolicy(self, actions):
        self.previousAction = actions[random.randint(0, len(actions) - 1)]
        return self.previousAction

    def Terminate(self):
        self.trials += 1
        return self.trials > 100


class QLearningAgent:
    def __init__(self, discountFactor=1):
        self.QTable = {}
        self.frequencyTable = {}
        self.previousState = ''
        self.previousAction = ''
        self.previousReward = 0
        self.dF = discountFactor
        self.epsilon = 0.2  # Chance to explore randomly
        self.trials = 0
        self.deadEnds = 0

    def runBestSolution(self, protein, possibleActions):
        currentState = protein.getCompactState()
        QValues = {}
        for action in possibleActions:
            try:
                QValues[action] = self.QTable[currentState + action]
            except KeyError:
                print(f'Never visited state {currentState + action}!')
        try:
            return max(QValues, key=lambda key: QValues[key])
        except ValueError:
            return None

    def PerceiveAndAct(self, protein, currentReward, possibleActions):
        if currentReward == -1: self.deadEnds += 1
        currentState = protein.getCompactState()

        # Determine best next action
        bestAction = None
        bestNextQValue = 0
        nextQValues = {}

        # Shuffle actions to ensure random picks when Q values are the same
        random.shuffle(possibleActions)
        for action in possibleActions:
            try:
                nextQValues[action] = self.QTable[currentState + action]
            except KeyError:
                pass
        try:
            bestAction = max(nextQValues, key=lambda key: nextQValues[key])
            # print(f"Best action {bestAction} in state {currentState} \t: {nextQValues}")
            bestNextQValue = self.QTable[currentState + bestAction]
        except:
            bestNextQValue = 0
            bestAction = None
            if len(possibleActions) > 0: bestAction = random.choice(possibleActions)

        # [s,a] = previous state + previous action = current state
        # [s',a'] = current state + current possible action = next possible states

        # Terminal
        if len(possibleActions) == 0:
            self.QTable[currentState] = currentReward

        # Q-Learning
        if currentState != '':
            try:
                self.frequencyTable[currentState] += 1
            except KeyError:
                self.frequencyTable[currentState] = 1
            try:  # Q[s,a] = Q[s,a] + (a(s,a) * (r + (y * argmax Q[s', a'] - Q[s,a]))
                self.QTable[currentState] += (self.StepSize(currentState) *
                                              (self.previousReward + (self.dF * bestNextQValue) - self.QTable[
                                                  currentState]))
            except KeyError:  # Q[s,a] = 0 + (a(s,a) * (r + (y * argmax Q[s', a'] - 0))
                self.QTable[currentState] = (self.StepSize(currentState) *
                                             (self.previousReward + (self.dF * bestNextQValue)))

        self.previousState = currentState
        if len(possibleActions) > 0:
            self.previousAction = self.exploration(bestAction, possibleActions)
        else:
            self.previousAction = None
        self.previousReward = currentReward

        return self.previousAction

    def exploration(self, bestAction, possibleActions):
        if random.random() < self.epsilon:
            return random.choice(possibleActions)
        else:
            return bestAction

    def StepSize(self, state):
        return 60 / (59 + self.frequencyTable[state])

    def reset(self):
        self.previousState = ''
        self.previousAction = ''
        self.previousReward = 0

    def terminate(self):
        self.trials += 1
        if self.trials >= 5000:
            print(f"Trials: {self.trials}")
            print(f"Dead ends met: {self.deadEnds}")
            return True
        else:
            if self.trials % 250 == 0:
                print(f"Trial counter: {self.trials}\r", end='')
            return False

    def export_values(self, sequence):
        with open(f'q_table/{sequence}.csv', 'w') as f:
            for key in self.QTable.keys():
                f.write("%s,%s\n" % (key, self.QTable[key]))
        with open(f'freq/{sequence}.csv', 'w') as f:
            for key in self.frequencyTable.keys():
                f.write("%s,%s\n" % (key, self.frequencyTable[key]))
        return True

    def import_values(self, sequence):
        try:
            with open(f'q_table/{sequence}.csv', 'r') as infile:
                reader = csv.reader(infile)
                self.QTable = {rows[0]:float(rows[1]) for rows in reader}

            with open(f'freq/{sequence}.csv', 'r') as infile:
                reader = csv.reader(infile)
                self.frequencyTable = {rows[0]:int(rows[1]) for rows in reader}
            return True
        except OSError:
            return False
