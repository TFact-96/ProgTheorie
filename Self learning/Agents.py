import random
import numpy as np


class PassiveTDAgent:
    def __init__(self, discountFactor=1):
        self.utilityTable = {}
        self.policyStates = {}
        self.frequencyTable = {}
        self.previousState = ''
        self.previousAction = ''
        self.previousReward = 0
        self.discountFactor = discountFactor
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
                    self.previousReward + self.discountFactor * self.utilityTable[stateAfterAction] -
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
        self.discountFactor = discountFactor
        self.stepReward = -0.25
        self.trials = 0
        self.Ne = 10
        self.exp = 0

    def PerceiveAndAct(self, protein, currentReward, moves):
        highest = {}
        PossibleActions = moves
        stateAfterAction = protein.getCompactState()
        episode_reward = currentReward

        # Check if we have a new state
        if stateAfterAction not in self.QTable.keys():
            self.QTable[stateAfterAction] = 0
            self.frequencyTable[stateAfterAction] = 0

        if len(moves) == 0:
            self.frequencyTable[stateAfterAction] += 1
            self.QTable[stateAfterAction] = episode_reward
            return False

        if self.previousState != '':
            self.frequencyTable[stateAfterAction] += 1
            self.QTable[self.previousState] += self.StepSize() * (self.frequencyTable[self.previousState]) * (
                    episode_reward + self.discountFactor * self.QTable[stateAfterAction] -
                    self.QTable[self.previousState])

        for choice in PossibleActions:
            try:
                q_val = self.QTable[stateAfterAction + choice]
                new_state = stateAfterAction + choice
            except KeyError as e:
                new_state = str(e)
                q_val = 0
            try:
                freq = self.frequencyTable[stateAfterAction + choice]
            except KeyError:
                freq = 0
                
            explor = self.exploration(new_state, q_val, freq)
            highest[choice] = self.exp + explor

        action = min(highest, key=highest.get)

        self.previousState = stateAfterAction
        self.previousReward = currentReward

        return action

    def exploration(self, new_state, q, freq):
        if freq < self.Ne:
            if len(new_state) < 5:
                return random.random()
            elif len(new_state) == 5:
                return 1
            else:
                return 2
        else:
            return q

    def epsilon_greedy_policy(self, Q, state, actions, epsilon=0.05):
        '''
        Create a policy in which epsilon dictates how likely it will
        take a random action.

        :param Q: links state -> action value (dictionary)
        :param state: state character is in (int)
        :param actions: actions that can be taken (list)
        :param epsilon: chance it will take a random move (float)
        :return: probability of each action to be taken (list)
        '''
        probs = np.ones(len(actions)) * epsilon / len(actions)
        best_action = np.argmax(Q[state])
        probs[best_action] += 1.0 - epsilon

        return probs

    def reset(self):
        self.previousState = ''
        self.previousAction = ''
        self.previousReward = 0
        self.exp = 0

    def StepSize(self):
        return 60 / (59 + self.frequencyTable[self.previousState])

    def terminate(self):
        self.trials += 1
        return self.trials > 500
