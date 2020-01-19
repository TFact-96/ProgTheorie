import random


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
