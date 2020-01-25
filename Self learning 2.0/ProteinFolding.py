# Import public code
import PySimpleGUI as sg  # The GUI
import matplotlib.pyplot as plt  # Visualizing proteins
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import traceback  # For debugging errors
import numpy as np

# Import personal code
import Environment
import Agents


def showProtein(protein, stability, progress=None):
    if progress is None:
        progress = []
    fig, ax = plt.subplots(figsize=(14, 7))

    # Turn grid on for both major and minor ticks and style minor slightly
    # differently.
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')

    x = []
    y = []

    if len(progress) > 0: plt.subplot(1, 2, 1)
    for amino in protein:
        x.append(amino.pos[0])
        y.append(amino.pos[1])
        # Draw each amino acid
        plt.scatter(amino.pos[0], amino.pos[1], color=amino.colour, zorder=2)
    plt.plot(x, y, "-", linewidth=3, color='black', zorder=1)

    # Plot the chain line between amino acids
    plt.title(f"Amino acid chain (Stability: {stability})")
    plt.grid(True)
    plt.ylabel('y-as')
    plt.xlabel('x-as')

    # Set axis ranges; by default this will put major ticks every 25.
    #    proteinLength = len(protein)
    #    ax.set_xlim(-proteinLength, proteinLength)
    #    ax.set_ylim(-proteinLength, proteinLength)

    # Turn grid on for both major and minor ticks and style minor slightly differently.
    ax.grid(which='major', color='#CCCCCC', linestyle='--')

    # Change major ticks to show every 20.
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))

    if len(progress) > 0: plt.subplot(1, 2, 2)

    plt.plot(list(range(1, len(stabilityProgress) + 1)), stabilityProgress, "-", linewidth=1, color='black', alpha=0.5,
             zorder=1)

    plt.show(block=False)


################ GUI layout ################################################################
sg.change_look_and_feel('Default1')  # Add a touch of color
frame_passive_td_agent = [
    [
        sg.Button('Run until converged', key='Run PTDA'), sg.Button('Run next trial', key='Run PTDA one trial'),
        sg.Button('Show Result', key='Show')
    ], [
        sg.Button('Reset', key='Reset PTDA')
    ]
]

layout = [
    [
        sg.Text('Protein sequence:'), sg.InputText(key='Protein Sequence'), sg.OK(button_text='OK/RESET', key='OK')
    ], [
        sg.Frame('Passive-TD-Agent', frame_passive_td_agent)
    ], [
        sg.Canvas(key='Canvas')
    ], [
        sg.Button('Exit')
    ]
]

################ Initializing ################################################################
window = sg.Window('Protein Folding', layout)  # Create window
protein = None
proteinSequence = ''
agent = None
reward = 0
stabilityProgress = []

################ Main ################################################################
try:
    while True:
        event, values = window.read()
        if event in (None, 'Exit'):  # if user closes window or clicks cancel
            break
        elif event in ('OK'):
            proteinSequence = window.FindElement('Protein Sequence').Get()
            protein = Environment.Protein(proteinSequence)
            # Initialize utilityTable
            agent = Agents.QLearningAgent()
            showProtein(protein.state, protein.getStability())
        elif event in 'Run PTDA':
            proteinSequence = window.FindElement('Protein Sequence').Get()
            # Create environment
            protein = Environment.Protein(proteinSequence)
            # Create agent
            agent = Agents.QLearningAgent()
            if agent.import_q_values(proteinSequence):
                print("Imported Q values")
            else:
                print("New Sequence")
            #            showProtein(protein.state, protein.getStability())
            while not agent.terminate():
                while True:
                    action = agent.PerceiveAndAct(protein, reward, protein.actions)
                    if action is None: break
                    reward = protein.ProcessAction(action)
                # Create new environment and reset agent
                stabilityProgress.append(protein.getStability())
                protein = Environment.Protein(proteinSequence)
                agent.reset()
            while True:
                action = agent.runBestSolution(protein, protein.actions)
                if action is None:
                    break
                reward = protein.ProcessAction(action)
                agent.export_q_values(proteinSequence)
            showProtein(protein.state, protein.getStability(), stabilityProgress)

        elif event in ('Run PTDA one trial'):
            proteinSequence = window.FindElement('Protein Sequence').Get()
            # Create environment and reset agent
            protein = Environment.Protein(proteinSequence)
            agent.reset()
            i = 0
            while True:
                i += 1
                action = agent.PerceiveAndAct(protein, reward, protein.actions)
                print(f"{i}. Actions for state {protein.getCompactState()}: {protein.actions} => {action}")
                if action == None: break
                reward = protein.ProcessAction(action)
            print(f"States visited: {agent.QTable}")
            print(f"Frequency: {agent.frequencyTable}")
            showProtein(protein.state, protein.getStability())
            # Reset environment
            protein = Environment.Protein(proteinSequence)
        elif event in ('Reset PTDA'):
            agent = Agents.PassiveTDAgent()
        elif event in ('Show'):
            best_movement = max(agent.QTable, key=agent.QTable.get)
            print(f"The best sequence is: {best_movement}")
except Exception:
    print(traceback.format_exc())
window.close()
