# Import public code
import PySimpleGUI as sg  # The GUI
import matplotlib.pyplot as plt  # Visualizing proteins
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import traceback  # For debugging errors
import numpy as np

# Import personal code
import Environment
import Agents


def showProtein(protein, stability):
    fig, ax = plt.subplots(figsize=(7, 7))

    # Turn grid on for both major and minor ticks and style minor slightly
    # differently.
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')

    x = []
    y = []
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
    proteinLength = len(protein)
    ax.set_xlim(-proteinLength, proteinLength)
    ax.set_ylim(-proteinLength, proteinLength)

    # Turn grid on for both major and minor ticks and style minor slightly differently.
    ax.grid(which='major', color='#CCCCCC', linestyle='--')

    # Change major ticks to show every 20.
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))

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
reward = 0 #initial value
moves = ['S', 'L', 'R']
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
        elif event in ('Run PTDA'):
            proteinSequence = window.FindElement('Protein Sequence').Get()
            protein = Environment.Protein(proteinSequence)
            # Initialize utilityTable
            agent = Agents.QLearningAgent()
            showProtein(protein.state, protein.getStability())

            while not agent.terminate():
                while True:
                    action = agent.PerceiveAndAct(protein, reward, moves)
                    reward, moves = protein.ProcessAction(action)
                    if len(moves) == 0:
                        action = agent.PerceiveAndAct(protein, reward, moves)
                        break

                # Reset environment and agent
                protein = Environment.Protein(proteinSequence)
                agent.reset()
                moves = ['S', 'L', 'R']
            print(f"QTable: {agent.QTable}")
            print(f"Frequency: {agent.frequencyTable}")
        elif event in ('Run PTDA one trial'):
            while len(protein.actions) > 0:
                action = agent.PerceiveAndAct(protein.actions)
                reward = protein.ProcessAction(action)
                agent.ProcessReward(reward, protein.getCompactState())
            print(f"States visited: {agent.utilityTable}")
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
