import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)


def showProtein(protein, stability, progress=None):
    if progress is None:
        progress = []
    fig, ax = plt.subplots(figsize=(14, 7))

    # Turn grid on for both major and minor ticks and style minor slightly
    # differently.
    ax.grid(which='major', color='#CCCCCC', linestyle='--')
    ax.grid(which='minor', color='#CCCCCC', linestyle=':')

    x, y = [], []

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
    plt.ylabel('y')
    plt.xlabel('x')

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

    plt.plot(list(range(1, len(progress) + 1)), progress, "-", linewidth=1, color='black', alpha=0.5,
             zorder=1)

    plt.show()
