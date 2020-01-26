import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D

def Plot3D(x, y, z, hh_bonds, ch_bonds, cc_bonds, stability, color):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # plot the bond lines
    for bond in hh_bonds:
        ax.plot3D(
            bond[0], bond[1], bond[2],
            "--",
            markersize=1,
            color='red',
            zorder=-1
        )

    for bond in ch_bonds:
        ax.plot3D(
            bond[0], bond[1], bond[2],
            "--",
            markersize=1,
            color='orange',
            zorder=-1
        )

    for bond in cc_bonds:
        ax.plot3D(
            bond[0], bond[1], bond[2],
            "--",
            markersize=1,
            color='yellow',
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
    ax.set_title(f"Protein chain (Stability: {stability})")
    plt.show()                
