import matplotlib.pyplot as plt

def data_plot_hillclimb(stability_over_time, amino):
    min_stability = min(stability_over_time)    
    plt.plot(stability_over_time, label=f"Min stability: {min_stability}")
    plt.xlabel("Chainlength pulls")
    plt.ylabel("Stability")
    plt.title(f"Restart Hill Climb Algorithm, Amino: {amino}")
    plt.legend()
    plt.show()
    return