import matplotlib.pyplot as plt

# Plotting stability over time for the Restart Hillclimb
def data_plot_hillclimb(stability_over_time, amino):
    min_stability = min(stability_over_time)    
    plt.plot(stability_over_time, label=f"Min stability: {min_stability}")
    plt.xlabel("Iteration")
    plt.ylabel("Stability")
    plt.title(f"Restart Hill Climb Algorithm, Amino: {amino}")
    plt.legend()
    plt.show()
    return
    
# Plotting stability over time for Simulated Annealing
def data_plot_annealing(stability_over_time, amino):
    min_stability = min(stability_over_time)    
    plt.plot(stability_over_time, label=f"Min stability: {min_stability}")
    plt.xlabel("Nodepull iteration")
    plt.ylabel("Stability")
    plt.title(f"Simulated annealing, Amino: {amino}")
    plt.legend()
    plt.show()
    return