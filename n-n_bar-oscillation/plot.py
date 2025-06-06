import numpy as np
import matplotlib.pyplot as plt

def plot_wilson_coefficients(list_of_energy_solutions, list_of_wc, labels=None, title=None, save_path=None):
    colors = ["blue", "green", "red", "orange", "black", "gray"]
    line_styles = ["solid", "solid", "solid", "solid"]
    fig, axs = plt.subplots(3, 3, figsize=(12, 10))
    axs = axs.flatten()

    # Iterate over all datasets
    for idx, (energy_solutions, wilson_coefficients) in enumerate(zip(list_of_energy_solutions, list_of_wc)):

        for i in range(9):
            flag = True
            for e, wc in zip(energy_solutions, wilson_coefficients):
                if flag:
                    label = labels[idx] if labels else None
                    flag = False
                else:
                    label = None
                axs[i].plot(e, wc[i], label=label, color=colors[idx], linestyle=line_styles[idx % len(line_styles)])

    for i in range(9):
        axs[i].set_xlabel("E (TeV)")
        axs[i].set_ylabel(f"$C^{(i+1)}$")
        axs[i].grid(True)
        #axs[i].set_yscale('log')
        if labels:
            axs[i].legend()

    if title:
        plt.suptitle(title)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    plt.show()

def plot_Lambda_bar_graph(Lambda, title=None, save_path=None):
    categories = [f"$C^{{({i+1})}}$" for i in range(9)]

    # Bar values (replace with your real data)
    with_rge = [Lambda[-1][i][-1]*1e-12 for i in range(9)]
    without_rge = [Lambda[0][i][0]*1e-12 for i in range(9)]

    print(title)
    for i in range(9):
        print(f"C{i+1}: with RGE = {with_rge[i]} TeV, without RGE = {without_rge[i]} TeV")
    # X locations
    x = np.arange(len(categories))  # [0, 1, 2, ..., 8]
    width = 0.35  # Width of the bars

    # Plot side-by-side bars
    plt.bar(x - width/2, with_rge, width, label="with RGE")
    plt.bar(x + width/2, without_rge, width, label="without RGE", alpha=0.5, color="blue")

    # Log scale on y-axis
    #plt.yscale('log')

    # Labels and formatting
    plt.xticks(x, categories)
    plt.ylabel(f"$\Lambda$ (TeV)")
    if title:
        plt.title(title)
    plt.legend(loc="lower right")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)

    plt.show()
