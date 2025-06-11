import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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

def plot_correlation_matrix(rho_matrix, title = None, save_path=None):
    mask = np.abs(rho_matrix) == 0.0
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(
        rho_matrix,
        cmap='coolwarm',
        mask=mask,
        xticklabels=[f'$C^{{({i+1})}}$' for i in range(9)],
        yticklabels=[f'$C^{{({i+1})}}$' for i in range(9)],
        annot=True,
        #fmt=".3f",   
        linewidths=0.5,  # Thickness of the lines
        linecolor='black',
        ax = ax
    )

    n = rho_matrix.shape[0]
    ax.add_patch(plt.Rectangle((0, 0), n, n, fill=False, edgecolor='black', lw=0.5, clip_on=False))

    plt.tight_layout()
    if title:
        plt.title(title)
    if save_path:
        plt.savefig(save_path)
    plt.show()