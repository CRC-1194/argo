import os
import pandas as pd
from math import pi
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from math import ceil,exp

rcParams["text.usetex"] = True
rcParams["figure.dpi"] = 300 
rcParams["font.size"] = 9
rcParams["axes.titlesize"] = 9
rcParams["axes.labelsize"] = 9
global_markers = ['o', 'x', '^', 'v', 'd']

def read_data(pattern, csv_file):
    param_dirs = [param_dir for param_dir in os.listdir(os.curdir) if os.path.isdir(param_dir) and \
                  param_dir.startswith(pattern)] 

    param_dirs.sort()

    csv_files = [os.path.join(param_dir,csv_file) for param_dir in param_dirs \
                 if os.path.exists(os.path.join(param_dir,csv_file))]

    if (len(csv_files) == 0):
        raise ValueError("No CSV files found, did the parameter variation run?")

    data = pd.concat((pd.read_csv(csv_file) for csv_file in csv_files))

    return data

def n_cells_sorted(data):
    n_cells = list(set(data["N_CELLS"]))
    n_cells.sort()

    return n_cells

def compute_dependent_data(data, exact_volume, pattern, alg_name, data_write_dir):
    n_cells = n_cells_sorted(data)
    n_triangles = list(set(data["N_TRIANGLES"]))
    
    data["N_TRIANGLES_ROOT"] = data["N_TRIANGLES"] ** 0.5
    data["N_C"] = np.ceil(data["N_CELLS"]**(1./3.))
    data["VOLUME_ERROR_FROM_EXACT_VOLUME"] = (data["VOLUME_FROM_VOLUME_FRACTION"] - exact_volume).abs() / exact_volume
    data["CPU_TIME_SECONDS"] = data["CPU_TIME_MICROSECONDS"] / 1e06 

    data.to_csv(os.path.join(data_write_dir,"Data_Pandas_Dframe-%s-%s.csv" % (pattern, alg_name)), index=False)
    data.to_latex(os.path.join(data_write_dir, "Data_Pandas_Dframe-%s-%s.tex" % (pattern, alg_name)), index=False)

def plot_convergence(data, pattern, alg_name, data_write_dir):
    n_cells = n_cells_sorted(data)
    fig_conv, ax_conv = plt.subplots()
    for i,n_cell in enumerate(n_cells):
        n_c = ceil(n_cell**(1./3.))
        n_cell_data = data[data["N_CELLS"] == n_cell] 
        ax_conv.plot(n_cell_data["N_TRIANGLES_ROOT"], 
                     n_cell_data["VOLUME_ERROR_FROM_EXACT_VOLUME"], 
                     label = r"$N_c = %d$" % n_c, 
                     marker=global_markers[i % len(n_cells)])

    ax_conv.set_ylabel(r"$E_v$")
    ax_conv.set_xlabel(r"$\sqrt{N_T}$")
    ax_conv.loglog()

    box = ax_conv.get_position()
    ax_conv.legend(bbox_to_anchor=(1.0, 1.), fancybox=False, frameon=False, framealpha=0)

    # Save the figure.
    fig_conv.set_size_inches(4,3)
    fig_conv.savefig(os.path.join(data_write_dir, "Ev-%s-%s.pdf" % (pattern,alg_name)),
                     bbox_inches='tight')

def plot_cpu_times(data, pattern, alg_name, data_write_dir):
    n_cells =n_cells_sorted(data)
    fig_cpu, ax_cpu = plt.subplots()
    fig_cpu.set_size_inches(4,3)
    for i, n_cell in enumerate(n_cells):
        n_c = ceil(n_cell**(1./3.))
        n_cell_data = data[data["N_CELLS"] == n_cell] 
        n_cell_data = n_cell_data[n_cell_data["N_TRIANGLES_ROOT"] < 400]
        ax_cpu.plot(n_cell_data["N_TRIANGLES_ROOT"], 
                    n_cell_data["CPU_TIME_SECONDS"], 
                    label = r"$N_c = %d$" % n_c, 
                    marker=global_markers[i % len(n_cells)])
    ax_cpu.set_ylabel(r"CPU time in seconds")
    ax_cpu.set_xlabel(r"$\sqrt{N_T}$")

    ax_cpu.legend(bbox_to_anchor=[0.01,0.99], loc="upper left", framealpha=0)
    # Save the figure.
    fig_cpu.savefig(os.path.join(data_write_dir, "CPUtime-%s-%s.pdf" % (pattern,alg_name)),
                    bbox_inches='tight')

def plot_smca_refinement_convergence(pattern, exact_volume, data_write_dir="", csv_file="smcaVofInit.csv"):
    data = read_data(pattern, csv_file)
    compute_dependent_data(data, exact_volume, pattern, "SMCA", data_write_dir)
    fig, axis = plt.subplots()
    fig.set_size_inches(4,3)
    axis.grid()
    axis.semilogy(data["MAX_REFINEMENT_LEVEL"], data["VOLUME_ERROR_FROM_EXACT_VOLUME"],
              linewidth=0, marker='x', color='k')
    axis.set_xlabel(r"$l_\mathrm{max}$")
    axis.set_ylabel(r"$E_v$")
    fig.savefig(os.path.join(data_write_dir, "refinement-convergence-%s-SCMA.pdf" % pattern),
                    bbox_inches='tight')

def plot_study(pattern, alg_name, exact_volume, data_write_dir="", csv_file="surfaceCellVofInit.csv"):
    data = read_data(pattern, csv_file)
    compute_dependent_data(data, exact_volume, pattern, alg_name, data_write_dir)
    plot_convergence(data, pattern, alg_name, data_write_dir)
    plot_cpu_times(data, pattern, alg_name, data_write_dir)

    return data
