"""
Some specific utility functions, such as plotting.

"""

import os

import matplotlib.pyplot as plt
import numpy as np

import shapefile

import matplotlib.ticker as tick
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter

from matplotlib.tri import Triangulation

import meshio

import geopandas as gpd

import vtk

import pyvista as pv

import json

from pathlib import Path

plt.rc('text', usetex=True)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def read_data_from_vtk(vtk_file_name):

    # read data from vtk file
    vtk_result = meshio.read(vtk_file_name)

    # Extract points and cells (triangles)
    points = vtk_result.points[:, :2]  # Take only the X and Y coordinates
    cells = vtk_result.cells_dict.get("triangle")

    if cells is None:
        raise ValueError("The VTK file does not contain triangular cells.")

    # number of triangles
    nTri = vtk_result.cells[0].data.shape[0]

    tris = np.zeros((nTri, 3))
    for i in range(nTri):
        tris[i, 0] = vtk_result.cells[0].data[i, 0]
        tris[i, 1] = vtk_result.cells[0].data[i, 1]
        tris[i, 2] = vtk_result.cells[0].data[i, 2]

    xCoord = vtk_result.points[:, 0]
    yCoord = vtk_result.points[:, 1]

    xl = xCoord.min()
    xh = xCoord.max()
    yl = yCoord.min()
    yh = yCoord.max()

    # get the data
    WSE = np.squeeze(vtk_result.point_data['WSE'])
    Vel = np.squeeze(vtk_result.point_data['U'])

    print("Vel: ", Vel.shape)

    #Vel is a vector. Get its components U and V
    U = Vel[:,0]
    V = Vel[:,1]

    return points, cells, tris, xCoord, yCoord, xl, xh, yl, yh, WSE, U, V

def plot_one_subplot(iRow, iCol, sensitivity_name, fig, axs, xCoord, yCoord, tris, sensitivity, local_levels,
                     min_value, max_value, xl, xh, yl, yh, bShow_x_title=False, bShow_y_title=False):
    #just plot a specified subplot

    cf_zb = axs.tricontourf(xCoord, yCoord, tris, sensitivity, local_levels,
                                  vmin=min_value, vmax=max_value, cmap=plt.cm.jet, extend='neither')
    axs.set_xlim([xl, xh])
    axs.set_ylim([yl, yh])
    #axs.set_aspect('equal')

    if iRow == 0 and iCol == 0:
        axs.set_title('Truth', fontsize=14)
    elif iRow == 0 and iCol == 1:
        axs.set_title('Inversion', fontsize=14)
    elif iRow == 0 and iCol == 2:
        axs.set_title('Difference', fontsize=14)

    axs.text(600, 340, sensitivity_name, ha='center', va='bottom', fontsize=14)

    if bShow_x_title:
        axs.set_xlabel('$x$ (m)', fontsize=16)
        axs.tick_params(axis='x', labelsize=14)
    else:
        axs.tick_params(axis='x',length=0, labelbottom=False, labelleft=False)

    if bShow_y_title:
        axs.set_ylabel('$y$ (m)', fontsize=16)
        axs.tick_params(axis='y', labelsize=14)
    else:
        axs.tick_params(axis='y', length=0, labelbottom=False, labelleft=False)

    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="3%", pad=0.2)
    clb = fig.colorbar(cf_zb, ticks=np.linspace(min_value, max_value, 7), cax=cax)
    if iCol==0 or iCol==1:
        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    else:
        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1e'))

    clb.ax.tick_params(labelsize=12)

    if iRow==0: #WSE
        clb.ax.set_title("(m)", loc='center', fontsize=12)
    else: #U or V
        clb.ax.set_title("(m/s)", loc='center', fontsize=12)

def plot_with_pyvista():
    # Load the VTK file
    file_path = 'forward_simulation_results_0101.vtk'
    grid = pv.read(file_path)

    # Check available scalar data
    print("Available data arrays:", grid.array_names)

    # Extract and visualize the 'zb_cell' data
    if 'zb_cell' in grid.array_names:
        # Add 'zb_cell' data to the active scalars
        grid.set_active_scalars('zb_cell')

        # Customize the plot
        plotter = pv.Plotter(off_screen=True)  # Use off_screen=True for saving without showing
        plotter.add_mesh(
            grid,
            cmap='viridis',  # Colormap for the contour
            scalar_bar_args={
                'title': 'zb_cell',  # Title of the scalar bar
                'title_font_size': 10,
                'label_font_size': 8,
            },
            show_edges=False,  # Disable mesh edges for cleaner look
        )
        plotter.view_xy()  # View along the XY plane

        # Add labels and title
        plotter.add_text(
            'Contour Plot of zb_cell',
            font_size=12,
            position='upper_left'
        )

        # Save the figure as a high-quality image
        save_path = 'zb_cell_contour.png'
        plotter.screenshot(save_path, window_size=[1920, 1080])  # Specify resolution
        print(f"Figure saved to {save_path}")
    else:
        print("The 'zb_cell' data was not found in the VTK file.")

def triangulate_vtk_files(vtk_file_names):
    #triangulate the vtk

    for vtk_file_name in vtk_file_names:
        print(f"Triangulating {vtk_file_name}")

        file_name_without_extension = os.path.splitext(vtk_file_name)[0]

        # Load the VTK file
        grid = pv.read(vtk_file_name)

        # Triangulate the grid
        triangulated = grid.triangulate()

        # Interpolate cell data to point data
        interpolated = triangulated.cell_data_to_point_data()

        # Save the triangulated grid to a new VTK file
        triangulated_file_path = file_name_without_extension + '_tri.vtk'
        interpolated.save(triangulated_file_path)

        print(f"Triangulated file saved to: {triangulated_file_path}")

        print("Available point data:", interpolated.point_data.keys())

def plot_loss_history(nIterations):
    """
    plot the inversion loss history

    """

    #read the loss data from the "inversion_results_losses_and_ODE_solution.json" file
    with open('inversion_results_losses_and_ODE_solution.json', 'r') as file:
        data = json.load(file)

    loss_preds = data['loss_preds']
    loss_pred_WSEs = data['loss_pred_WSEs']
    loss_bounds = data['loss_bounds']
    loss_pred_uvs = data['loss_pred_uvs']

    iter_numbers = np.arange(len(loss_preds))

    # plot losses
    plt.plot(iter_numbers, loss_preds, 'k', label=r'$l_{pred, total}$')
    plt.plot(iter_numbers, loss_pred_WSEs, 'r--', label=r'$l_{pred, WSE}$')
    plt.plot(iter_numbers, loss_pred_uvs, 'b-.', label=r'$l_{pred, uv}$')
    #plt.plot(iter_numbers, loss_bounds, 'r--', label=r'$l_{bound}$')

    ax = plt.gca()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    #ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    plt.tick_params(axis='both', which='major', labelsize=14)

    plt.xlabel('Iteration', fontsize=16)
    plt.ylabel('Loss', fontsize=16)
    #plt.xlim([0, 100])
    # plt.ylim([1e-5, 1e-2])
    # plt.ylim([-0.0001,0.0015])
    plt.yscale('log')
    # plt.ylim([0, 0.01])
    plt.legend(loc='upper right', fontsize=14, frameon=False)
    plt.savefig("inversion_losses.png", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()


def plot_n_history_trajectories(nIterations):
    """
    Plot the n history and trajectories in one figure

    Parameters
    ----------
    nIterations

    Returns
    -------

    """

    n_truth_values = [0.03, 0.04, 0.05, 0.03, 0.045, 0.05]

    n_0 = []
    n_1 = []
    n_2 = []
    n_3 = []
    n_4 = []
    n_5 = []

    iter_numbers = np.arange(nIterations - 1)

    # loop over all inversion_callback_save_iter_*.json files to read the n value histories
    for iter in range(1, nIterations):
        print("reading n values from iteration ", iter)
        file_name = "inversion_callback_save_iter_" + str(iter) + ".json"
        with open(file_name, 'r') as file:
            data = json.load(file)
            n_0.append(data['theta'][0])
            n_1.append(data['theta'][1])
            n_2.append(data['theta'][2])
            n_3.append(data['theta'][3])
            n_4.append(data['theta'][4])
            n_5.append(data['theta'][5])

    fig, axs = plt.subplots(1, 2, figsize=(16, 6), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    # fig.subplots_adjust(hspace=.15, wspace=.01)

    # plot n histories
    axs[0].plot(iter_numbers, n_1, 'k', label=r'$n_1$')
    axs[0].plot(iter_numbers, n_2, 'r--', label=r'$n_2$')
    axs[0].plot(iter_numbers, n_3, 'b-.', label=r'$n_3$')
    axs[0].plot(iter_numbers, n_4, 'g:', label=r'$n_4$')
    axs[0].plot(iter_numbers[::4], n_5[::4], 'c-o', markersize=2, label=r'$n_5$',
             markerfacecoloralt='c', markerfacecolor='none', markeredgecolor='c', alpha=0.5)

    # plot truth values
    x_truth_text = 200
    axs[0].plot([x_truth_text, nIterations + 10], [n_truth_values[1], n_truth_values[1]], 'k', alpha=0.5)
    axs[0].plot([x_truth_text, nIterations + 10], [n_truth_values[2], n_truth_values[2]], 'r', alpha=0.5)
    axs[0].plot([x_truth_text, nIterations + 10], [n_truth_values[3], n_truth_values[3]], 'b', alpha=0.5)
    axs[0].plot([x_truth_text, nIterations + 10], [n_truth_values[4], n_truth_values[4]], 'g', alpha=0.5)
    axs[0].plot([x_truth_text, nIterations + 10], [n_truth_values[5], n_truth_values[5]], 'c', alpha=0.5)

    # add truth text
    axs[0].text(x_truth_text, n_truth_values[1] + 0.0005, r'$n_1$ truth', fontsize=14)
    axs[0].text(x_truth_text, n_truth_values[2] + 0.0005, r'$n_2$ and $n_5$ truth', fontsize=14)
    axs[0].text(x_truth_text, n_truth_values[3] + 0.0005, r'$n_3$ truth', fontsize=14)
    axs[0].text(x_truth_text, n_truth_values[4] + 0.0005, r'$n_4$ truth', fontsize=14)

    axs[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    axs[0].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    axs[0].tick_params(axis='both', which='major', labelsize=14)

    axs[0].set_xlabel('Iteration', fontsize=16)
    axs[0].set_ylabel('Manning\'s $n$ values', fontsize=16)
    # axs[0].set_xlim([0, 100])
    # axs[0].set_ylim([0, 0.01])
    # axs[0].set_yscale('log')
    legend = axs[0].legend(loc='lower right', fontsize=14, frameon=True, facecolor='white')

    #plot n trajectories
    # Create color array from indices
    colors = np.arange(len(n_1))

    print("len(n_1)", len(n_1))

    # Plot line
    plot_trajectory = axs[1].plot(n_1, n_5, 'b-', alpha=0.25, label='Inversion')

    # Add colored markers
    scatter = axs[1].scatter(n_1, n_5, c=colors, cmap='jet')
    clb = fig.colorbar(scatter, ticks=np.linspace(0, len(n_1) + 1, 11))  # Optional: add colorbar
    clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.0f'))
    clb.ax.tick_params(labelsize=14)
    clb.ax.set_title("Iteration", loc='center', fontsize=14)

    # plot the truth (target)
    axs[1].plot(n_truth_values[1], n_truth_values[5], 'x',
             markersize=15,
             markeredgewidth=2,  # Line width of the cross
             color='red',
             label='Truth')

    axs[1].text(0.04, 0.047, 'Truth', fontsize=14, ha='center')

    axs[1].text(0.033, 0.029, 'Inversion trajectory', rotation=40, fontsize=14, ha='left')
    # Draw an arrow
    axs[1].annotate('', xy=(0.042, 0.0402), xytext=(0.032, 0.030),
                 arrowprops=dict(facecolor='blue', edgecolor='blue', linewidth=0.1))

    axs[1].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    axs[1].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    axs[1].tick_params(axis='both', which='major', labelsize=14)

    axs[1].set_xlabel('Manning\'s $n_1$', fontsize=16)
    axs[1].set_ylabel('Manning\'s $n_5$', fontsize=16)
    axs[1].set_xlim([0.025, 0.05])
    axs[1].set_ylim([0.025, 0.055])

    # add caption
    axs[0].text(-0.1, 1.05, "(a)", size=16, ha="center", transform=axs[0].transAxes)  # upper left
    axs[1].text(-0.1, 1.05, "(b)", size=16, ha="center", transform=axs[1].transAxes)

    plt.savefig("inversion_ManningN_histories_trajectories.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.show()

def plot_n_history(nIterations):
    """
    Plot this history of n values vs iterations
    """

    n_truth_values = [0.03, 0.04, 0.05, 0.03, 0.045, 0.05]

    n_0 = []
    n_1 = []
    n_2 = []
    n_3 = []
    n_4 = []
    n_5 = []

    iter_numbers = np.arange(nIterations-1)

    #loop over all inversion_callback_save_iter_*.json files to read the n value histories
    for iter in range(1, nIterations):
        print("reading n values from iteration ", iter)
        file_name = "inversion_callback_save_iter_"+ str(iter)+ ".json"
        with open(file_name, 'r') as file:
            data = json.load(file)
            n_0.append(data['theta'][0])
            n_1.append(data['theta'][1])
            n_2.append(data['theta'][2])
            n_3.append(data['theta'][3])
            n_4.append(data['theta'][4])
            n_5.append(data['theta'][5])


    # plot n histories
    plt.plot(iter_numbers, n_1, 'k', label=r'$n_1$')
    plt.plot(iter_numbers, n_2, 'r--', label=r'$n_2$')
    plt.plot(iter_numbers, n_3, 'b-.', label=r'$n_3$')
    plt.plot(iter_numbers, n_4, 'g:', label=r'$n_4$')
    plt.plot(iter_numbers[::4], n_5[::4], 'c-o', markersize=2, label=r'$n_5$',
             markerfacecoloralt='c', markerfacecolor='none', markeredgecolor='c', alpha=0.5)

    # plot truth values
    x_truth_text = 200
    plt.plot([x_truth_text, nIterations+10], [n_truth_values[1], n_truth_values[1]], 'k', alpha=0.5)
    plt.plot([x_truth_text, nIterations+10], [n_truth_values[2], n_truth_values[2]], 'r', alpha=0.5)
    plt.plot([x_truth_text, nIterations+10], [n_truth_values[3], n_truth_values[3]], 'b', alpha=0.5)
    plt.plot([x_truth_text, nIterations+10], [n_truth_values[4], n_truth_values[4]], 'g', alpha=0.5)
    plt.plot([x_truth_text, nIterations+10], [n_truth_values[5], n_truth_values[5]], 'c', alpha=0.5)

    # add truth text

    plt.text(x_truth_text, n_truth_values[1] + 0.0005, r'$n_1$ truth', fontsize=14)
    plt.text(x_truth_text, n_truth_values[2] + 0.0005, r'$n_2$ and $n_5$ truth', fontsize=14)
    plt.text(x_truth_text, n_truth_values[3] + 0.0005, r'$n_3$ truth', fontsize=14)
    plt.text(x_truth_text, n_truth_values[4] + 0.0005, r'$n_4$ truth', fontsize=14)


    ax = plt.gca()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    plt.tick_params(axis='both', which='major', labelsize=14)

    plt.xlabel('Iteration', fontsize=16)
    plt.ylabel('Manning\'s $n$ values', fontsize=16)
    #plt.xlim([0, 100])
    # plt.ylim([0, 0.01])
    #plt.yscale('log')
    legend = plt.legend(loc='lower right', fontsize=14, frameon=True, facecolor='white')

    plt.savefig("inversion_ManningN_histories.png", dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.show()

def plot_n_trajectories(nIterations):
    """
    Plot this history of n values vs iterations
    """

    n_truth_values = [0.03, 0.04, 0.05, 0.03, 0.045, 0.05]

    n_0 = []
    n_1 = []
    n_2 = []
    n_3 = []
    n_4 = []
    n_5 = []

    iter_numbers = np.arange(nIterations-1)

    #loop over all inversion_callback_save_iter_*.json files to read the n value histories
    for iter in range(1, nIterations):
        print("reading n values from iteration ", iter)
        file_name = "inversion_callback_save_iter_"+ str(iter)+ ".json"
        with open(file_name, 'r') as file:
            data = json.load(file)
            n_0.append(data['theta'][0])
            n_1.append(data['theta'][1])
            n_2.append(data['theta'][2])
            n_3.append(data['theta'][3])
            n_4.append(data['theta'][4])
            n_5.append(data['theta'][5])


    # plot n histories
    #plt.plot(n_1, n_5, 'b-o', alpha=0.25, label='Iterations')

    # Create color array from indices
    colors = np.arange(len(n_1))

    print("len(n_1)", len(n_1))

    # Plot line
    plt.plot(n_1, n_5, 'b-', alpha=0.25, label='Inversion')

    # Add colored markers
    scatter = plt.scatter(n_1, n_5, c=colors, cmap='jet')
    clb = plt.colorbar(scatter, ticks=np.linspace(0, len(n_1)+1, 11))  # Optional: add colorbar
    clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.0f'))
    clb.ax.tick_params(labelsize=14)
    clb.ax.set_title("Iteration", loc='center', fontsize=14)

    #plot the truth (target)
    plt.plot(n_truth_values[1], n_truth_values[5], 'x',
             markersize=15,
             markeredgewidth=2,  # Line width of the cross
             color='red',
             label='Truth')

    plt.text(0.04, 0.047, 'Truth', fontsize=14, ha='center')

    plt.text(0.033, 0.029, 'Inversion trajectory', rotation=40, fontsize=14, ha='left')
    # Draw an arrow
    plt.annotate('', xy=(0.042, 0.0402), xytext=(0.032, 0.030),
                 arrowprops=dict(facecolor='blue', edgecolor='blue', linewidth=0.1))

    ax = plt.gca()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    plt.tick_params(axis='both', which='major', labelsize=14)

    plt.xlabel('Manning\'s $n_1$', fontsize=16)
    plt.ylabel('Manning\'s $n_5$', fontsize=16)
    plt.xlim([0.025, 0.05])
    plt.ylim([0.025, 0.055])
    #legend = plt.legend(loc='upper left', fontsize=14, frameon=True, facecolor='white')

    plt.savefig("inversion_ManningN_trajectory.png", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()

def plot_flow_field_comparison(inversion_result_vtk_file_name):
    """
    plot flow field comparison between inversion results and inversion
    Parameters
    ----------
    inversion_result_vtk_file_name: inversion result vtk file name

    Returns
    -------

    """

    # triangulate the inversion result vtk file first
    triangulate_vtk_files([inversion_result_vtk_file_name])

    # get the new tri file name for vtk
    file_name_without_extension = os.path.splitext(inversion_result_vtk_file_name)[0]
    inversion_result_tri_vtk_file_name = file_name_without_extension + '_tri.vtk'

    # this should be from forward simulation, which is the truth.
    truth_result_vtk_file_name = 'flow_filed_truth_results_tri.vtk'

    #read data from vtk files
    points_inv, cells_inv, tris_inv, xCoord_inv, yCoord_inv, xl_inv, xh_inv, yl_inv, yh_inv, WSE_inv, U_inv, V_inv = read_data_from_vtk(inversion_result_tri_vtk_file_name)
    points_truth, cells_truth, tris_truth, xCoord_truth, yCoord_truth, xl_truth, xh_truth, yl_truth, yh_truth, WSE_truth, U_truth, V_truth = read_data_from_vtk(truth_result_vtk_file_name)

    fig, axs = plt.subplots(3, 3, figsize=(18, 8), sharex=True, sharey=False, facecolor='w', edgecolor='k')

    #plot truth
    # plot WSE_truth
    min_value = np.min(WSE_truth)
    max_value = np.max(WSE_truth)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(0, 0, "$WSE_{truth}$", fig, axs[0, 0], xCoord_truth, yCoord_truth, tris_truth, WSE_truth,
                     local_levels, min_value, max_value, xl_truth, xh_truth, yl_truth, yh_truth, bShow_x_title=False, bShow_y_title=True)

    # plot U_truth
    min_value = np.min(U_truth)
    max_value = np.max(U_truth)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(1, 0, "$u_{truth}$", fig, axs[1, 0], xCoord_truth, yCoord_truth, tris_truth, U_truth,
                     local_levels, min_value, max_value, xl_truth, xh_truth, yl_truth, yh_truth, bShow_x_title=False,
                     bShow_y_title=True)

    # plot V_truth
    min_value = np.min(V_truth)
    max_value = np.max(V_truth)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(2, 0, "$v_{truth}$", fig, axs[2, 0], xCoord_truth, yCoord_truth, tris_truth, V_truth,
                     local_levels, min_value, max_value, xl_truth, xh_truth, yl_truth, yh_truth, bShow_x_title=True,
                     bShow_y_title=True)

    # plot inversion
    # plot WSE_inv
    min_value = np.min(WSE_inv)
    max_value = np.max(WSE_inv)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(0, 1, "$WSE_{inversion}$", fig, axs[0, 1], xCoord_inv, yCoord_inv, tris_inv, WSE_inv,
                     local_levels, min_value, max_value, xl_inv, xh_inv, yl_inv, yh_inv, bShow_x_title=False,
                     bShow_y_title=False)

    # plot U_inv
    min_value = np.min(U_inv)
    max_value = np.max(U_inv)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(1, 1, "$u_{inversion}$", fig, axs[1, 1], xCoord_inv, yCoord_inv, tris_inv, U_inv,
                     local_levels, min_value, max_value, xl_inv, xh_inv, yl_inv, yh_inv, bShow_x_title=False,
                     bShow_y_title=False)

    # plot V_inv
    min_value = np.min(V_inv)
    max_value = np.max(V_inv)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(2, 1, "$v_{inversion}$", fig, axs[2, 1], xCoord_inv, yCoord_inv, tris_inv, V_inv,
                     local_levels, min_value, max_value, xl_inv, xh_inv, yl_inv, yh_inv, bShow_x_title=True,
                     bShow_y_title=False)

    # plot diff
    WSE_diff = WSE_truth - WSE_inv
    U_diff = U_truth - U_inv
    V_diff = V_truth - V_inv
    
    # plot WSE_diff
    min_value = np.min(WSE_diff)
    max_value = np.max(WSE_diff)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(0, 2, "$WSE_{diff}$", fig, axs[0, 2], xCoord_inv, yCoord_inv, tris_inv, WSE_diff,
                     local_levels, min_value, max_value, xl_inv, xh_inv, yl_inv, yh_inv, bShow_x_title=False,
                     bShow_y_title=False)

    # plot U_diff
    min_value = np.min(U_diff)
    max_value = np.max(U_diff)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(1, 2, "$u_{diff}$", fig, axs[1, 2], xCoord_inv, yCoord_inv, tris_inv, U_diff,
                     local_levels, min_value, max_value, xl_inv, xh_inv, yl_inv, yh_inv, bShow_x_title=False,
                     bShow_y_title=False)

    # plot V_diff
    min_value = np.min(V_diff)
    max_value = np.max(V_diff)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(2, 2, "$v_{diff}$", fig, axs[2, 2], xCoord_inv, yCoord_inv, tris_inv, V_diff,
                     local_levels, min_value, max_value, xl_inv, xh_inv, yl_inv, yh_inv, bShow_x_title=True,
                     bShow_y_title=False)

    #add caption
    axs[0, 0].text(-0.05, 1.05, "(a)", size=16, ha="center", transform=axs[0, 0].transAxes)  # upper left
    axs[0, 1].text(-0.05, 1.05, "(b)", size=16, ha="center", transform=axs[0, 1].transAxes)
    axs[0, 2].text(-0.05, 1.05, "(c)", size=16, ha="center", transform=axs[0, 2].transAxes)
    axs[1, 0].text(-0.05, 1.05, "(d)", size=16, ha="center", transform=axs[1, 0].transAxes)  # upper left
    axs[1, 1].text(-0.05, 1.05, "(e)", size=16, ha="center", transform=axs[1, 1].transAxes)
    axs[1, 2].text(-0.05, 1.05, "(f)", size=16, ha="center", transform=axs[1, 2].transAxes)
    axs[2, 0].text(-0.05, 1.05, "(g)", size=16, ha="center", transform=axs[2, 0].transAxes)  # upper left
    axs[2, 1].text(-0.05, 1.05, "(h)", size=16, ha="center", transform=axs[2, 1].transAxes)
    axs[2, 2].text(-0.05, 1.05, "(i)", size=16, ha="center", transform=axs[2, 2].transAxes)

    plt.savefig("inversion_Savana_flow_filed_comparison.png", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()

if __name__ == '__main__':

    nIterations = 300

    plot_n_history_trajectories(nIterations)

    #plot comparison of flow field (WSE, u, v) between truth and inversion result
    inversion_result_vtk_file_name = 'forward_simulation_results_0300.vtk'
    #plot_flow_field_comparison(inversion_result_vtk_file_name)

    print("Done!")