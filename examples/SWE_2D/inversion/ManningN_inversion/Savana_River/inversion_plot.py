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
    dhdn = np.squeeze(vtk_result.point_data['dh_dparam'])
    dhudn = np.squeeze(vtk_result.point_data['dhu_dparam'])
    dhvdn = np.squeeze(vtk_result.point_data['dhv_dparam'])

    return points, cells, tris, xCoord, yCoord, xl, xh, yl, yh, dhdn, dhudn, dhvdn

def plot_one_subplot(iRow, iCol, sensitivity_name, fig, axs, xCoord, yCoord, tris, sensitivity, local_levels,
                     min_value, max_value, xl, xh, yl, yh, bShow_x_title=False, bShow_y_title=False):
    #just plot a specified subplot

    cf_zb = axs.tricontourf(xCoord, yCoord, tris, sensitivity, local_levels,
                                  vmin=min_value, vmax=max_value, cmap=plt.cm.jet, extend='neither')
    axs.set_xlim([xl, xh])
    axs.set_ylim([yl, yh])
    axs.set_aspect('equal')
    #axs.set_title(sensitivity_name, fontsize=14)

    axs.text(600, 340, sensitivity_name, ha='center', va='bottom', fontsize=14)

    if bShow_x_title: axs.set_xlabel('$x$ (m)', fontsize=16)
    axs.tick_params(axis='x', labelsize=14)
    if bShow_y_title: axs.set_ylabel('$y$ (m)', fontsize=16)
    axs.tick_params(axis='y', labelsize=14)

    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="3%", pad=0.2)
    clb = fig.colorbar(cf_zb, ticks=np.linspace(min_value, max_value, 7), cax=cax)
    if iCol==0:
        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    else:
        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))
    clb.ax.tick_params(labelsize=12)
    #clb.ax.set_title("(-)", loc='center', fontsize=12)

def plot_sensitivity_result():
    """
    Make plots for the sensitivity analysis:
    1. dWSE/dn
    2. dhu/dn
    3. dhv/dn

    :return:
    """

    #read data from vtk files
    #Note: 1 is the default Manning's n; 2 is for zone 1 n, etc.

    points_1, cells_1, tris_1, xCoord_1, yCoord_1, xl_1, xh_1, yl_1, yh_1, dhdn_1, dhudn_1, dhvdn_1 = read_data_from_vtk(
        'sensitivity_results_ManningN_2_tri.vtk')
    points_2, cells_2, tris_2, xCoord_2, yCoord_2, xl_2, xh_2, yl_2, yh_2, dhdn_2, dhudn_2, dhvdn_2 = read_data_from_vtk(
        'sensitivity_results_ManningN_3_tri.vtk')
    points_3, cells_3, tris_3, xCoord_3, yCoord_3, xl_3, xh_3, yl_3, yh_3, dhdn_3, dhudn_3, dhvdn_3 = read_data_from_vtk(
        'sensitivity_results_ManningN_4_tri.vtk')
    points_4, cells_4, tris_4, xCoord_4, yCoord_4, xl_4, xh_4, yl_4, yh_4, dhdn_4, dhudn_4, dhvdn_4 = read_data_from_vtk(
        'sensitivity_results_ManningN_5_tri.vtk')
    points_5, cells_5, tris_5, xCoord_5, yCoord_5, xl_5, xh_5, yl_5, yh_5, dhdn_5, dhudn_5, dhvdn_5 = read_data_from_vtk(
        'sensitivity_results_ManningN_6_tri.vtk')

    fig, axs = plt.subplots(5, 3, figsize=(13, 9), sharex=True, sharey=False, facecolor='w', edgecolor='k')

    #fig.subplots_adjust(hspace=.15, wspace=.01)

    #plot sensitivity (dWSE/dn)
    #dWSE/dn1
    min_value = np.min(dhdn_1)
    max_value = np.max(dhdn_1)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(0, 0,"$\partial WSE/\partial n_1$", fig, axs[0, 0], xCoord_1, yCoord_1, tris_1, dhdn_1, local_levels,
                     min_value, max_value, xl_1, xh_1, yl_1, yh_1, bShow_x_title=False, bShow_y_title=True)

    # dWSE/dn2
    min_value = np.min(dhdn_2)
    max_value = np.max(dhdn_2)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(1, 0,"$\partial WSE/\partial n_2$", fig, axs[1, 0], xCoord_2, yCoord_2, tris_2, dhdn_2, local_levels,
                     min_value, max_value, xl_2, xh_2, yl_2, yh_2, bShow_x_title=False, bShow_y_title=True)

    # dWSE/dn3
    min_value = np.min(dhdn_3)
    max_value = np.max(dhdn_3)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(2, 0,"$\partial WSE/\partial n_3$", fig, axs[2, 0], xCoord_3, yCoord_3, tris_3, dhdn_3, local_levels,
                     min_value, max_value, xl_3, xh_3, yl_3, yh_3, bShow_x_title=False, bShow_y_title=True)

    # dWSE/dn4
    min_value = np.min(dhdn_4)
    max_value = np.max(dhdn_4)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(3, 0,"$\partial WSE/\partial n_4$", fig, axs[3, 0], xCoord_4, yCoord_4, tris_4, dhdn_4, local_levels,
                     min_value, max_value, xl_4, xh_4, yl_4, yh_4, bShow_x_title=False, bShow_y_title=True)

    # dWSE/dn5
    min_value = np.min(dhdn_5)
    max_value = np.max(dhdn_5)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(4, 0,"$\partial WSE/\partial n_5$", fig, axs[4, 0], xCoord_5, yCoord_5, tris_5, dhdn_5, local_levels,
                     min_value, max_value, xl_5, xh_5, yl_5, yh_5, bShow_x_title=True, bShow_y_title=True)

    # plot sensitivity (dhu/dn)
    # dhu/dn1
    min_value = np.min(dhudn_1)
    max_value = np.max(dhudn_1)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(0, 1,"$\partial hu/\partial n_1$", fig, axs[0, 1], xCoord_1, yCoord_1, tris_1, dhudn_1, local_levels,
                     min_value, max_value, xl_1, xh_1, yl_1, yh_1, bShow_x_title=False, bShow_y_title=False)

    # dhu/dn2
    min_value = np.min(dhudn_2)
    max_value = np.max(dhudn_2)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(1, 1,"$\partial hu/\partial n_2$", fig, axs[1, 1], xCoord_2, yCoord_2, tris_2, dhudn_2, local_levels,
                     min_value, max_value, xl_2, xh_2, yl_2, yh_2, bShow_x_title=False, bShow_y_title=False)

    # dhu/dn3
    min_value = np.min(dhudn_3)
    max_value = np.max(dhudn_3)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(2, 1,"$\partial hu/\partial n_3$", fig, axs[2, 1], xCoord_3, yCoord_3, tris_3, dhudn_3, local_levels,
                     min_value, max_value, xl_3, xh_3, yl_3, yh_3, bShow_x_title=False, bShow_y_title=False)

    # dhu/dn4
    min_value = np.min(dhudn_4)
    max_value = np.max(dhudn_4)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(3, 1,"$\partial hu/\partial n_4$", fig, axs[3, 1], xCoord_4, yCoord_4, tris_4, dhudn_4, local_levels,
                     min_value, max_value, xl_4, xh_4, yl_4, yh_4, bShow_x_title=False, bShow_y_title=False)

    # dhu/dn5
    min_value = np.min(dhudn_5)
    max_value = np.max(dhudn_5)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(4, 1, "$\partial hu/\partial n_5$", fig, axs[4, 1], xCoord_5, yCoord_5, tris_5, dhudn_5, local_levels,
                     min_value, max_value, xl_5, xh_5, yl_5, yh_5, bShow_x_title=True, bShow_y_title=False)

    # plot sensitivity (dhv/dn)
    # dhv/dn1
    min_value = np.min(dhvdn_1)
    max_value = np.max(dhvdn_1)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(0, 2,"$\partial hv/\partial n_1$", fig, axs[0, 2], xCoord_1, yCoord_1, tris_1, dhvdn_1, local_levels,
                     min_value, max_value, xl_1, xh_1, yl_1, yh_1, bShow_x_title=False, bShow_y_title=False)

    # dhv/dn2
    min_value = np.min(dhvdn_2)
    max_value = np.max(dhvdn_2)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(1,2,"$\partial hv/\partial n_2$", fig, axs[1, 2], xCoord_2, yCoord_2, tris_2, dhvdn_2, local_levels,
                     min_value, max_value, xl_2, xh_2, yl_2, yh_2, bShow_x_title=False, bShow_y_title=False)

    # dhv/dn3
    min_value = np.min(dhvdn_3)
    max_value = np.max(dhvdn_3)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(2, 2,"$\partial hv/\partial n_3$", fig, axs[2, 2], xCoord_3, yCoord_3, tris_3, dhvdn_3, local_levels,
                     min_value, max_value, xl_3, xh_3, yl_3, yh_3, bShow_x_title=False, bShow_y_title=False)

    # dhv/dn4
    min_value = np.min(dhvdn_4)
    max_value = np.max(dhvdn_4)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(3, 2,"$\partial hv/\partial n_4$", fig, axs[3, 2], xCoord_4, yCoord_4, tris_4, dhvdn_4, local_levels,
                     min_value, max_value, xl_4, xh_4, yl_4, yh_4, bShow_x_title=False, bShow_y_title=False)

    # dhv/dn5
    min_value = np.min(dhvdn_5)
    max_value = np.max(dhvdn_5)
    print("min and max values: ", min_value, max_value)

    local_levels = np.linspace(min_value, max_value, 51)

    plot_one_subplot(4,2,"$\partial hv/\partial n_5$", fig, axs[4, 2], xCoord_5, yCoord_5, tris_5, dhvdn_5, local_levels,
                     min_value, max_value, xl_5, xh_5, yl_5, yh_5, bShow_x_title=True, bShow_y_title=False)

    #add caption
    #axs[0,0].text(-0.1, 1.05, "(a)", size=16, ha="center", transform=axs[0,0].transAxes)  # upper left
    #axs[0,1].text(-0.1, 1.05, "(b)", size=16, ha="center", transform=axs[0,1].transAxes)
    #axs[1, 0].text(-0.1, 1.05, "(c)", size=16, ha="center", transform=axs[1, 0].transAxes)  # upper left
    #axs[1, 1].text(-0.1, 1.05, "(d)", size=16, ha="center", transform=axs[1, 1].transAxes)

    plt.savefig("Savana_sensitivity_result.png", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()

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
    plt.plot(iter_numbers[::4], n_5[::4], 'c-o', markersize=3, label=r'$n_5$')

    # plot truth values
    plt.plot([nIterations-50, nIterations], [n_truth_values[1], n_truth_values[1]], 'k')
    plt.plot([nIterations-50, nIterations], [n_truth_values[2], n_truth_values[2]], 'r')
    plt.plot([nIterations-50, nIterations], [n_truth_values[3], n_truth_values[3]], 'b')
    plt.plot([nIterations-50, nIterations], [n_truth_values[4], n_truth_values[4]], 'g')
    plt.plot([nIterations-50, nIterations], [n_truth_values[5], n_truth_values[5]], 'c')


    ax = plt.gca()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    plt.tick_params(axis='both', which='major', labelsize=14)

    plt.xlabel('Iteration', fontsize=16)
    plt.ylabel('Manning\'s $n$ values', fontsize=16)
    #plt.xlim([0, 100])
    # plt.ylim([0, 0.01])
    legend = plt.legend(loc='upper left', fontsize=14, frameon=True, facecolor='white')

    plt.savefig("ManningN_histories.png", dpi=300, bbox_inches='tight', pad_inches=0)
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

    ax = plt.gca()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    plt.tick_params(axis='both', which='major', labelsize=14)

    plt.xlabel('Manning\'s $n_1$', fontsize=16)
    plt.ylabel('Manning\'s $n_5$', fontsize=16)
    plt.xlim([0.025, 0.05])
    plt.ylim([0.025, 0.055])
    legend = plt.legend(loc='upper left', fontsize=14, frameon=True, facecolor='white')

    plt.savefig("ManningN_trajectory.png", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()


if __name__ == '__main__':

    #triangulate all
    #Note: 1 is for default zone, 2 is for zone 1, etc.
    #vtk_file_names = ['sensitivity_results_ManningN_1.vtk',
    #                  'sensitivity_results_ManningN_2.vtk']
    #triangulate_vtk_files(vtk_file_names)

    nIterations = 100

    plot_loss_history(nIterations)

    #plot_n_history(nIterations)

    #plot_n_trajectories(nIterations)



    print("Done!")