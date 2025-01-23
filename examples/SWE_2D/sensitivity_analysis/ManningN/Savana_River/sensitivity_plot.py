"""
Some specific utility functions, such as plotting.

"""

import os

import matplotlib.pyplot as plt
import numpy as np

import shapefile

import matplotlib.ticker as tick
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

    if iRow == 0 and iCol == 0:
        axs.set_title('$WSE$\nsensitivity', fontsize=14)
    elif iRow == 0 and iCol == 1:
        axs.set_title('$hu$\nsensitivity', fontsize=14)
    elif iRow == 0 and iCol == 2:
        axs.set_title('$hv$\nsensitivity', fontsize=14)

    axs.text(600, 340, sensitivity_name, ha='center', va='bottom', fontsize=14)

    if bShow_x_title:
        axs.set_xlabel('$x$ (m)', fontsize=16)
        axs.tick_params(axis='x', labelsize=14)
    else:
        axs.tick_params(axis='x', length=0, labelbottom=False, labelleft=False)

    if bShow_y_title:
        axs.set_ylabel('$y$ (m)', fontsize=16)
        axs.tick_params(axis='y', labelsize=14)
    else:
        axs.tick_params(axis='y', length=0, labelbottom=False, labelleft=False)

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
    axs[0, 0].text(-0.05, 1.05, "(a)", size=16, ha="center", transform=axs[0, 0].transAxes)  # upper left
    axs[0, 1].text(-0.05, 1.05, "(b)", size=16, ha="center", transform=axs[0, 1].transAxes)
    axs[0, 2].text(-0.05, 1.05, "(c)", size=16, ha="center", transform=axs[0, 2].transAxes)
    axs[1, 0].text(-0.05, 1.05, "(d)", size=16, ha="center", transform=axs[1, 0].transAxes)  # upper left
    axs[1, 1].text(-0.05, 1.05, "(e)", size=16, ha="center", transform=axs[1, 1].transAxes)
    axs[1, 2].text(-0.05, 1.05, "(f)", size=16, ha="center", transform=axs[1, 2].transAxes)
    axs[2, 0].text(-0.05, 1.05, "(g)", size=16, ha="center", transform=axs[2, 0].transAxes)  # upper left
    axs[2, 1].text(-0.05, 1.05, "(h)", size=16, ha="center", transform=axs[2, 1].transAxes)
    axs[2, 2].text(-0.05, 1.05, "(i)", size=16, ha="center", transform=axs[2, 2].transAxes)
    axs[3, 0].text(-0.05, 1.05, "(j)", size=16, ha="center", transform=axs[3, 0].transAxes)  # upper left
    axs[3, 1].text(-0.05, 1.05, "(k)", size=16, ha="center", transform=axs[3, 1].transAxes)
    axs[3, 2].text(-0.05, 1.05, "(l)", size=16, ha="center", transform=axs[3, 2].transAxes)
    axs[4, 0].text(-0.05, 1.05, "(m)", size=16, ha="center", transform=axs[4, 0].transAxes)  # upper left
    axs[4, 1].text(-0.05, 1.05, "(n)", size=16, ha="center", transform=axs[4, 1].transAxes)
    axs[4, 2].text(-0.05, 1.05, "(o)", size=16, ha="center", transform=axs[4, 2].transAxes)

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


if __name__ == '__main__':

    #triangulate all
    #Note: 1 is for default zone, 2 is for zone 1, etc.
    vtk_file_names = ['sensitivity_results_ManningN_1.vtk',
                      'sensitivity_results_ManningN_2.vtk',
                      'sensitivity_results_ManningN_3.vtk',
                      'sensitivity_results_ManningN_4.vtk',
                      'sensitivity_results_ManningN_5.vtk',
                      'sensitivity_results_ManningN_6.vtk']
    #triangulate_vtk_files(vtk_file_names)

    plot_sensitivity_result()

    print("Done!")