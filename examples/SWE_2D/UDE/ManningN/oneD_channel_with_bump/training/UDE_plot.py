"""
Some specific utility functions, such as plotting.

"""

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

import pyHMT2D
from pyHMT2D.Misc.tools import setNumpyArrayValueToNaN

plt.rc('text', usetex=True)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def plot_UDE_result():
    """
    Make plots for the UDE result:
    1. loss history
    2. WSE over zb
    3. Manning n truth vs NN
    4. sigmoid vs NN

    :return:
    """

    fig, axs = plt.subplots(2, 2, figsize=(14, 9), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    fig.subplots_adjust(hspace=0.3, wspace=.2)

    #plot WSE: simulated vs truth
    # Define the start and end points of the line (example: from (0, 0) to (1, 1))
    p1 = [0.0, 0.5, 0.0]
    p2 = [25.0, 0.5, 0.0]

    # read the sampling input file
    vtk_handler = pyHMT2D.Misc.vtkHandler()

    # Define the numer of interpolation points
    numPoints = 500

    # for truth
    reader = vtk_handler.readVTK_UnstructuredGrid('forward_simulation_truth.vtk')  # read the VTKfile
    line = vtk_handler.createVtkLine(p1, p2, numPoints)  # Create the line
    points_truth, WSE_truth, zb_truth = vtk_handler.probeUnstructuredGridVTKOverLine(line, reader, 'WSE')  # interpolate the data over the line

    WSE_truth = setNumpyArrayValueToNaN(WSE_truth, 0)  # Set the zero's to NaN's

    # for UDE
    reader = vtk_handler.readVTK_UnstructuredGrid('forward_simulation_results_0100.vtk')  # read the VTKfile
    line = vtk_handler.createVtkLine(p1, p2, numPoints)  # Create the line
    points_UDE, WSE_UDE, zb_UDE = vtk_handler.probeUnstructuredGridVTKOverLine(line, reader, 'WSE')
    points_UDE, h_UDE, zb_UDE = vtk_handler.probeUnstructuredGridVTKOverLine(line, reader, 'h')
    points_UDE, ManningN_UDE, zb_UDE = vtk_handler.probeUnstructuredGridVTKOverLine(line, reader, 'ManningN_cells')
    points_UDE, ManningN_truth, zb_UDE = vtk_handler.probeUnstructuredGridVTKOverLine(line, reader, 'ManningN_cell_truth')

    WSE_truth = setNumpyArrayValueToNaN(WSE_truth, 0)  # Set the zero's to NaN's
    ManningN_UDE = setNumpyArrayValueToNaN(ManningN_UDE, 0)  # Set the zero's to NaN's
    ManningN_truth = setNumpyArrayValueToNaN(ManningN_truth, 0)  # Set the zero's to NaN's

    # plot WSE
    axs[0, 0].plot(points_UDE[:,0], WSE_UDE, color='blue', linestyle='--', label='From UDE')
    axs[0, 0].plot(points_truth[:,0], WSE_truth, color='red', linestyle='-.', label='Truth')

    # plot the bed
    axs[0, 0].plot(points_UDE[:,0], zb_UDE, color='gray', label='Bed')

    # set the limit for the x and y axes
    axs[0, 0].set_xlim([min(points_UDE[:,0]), max(points_UDE[:,0])])
    axs[0, 0].set_ylim([0, 0.6])

    # set x and y axes label and font size
    axs[0, 0].set_xlabel('x (m)', fontsize=16)
    axs[0, 0].set_ylabel('z (m)', fontsize=16)

    # show the ticks on both axes and set the font size
    axs[0, 0].tick_params(axis='both', which='major', labelsize=14)

    # show legend, set its location, font size, and turn off the frame
    axs[0, 0].legend(loc='upper right', fontsize=14, frameon=False)


    #plot loss history
    # read the loss data from the "inversion_results_losses_and_ODE_solution.json" file
    with open('UDE_results_losses_and_ODE_solution.json', 'r') as file:
        data = json.load(file)

    loss_preds = data['loss_preds']
    loss_pred_WSEs = data['loss_pred_WSEs']
    loss_pred_uvs = data['loss_pred_uvs']

    iter_numbers = np.arange(len(loss_preds))

    # plot losses
    axs[0, 1].plot(iter_numbers, loss_preds, 'k', label=r'$l_{pred, total}$')
    axs[0, 1].plot(iter_numbers, loss_pred_WSEs, 'r--', label=r'$l_{pred, WSE}$')
    axs[0, 1].plot(iter_numbers, loss_pred_uvs, 'b-.', label=r'$l_{pred, uv}$')
    # plt.plot(iter_numbers, loss_bounds, 'r--', label=r'$l_{bound}$')

    axs[0, 1].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    # ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    axs[0, 1].tick_params(axis='both', which='major', labelsize=14)

    axs[0, 1].set_xlabel('Iteration', fontsize=16)
    axs[0, 1].set_ylabel('Loss', fontsize=16)
    # axs[0, 1].xlim([0, 100])
    # axs[0, 1].ylim([1e-5, 1e-2])
    axs[0, 1].set_yscale('log')

    axs[0, 1].legend(loc='upper right', fontsize=14, frameon=False)

    #plot scatter n values from NN agains truth
    axs[1, 0].scatter(ManningN_UDE, ManningN_truth,  edgecolors='blue', facecolors='none', s=50)
    #draw the diagonal line
    axs[1, 0].plot([0.02,0.07],[0.02, 0.07], color='black')

    axs[1, 0].set_xlim([0.02, 0.07])
    axs[1, 0].set_ylim([0.02, 0.07])

    axs[1, 0].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axs[1, 0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    axs[1, 0].tick_params(axis='both', which='major', labelsize=14)

    axs[1, 0].set_xlabel('Manning\'s $n$ from NN', fontsize=16)
    axs[1, 0].set_ylabel('Manning\'s $n$ truth', fontsize=16)
    axs[1, 0].set_aspect('equal')

    #plot h vs n
    h_truth = np.linspace(0.0, 0.6, 100)
    n_lower = 0.03
    n_upper = 0.06
    k = 100.0
    h_mid = 0.3
    n_truth = n_lower + (n_upper - n_lower)/(1.0 + np.exp(k*(h_truth-h_mid)))

    axs[1, 1].plot(n_truth, h_truth, color='black', label='Truth (Sigmoid function)')
    axs[1, 1].scatter(ManningN_UDE, h_UDE, edgecolors='blue', facecolors='none', s=50, label='NN approximation')

    axs[1, 1].set_xlim([0.02, 0.07])
    axs[1, 1].set_ylim([0.0, 0.6])

    #axs[1, 1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    #axs[1, 1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    axs[1, 1].tick_params(axis='both', which='major', labelsize=14)

    axs[1, 1].set_xlabel('Manning\'s $n$', fontsize=16)
    axs[1, 1].set_ylabel('Water depth $h$', fontsize=16)

    axs[1, 1].legend(loc='upper right', fontsize=14, frameon=False)

    #add caption
    axs[0, 0].text(-0.1, 1.05, "(a)", size=16, ha="center", transform=axs[0,0].transAxes)  # upper left
    axs[0, 1].text(-0.1, 1.05, "(b)", size=16, ha="center", transform=axs[0,1].transAxes)
    axs[1, 0].text(-0.2, 1.05, "(c)", size=16, ha="center", transform=axs[1, 0].transAxes)  # upper left
    axs[1, 1].text(-0.1, 1.05, "(d)", size=16, ha="center", transform=axs[1, 1].transAxes)

    plt.savefig("oneD_channel_with_bump_UDE_result.png", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()


if __name__ == '__main__':

    #plot UDE result
    plot_UDE_result()

    print("Done!")