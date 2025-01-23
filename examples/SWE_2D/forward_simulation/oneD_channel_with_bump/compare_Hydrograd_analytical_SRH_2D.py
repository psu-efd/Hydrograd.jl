import numpy as np
import meshio
import vtk
from vtk.util import numpy_support as VN
import h5py
import os

import matplotlib.pyplot as plt

import pyHMT2D
from pyHMT2D.Misc.tools import setNumpyArrayValueToNaN

plt.rc('text', usetex=True)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def plot_surface(x_Hydrograd, zb_Hydrograd, WSE_Hydrograd,
                 x_srh_2d, WSE_srh_2d,
                 x_analytical, wse_analytical, zb_analytical):
    """

    Parameters
    ----------
    surface
    sampled_points
    sampled_scalars

    Returns
    -------

    """

    # plot WSE
    plt.plot(x_analytical, wse_analytical, color='black', label='Analytical')
    plt.plot(x_Hydrograd, WSE_Hydrograd, color='blue', linestyle='--', label='Hydrograd')
    plt.plot(x_srh_2d, WSE_srh_2d, color='red', linestyle='-.', label='SRH-2D')

    # plot the bed
    plt.plot(x_Hydrograd, zb_Hydrograd, color='gray', label='Bed')
    # plt.plot(x_analytical, zb_analytical, color='blue', linestyle='--', label='Bed(ana)')

    # set the limit for the x and y axes
    plt.xlim([min(x_Hydrograd), max(x_Hydrograd)])
    plt.ylim([0, 0.6])

    # set x and y axes label and font size
    plt.xlabel('x', fontsize=16)
    plt.ylabel('z', fontsize=16)

    # show the ticks on both axes and set the font size
    plt.tick_params(axis='both', which='major', labelsize=14)

    # show legend, set its location, font size, and turn off the frame
    plt.legend(loc='upper right', fontsize=14, frameon=False)

    #save to file
    plt.savefig('oneD_channel_with_bump_forward_validation.png', dpi=300, bbox_inches='tight', pad_inches=0)

    plt.show()



if __name__ == "__main__":

    #Load the analytical solution
    data_analytical = np.loadtxt("analytical_solution_from_FullSWOF.dat", skiprows=1)
    x_analytical = data_analytical[:, 0]
    h_analytical = data_analytical[:, 1]
    zb_analytical = data_analytical[:, 2]
    wse_analytical = data_analytical[:, 3]
    #print(x_analytical)
    #print(wse_analytical)

    #exit()

    # Define the filename of VTK file
    srh_filename = os.path.join('SRH_2D', 'SRH2D_oneD_channel_with_bump_refined_C_0010.vtk')

    Hydrograd_filename = 'forward_simulation_results_0101.vtk'

    num_samples = 100

    #print(vtk_data_AdHydraulics)

    # Define the start and end points of the line (example: from (0, 0) to (1, 1))
    p1 = [0.0, 0.5, 0.0]
    p2 = [25.0, 0.5, 0.0]

    # read the sampling input file
    vtk_handler = pyHMT2D.Misc.vtkHandler()

    # Define the numer of interpolation points
    numPoints = 500

    # for SRH-2D
    reader = vtk_handler.readVTK_UnstructuredGrid(srh_filename)  # read the VTKfile
    line = vtk_handler.createVtkLine(p1, p2, numPoints)  # Create the line
    points_srh_2d, WSE_srh_2d, zb_srh_2d = vtk_handler.probeUnstructuredGridVTKOverLine(line, reader, 'Water_Elev_m')  # interpolate the data over the line
    # print(points)
    # print(U)

    #max_x = np.max(points_srh_2d[:, 0])
    #min_x = np.min(points_srh_2d[:, 0])

    WSE_srh_2d = setNumpyArrayValueToNaN(WSE_srh_2d, 0)  # Set the zero's to NaN's

    # for Hydrograd
    reader = vtk_handler.readVTK_UnstructuredGrid(Hydrograd_filename)  # read the VTKfile
    line = vtk_handler.createVtkLine(p1, p2, numPoints)  # Create the line
    points_Hydrograd, WSE_Hydrograd, zb_Hydrograd = vtk_handler.probeUnstructuredGridVTKOverLine(line, reader, 'WSE')  # interpolate the data over the line
    # print(points)
    # print(U)

    max_x = np.max(points_Hydrograd[:, 0])
    min_x = np.min(points_Hydrograd[:, 0])

    WSE_Hydrograd = setNumpyArrayValueToNaN(WSE_Hydrograd, 0)  # Set the zero's to NaN's

    plot_surface(points_Hydrograd[:,0], zb_Hydrograd, WSE_Hydrograd,
                 points_srh_2d[:,0], WSE_srh_2d,
                 x_analytical, wse_analytical, zb_analytical)

    print("Done!")