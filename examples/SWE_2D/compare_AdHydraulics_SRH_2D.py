import numpy as np
import meshio
import vtk
from vtk.util import numpy_support as VN

import matplotlib.pyplot as plt

import pyHMT2D
from pyHMT2D.Misc.tools import setNumpyArrayValueToNaN

plt.rc('text', usetex=True)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def flatten_surface(vtk_file):
    """
    Flattens a 3D VTK surface to the XY plane by setting Z-coordinates to zero.

    Parameters:
        vtk_file (str): Path to the input VTK file.

    Returns:
        vtkPolyData: Flattened VTK surface.
    """
    # Read the VTK file
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtk_file)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    vtk_surface = reader.GetOutput()

    # Create a transform to flatten the surface
    transform = vtk.vtkTransform()
    transform.Scale(1, 1, 0)  # Flatten by setting Z-scale to 0

    # Apply the transform to the surface
    transform_filter = vtk.vtkTransformFilter()
    transform_filter.SetInputData(vtk_surface)
    transform_filter.SetTransform(transform)
    transform_filter.Update()

    return transform_filter.GetOutput()


def sample_on_line(flattened_surface, p1, p2, num_samples):
    """
    Samples scalar data on a line from a flattened VTK surface.

    Parameters:
        flattened_surface (vtkPolyData): Flattened VTK surface.
        p1 (tuple): Starting point of the line (x, y, 0).
        p2 (tuple): Ending point of the line (x, y, 0).
        num_samples (int): Number of points to sample along the line.

    Returns:
        sampled_points (list): Sampled (x, y, z) coordinates.
        sampled_scalars (list): Scalar values sampled along the line.
    """
    # Define the sampling line
    line = vtk.vtkLineSource()
    line.SetPoint1(p1)
    line.SetPoint2(p2)
    line.SetResolution(num_samples - 1)
    line.Update()

    # Probe the data using the line
    probe_filter = vtk.vtkProbeFilter()
    probe_filter.SetInputConnection(line.GetOutputPort())
    probe_filter.SetSourceData(flattened_surface)
    probe_filter.Update()

    # Extract the probed output
    probed_data = probe_filter.GetOutput()

    print(probed_data)

    # Retrieve sampled scalar values
    scalar_array = probed_data.GetPointData().GetScalars()
    if scalar_array is None:
        raise ValueError("No scalar data found in the flattened surface.")

    sampled_scalars = [scalar_array.GetValue(i) for i in range(scalar_array.GetNumberOfTuples())]

    # Retrieve the sampled points
    sampled_points = [
        probed_data.GetPoint(i) for i in range(probed_data.GetNumberOfPoints())
    ]

    return sampled_points, sampled_scalars

def plot_surface(x, zb, WSE_AdHydraulics, WSE_srh_2dD):
    """

    Parameters
    ----------
    surface
    sampled_points
    sampled_scalars

    Returns
    -------

    """

    # plot the bed
    plt.plot(x, zb, color='black', label='Bed')

    # plot WSE
    plt.plot(x, WSE_AdHydraulics, color='blue', label='AdHydraulics')
    plt.plot(x, WSE_srh_2dD, color='red', label='SRH')

    # set the limit for the x and y axes
    plt.xlim([min(x), max(x)])
    plt.ylim([0, 1])

    # set x and y axes label and font size
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)

    # show the ticks on both axes and set the font size
    plt.tick_params(axis='both', which='major', labelsize=12)

    # show title and set font size
    #plt.title('Plot of data read with Numpy loadtxt function', fontsize=16)

    # show legend, set its location, font size, and turn off the frame
    plt.legend(loc='upper right', fontsize=14, frameon=False)
    plt.show()



if __name__ == "__main__":

    # Define the filename of VTK file
    srh_filename = 'SRH2D_twoD_channel_with_bump_C_0001.vtk'
    #AdHydraulics_filename = 'solution_final.vtk'
    AdHydraulics_filename = 'solution_3800.vtk'

    num_samples = 100

    #print(vtk_data_AdHydraulics)

    # Define the start and end points of the line (example: from (0, 0) to (1, 1))
    p1 = [0.0, 0.5, 0.0]
    p2 = [25.0, 0.5, 0.0]

    # read the sampling input file
    vtk_handler = pyHMT2D.Misc.vtkHandler()

    # Define the numer of interpolation points
    numPoints = 100

    # for SRH-2D
    reader = vtk_handler.readVTK_UnstructuredGrid(srh_filename)  # read the VTKfile
    line = vtk_handler.createVtkLine(p1, p2, numPoints)  # Create the line
    points, WSE_srh_2d, zb_srh_2d = vtk_handler.probeUnstructuredGridVTKOverLine(line, reader, 'Water_Elev_m')  # interpolate the data over the line
    # print(points)
    # print(U)

    max_x = np.max(points[:, 0])
    min_x = np.min(points[:, 0])

    WSE_srh_2d = setNumpyArrayValueToNaN(WSE_srh_2d, 0)  # Set the zero's to NaN's

    # for AdHydraulics
    reader = vtk_handler.readVTK_UnstructuredGrid(AdHydraulics_filename)  # read the VTKfile
    line = vtk_handler.createVtkLine(p1, p2, numPoints)  # Create the line
    points, WSE_AdHydraulics, zb_AdHydraulics = vtk_handler.probeUnstructuredGridVTKOverLine(line, reader, 'WSE')  # interpolate the data over the line
    # print(points)
    # print(U)

    max_x = np.max(points[:, 0])
    min_x = np.min(points[:, 0])

    WSE_AdHydraulics = setNumpyArrayValueToNaN(WSE_AdHydraulics, 0)  # Set the zero's to NaN's

    plot_surface(points[:,0],zb_AdHydraulics, WSE_AdHydraulics, WSE_srh_2d)

    print("Done!")