"""
Some specific utility functions, such as plotting.

"""

import matplotlib.pyplot as plt
import numpy as np

import shapefile
import cv2

import matplotlib.ticker as tick
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter

from matplotlib.tri import Triangulation

from matplotlib.colors import Normalize

import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np

from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

import meshio

import geopandas as gpd

import vtk

import pyvista as pv

import json

from pathlib import Path

#import pyHMT2D
#from pyHMT2D.Misc.tools import setNumpyArrayValueToNaN

plt.rc('text', usetex=True)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def plot_UDE_training_history_ManningN(iter_number, bPlotTitle=False):
    """
    Make plots for the UDE result for each iteration:
    1. loss history
    2. Manning n truth vs NN
    3. Manning n on f-Re plot (truth)
    4. Manning n on f-Re plot (NN)    

    :param iter_number: int
    :param bPlotTitle: bool, whether to plot the title
    :return:
    """

    fig, axs = plt.subplots(2, 2, figsize=(14, 9), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    fig.subplots_adjust(hspace=0.3, wspace=.2)

    #1. plot loss history
    # read the loss data from the "inversion_results_losses_and_ODE_solution.json" file
    with open('UDE_results_losses_and_ODE_solution.json', 'r') as file:
        data = json.load(file)

    loss_preds = data['loss_preds']
    loss_pred_WSEs = data['loss_pred_WSEs']
    loss_pred_uvs = data['loss_pred_uvs']

    iter_numbers = np.arange(len(loss_preds))

    # plot losses (only the first iter_number iterations)
    axs[0, 0].plot(iter_numbers[0:iter_number], loss_preds[0:iter_number], 'k', label=r'$l_{pred, total}$')
    axs[0, 0].plot(iter_numbers[0:iter_number], loss_pred_WSEs[0:iter_number], 'r--', label=r'$l_{pred, WSE}$')
    axs[0, 0].plot(iter_numbers[0:iter_number], loss_pred_uvs[0:iter_number], 'b-.', label=r'$l_{pred, uv}$')    

    axs[0, 0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    axs[0, 0].tick_params(axis='both', which='major', labelsize=14)

    axs[0, 0].set_xlabel('Iteration', fontsize=16)
    axs[0, 0].set_ylabel('Loss', fontsize=16)
    axs[0, 0].set_xlim([0, 150])
    axs[0, 0].set_ylim([5e-8, 4e-3])
    axs[0, 0].set_yscale('log')

    axs[0, 0].legend(loc='upper right', fontsize=14, frameon=False)

    #2. plot Manning's n truth vs Manning's n UDE 
    # read the Manning's n truth and UDE data from the "forward_simulation_results_iter_number.vtk" file
    vtk_file_name = f"UDE_iteration_results_{iter_number:04d}.vtk"
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtk_file_name)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    vtk_data = reader.GetOutput()

    # Get Manning's n truth and UDE data
    ManningN_truth = vtk_to_numpy(vtk_data.GetCellData().GetArray('ManningN_cell_truth'))
    ManningN_UDE = vtk_to_numpy(vtk_data.GetCellData().GetArray('ManningN_cells'))

    # compute RMSE 
    RMSE = np.sqrt(np.mean((ManningN_UDE - ManningN_truth)**2))
    print(f"RMSE = {RMSE:.4f}")

    # Plot Manning's n truth vs Manning's n UDE
    axs[0, 1].scatter(ManningN_UDE, ManningN_truth, edgecolors='blue', facecolors='none', s=50)

    #draw the diagonal line
    axs[0, 1].plot([0.02,0.07],[0.02, 0.07], color='black')

    #add RMSE text
    axs[0, 1].text(0.022, 0.037, f"RMSE = {RMSE:.4f}", fontsize=14)

    axs[0, 1].set_xlim([0.02, 0.04])
    axs[0, 1].set_ylim([0.02, 0.04])

    axs[0, 1].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    axs[0, 1].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    axs[0, 1].set_xlabel('Manning\'s $n$ from USWEs', fontsize=16)
    axs[0, 1].set_ylabel('Manning\'s $n$ truth', fontsize=16)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=14)

    #add caption
    axs[0, 0].text(-0.1, 1.05, "(a)", size=16, ha="center", transform=axs[0, 0].transAxes)  # upper left
    axs[0, 1].text(-0.1, 1.05, "(b)", size=16, ha="center", transform=axs[0, 1].transAxes)
    axs[1, 0].text(-0.2, 1.05, "(c)", size=16, ha="center", transform=axs[1, 0].transAxes)  # upper left
    axs[1, 1].text(-0.1, 1.05, "(d)", size=16, ha="center", transform=axs[1, 1].transAxes)

    #3. plot the f-Re plot for the truth data
    plot_f_h_ks_Re(plt, fig, axs[1, 0], iter_number, bPlotTruth=True)

    #4. plot the f-Re plot for the UDE data
    plot_f_h_ks_Re(plt, fig, axs[1, 1], iter_number, bPlotTruth=False)

    if bPlotTitle:
        plt.suptitle(f"USWEs training history, iter_number = {iter_number:04d}", fontsize=16)
        plt.savefig(f"Savana_River_UDE_result_training_history_ManningN_{iter_number:04d}_with_title.png", dpi=300, bbox_inches='tight', pad_inches=0)
    else:
        plt.savefig(f"Savana_River_UDE_result_training_history_ManningN_{iter_number:04d}.png", dpi=300, bbox_inches='tight', pad_inches=0)

    #plt.show()


def triangulate_vtk(vtk_file_name):
    """
    Triangulate the VTK file.
    """

    #extract the file name from the vtk_file_name
    file_name = vtk_file_name.split('.')[0]

    # Load the VTK file
    file_path = vtk_file_name
    grid = pv.read(file_path)


    # Triangulate the grid
    triangulated = grid.triangulate()

    # Interpolate cell data to point data
    interpolated = triangulated.cell_data_to_point_data()

    # Save the triangulated grid to a new VTK file
    triangulated_file_path = f'{file_name}_tri.vtk'
    interpolated.save(triangulated_file_path)


def plot_f_h_ks_Re(plt, fig, ax, iter_number, bPlotTruth=False):
    """
    Plot the function of f(h/ks, Re) and the scatter for simulation data.


    Reference: Cheng (2008) JHE, "Formulas for Friction Factor in Transitional Regimes", Eq. (17)

    :param plt: matplotlib.pyplot
    :param fig: matplotlib.figure
    :param ax: matplotlib.axes
    :param iter_number: int
    :param bPlotTrue: bool, whether it is for the truth data or the UDE data
    """

    #get simulation data from the VTK file
    vtk_fileName_truth = 'truth_results_0100.vtk'
    vtk_fileName_UDE = f'UDE_iteration_results_{iter_number:04d}.vtk'

     # Read truth data
    reader_truth = vtk.vtkUnstructuredGridReader()
    reader_truth.SetFileName(vtk_fileName_truth)
    reader_truth.ReadAllScalarsOn()
    reader_truth.ReadAllVectorsOn()
    reader_truth.Update()
    vtk_data_truth = reader_truth.GetOutput()

    # Read UDE data
    reader_UDE = vtk.vtkUnstructuredGridReader()
    reader_UDE.SetFileName(vtk_fileName_UDE)
    reader_UDE.ReadAllScalarsOn()
    reader_UDE.ReadAllVectorsOn()
    reader_UDE.Update()
    vtk_data_UDE = reader_UDE.GetOutput()

    # Print array names correctly
    print("\nArrays in truth data:")
    cell_data_truth = vtk_data_truth.GetCellData()
    for i in range(cell_data_truth.GetNumberOfArrays()):
        print(f"  {cell_data_truth.GetArrayName(i)}")

    print("\nArrays in UDE data:")
    cell_data_UDE = vtk_data_UDE.GetCellData()
    for i in range(cell_data_UDE.GetNumberOfArrays()):
        print(f"  {cell_data_UDE.GetArrayName(i)}")

    #get the data
    h_ks_simulation_truth = np.squeeze(vtk_data_truth.GetCellData().GetArray('h_ks'))
    Re_simulation_truth = np.squeeze(vtk_data_truth.GetCellData().GetArray('Re'))
    friction_factor_simulation_truth = np.squeeze(vtk_data_truth.GetCellData().GetArray('friction_factor'))
    n_simulation_truth = np.squeeze(vtk_data_truth.GetCellData().GetArray('ManningN_cells'))

    h_ks_simulation_UDE = np.squeeze(vtk_data_UDE.GetCellData().GetArray('h_ks'))
    Re_simulation_UDE = np.squeeze(vtk_data_UDE.GetCellData().GetArray('Re'))
    friction_factor_simulation_UDE = np.squeeze(vtk_data_UDE.GetCellData().GetArray('friction_factor_from_formula'))
    n_simulation_UDE = np.squeeze(vtk_data_UDE.GetCellData().GetArray('ManningN_cells'))

    #print("h_ks_simulation_truth = ", h_ks_simulation_truth)
    #print("Re_simulation_truth = ", Re_simulation_truth)
    #print("friction_factor_simulation_truth = ", friction_factor_simulation_truth)
    #print("n_simulation_truth = ", n_simulation_truth)

    #print("Re_simulation_truth.shape = ", Re_simulation_truth.shape)
    #print("friction_factor_simulation_truth.shape = ", friction_factor_simulation_truth.shape)
    #print("n_simulation_truth.shape = ", n_simulation_truth.shape)

    #print("Re_simulation_UDE.shape = ", Re_simulation_UDE.shape)
    #print("friction_factor_simulation_UDE.shape = ", friction_factor_simulation_UDE.shape)
    #print("n_simulation_UDE.shape = ", n_simulation_UDE.shape)

    #define the range of h/ks and Re
    h_ks = [2.0, 5.0, 10.0, 20.0, 50.0, 100.0]

    #number of points on each line
    n_points = 200

    #Create logarithmically spaced points for Re
    Re = np.logspace(2, 8, n_points)  # 100 points from 10^2 to 10^8

    alpha = 1.0/(1.0  + (Re/850.0)**9)    

    #compute f values for each h/ks
    f = []
    for i in range(len(h_ks)):
        h_ks_i = h_ks[i]

        f_i = np.zeros(n_points)
        part1 = np.zeros(n_points)
        part2 = np.zeros(n_points)
        part3 = np.zeros(n_points)

        #compute in parts
        for iPoint in range(n_points):
            beta_iPoint = 1.0/(1.0 + (Re[iPoint]/(160*h_ks_i))**2)

            part1[iPoint] = (Re[iPoint]/24)**alpha[iPoint]
            part2[iPoint] = (1.8*np.log10(Re[iPoint]/2.1))**(2*(1-alpha[iPoint])*beta_iPoint)
            part3[iPoint] = (2*np.log10(11.8*h_ks_i))**(2*(1-alpha[iPoint])*(1-beta_iPoint))

            f_i[iPoint] = 1.0/(part1[iPoint] * part2[iPoint] * part3[iPoint])

        f.append(f_i)

    # Create normalizer for the color mapping
    h_ks_min = 2
    h_ks_max = 100
    norm = Normalize(vmin=h_ks_min, vmax=h_ks_max)

    #loop over the h_ks
    for i in range(len(h_ks)):
        color = plt.cm.jet(norm(h_ks[i]))  # Get color from colormap
        ax.plot(Re, f[i], color=color, label=f'h/ks = {h_ks[i]}')

    #plot the scatter for simulation data and colored by h/ks
    if bPlotTruth:
        scatter = ax.scatter(Re_simulation_truth, friction_factor_simulation_truth, c=h_ks_simulation_truth, 
                             alpha=0.5, edgecolors='black', linewidths=0.5, norm=norm, cmap='jet', label='Truth data')
        
        #add text annotation to the scatter
        ax.text(1E3, 0.17, 'Truth', fontsize=14)
    else:
        scatter = ax.scatter(Re_simulation_UDE, friction_factor_simulation_UDE, c=h_ks_simulation_UDE, 
                             alpha=0.5, edgecolors='black', linewidths=0.5, norm=norm, cmap='jet', label='UDE data')
        
        #add text annotation to the scatter
        ax.text(1E3, 0.17, 'USWEs', fontsize=14)

    #add a color bar and customize the color bar
    cbar = plt.colorbar(scatter, location='right', pad=0.05)
    cbar.ax.set_title("$h/k_s$", loc='center', fontsize=16)
    cbar.ax.tick_params(labelsize=14)

    #axs.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')

    #set axis limits
    ax.set_xlim(100, 2E7)
    ax.set_ylim(0.02, 0.2)

    #set aspect ratio 
    #ax.set_aspect('equal')
    ax.set_box_aspect(0.8)

    ax.set_xlabel('$Re$', fontsize=16)
    ax.set_ylabel('$f$', fontsize=16)

    ax.tick_params(axis='x', which='both', labelsize=14)
    ax.tick_params(axis='y', which='both', labelsize=14)

    #set y tick format
    ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    ax.yaxis.set_minor_formatter(tick.FormatStrFormatter('%.2f'))

    #add text to annotate the h/ks values for each line
    for i in range(len(h_ks)):
        ax.text(2E4, f[i][-1], f'$h/k_s$ = {h_ks[i]:.0f}', fontsize=12)

def create_animation_training_history():
    """
    Create an animation of the training history ManningN by combining all the plots.
    """

    # Read the first image to get the dimensions
    first_image_path = "Savana_River_UDE_result_training_history_ManningN_0001_with_title.png"
    frame = cv2.imread(first_image_path)
    height, width, layers = frame.shape

    # Define the video codec and create a VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Codec for .mp4 files
    video = cv2.VideoWriter('Savana_River_UDE_result_training_history_ManningN.mp4', fourcc, 10.0, (width, height))  

    #create a list to store the images
    images = []
    for iImage in range(1, 151):
        image_path = f"Savana_River_UDE_result_training_history_ManningN_{iImage:04d}_with_title.png"
        image = cv2.imread(image_path)
        video.write(image)

    # Release the VideoWriter
    video.release()

    print(f"Video saved as Savana_River_UDE_result_training_history_ManningN.mp4")
    

def create_animation_ManningN_iterations():
    """
    Create an animation of the Manning's n iterations.
    """

    # Read the first image to get the dimensions
    first_image_path = "Savana_UDE_ManningN_iteration_0001.png"
    frame = cv2.imread(first_image_path)
    height, width, layers = frame.shape

    # Define the video codec and create a VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Codec for .mp4 files
    video = cv2.VideoWriter('Savana_UDE_ManningN_iterations.mp4', fourcc, 10.0, (width, height))  

    #create a list to store the images
    images = []
    for iImage in range(1, 151):
        image_path = f"Savana_UDE_ManningN_iteration_{iImage:04d}.png"
        image = cv2.imread(image_path)
        video.write(image)

    # Release the VideoWriter
    video.release()

    print(f"Video saved as Savana_UDE_ManningN_iterations.mp4")
    
def plot_ks_contour():

    """
    Plot the ks contour
    """

    fig, ax = plt.subplots(1, 1, figsize=(8, 6), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    # Read VTK file
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName("UDE_iteration_results_0100.vtk")
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    grid = reader.GetOutput()

    # Get points/coordinates
    points = vtk_to_numpy(grid.GetPoints().GetData())

    xCoord = points[:, 0]
    yCoord = points[:, 1]

    xl = xCoord.min()
    xh = xCoord.max()
    yl = yCoord.min()
    yh = yCoord.max()

    ks_min = 0.1
    ks_max = 0.3

    local_levels = [0.1, 0.15, 0.2, 0.25, 0.3]
    
    # Get cell data (ks field)
    ks = vtk_to_numpy(grid.GetCellData().GetArray('ks'))
    
    # Create figure    
    polygons = []
    
    # Loop through all cells
    for i in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(i)
        n_points = cell.GetNumberOfPoints()
        
        # Get vertex indices for this cell
        vertex_indices = [cell.GetPointId(j) for j in range(n_points)]
        
        # Get vertex coordinates
        vertices = points[vertex_indices]
        
        # Create polygon
        polygon = Polygon(vertices[:, :2])  # Only use x,y coordinates
        polygons.append(polygon)
    
    # Create PatchCollection with the cell data
    p = PatchCollection(polygons, cmap='jet', alpha=1.0)
    p.set_array(ks)
    
    # Add to plot
    ax.add_collection(p)

    ax.set_xlim([xl, xh])
    ax.set_ylim([yl, yh])

    ax.set_aspect('equal')

    ax.set_xlabel('$x$ (m)', fontsize=16)
    ax.tick_params(axis='x', labelsize=14)
    ax.set_ylabel('$y$ (m)', fontsize=16)
    ax.tick_params(axis='y', labelsize=14)


    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.2)
    clb = fig.colorbar(p, ticks=local_levels, cax=cax)
    clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    clb.ax.tick_params(labelsize=14)
    clb.ax.set_title("(m)", loc='center', fontsize=14)

    #add annotation to ks zones
    # Add annotation
    ax.annotate(
        '$k_{s,1}$ = 0.25',  # Text to display
        xy=(380, 300),  # Point to annotate
        xytext=(400, 350),  # Location of the text
        arrowprops=dict(
            facecolor='black',  # Arrow color
            arrowstyle='->',  # Arrow style
            lw=1  # Line width
        ),
        fontsize=14,  # Text size
        color='black'  # Text color
    )

    ax.annotate(
        '$k_{s,2}$ = 0.3',  # Text to display
        xy=(800, 210),  # Point to annotate
        xytext=(670, 300),  # Location of the text
        arrowprops=dict(

            facecolor='black',  # Arrow color
            arrowstyle='->',  # Arrow style
            lw=1  # Line width
        ),
        fontsize=14,  # Text size
        color='black'  # Text color
    )

    ax.annotate(
        '$k_{s,3}$ = 0.1',  # Text to display
        xy=(700, 70),  # Point to annotate
        xytext=(780, 20),  # Location of the text
        arrowprops=dict(

            facecolor='black',  # Arrow color
            arrowstyle='->',  # Arrow style
            lw=1  # Line width
        ),
        fontsize=14,  # Text size
        color='black'  # Text color
    )

    ax.annotate(
        '$k_{s,4}$ = 0.2',  # Text to display
        xy=(900, 200),  # Point to annotate
        xytext=(820, 100),  # Location of the text
        arrowprops=dict(

            facecolor='black',  # Arrow color
            arrowstyle='->',  # Arrow style
            lw=1  # Line width
        ),
        fontsize=14,  # Text size
        color='black'  # Text color
    )

    ax.annotate(
        '$k_{s,5}$ = 0.3',  # Text to display
        xy=(300, 280),  # Point to annotate
        xytext=(150, 200),  # Location of the text
        arrowprops=dict(

            facecolor='black',  # Arrow color
            arrowstyle='->',  # Arrow style
            lw=1  # Line width
        ),
        fontsize=14,  # Text size
        color='black'  # Text color
    )
   
    #read the domain outer boundary
    # Path to the shapefile
    shapefile_path = 'domain_outer_boundary.shp'

    # Read the shapefile
    gdf = gpd.read_file(shapefile_path)

    #plot the outer boundary
    for geometry in gdf.geometry:
        if geometry.geom_type == 'Polygon':
            # For a single polygon
            x, y = geometry.exterior.xy
            ax.plot(x, y, color='blue', linewidth=0.5)
            ax.fill(x, y, color='lightblue', alpha=0.3)

    plt.savefig("Savana_UDE_ks_contour.png", dpi=300, bbox_inches='tight', pad_inches=0)
    
    #plt.show()


def plot_ManningN_one_iteration(iter_number):
    """
    Plot the Manning's n at one iteration.

    It plots the UDE Manning's n, the truth Manning's n, and the difference between the UDE and the truth Manning's n.
    """

    fig, axs = plt.subplots(1, 3, figsize=(16, 3), sharex=False, sharey=True, facecolor='w', edgecolor='k')

    #set the spacing between the subplots
    plt.subplots_adjust(wspace=0.2, hspace=0.2)

   
    print(f"Plotting Manning's n, iter_number = {iter_number:04d}")

    # read the Manning's n values from the vtk file
    vtk_file_name = f"UDE_iteration_results_{iter_number:04d}.vtk"
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtk_file_name)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.Update()
    vtk_data = reader.GetOutput()

    # Get points/coordinates
    points = vtk_to_numpy(vtk_data.GetPoints().GetData())

    xCoord = points[:, 0]
    yCoord = points[:, 1]

    xl = xCoord.min()
    xh = xCoord.max()
    yl = yCoord.min()
    yh = yCoord.max()

    # Get Manning's n values
    ManningN_truth = vtk_to_numpy(vtk_data.GetCellData().GetArray('ManningN_cell_truth'))
    ManningN_UDE = vtk_to_numpy(vtk_data.GetCellData().GetArray('ManningN_cells'))

    #ManningN_min = 0.01
    #ManningN_max = 0.05
    #local_levels = [0.01, 0.02, 0.03, 0.04, 0.05]
    ManningN_min = ManningN_truth.min()
    ManningN_max = ManningN_truth.max()
    local_levels = np.linspace(ManningN_min, ManningN_max, 6)
    norm_ManningN = Normalize(vmin=ManningN_min, vmax=ManningN_max)

    #compute the difference
    diff = ManningN_truth - ManningN_UDE

    print("diff min and max = ", diff.min(), diff.max())

    diff_min = -0.006
    diff_max = 0.006
    local_diff_levels = np.linspace(diff_min, diff_max, 7)
    norm_diff = Normalize(vmin=diff_min, vmax=diff_max)

    #compute rmse
    rmse = np.sqrt(np.mean(diff**2))

    # Create figure for ManningN_truth
    polygons = []        

    # Loop through all cells
    for i in range(vtk_data.GetNumberOfCells()):
        cell = vtk_data.GetCell(i)
        n_points = cell.GetNumberOfPoints()            

        # Get vertex indices for this cell
        vertex_indices = [cell.GetPointId(j) for j in range(n_points)]
            
        # Get vertex coordinates
        vertices = points[vertex_indices]
            
        # Create polygon
        polygon = Polygon(vertices[:, :2])  # Only use x,y coordinates
        polygons.append(polygon)
    
    # Create PatchCollection with the cell data
    p_truth = PatchCollection(polygons, cmap='jet', norm=norm_ManningN, alpha=1.0)
    p_UDE = PatchCollection(polygons, cmap='jet', norm=norm_ManningN, alpha=1.0)
    p_diff = PatchCollection(polygons, cmap='jet', norm=norm_diff, alpha=1.0)        

    p_truth.set_array(ManningN_truth)
    p_UDE.set_array(ManningN_UDE)
    p_diff.set_array(diff)    

    # Add to plot
    axs[0].add_collection(p_truth)
    axs[1].add_collection(p_UDE)
    axs[2].add_collection(p_diff)

    axs[0].set_xlim([xl, xh])
    axs[0].set_ylim([yl, yh])
    axs[0].set_aspect('equal')

    axs[1].set_xlim([xl, xh])
    axs[1].set_ylim([yl, yh])
    axs[1].set_aspect('equal')

    axs[2].set_xlim([xl, xh])
    axs[2].set_ylim([yl, yh])
    axs[2].set_aspect('equal')

    axs[0].set_xlabel('$x$ (m)', fontsize=16)        
    axs[0].tick_params(axis='x', labelsize=14)
    axs[0].set_ylabel('$y$ (m)', fontsize=16)
    axs[0].tick_params(axis='y', labelsize=14)

    axs[1].set_xlabel('$x$ (m)', fontsize=16)
    axs[1].tick_params(axis='x', labelsize=14)
    #axs[1].set_ylabel('$y$ (m)', fontsize=16)
    axs[1].tick_params(axis='y', labelsize=14)

    axs[2].set_xlabel('$x$ (m)', fontsize=16)
    axs[2].tick_params(axis='x', labelsize=14)
    #axs[2].set_ylabel('$y$ (m)', fontsize=16)
    axs[2].tick_params(axis='y', labelsize=14)

    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="3%", pad=0.2)
    clb = fig.colorbar(p_truth, ticks=local_levels, cax=cax)
    clb.set_ticks(local_levels)
    clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
    clb.ax.tick_params(labelsize=14)
    clb.ax.set_title("(s/m$^{1/3}$)", loc='center', fontsize=14)        
    
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("right", size="3%", pad=0.2)
    clb = fig.colorbar(p_UDE, ticks=local_levels, cax=cax)
    clb.set_ticks(local_levels)
    clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
    clb.ax.tick_params(labelsize=14)
    clb.ax.set_title("(s/m$^{1/3}$)", loc='center', fontsize=14)      

    divider = make_axes_locatable(axs[2])
    cax = divider.append_axes("right", size="3%", pad=0.2)
    clb = fig.colorbar(p_diff, ticks=local_diff_levels, cax=cax)
    clb.set_ticks(local_diff_levels)
    clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
    clb.ax.tick_params(labelsize=14)

    clb.ax.set_title("(s/m$^{1/3}$)", loc='center', fontsize=14)  

    axs[2].text(400, 350, f"RMSE = {rmse:.5f}", fontsize=12, ha='left', va='center')

    #add title above the first row
    axs[0].set_title("Truth", fontsize=16)
    axs[1].set_title("From USWEs", fontsize=16)
    axs[2].set_title("Difference", fontsize=16)
    
    #add text for iter_number
    fig.suptitle(f"Manning's n, Iteration = {iter_number}", fontsize=16)
    
    
    plt.savefig(f"Savana_UDE_ManningN_iteration_{iter_number:04d}.png", dpi=300, bbox_inches='tight', pad_inches=0)



def plot_ManningN_iterations(iter_numbers):
    """
    Plot the Manning's n at selected iterations.

    For each row, it plots the UDE Manning's n, the truth Manning's n, and the difference between the UDE and the truth Manning's n.
    """

    fig, axs = plt.subplots(4, 3, figsize=(16, 9), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    #set the spacing between the subplots
    plt.subplots_adjust(wspace=0.2, hspace=0.2)

    # loop over the iter_numbers
    for index, iter_number in enumerate(iter_numbers):
        print(f"Plotting Manning's n, iter_number = {iter_number:04d}")

        # read the Manning's n values from the vtk file
        vtk_file_name = f"UDE_iteration_results_{iter_number:04d}.vtk"
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(vtk_file_name)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        vtk_data = reader.GetOutput()

        # Get points/coordinates
        points = vtk_to_numpy(vtk_data.GetPoints().GetData())

        xCoord = points[:, 0]
        yCoord = points[:, 1]

        xl = xCoord.min()
        xh = xCoord.max()
        yl = yCoord.min()
        yh = yCoord.max()

        # Get Manning's n values
        ManningN_truth = vtk_to_numpy(vtk_data.GetCellData().GetArray('ManningN_cell_truth'))
        ManningN_UDE = vtk_to_numpy(vtk_data.GetCellData().GetArray('ManningN_cells'))

        #ManningN_min = 0.01
        #ManningN_max = 0.05
        #local_levels = [0.01, 0.02, 0.03, 0.04, 0.05]
        ManningN_min = ManningN_truth.min()
        ManningN_max = ManningN_truth.max()
        local_levels = np.linspace(ManningN_min, ManningN_max, 6)
        norm_ManningN = Normalize(vmin=ManningN_min, vmax=ManningN_max)

        #compute the difference
        diff = ManningN_truth - ManningN_UDE

        print("diff min and max = ", diff.min(), diff.max())

        diff_min = -0.006
        diff_max = 0.006
        local_diff_levels = np.linspace(diff_min, diff_max, 7)
        norm_diff = Normalize(vmin=diff_min, vmax=diff_max)

        #compute rmse
        rmse = np.sqrt(np.mean(diff**2))

        # Create figure for ManningN_truth
        polygons = []        

        # Loop through all cells
        for i in range(vtk_data.GetNumberOfCells()):
            cell = vtk_data.GetCell(i)
            n_points = cell.GetNumberOfPoints()            

            # Get vertex indices for this cell
            vertex_indices = [cell.GetPointId(j) for j in range(n_points)]
            
            # Get vertex coordinates
            vertices = points[vertex_indices]
            
            # Create polygon
            polygon = Polygon(vertices[:, :2])  # Only use x,y coordinates
            polygons.append(polygon)
    
        # Create PatchCollection with the cell data
        p_truth = PatchCollection(polygons, cmap='jet', norm=norm_ManningN, alpha=1.0)
        p_UDE = PatchCollection(polygons, cmap='jet', norm=norm_ManningN, alpha=1.0)
        p_diff = PatchCollection(polygons, cmap='jet', norm=norm_diff, alpha=1.0)        

        p_truth.set_array(ManningN_truth)
        p_UDE.set_array(ManningN_UDE)
        p_diff.set_array(diff)    

        # Add to plot
        axs[index, 0].add_collection(p_truth)
        axs[index, 1].add_collection(p_UDE)
        axs[index, 2].add_collection(p_diff)

        axs[index, 0].set_xlim([xl, xh])
        axs[index, 0].set_ylim([yl, yh])
        axs[index, 0].set_aspect('equal')

        axs[index, 1].set_xlim([xl, xh])
        axs[index, 1].set_ylim([yl, yh])
        axs[index, 1].set_aspect('equal')

        axs[index, 2].set_xlim([xl, xh])
        axs[index, 2].set_ylim([yl, yh])
        axs[index, 2].set_aspect('equal')

        if index == 3: axs[index, 0].set_xlabel('$x$ (m)', fontsize=16)        
        axs[index, 0].tick_params(axis='x', labelsize=14)
        axs[index, 0].set_ylabel('$y$ (m)', fontsize=16)
        axs[index, 0].tick_params(axis='y', labelsize=14)

        if index == 3: axs[index, 1].set_xlabel('$x$ (m)', fontsize=16)
        axs[index, 1].tick_params(axis='x', labelsize=14)
        #axs[index, 1].set_ylabel('$y$ (m)', fontsize=16)
        axs[index, 1].tick_params(axis='y', labelsize=14)

        if index == 3: axs[index, 2].set_xlabel('$x$ (m)', fontsize=16)
        axs[index, 2].tick_params(axis='x', labelsize=14)
        #axs[index, 2].set_ylabel('$y$ (m)', fontsize=16)
        axs[index, 2].tick_params(axis='y', labelsize=14)

        divider = make_axes_locatable(axs[index, 0])
        cax = divider.append_axes("right", size="3%", pad=0.2)
        clb = fig.colorbar(p_truth, ticks=local_levels, cax=cax)
        clb.set_ticks(local_levels)
        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
        clb.ax.tick_params(labelsize=14)
        clb.ax.set_title("(s/m$^{1/3}$)", loc='center', fontsize=14)        
    
        divider = make_axes_locatable(axs[index, 1])
        cax = divider.append_axes("right", size="3%", pad=0.2)
        clb = fig.colorbar(p_UDE, ticks=local_levels, cax=cax)
        clb.set_ticks(local_levels)
        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
        clb.ax.tick_params(labelsize=14)
        clb.ax.set_title("(s/m$^{1/3}$)", loc='center', fontsize=14)      

        divider = make_axes_locatable(axs[index, 2])
        cax = divider.append_axes("right", size="3%", pad=0.2)
        clb = fig.colorbar(p_diff, ticks=local_diff_levels, cax=cax)
        clb.set_ticks(local_diff_levels)
        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
        clb.ax.tick_params(labelsize=14)

        clb.ax.set_title("(s/m$^{1/3}$)", loc='center', fontsize=14)  

        axs[index, 2].text(400, 350, f"RMSE = {rmse:.5f}", fontsize=12, ha='left', va='center')

    #add title above the first row
    axs[0, 0].set_title("Truth", fontsize=16)
    axs[0, 1].set_title("From USWEs", fontsize=16)
    axs[0, 2].set_title("Difference", fontsize=16)

    #add text for iter_number
    axs[0, 1].text(400, 420, f"iter = {iter_numbers[0]}", fontsize=14, ha='left', va='center')
    axs[1, 1].text(400, 420, f"iter = {iter_numbers[1]}", fontsize=14, ha='left', va='center')
    axs[2, 1].text(400, 420, f"iter = {iter_numbers[2]}", fontsize=14, ha='left', va='center')
    axs[3, 1].text(400, 420, f"iter = {iter_numbers[3]}", fontsize=14, ha='left', va='center')

    axs[0, 2].text(400, 420, f"iter = {iter_numbers[0]}", fontsize=14, ha='left', va='center')
    axs[1, 2].text(400, 420, f"iter = {iter_numbers[1]}", fontsize=14, ha='left', va='center')
    axs[2, 2].text(400, 420, f"iter = {iter_numbers[2]}", fontsize=14, ha='left', va='center')
    axs[3, 2].text(400, 420, f"iter = {iter_numbers[3]}", fontsize=14, ha='left', va='center')

    plt.savefig("Savana_UDE_ManningN_iterations.png", dpi=300, bbox_inches='tight', pad_inches=0)


def plot_WSE_iterations(iter_numbers):
    """
    Plot the WSE at selected iterations.

    For each row, it plots the UDE WSE, the truth WSE, and the difference between the UDE and the truth WSE.
    """

    fig, axs = plt.subplots(4, 3, figsize=(16, 9), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    #set the spacing between the subplots
    plt.subplots_adjust(wspace=0.2, hspace=0.2)    

    # loop over the iter_numbers
    for index, iter_number in enumerate(iter_numbers):
        print(f"Plotting Manning's n, iter_number = {iter_number:04d}")

        # read the Manning's n values from the vtk files
        vtk_file_name_UDE = f"UDE_iteration_results_{iter_number:04d}.vtk"
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(vtk_file_name_UDE)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()
        vtk_data_UDE = reader.GetOutput()

        # Get points/coordinates
        points = vtk_to_numpy(vtk_data_UDE.GetPoints().GetData())

        xCoord = points[:, 0]
        yCoord = points[:, 1]

        xl = xCoord.min()
        xh = xCoord.max()
        yl = yCoord.min()
        yh = yCoord.max()

        # Get Manning's n values
        WSE_UDE = vtk_to_numpy(vtk_data_UDE.GetCellData().GetArray('WSE'))
        WSE_truth = vtk_to_numpy(vtk_data_UDE.GetCellData().GetArray('WSE_truth'))
      
        WSE_min = WSE_truth.min()
        WSE_max = WSE_truth.max()
        local_levels = np.linspace(WSE_min, WSE_max, 6)

        norm_WSE = Normalize(vmin=WSE_min, vmax=WSE_max)

        #compute the difference
        diff = WSE_truth - WSE_UDE

        print("diff min and max = ", diff.min(), diff.max())

        diff_min = -0.03
        diff_max = 0.03
        local_diff_levels = np.linspace(diff_min, diff_max, 7)
        norm_diff = Normalize(vmin=diff_min, vmax=diff_max)

        #compute rmse
        rmse = np.sqrt(np.mean(diff**2))

        # Create figure for ManningN_truth
        polygons = []        

        # Loop through all cells
        for i in range(vtk_data_UDE.GetNumberOfCells()):
            cell = vtk_data_UDE.GetCell(i)
            n_points = cell.GetNumberOfPoints()     

            # Get vertex indices for this cell
            vertex_indices = [cell.GetPointId(j) for j in range(n_points)]
            
            # Get vertex coordinates
            vertices = points[vertex_indices]
            
            # Create polygon
            polygon = Polygon(vertices[:, :2])  # Only use x,y coordinates
            polygons.append(polygon)
    
        # Create PatchCollection with the cell data
        p_truth = PatchCollection(polygons, cmap='jet', norm=norm_WSE, alpha=1.0)
        p_UDE = PatchCollection(polygons, cmap='jet', norm=norm_WSE, alpha=1.0)
        p_diff = PatchCollection(polygons, cmap='jet', norm=norm_diff, alpha=1.0)        

        p_truth.set_array(WSE_truth)
        p_UDE.set_array(WSE_UDE)
        p_diff.set_array(diff)    

        # Add to plot
        axs[index, 0].add_collection(p_truth)
        axs[index, 1].add_collection(p_UDE)
        axs[index, 2].add_collection(p_diff)

        axs[index, 0].set_xlim([xl, xh])
        axs[index, 0].set_ylim([yl, yh])
        axs[index, 0].set_aspect('equal')

        axs[index, 1].set_xlim([xl, xh])
        axs[index, 1].set_ylim([yl, yh])
        axs[index, 1].set_aspect('equal')

        axs[index, 2].set_xlim([xl, xh])
        axs[index, 2].set_ylim([yl, yh])
        axs[index, 2].set_aspect('equal')

        if index == 3: axs[index, 0].set_xlabel('$x$ (m)', fontsize=16)        
        axs[index, 0].tick_params(axis='x', labelsize=14)
        axs[index, 0].set_ylabel('$y$ (m)', fontsize=16)
        axs[index, 0].tick_params(axis='y', labelsize=14)

        if index == 3: axs[index, 1].set_xlabel('$x$ (m)', fontsize=16)
        axs[index, 1].tick_params(axis='x', labelsize=14)
        #axs[index, 1].set_ylabel('$y$ (m)', fontsize=16)
        axs[index, 1].tick_params(axis='y', labelsize=14)

        if index == 3: axs[index, 2].set_xlabel('$x$ (m)', fontsize=16)
        axs[index, 2].tick_params(axis='x', labelsize=14)
        #axs[index, 2].set_ylabel('$y$ (m)', fontsize=16)
        axs[index, 2].tick_params(axis='y', labelsize=14)

        divider = make_axes_locatable(axs[index, 0])
        cax = divider.append_axes("right", size="3%", pad=0.2)
        clb = fig.colorbar(p_truth, ticks=local_levels, cax=cax)
        clb.set_ticks(local_levels)
        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
        clb.ax.tick_params(labelsize=14)
        clb.ax.set_title("WSE (m)", loc='center', fontsize=14)            

        divider = make_axes_locatable(axs[index, 1])
        cax = divider.append_axes("right", size="3%", pad=0.2)
        clb = fig.colorbar(p_UDE, ticks=local_levels, cax=cax)
        clb.set_ticks(local_levels)
        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
        clb.ax.tick_params(labelsize=14)
        clb.ax.set_title("WSE (m)", loc='center', fontsize=14)      

        divider = make_axes_locatable(axs[index, 2])
        cax = divider.append_axes("right", size="3%", pad=0.2)
        clb = fig.colorbar(p_diff, ticks=local_diff_levels, cax=cax)
        clb.set_ticks(local_diff_levels)
        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
        clb.ax.tick_params(labelsize=14)

        clb.ax.set_title("Diff (m)", loc='center', fontsize=14)  

        axs[index, 2].text(400, 350, f"RMSE = {rmse:.5f} m", fontsize=12, ha='left', va='center')

    #add title above the first row
    axs[0, 0].set_title("Truth", fontsize=16)
    axs[0, 1].set_title("From USWEs", fontsize=16)
    axs[0, 2].set_title("Difference", fontsize=16)

    #add text for iter_number
    axs[0, 1].text(400, 420, f"iter = {iter_numbers[0]}", fontsize=14, ha='left', va='center')
    axs[1, 1].text(400, 420, f"iter = {iter_numbers[1]}", fontsize=14, ha='left', va='center')
    axs[2, 1].text(400, 420, f"iter = {iter_numbers[2]}", fontsize=14, ha='left', va='center')
    axs[3, 1].text(400, 420, f"iter = {iter_numbers[3]}", fontsize=14, ha='left', va='center')

    axs[0, 2].text(400, 420, f"iter = {iter_numbers[0]}", fontsize=14, ha='left', va='center')
    axs[1, 2].text(400, 420, f"iter = {iter_numbers[1]}", fontsize=14, ha='left', va='center')
    axs[2, 2].text(400, 420, f"iter = {iter_numbers[2]}", fontsize=14, ha='left', va='center')
    axs[3, 2].text(400, 420, f"iter = {iter_numbers[3]}", fontsize=14, ha='left', va='center')

    plt.savefig("Savana_UDE_WSE_iterations.png", dpi=300, bbox_inches='tight', pad_inches=0)


def plot_Velocity_iterations(iter_numbers, bPlotU=True):
    """
    Plot the velocity at selected iterations.


    For each row, it plots the UDE velocity, the truth velocity, and the difference between the UDE and the truth velocity.
    """

    fig, axs = plt.subplots(4, 3, figsize=(16, 9), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    #set the spacing between the subplots
    plt.subplots_adjust(wspace=0.2, hspace=0.2)    

    # loop over the iter_numbers
    for index, iter_number in enumerate(iter_numbers):
        print(f"Plotting velocity, iter_number = {iter_number:04d}")

        # read the velocity values from the vtk files
        vtk_file_name_UDE = f"UDE_iteration_results_{iter_number:04d}.vtk"
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(vtk_file_name_UDE)
        reader.ReadAllScalarsOn()
        reader.ReadAllVectorsOn()
        reader.Update()

        vtk_data_UDE = reader.GetOutput()

        # Get points/coordinates
        points = vtk_to_numpy(vtk_data_UDE.GetPoints().GetData())

        xCoord = points[:, 0]
        yCoord = points[:, 1]

        xl = xCoord.min()
        xh = xCoord.max()
        yl = yCoord.min()
        yh = yCoord.max()

        # Get velocity values
        u_UDE = vtk_to_numpy(vtk_data_UDE.GetCellData().GetArray('U'))[:, 0]
        v_UDE = vtk_to_numpy(vtk_data_UDE.GetCellData().GetArray('U'))[:, 1]

        u_truth = vtk_to_numpy(vtk_data_UDE.GetCellData().GetArray('u_truth'))
        v_truth = vtk_to_numpy(vtk_data_UDE.GetCellData().GetArray('v_truth'))      

        u_min = u_truth.min()
        u_max = u_truth.max()
        local_levels_u = np.linspace(u_min, u_max, 6)

        v_min = v_truth.min()
        v_max = v_truth.max()
        local_levels_v = np.linspace(v_min, v_max, 6)

        norm_u = Normalize(vmin=u_min, vmax=u_max)
        norm_v = Normalize(vmin=v_min, vmax=v_max)

        #compute the difference
        diff_u = u_truth - u_UDE
        diff_v = v_truth - v_UDE

        print("diff_u min and max = ", diff_u.min(), diff_u.max())
        print("diff_v min and max = ", diff_v.min(), diff_v.max())

        diff_min_u = -0.1
        diff_max_u = 0.1
        local_diff_levels_u = np.linspace(diff_min_u, diff_max_u, 7)
        norm_diff_u = Normalize(vmin=diff_min_u, vmax=diff_max_u)

        diff_min_v = -0.1
        diff_max_v = 0.1
        local_diff_levels_v = np.linspace(diff_min_v, diff_max_v, 7)
        norm_diff_v = Normalize(vmin=diff_min_v, vmax=diff_max_v)

        #compute rmse
        rmse_u = np.sqrt(np.mean(diff_u**2))
        rmse_v = np.sqrt(np.mean(diff_v**2))

        # Create figure for ManningN_truth
        polygons = []        

        # Loop through all cells
        for i in range(vtk_data_UDE.GetNumberOfCells()):
            cell = vtk_data_UDE.GetCell(i)
            n_points = cell.GetNumberOfPoints()     

            # Get vertex indices for this cell
            vertex_indices = [cell.GetPointId(j) for j in range(n_points)]
            
            # Get vertex coordinates
            vertices = points[vertex_indices]
            
            # Create polygon
            polygon = Polygon(vertices[:, :2])  # Only use x,y coordinates
            polygons.append(polygon)
    
        # Create PatchCollection with the cell data
        p_truth_u = PatchCollection(polygons, cmap='jet', norm=norm_u, alpha=1.0)
        p_UDE_u = PatchCollection(polygons, cmap='jet', norm=norm_u, alpha=1.0)
        p_diff_u = PatchCollection(polygons, cmap='jet', norm=norm_diff_u, alpha=1.0)        

        p_truth_v = PatchCollection(polygons, cmap='jet', norm=norm_v, alpha=1.0)
        p_UDE_v = PatchCollection(polygons, cmap='jet', norm=norm_v, alpha=1.0)
        p_diff_v = PatchCollection(polygons, cmap='jet', norm=norm_diff_v, alpha=1.0)        

        p_truth_u.set_array(u_truth)
        p_UDE_u.set_array(u_UDE)
        p_diff_u.set_array(diff_u)    

        p_truth_v.set_array(v_truth)
        p_UDE_v.set_array(v_UDE)
        p_diff_v.set_array(diff_v)    

        # Add to plot
        if bPlotU:
            axs[index, 0].add_collection(p_truth_u)
            axs[index, 1].add_collection(p_UDE_u)
            axs[index, 2].add_collection(p_diff_u)
        else:
            axs[index, 0].add_collection(p_truth_v)
            axs[index, 1].add_collection(p_UDE_v)
            axs[index, 2].add_collection(p_diff_v)

        axs[index, 0].set_xlim([xl, xh])
        axs[index, 0].set_ylim([yl, yh])
        axs[index, 0].set_aspect('equal')

        axs[index, 1].set_xlim([xl, xh])
        axs[index, 1].set_ylim([yl, yh])
        axs[index, 1].set_aspect('equal')

        axs[index, 2].set_xlim([xl, xh])
        axs[index, 2].set_ylim([yl, yh])
        axs[index, 2].set_aspect('equal')

        if index == 3: axs[index, 0].set_xlabel('$x$ (m)', fontsize=16)        
        axs[index, 0].tick_params(axis='x', labelsize=14)
        axs[index, 0].set_ylabel('$y$ (m)', fontsize=16)
        axs[index, 0].tick_params(axis='y', labelsize=14)

        if index == 3: axs[index, 1].set_xlabel('$x$ (m)', fontsize=16)
        axs[index, 1].tick_params(axis='x', labelsize=14)
        #axs[index, 1].set_ylabel('$y$ (m)', fontsize=16)
        axs[index, 1].tick_params(axis='y', labelsize=14)

        if index == 3: axs[index, 2].set_xlabel('$x$ (m)', fontsize=16)
        axs[index, 2].tick_params(axis='x', labelsize=14)
        #axs[index, 2].set_ylabel('$y$ (m)', fontsize=16)
        axs[index, 2].tick_params(axis='y', labelsize=14)

        divider = make_axes_locatable(axs[index, 0])
        cax = divider.append_axes("right", size="3%", pad=0.2)

        if bPlotU:
            clb = fig.colorbar(p_truth_u, ticks=local_levels_u, cax=cax)
            clb.set_ticks(local_levels_u)
        else:
            clb = fig.colorbar(p_truth_v, ticks=local_levels_v, cax=cax)
            clb.set_ticks(local_levels_v)

        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
        clb.ax.tick_params(labelsize=14)

        if bPlotU:
            clb.ax.set_title("u (m/s)", loc='center', fontsize=14)            
        else:
            clb.ax.set_title("v (m/s)", loc='center', fontsize=14)            


        divider = make_axes_locatable(axs[index, 1])
        cax = divider.append_axes("right", size="3%", pad=0.2)

        if bPlotU:
            clb = fig.colorbar(p_UDE_u, ticks=local_levels_u, cax=cax)
            clb.set_ticks(local_levels_u)
        else:
            clb = fig.colorbar(p_UDE_v, ticks=local_levels_v, cax=cax)
            clb.set_ticks(local_levels_v)

        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
        clb.ax.tick_params(labelsize=14)

        if bPlotU:
            clb.ax.set_title("u (m/s)", loc='center', fontsize=14)      
        else:
            clb.ax.set_title("v (m/s)", loc='center', fontsize=14)      


        divider = make_axes_locatable(axs[index, 2])
        cax = divider.append_axes("right", size="3%", pad=0.2)

        if bPlotU:
            clb = fig.colorbar(p_diff_u, ticks=local_diff_levels_u, cax=cax)
            clb.set_ticks(local_diff_levels_u)
        else:
            clb = fig.colorbar(p_diff_v, ticks=local_diff_levels_v, cax=cax)
            clb.set_ticks(local_diff_levels_v)

        clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
        clb.ax.tick_params(labelsize=14)


        clb.ax.set_title("(m/s)", loc='center', fontsize=14)  

        if bPlotU:
            axs[index, 2].text(400, 350, f"RMSE = {rmse_u:.3f} m/s", fontsize=12, ha='left', va='center')
        else:
            axs[index, 2].text(400, 350, f"RMSE = {rmse_v:.3f} m/s", fontsize=12, ha='left', va='center')

    #add title above the first row
    axs[0, 0].set_title("Truth", fontsize=16)
    axs[0, 1].set_title("From USWEs", fontsize=16)
    axs[0, 2].set_title("Difference", fontsize=16)

    #add text for iter_number
    axs[0, 1].text(400, 420, f"iter = {iter_numbers[0]}", fontsize=14, ha='left', va='center')
    axs[1, 1].text(400, 420, f"iter = {iter_numbers[1]}", fontsize=14, ha='left', va='center')
    axs[2, 1].text(400, 420, f"iter = {iter_numbers[2]}", fontsize=14, ha='left', va='center')
    axs[3, 1].text(400, 420, f"iter = {iter_numbers[3]}", fontsize=14, ha='left', va='center')

    axs[0, 2].text(400, 420, f"iter = {iter_numbers[0]}", fontsize=14, ha='left', va='center')
    axs[1, 2].text(400, 420, f"iter = {iter_numbers[1]}", fontsize=14, ha='left', va='center')
    axs[2, 2].text(400, 420, f"iter = {iter_numbers[2]}", fontsize=14, ha='left', va='center')
    axs[3, 2].text(400, 420, f"iter = {iter_numbers[3]}", fontsize=14, ha='left', va='center')

    if bPlotU:
        plt.savefig("Savana_UDE_U_iterations.png", dpi=300, bbox_inches='tight', pad_inches=0)
    else:
        plt.savefig("Savana_UDE_V_iterations.png", dpi=300, bbox_inches='tight', pad_inches=0)


if __name__ == '__main__':

    #plot ks contour 
    #plot_ks_contour()

    #plot UDE training history, Manning's n vs Manning's n truth on the f-Re plot
    #for iter_number in range(1, 151):
        #print(f"Plotting UDE training history, iter_number = {iter_number:04d}")
        #plot_UDE_training_history_ManningN(iter_number=iter_number, bPlotTitle=True)

    #plot the last iteration
    #plot_UDE_training_history_ManningN(iter_number=150, bPlotTitle=False)

    #create animation
    #create_animation_training_history()

    #plot Manning's n iterations
    #iter_numbers = [1, 10, 50, 150]
    #plot_ManningN_iterations(iter_numbers)

    #plot WSE iterations
    #plot_WSE_iterations(iter_numbers)

    #plot velocity iterations
    #plot_Velocity_iterations(iter_numbers, bPlotU=True)
    #plot_Velocity_iterations(iter_numbers, bPlotU=False)

    #plot Manning's n at all iterations
    for iter_number in range(1, 151):
        print(f"Plotting Manning's n, iter_number = {iter_number:04d}")
        plot_ManningN_one_iteration(iter_number)

    #create animation of Manning's n iterations
    create_animation_ManningN_iterations()

    print("Done!")