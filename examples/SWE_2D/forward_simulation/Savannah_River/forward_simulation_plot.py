"""
Some specific utility functions, such as plotting.

"""

import matplotlib.pyplot as plt
import numpy as np

#import shapefile

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


def plot_forward_simulation_result():
    """
    Make plots for the forward simulation:
    1. bathy
    2. Manning's n zones
    3. WSE
    4. Velocity vector

    :return:
    """

    fig, axs = plt.subplots(2, 2, figsize=(16, 9), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    #fig.subplots_adjust(hspace=.15, wspace=.01)

    zb_min = 21
    zb_max = 29

    local_levels = np.linspace(zb_min, zb_max, 51)

    #1. plot bathy in xy and center line
    vtk_fileName = 'forward_simulation_results_0101_tri.vtk'


    #read data from vtk file
    vtk_result = meshio.read(vtk_fileName)

    # Extract points and cells (triangles)
    points = vtk_result.points[:, :2]  # Take only the X and Y coordinates
    cells = vtk_result.cells_dict.get("triangle")

    if cells is None:
        raise ValueError("The VTK file does not contain triangular cells.")

    # Create a Triangulation object for plotting
    triangulation = Triangulation(points[:, 0], points[:, 1], cells)

    # number of triangles
    nTri = vtk_result.cells[0].data.shape[0]

    tri = np.zeros((nTri, 3))
    for i in range(nTri):
        tri[i, 0] = vtk_result.cells[0].data[i, 0]
        tri[i, 1] = vtk_result.cells[0].data[i, 1]
        tri[i, 2] = vtk_result.cells[0].data[i, 2]

    xCoord = vtk_result.points[:, 0]
    yCoord = vtk_result.points[:, 1]

    xl = xCoord.min()
    xh = xCoord.max()
    yl = yCoord.min()
    yh = yCoord.max()

    #get the data
    ManningN = np.squeeze(vtk_result.point_data['ManningN'])
    zb = np.squeeze(vtk_result.point_data['zb_cell'])
    WSE = np.squeeze(vtk_result.point_data['WSE'])
    U = np.squeeze(vtk_result.point_data['U'])
    Ux = U[:,0]
    Uy = U[:,1]
    Umag = np.sqrt(Ux**2 + Uy**2)

    print("WSE = ", WSE)
    print("U = ", U)
    print("Ux = ", Ux)
    print("Uy = ", Uy)
    print("Umag = ", Umag)

    #plot bathymetry
    cf_zb = axs[0, 0].tricontourf(xCoord, yCoord, tri, zb, local_levels,
                                        vmin=zb_min, vmax=zb_max, cmap=plt.cm.terrain, extend='neither')
    axs[0, 0].set_xlim([xl, xh])
    axs[0, 0].set_ylim([yl, yh])
    axs[0, 0].set_aspect('equal')
    axs[0, 0].set_title("Bathymetry", fontsize=14)

    axs[0, 0].set_xlabel('$x$ (m)', fontsize=16)
    axs[0, 0].tick_params(axis='x', labelsize=14)
    axs[0, 0].set_ylabel('$y$ (m)', fontsize=16)
    axs[0, 0].tick_params(axis='y', labelsize=14)

    #add an arrow for the flow direction
    axs[0, 0].annotate(
    'Flow',  # No text
    xy=(250, 200),  # End point
    xytext=(50, 300),  # Start point
    fontsize=16,  # Text size
    arrowprops=dict(
        facecolor='blue',  # Arrow color
        arrowstyle='->',  # Arrow style
        lw=2  # Line width
        )
    )

    divider = make_axes_locatable(axs[0, 0])
    cax = divider.append_axes("right", size="3%", pad=0.2)
    clb = fig.colorbar(cf_zb, ticks=np.linspace(21, 29, 7), cax=cax)
    clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.0f'))
    clb.ax.tick_params(labelsize=12)
    clb.ax.set_title("(m)", loc='center', fontsize=12)

    # plot Manning's n
    ManningN_min = 0.03
    ManningN_max = 0.05
    local_levels = np.linspace(ManningN_min, ManningN_max, 6)

    cf_zb = axs[0, 1].tricontourf(xCoord, yCoord, tri, ManningN, local_levels,
                                  vmin=ManningN_min, vmax=ManningN_max, cmap=plt.cm.jet, extend='neither')

    #add annotation to Manning's n zones
    # Add annotation
    axs[0,1].annotate(
        '$n_1$ = 0.04',  # Text to display
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

    axs[0, 1].annotate(
        '$n_2$ = 0.05',  # Text to display
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

    axs[0, 1].annotate(
        '$n_3$ = 0.03',  # Text to display
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

    axs[0, 1].annotate(
        '$n_4$ = 0.045',  # Text to display
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

    axs[0, 1].annotate(
        '$n_5$ = 0.05',  # Text to display
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

    axs[0, 1].set_xlim([xl, xh])
    axs[0, 1].set_ylim([yl, yh])
    #axs[0, 1].set_aspect('equal')
    axs[0, 1].set_aspect(1)
    axs[0, 1].set_title("Manning's $n$", fontsize=14)

    divider = make_axes_locatable(axs[0, 1])
    cax = divider.append_axes("right", size="3%", pad=0.2)
    clb = fig.colorbar(cf_zb, ticks=[0.03,0.035, 0.04,0.045,0.05], cax=cax)
    clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.3f'))
    clb.ax.tick_params(labelsize=12)
    clb.ax.set_title("(s/m$^{1/3}$)", loc='center', fontsize=12)

    axs[0, 1].set_xlabel('$x$ (m)', fontsize=16)
    axs[0, 1].tick_params(axis='x', labelsize=14)
    axs[0, 1].set_ylabel('$y$ (m)', fontsize=16)
    axs[0, 1].tick_params(axis='y', labelsize=14)

    #plot WSE
    WSE_min = min(WSE)
    WSE_max = max(WSE)
    local_levels = np.linspace(WSE_min, WSE_max, 51)

    cf_zb = axs[1, 0].tricontourf(xCoord, yCoord, tri, WSE, local_levels,
                                  vmin=WSE_min, vmax=WSE_max, cmap=plt.cm.jet, extend='neither')
    axs[1, 0].set_xlim([xl, xh])
    axs[1, 0].set_ylim([yl, yh])
    axs[1, 0].set_aspect('equal')
    axs[1, 0].set_title("WSE", fontsize=14)
    axs[1, 0].set_xlabel('$x$ (m)', fontsize=16)
    axs[1, 0].set_ylabel('$y$ (m)', fontsize=16)
    axs[1, 0].tick_params(axis='x', labelsize=14)
    axs[1, 0].tick_params(axis='y', labelsize=14)

    divider = make_axes_locatable(axs[1, 0])
    cax = divider.append_axes("right", size="3%", pad=0.2)
    clb = fig.colorbar(cf_zb, ticks=np.linspace(WSE_min, WSE_max, 7), cax=cax)
    clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    clb.ax.tick_params(labelsize=12)
    clb.ax.set_title("(m)", loc='center', fontsize=12)

    #Velocity
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
            axs[1,1].plot(x, y, color='blue', linewidth=0.5)
            axs[1,1].fill(x, y, color='lightblue', alpha=0.3)

    Umag_min = min(Umag)
    Umag_max = max(Umag)
    local_levels = np.linspace(Umag_min, Umag_max, 51)

    #velocity magnitude contour
    #axs[1, 1].triplot(triangulation, color="black", linewidth=0.5)
    #cf_zb = axs[1, 1].tricontourf(xCoord, yCoord, tri, Umag, local_levels, vmin=Umag_min, vmax=Umag_max, cmap=plt.cm.jet, extend='neither')

    #velocity vector plot
    #axs[1, 1].quiver(xCoord, yCoord, Ux, Uy, angles='xy', scale_units='xy', scale=0.05, headwidth=1, headlength=1, color='black')
    cf_zb = axs[1, 1].quiver(xCoord, yCoord, Ux, Uy, Umag, cmap=plt.cm.jet)

    axs[1, 1].set_xlim([xl, xh])
    axs[1, 1].set_ylim([yl, yh])
    axs[1, 1].set_aspect('equal')
    axs[1, 1].set_title("Velocity", fontsize=14)
    axs[1, 1].set_xlabel('$x$ (m)', fontsize=16)
    axs[1, 1].set_ylabel('$y$ (m)', fontsize=16)
    axs[1, 1].tick_params(axis='x', labelsize=14)
    axs[1, 1].tick_params(axis='y', labelsize=14)

    divider = make_axes_locatable(axs[1, 1])
    cax = divider.append_axes("right", size="3%", pad=0.2)
    clb = fig.colorbar(cf_zb, ticks=np.linspace(Umag_min, Umag_max, 7), cax=cax)
    clb.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.1f'))
    clb.ax.tick_params(labelsize=12)
    clb.ax.set_title("(m/s)", loc='center', fontsize=12)

    #add caption
    axs[0,0].text(-0.1, 1.05, "(a)", size=16, ha="center", transform=axs[0,0].transAxes)  # upper left
    axs[0,1].text(-0.1, 1.05, "(b)", size=16, ha="center", transform=axs[0,1].transAxes)
    axs[1, 0].text(-0.1, 1.05, "(c)", size=16, ha="center", transform=axs[1, 0].transAxes)  # upper left
    axs[1, 1].text(-0.1, 1.05, "(d)", size=16, ha="center", transform=axs[1, 1].transAxes)

    plt.savefig("Savana_forward_simulation_result.png", dpi=300, bbox_inches='tight', pad_inches=0)
    #plt.show()

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

def triangulate_vtk():
    #triangulate the vtk

    # Load the VTK file
    file_path = 'forward_simulation_results_0101.vtk'
    grid = pv.read(file_path)

    # Triangulate the grid
    triangulated = grid.triangulate()

    # Interpolate cell data to point data
    interpolated = triangulated.cell_data_to_point_data()

    # Save the triangulated grid to a new VTK file
    triangulated_file_path = 'forward_simulation_results_0101_tri.vtk'
    interpolated.save(triangulated_file_path)

    print(f"Triangulated file saved to: {triangulated_file_path}")

    print("Available point data:", interpolated.point_data.keys())


if __name__ == '__main__':

    #triangulate_vtk()

    plot_forward_simulation_result()

    print("Done!")