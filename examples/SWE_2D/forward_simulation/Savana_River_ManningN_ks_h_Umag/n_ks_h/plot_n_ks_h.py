"""
Plot the functions of n(ks, h).

"""

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.ticker as tick
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib.cm import jet
from matplotlib.colors import Normalize
import json

import meshio
import pyvista as pv



from pathlib import Path

plt.rc('text', usetex=True)  #allow the use of Latex for math expressions and equations
plt.rc('font', family='serif') #specify the default font family to be "serif"

def triangulate_vtk():
    #triangulate the vtk

    # Load the VTK file
    file_path = '../forward_simulation_results_0101.vtk'
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
    h_ks = np.squeeze(vtk_result.point_data['h_ks'])
    Re = np.squeeze(vtk_result.point_data['Re'])
    friction_factor = np.squeeze(vtk_result.point_data['friction_factor'])
    n = np.squeeze(vtk_result.point_data['ManningN'])

    print("h_ks = ", h_ks)
    print("Re = ", Re)
    print("friction_factor = ", friction_factor)
    print("n = ", n)

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
    plt.show()

def plot_n_ks_h():
    """
    Plot an example of the function n(ks, h) as a 2D contour plot.
    """

    fig, axs = plt.subplots(1, 1, figsize=(16, 9), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    #define the range of ks and h 
    ks = np.linspace(0, 1, 100)
    h = np.linspace(0, 1, 100)

    #generate the mesh grid
    KS, H = np.meshgrid(ks, h)


    #plot the function of n(ks, h)
    #read the data from the file
    data = np.loadtxt("n_ks_h.txt")
    
    plt.savefig("Savana_forward_simulation_result.png", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()


def plot_f_h_ks_Re():
    """
    Plot the function of f(h/ks, Re) and the scatter for simulation data.

    Reference: Cheng (2008) JHE, "Formulas for Friction Factor in Transitional Regimes", Eq. (17)
    """

    #get simulatio data from the VTK file
    vtk_fileName = 'forward_simulation_results_0101_tri.vtk'

    #read data from vtk file
    vtk_result = meshio.read(vtk_fileName)

    # Extract points and cells (triangles)
    points = vtk_result.points[:, :2]  # Take only the X and Y coordinates
    cells = vtk_result.cells_dict.get("triangle")

    if cells is None:
        raise ValueError("The VTK file does not contain triangular cells.")

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
    h_ks_simulation = np.squeeze(vtk_result.point_data['h_ks'])
    Re_simulation = np.squeeze(vtk_result.point_data['Re'])
    friction_factor_simulation = np.squeeze(vtk_result.point_data['friction_factor'])
    n_simulation = np.squeeze(vtk_result.point_data['ManningN'])


    print("h_ks_simulation = ", h_ks_simulation)
    print("Re_simulation = ", Re_simulation)
    print("friction_factor_simulation = ", friction_factor_simulation)
    print("n_simulation = ", n_simulation)

    print("Re_simulation.shape = ", Re_simulation.shape)
    print("friction_factor_simulation.shape = ", friction_factor_simulation.shape)
    print("n_simulation.shape = ", n_simulation.shape)


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


    #plot the function of f(h/ks, Re)
    fig, axs = plt.subplots(1, 1, figsize=(16, 9), sharex=False, sharey=False, facecolor='w', edgecolor='k')

    # Create normalizer for the color mapping
    norm = Normalize(vmin=min(h_ks), vmax=max(h_ks))

    #loop over the h_ks
    for i in range(len(h_ks)):
        color = jet(norm(h_ks[i]))  # Get color from colormap
        axs.plot(Re, f[i], color=color, label=f'h/ks = {h_ks[i]}')

    #plot the scatter for simulation data and colored by h/ks
    scatter = axs.scatter(Re_simulation, friction_factor_simulation, c=h_ks_simulation, alpha=0.5, edgecolors='black', linewidths=0.5, norm=norm, cmap='jet', label='Simulation data')

    #add a color bar and customize the color bar
    cbar = plt.colorbar(scatter, location='right', pad=0.02)
    cbar.ax.set_title("$h/k_s$", loc='center', fontsize=16)
    cbar.ax.tick_params(labelsize=14)

    #axs.legend()
    axs.set_xscale('log')
    axs.set_yscale('log')

    #set axis limits
    axs.set_xlim(100, 3E7)
    axs.set_ylim(0.02, 0.2)

    axs.set_xlabel('$Re$', fontsize=16)
    axs.set_ylabel('$f$', fontsize=16)

    axs.tick_params(axis='x', which='both', labelsize=14)
    axs.tick_params(axis='y', which='both', labelsize=14)

    #set y tick format
    axs.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    axs.yaxis.set_minor_formatter(tick.FormatStrFormatter('%.2f'))

    #add text to annotate the h/ks values for each line
    for i in range(len(h_ks)):
        axs.text(6E6, f[i][-1], f'$h/k_s$ = {h_ks[i]:.0f}', fontsize=14)

    plt.savefig("f_h_ks_Re.png", dpi=300, bbox_inches='tight', pad_inches=0.01)
    #plt.show()


if __name__ == '__main__':

    #triangulate the vtk file
    triangulate_vtk()

    #plot the function of f(h/ks, Re)
    plot_f_h_ks_Re()


    print("Done!")