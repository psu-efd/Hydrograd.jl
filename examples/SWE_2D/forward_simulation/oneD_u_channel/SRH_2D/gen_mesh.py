import numpy as np
import meshio

from pyHMT2D.Misc import gmsh2d_to_srh

def gmsh2d_to_3d(gmsh2d_fileName, gmsh3d_fileName):
    """
    Assign the elevation of points so the 2D mesh becomes 3D (still surface though)
    :return:
    """

    print("Assign elevation to points in Gmsh 2D mesh to make it 3D ...")

    #read in the Gmsh MSH with meshio
    mesh = meshio.read(gmsh2d_fileName)

    #output the 2D mesh as vtk for checking
    mesh.write("check_mesh_before.vtk")

    zb_min = - 0.5 #mininum zb in the middle
    Lx = 10.0 #channel length
    a = 0.1  #parameter to control the curvature

    #loop over all points in the mesh
    for i in range(mesh.points.shape[0]):
        x = mesh.points[i, 0]
        y = mesh.points[i, 1]

        #mesh.points[i,2] = max(0.0, 0.2 - 0.05*((x-10.0)**2.0 + (y-2.5)**2.0))
        mesh.points[i, 2] = min(-0.1, a*(x-Lx/2)**2 + zb_min)

    #output the 3D mesh as vtk for checking
    mesh.write("check_mesh_after.vtk")

    #write out the new mesh
    mesh.write(gmsh3d_fileName, file_format="gmsh22", binary=False)


def convert_gmsh_to_srh(gmsh3d_fileName, srh2d_caseName):

    monitoringLine1 = np.zeros((2,2))
    monitoringLine1[0, 0] = 2.0   #xML_start
    monitoringLine1[0, 1] = 0.0   #yML_start
    monitoringLine1[1, 0] = 2.0   #xML_end
    monitoringLine1[1, 1] = 1.0   #yML_end

    monitoringLine2 = np.zeros((2,2))
    monitoringLine2[0, 0] = 2.0   #xML_start
    monitoringLine2[0, 1] = 1.0   #yML_start
    monitoringLine2[1, 0] = 2.0   #xML_end
    monitoringLine2[1, 1] = 3.0  #yML_end

    monitoringLines = []
    monitoringLines.append(monitoringLine1)
    monitoringLines.append(monitoringLine2)

    gmsh2d_to_srh(gmsh3d_fileName, srh2d_caseName, units="Meters",
                  bAddMonitoringLines=False, monitoringLines=monitoringLines)

if __name__ == "__main__":
    gmsh2d_fileName = "rect_structured.msh"
    gmsh3d_fileName = "rect_structured_with_zb.msh"
    srh2d_caseName = "oneD_u_channel"

    gmsh2d_to_3d(gmsh2d_fileName, gmsh3d_fileName)

    convert_gmsh_to_srh(gmsh3d_fileName, srh2d_caseName)

    print("Done!")