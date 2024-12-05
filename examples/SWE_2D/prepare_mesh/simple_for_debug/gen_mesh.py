import numpy as np

from pyHMT2D.Misc import gmsh2d_to_srh

def convert_gmsh_to_srh():

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

    gmsh2d_to_srh("rect_structured.msh", "simple", units="Meters",
                  bAddMonitoringLines=True, monitoringLines=monitoringLines)

    #gmsh2d_to_srh("rectangular_domain.msh", "test", units="Meters")
    #gmsh2d_to_srh("rectangular_domain_with_embedded_line.msh", "rectangular_domain_with_embedded_line", units="Meters")


if __name__ == "__main__":
    convert_gmsh_to_srh()

    print("Done!")