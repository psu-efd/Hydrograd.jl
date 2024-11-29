using Revise

using IterTools

using PyCall

using AdHydraulics

include("process_bed_2D.jl")
include("process_initial_condition_2D.jl")
include("semi_discretize_2D.jl")
include("misc_utilities_2D.jl")


println("Python version: ", PyCall.pyversion)
println("Python path: ", PyCall.pyprogramname)

pyHMT2D = pyimport("pyHMT2D")
pyHMT2D.gVerbose = true

# Call a function from the custom package
#file_path = joinpath(@__DIR__, "backwater.srhhydro" ) 
file_path = joinpath(@__DIR__, "simple.srhhydro" ) 
my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data(file_path)

ManningN_cell = my_srh_2d_data.ManningN_cell
ManningN_node = my_srh_2d_data.ManningN_node

#println("ManningN_cell: ", ManningN_cell)
#println("ManningN_node: ", ManningN_node)

println("Result from Python function: ", my_srh_2d_data)
#println("srhhydro file: ", my_srh_2d_data.srhhydro_obj.srhhydro_content)
for (k, v) in my_srh_2d_data.srhhydro_obj.srhhydro_content
    println("$k => $v")
end

#println("srhgeom file: ", my_srh_2d_data.srhgeom_obj.nodeCoordinates)
println("srhgeom file: ", typeof(my_srh_2d_data.srhgeom_obj.nodeCoordinates))
println("srhgeom file: ", size(my_srh_2d_data.srhgeom_obj.nodeCoordinates))
println("getNumOfElementsNodes: ", my_srh_2d_data.srhgeom_obj.getNumOfElementsNodes())

# Access the srhhydro data
srhhydro_obj = my_srh_2d_data.srhhydro_obj
srhhydro_ManningsN = srhhydro_obj.srhhydro_content["ManningsN"]
srhhydro_BC = srhhydro_obj.srhhydro_content["BC"]

srhhydro_MONITORINGLines = Dict()
if haskey(srhhydro_obj.srhhydro_content, "MONITORING")
    println("Key MONITORING exists in the Python dictionary srhhydro_content.")
    srhhydro_MONITORINGLines = srhhydro_obj.srhhydro_content["MONITORING"]
else
    println("Key MONITORING does not exist in the Python dictionary srhhydro_content.")
end

srhhydro_Grid = srhhydro_obj.srhhydro_content["Grid"]
srhhydro_RunType = srhhydro_obj.srhhydro_content["RunType"]
srhhydro_TurbulenceModel = srhhydro_obj.srhhydro_content["TurbulenceModel"]
srhhydro_OutputFormat = srhhydro_obj.srhhydro_content["OutputFormat"]
srhhydro_OutputInterval = srhhydro_obj.srhhydro_content["OutputInterval"]
srhhydro_Case = srhhydro_obj.srhhydro_content["Case"]
srhhydro_ParabolicTurbulence = srhhydro_obj.srhhydro_content["ParabolicTurbulence"]
srhhydro_Description = srhhydro_obj.srhhydro_content["Description"]
srhhydro_OutputOption  = srhhydro_obj.srhhydro_content["OutputOption"]
srhhydro_HydroMat = srhhydro_obj.srhhydro_content["HydroMat"]
srhhydro_InitCondOption = srhhydro_obj.srhhydro_content["InitCondOption"]
srhhydro_SimTime = srhhydro_obj.srhhydro_content["SimTime"]
srhhydro_SRHHYDRO = srhhydro_obj.srhhydro_content["SRHHYDRO"]
srhhydro_ModelTemp = srhhydro_obj.srhhydro_content["ModelTemp"]
srhhydro_UnsteadyOutput = srhhydro_obj.srhhydro_content["UnsteadyOutput"]

srhhydro_IQParams = Dict()
if haskey(srhhydro_obj.srhhydro_content, "IQParams")
    println("Key IQParams exists in the Python dictionary srhhydro_content.")
    srhhydro_IQParams = srhhydro_obj.srhhydro_content["IQParams"]
else
    println("Key IQParams does not exist in the Python dictionary srhhydro_content.")
end

#number of INLET-Q boundary conditions
nInletQ_BCs = length(srhhydro_IQParams)

srhhydro_EWSParamsC = Dict()
if haskey(srhhydro_obj.srhhydro_content, "EWSParamsC")
    println("Key EWSParamsC exists in the Python dictionary srhhydro_content.")
    srhhydro_EWSParamsC = srhhydro_obj.srhhydro_content["EWSParamsC"]
else
    println("Key EWSParamsC does not exist in the Python dictionary srhhydro_content.")
end

#number of EXIT-H boundary conditions
nExitH_BCs = length(srhhydro_EWSParamsC)

#srhhydro_ISupCrParams = srhhydro_obj.srhhydro_content["ISupCrParams"]
#srhhydro_EWSParamsRC = srhhydro_obj.srhhydro_content["EWSParamsRC"]
#srhhydro_EQParams = srhhydro_obj.srhhydro_content["EQParams"]
#srhhydro_NDParams = srhhydro_obj.srhhydro_content["NDParams"]

println("ManningsN: ", srhhydro_ManningsN)
println("BC: ", srhhydro_BC)
println("MONITORING: ", srhhydro_MONITORINGLines)
println("IQParams: ", srhhydro_IQParams)
println("EWSParamsC: ", srhhydro_EWSParamsC)
#println("ISupCrParams: ", srhhydro_ISupCrParams)
#println("EWSParamsRC: ", srhhydro_EWSParamsRC)
#println("EQParams: ", srhhydro_EQParams)
#println("NDParams: ", srhhydro_NDParams)

# Aceess the mesh geometry data 
srhgeom_obj = my_srh_2d_data.srhgeom_obj
numOfCells = srhgeom_obj.numOfElements    # Number of cells
numOfNodes = srhgeom_obj.numOfNodes       # Number of nodes

numOfNodeStrings = srhgeom_obj.numOfNodeStrings    # Number of node strings (number of boundaries; it does not include the default boundary)
numOfNonDefaultBoundaries = length(srhhydro_BC)    # Number of non-default boundaries (excluding the default boundary and monitoring lines)
numOfAllBoundaries = numOfNonDefaultBoundaries + 1  # Number of all boundaries (including the default boundary)

cellNodesList = srhgeom_obj.elementNodesList    # List of cell's nodes  (2D Int array: [cellID, gMax_Nodes_per_Cell])
cellNodesCount = srhgeom_obj.elementNodesCount  # Count of cell's nodes: how many nodes for each cell (1D Int array: [numOfCells])

cellFacesList = srhgeom_obj.elementEdgesList    # List of cell's faces  (2D Int array: [cellID, gMax_Nodes_per_Cell]). The numbers of nodes and faces are the same for each cell.

println("cellNodesCount: ", cellNodesCount)
println("cellNodesList: ", cellNodesList)
println("cellFacesList: ", cellFacesList)

nodeCoordinates = srhgeom_obj.nodeCoordinates      # Node coordinates: Float64 2D array [numOfNodes, 3]
twoDMeshBoundingbox = srhgeom_obj.twoDMeshBoundingbox  # 2D mesh bounding box: array [xmin, ymin, zmin, xmax, ymax, zmax]
cellBedElevation = srhgeom_obj.elementBedElevation  # Element bed elevation: Float64 1D array [numOfCells]
nodeStringsDict = srhgeom_obj.nodeStringsDict        #Dictionary of node strings: key: node string ID, value: node string nodes
nodeCellsList = srhgeom_obj.nodeElementsList      # List of node's cells: list of list. Each list contains the cells for each node. 
nodeCellsCount = srhgeom_obj.nodeElementsCount    # Count of node's cells: how many cells for each node. 1D Int array: [numOfNodes]

# Convert srhgeom_obj.edges to Julia dictionary
#face_Dict: key: (node1, node2), value: face ID
faces_Dict = Dict{Tuple{Int, Int}, Int}(tuple(convert(Int64,k[1]), convert(Int64,k[2])) => convert(Int64,v) for (k, v) in srhgeom_obj.edges)
#faces_r_Dict (reverse of faces_Dict): key: face ID, value: (node1, node2)
faces_r_Dict = Dict{Int, Tuple{Int, Int}}(convert(Int64,k) => tuple(convert(Int64,v[1]), convert(Int64,v[2])) for (k, v) in srhgeom_obj.edges_r)

println("faces_Dict: ")
for key in keys(faces_Dict)
    println("Key: $key, Value: $(faces_Dict[key])")
end

println("faces_r_Dict: ")
for key in sort(collect(keys(faces_r_Dict)))
    println("Key: $key, Value: $(faces_r_Dict[key])")
end

faceCells_Dict = srhgeom_obj.edgeElements       # Dictionary for the List of face's cells: dictionary: {faceID: [cell list]}
numOfFaces = length(faceCells_Dict)             # Number of faces

println("numOfFaces: ", numOfFaces)

println("faceCells_Dict: ")
for key in sort(collect(keys(faceCells_Dict)))
    println("Key: $key, Value: $(faceCells_Dict[key])")
end

boundaryFaces_Dict = srhgeom_obj.boundaryEdges    # Dictionary for the List of boundary faces: dictionary: {boundaryID: [list of face IDs]}. It also has the default boundary (default: wall)
allBoundaryFacesIDs_List = srhgeom_obj.allBoundaryEdgeIDs  # List of all boundary faces IDs: all lumped to one list
numOfAllBounaryFaces = length(allBoundaryFacesIDs_List)  # Number of all boundary faces

println("boundaryFaces_Dict: ", boundaryFaces_Dict)
println("allBoundaryFacesIDs_List: ", allBoundaryFacesIDs_List)

#ghost cells for boundary faces: ghost cell's index is the same as their index in allBoundaryFacesIDs_List, e.g., the first ghost cell is for the first boundary face in allBoundaryFacesIDs_List.
#boundaryFacesGhostCells_Dict = Dict{Int, Tuple{Int, Int}}()  #Dictionary for the ghost cells of boundary faces: {boundaryFaceID: (left ghost cell, right ghost cell)}

#total number of cells (including ghost cells)
numOfTotalCells = numOfCells + numOfAllBounaryFaces

#ghost cell IDs: ghost cells are added to the end of the cell list
ghostCellIDs = [numOfCells + i for i in 1:numOfAllBounaryFaces]

#boundary face ID to ghost cell ID: key: boundary face ID, value: ghost cell ID
boundaryFaceID_to_ghostCellID_Dict = Dict{Int, Int}()
#loop over all boundary faces to find the ghost cell ID for each boundary face
for iBoundaryFace in 1:numOfAllBounaryFaces
    boundaryFaceID_to_ghostCellID_Dict[allBoundaryFacesIDs_List[iBoundaryFace]] = ghostCellIDs[iBoundaryFace]
end

#ghost cell ID to boundary face ID: key: ghost cell ID, value: boundary face ID
ghostCellID_to_boundaryFaceID_Dict = Dict{Int, Int}()
#loop over all ghost cells to find the boundary face ID for each ghost cell
for iGhostCell in 1:numOfAllBounaryFaces
    ghostCellID_to_boundaryFaceID_Dict[ghostCellIDs[iGhostCell]] = allBoundaryFacesIDs_List[iGhostCell]
end

#cell's neighbors: key: cell ID, value: list of neighbor cell IDs
cellNeighbors_Dict = Dict{Int, Vector{Int}}()

#loop over all cells to find the neighbors of each cell
for iCell in 1:numOfCells
    cellNeighbors_Dict[iCell] = []
    
    nNodes = cellNodesCount[iCell]
    
    #loop over all the faces of the current cell
    for iFace in 1:nNodes
        faceID = cellFacesList[iCell, iFace]
        
        faceCells = faceCells_Dict[abs(faceID)]
        
        if length(faceCells) == 2   #internal face
            if faceCells[1] == iCell
                push!(cellNeighbors_Dict[iCell], faceCells[2])
            else
                push!(cellNeighbors_Dict[iCell], faceCells[1])
            end
        else  #boundary face
            ghostCellID = boundaryFaceID_to_ghostCellID_Dict[abs(faceID)]
            push!(cellNeighbors_Dict[iCell], -ghostCellID)  #negative ghost cell ID is used to indicate the ghost cell
        end
    end
end

println("cellNeighbors_Dict: ")
for key in sort(collect(keys(cellNeighbors_Dict)))
    println("Key: $key, Value: $(cellNeighbors_Dict[key])")
end

#check cell's nodes counter-clockwise
check_cell_nodes_counter_clockwise_srhgeom(numOfCells, cellNodesList, nodeCoordinates, cellNodesCount)

#compute mesh properties
cell_areas, cell_centroids, cell_normals, face_normals, face_lengths = compute_mesh_properties_srhgeom(numOfCells, numOfFaces, numOfNodes, nodeCoordinates, cellNodesList, cellNodesCount, faces_r_Dict)

let counter_temp = 0
    println("cell_areas: ")
    for area in cell_areas
        println(area)
        counter_temp += 1
        if counter_temp == 5
            break
        end
    end
end 

let counter_temp = 0
    println("cell_centroids: ")
    for i in 1:size(cell_centroids, 1)
        println(cell_centroids[i,:])
        counter_temp += 1
        if counter_temp == 5
            break
        end
    end
end 

let counter_temp = 0
    println("cell_normals: ")
    for normals in cell_normals
        println(counter_temp+1, ": ", normals)
        counter_temp += 1
        if counter_temp == 5
            #break
        end
    end
end 

let counter_temp = 0
    println("face_normals: ")
    for normals in face_normals
        println(counter_temp+1, ": ", normals)
        counter_temp += 1
        if counter_temp == 5
            #break
        end
    end
end 

let counter_temp = 0
    println("face_lengths: ")
    for face_length in face_lengths
        println(counter_temp+1, ": ", face_length)
        counter_temp += 1
        if counter_temp == 5
            #break
        end
    end
end 

# Access srhmat file data
srhmat_obj = my_srh_2d_data.srhmat_obj
srhmat_numOfMaterials = srhmat_obj.numOfMaterials
srhmat_matNameList = srhmat_obj.matNameList
srhmat_matZoneCells = srhmat_obj.matZoneCells

println("numOfMaterials: ", srhmat_numOfMaterials)
#println("matNameList: ", srhmat_matNameList)
#println("matZoneCells: ", srhmat_matZoneCells)

#  setup bed 
zb_face = zeros(Float64, numOfFaces)      #zb at faces 
zb_cell = zeros(Float64, numOfCells)      #zb at cell centers 
zb_ghostCells = zeros(Float64, numOfAllBounaryFaces)   #zb at ghost cell centers 
S0 = zeros(Float64, numOfCells, 2)          #bed slope at cell centers 

setup_bed!(numOfCells, numOfNodes, nodeCoordinates, cellNodesList, cellNodesCount, cell_centroids, zb_cell, true)

#setup zb_ghostCells: zb_ghostCells[i] has the same bed elevation at the neighbor cell 
update_ghost_cells_bed!(numOfAllBounaryFaces, allBoundaryFacesIDs_List, faceCells_Dict, zb_cell, zb_ghostCells)


# computer gradient of a scalar field
function compute_scalar_gradients!(numOfCells, cell_areas, cell_normals, face_lengths, cellNodesCount, cellFacesList, cellNeighbors_Dict, scalar_variable, grad_scalar_variable)
    #check the size of the grad_scalar_variable
    if size(grad_scalar_variable) != (numOfCells, 2)
        println("Error: grad_scalar_variable size is not correct.")
        exit(-1)
    end

    fill!(grad_scalar_variable, 0.0)
    
    #loop over all cells to compute the gradient of the scalar field
    for iCell in 1:numOfCells
        cell_gradient = [0.0, 0.0]  # Gradient accumulator for this cell
        
        #neighbor cells of the current cell
        cellNeighbors = cellNeighbors_Dict[iCell]
        
        #number of nodes for the current cell
        nNodes = cellNodesCount[iCell]
        
        cell_faces = cellFacesList[iCell,:]
        
        #loop over all faces of the current cell
        for iFace in 1:nNodes
            faceID = cell_faces[iFace]
            neighbor_cellID = cellNeighbors[iFace]
            
            # Value of the variable at the current cell and neighbor cell
            variable_c = scalar_variable[iCell]
            if neighbor_cellID < 0  # Boundary face
                variable_n = variable_c  # Assume zero gradient at boundary
            else
                variable_n = scalar_variable[neighbor_cellID]
            end
            
            # Compute the variable value on face (average of cell and neighbor for now; can be improved with interpolation)
            variable_f = (variable_c + variable_n) / 2.0
            
            # Compute flux contribution
            flux = cell_normals[iCell][iFace] * variable_f * face_lengths[abs(faceID)]
            cell_gradient .+= flux
        end
        
        # Finalize the gradient by dividing by the cell area
        grad_scalar_variable[iCell, :] = cell_gradient / cell_areas[iCell]
        
    end
end

# interpolate a scalar field from cell centers to face centers
function cells_to_faces_scalar!(numOfFaces, faceCells_Dict, scalar_variable_c, scalar_variable_f)
    
    #loop through faces  
    @inbounds for iFace in 1:numOfFaces
        #get current face's cells
        facecells = faceCells_Dict[iFace]
        
        if length(facecells) == 2   #internal face
            scalar_variable_f[iFace] = (scalar_variable_c[facecells[1]] + scalar_variable_c[facecells[2]]) / 2.0    
        else  #boundary face
            scalar_variable_f[iFace] = scalar_variable_c[facecells[1]]  
        end
    end
end

#interpolate zb from cell to face 
cells_to_faces_scalar!(numOfFaces, faceCells_Dict, zb_cell, zb_face)

#compute bed slope at cell centers
compute_scalar_gradients!(numOfCells, cell_areas, cell_normals, face_lengths, cellNodesCount, cellFacesList, cellNeighbors_Dict, zb_cell, S0)
#bed slope is negative of zb gradient 
S0 = -S0

println("zb_cell: ", zb_cell)   
println("zb_face: ", zb_face)
println("S0: ", S0) 

if false
    vector_data = [S0] 
    vector_names = ["S0"]

    scalar_data = [zb_cell]
    scalar_names = ["zb_cell"]

    file_path = joinpath(@__DIR__, "zb_S0.vtk" ) 
    export_to_vtk(file_path, nodeCoordinates, cellNodesList, cellNodesCount, scalar_data, scalar_names, vector_data, vector_names)    
    println("zb and S0 are saved to ", file_path)
    #exit(0)
end

#make a copy of the bathymetry truth: zb_cell_truth
zb_cell_truth = deepcopy(zb_cell)

# setup initial conditions 
eta = zeros(numOfCells)          #free surface elevation at cells 
h = zeros(numOfCells)            #water depth at cells 
q_x = zeros(numOfCells)          #q_x=hu at cells 
q_y = zeros(numOfCells)          #q_y=hv at cells 

eta_ghostCells = zeros(numOfAllBounaryFaces)          #free surface elevation at ghost cells 
h_ghostCells = zeros(numOfAllBounaryFaces)            #water depth at ghost cells 
q_x_ghostCells = zeros(numOfAllBounaryFaces)          #q_x=hu at ghost cells 
q_y_ghostCells = zeros(numOfAllBounaryFaces)          #q_y=hv at ghost cells 

total_water_volume = Float64[]   #total volume of water in the domain 

#setup initial condition for eta, h, q_x, q_y
setup_initial_eta!(numOfCells, nodeCoordinates, cell_centroids, eta, zb_cell, h, true)

#setup ghost cells for initial condition
#update_ghost_cells_eta_h!(numOfAllBounaryFaces, allBoundaryFacesIDs_List, faceCells_Dict, eta, h, eta_ghostCells, h_ghostCells)
for iBoundaryFace in 1:numOfAllBounaryFaces
    cellID_neighbor = faceCells_Dict[allBoundaryFacesIDs_List[iBoundaryFace]][1]
    
    eta_ghostCells[iBoundaryFace] = eta[cellID_neighbor]
    h_ghostCells[iBoundaryFace] = h[cellID_neighbor]
    q_x_ghostCells[iBoundaryFace] = q_x[cellID_neighbor]
    q_y_ghostCells[iBoundaryFace] = q_y[cellID_neighbor]
end
println("h_ghostCells: ", h_ghostCells)


#process boundary conditions 
#loop over all boundaries
for iBoundary in 1:numOfNonDefaultBoundaries
    if haskey(srhhydro_BC, iBoundary)
        println("Key IQParams exists in the Python dictionary srhhydro_content.")
        println("Inlet specific discharge: ", srhhydro_IQParams[iBoundary])
    end

    boundaryType = srhhydro_BC[iBoundary]   #boundary type: wall, inletQ, exitH
    
    if lowercase(boundaryType) == "inlet-q"
        println("INLET-Q boundary condition is set for boundary ", iBoundary)
        
        #find the corresponding Q in IQParams
        if haskey(srhhydro_IQParams, iBoundary)
            println("Key IQParams exists in the Python dictionary srhhydro_content.")
            println("Inlet specific discharge: ", srhhydro_IQParams[iBoundary])
            
            #update the boundary value: currenlty only support constant discharge
            Q_value = parse(Float64, srhhydro_IQParams[iBoundary][1])
                        
        else
            println("Key IQParams does not exist in the Python dictionary srhhydro_content.")
        end
        
    elseif lowercase(boundaryType) == "exit-h"
        println("EXIT-H boundary condition is set for boundary ", iBoundary)
        
    elseif lowercase(boundaryType) == "wall"
        println("WALL boundary condition is set for boundary ", iBoundary)
    end
    
end


#fluxes on faces 
#fluxes = zeros(Number, 2, mesh.nFaces)
#temp flux
#flux = zeros(Number, 2)

#set up ODE parameters. In this problem, the "parameter" is para=zb_cell.
para = Float64.(zb_cell)
Q0 = hcat(h, q_x, q_y)   #initial condition: Q = [h q_x q_y]

# time information 
tspan = (0.0, 1.0)
dt = 0.1
t = tspan[1]:dt:tspan[2]

dt_save = (tspan[2] - tspan[1])/10.0
t_save = tspan[1]:dt_save:tspan[2]

# define the ODE
ode_f = ODEFunction((dQdt, Q, para, t) -> swe_2D_rhs!(dQdt, Q, para, t,
mesh, swe_1d_constants, left_bcType, right_bcType, left_bcValue, right_bcValue, ManningN), 
jac_prototype=nothing)

prob = ODEProblem(ode_f, Q0, tspan, p)



println("All done!")
