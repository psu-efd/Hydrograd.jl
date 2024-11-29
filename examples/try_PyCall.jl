using Revise

using IterTools

using PyCall

using AdHydraulics

println("Python version: ", PyCall.pyversion)
println("Python path: ", PyCall.pyprogramname)

pyHMT2D = pyimport("pyHMT2D")
pyHMT2D.gVerbose = true

# Call a function from the custom package
file_path = joinpath(@__DIR__, "backwater.srhhydro" ) 
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

srhhydro_IQParams = srhhydro_obj.srhhydro_content["IQParams"]
#srhhydro_ISupCrParams = srhhydro_obj.srhhydro_content["ISupCrParams"]
srhhydro_EWSParamsC = srhhydro_obj.srhhydro_content["EWSParamsC"]
#srhhydro_EWSParamsRC = srhhydro_obj.srhhydro_content["EWSParamsRC"]
#srhhydro_EQParams = srhhydro_obj.srhhydro_content["EQParams"]
#srhhydro_NDParams = srhhydro_obj.srhhydro_content["NDParams"]

println("ManningsN: ", srhhydro_ManningsN)
println("BC: ", srhhydro_BC)
#println("MONITORING: ", res_MONITORING)
println("IQParams: ", srhhydro_IQParams)
#println("ISupCrParams: ", srhhydro_ISupCrParams)
println("EWSParamsC: ", srhhydro_EWSParamsC)
#println("EWSParamsRC: ", srhhydro_EWSParamsRC)
#println("EQParams: ", srhhydro_EQParams)
#println("NDParams: ", srhhydro_NDParams)

# Aceess the mesh geometry data 
srhgeom_obj = my_srh_2d_data.srhgeom_obj
numOfCells = srhgeom_obj.numOfElements    # Number of cells
numOfNodes = srhgeom_obj.numOfNodes       # Number of nodes
numOfNodeStrings = srhgeom_obj.numOfNodeStrings    # Number of node strings
cellNodesList = srhgeom_obj.elementNodesList    # List of cell's nodes  (2D Int array: [cellID, gMax_Nodes_per_Cell])
cellNodesCount = srhgeom_obj.elementNodesCount  # Count of cell's nodes: how many nodes for each cell (1D Int array: [numOfCells])
nodeCoordinates = srhgeom_obj.nodeCoordinates      # Node coordinates: Float64 2D array [numOfNodes, 3]
twoDMeshBoundingbox = srhgeom_obj.twoDMeshBoundingbox  # 2D mesh bounding box: array [xmin, ymin, zmin, xmax, ymax, zmax]
cellBedElevation = srhgeom_obj.elementBedElevation  # Element bed elevation: Float64 1D array [numOfCells]
nodeStringsDict = srhgeom_obj.nodeStringsDict        #Dictionary of node strings: key: node string ID, value: node string nodes
nodeCellsList = srhgeom_obj.nodeElementsList      # List of node's cells: list of list. Each list contains the cells for each node. 
nodeCellsCount = srhgeom_obj.nodeElementsCount    # Count of node's cells: how many cells for each node. 1D Int array: [numOfNodes]

# Convert srhgeom_obj.edges to Julia dictionary
#face_Dict: key: (node1, node2), value: face ID
faces_Dict = Dict{Tuple{Int, Int}, Int}(tuple(convert(Int64,k[1]), convert(Int64,k[2])) => convert(Int64,v) for (k, v) in srhgeom_obj.edges)
#faces_r_Dict (reverse of faces_Dict): key: face ID, value: (cell1, cell2)
faces_r_Dict = Dict{Int, Tuple{Int, Int}}(convert(Int64,k) => tuple(convert(Int64,v[1]), convert(Int64,v[2])) for (k, v) in srhgeom_obj.edges_r)

#println("edges_Dict: ", edges_Dict)
#println("edges_r_Dict: ", edges_r_Dict)

# let counter_temp = 0
#     for (k, v) in edges_Dict
#         println("$k => $v")
#         counter_temp += 1
#         if counter_temp == 5
#             break
#         end
#     end
# end

# let counter_temp = 0
#     for (k, v) in edges_r_Dict
#         println("$k => $v")
#         counter_temp += 1
#         if counter_temp == 5
#             break
#         end
#     end
# end 

faceCells_Dict = srhgeom_obj.edgeElements       # Dictionary for the List of face's cells: dictionary: {faceID: [cell list]}
numOfFaces = length(faceCells_Dict)             # Number of faces

boundaryFaces_Dict = srhgeom_obj.boundaryEdges    # Dictionary for the List of boundary faces: dictionary: {boundaryID: [list of face IDs]}
allBoundaryFacesIDs_List = srhgeom_obj.allBoundaryEdgeIDs  # List of all boundary faces IDs: all lumped to one list

#check cell's nodes counter-clockwise
check_cell_nodes_counter_clockwise_srhgeom(numOfCells, cellNodesList, nodeCoordinates, cellNodesCount)

#compute mesh properties
cell_areas, cell_centroids, cell_normals = compute_mesh_properties_srhgeom(numOfCells, numOfFaces, numOfNodes, nodeCoordinates, cellNodesList, cellNodesCount)

let counter_temp = 0
    for area in cell_areas
        println(area)
        counter_temp += 1
        if counter_temp == 5
            break
        end
    end
end 

let counter_temp = 0
    for i in 1:size(cell_centroids, 1)
        println(cell_centroids[i,:])
        counter_temp += 1
        if counter_temp == 5
            break
        end
    end
end 

let counter_temp = 0
    for normals in cell_normals
        println(normals)
        counter_temp += 1
        if counter_temp == 5
            break
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
#Functions to setup the bed profile (this is just for some examples)
function my_gauss(x::Float64; sigma::Float64=1.0, h::Float64=1.0, mid::Float64=0.0)
    variance = sigma^2
    return h * exp(-(x - mid)^2 / (2 * variance))
end

function setup_bed!(numOfCells, cell_centroids, zb_cell)

    # parameters for bed setup
    b_bump_height = 0.3

    #loop through cells
    @inbounds for i in 1:numOfCells
        zb_cell[i] = my_gauss(cell_centroids[i,1], sigma=1.0, h=b_bump_height, 
                              mid=maximum(mesh.xFaces) / 2.0)
    end

    #interploate_zb_from_cell_to_face_and_compute_S0!(mesh, zb_face, zb_cell, S0)

    # optionally plot the bed for checking
    #plot_bed(mesh, zb_cell, zb_face, S0)
end

# setup initial condition: free surface 
function setup_initial_eta!(numOfCells, coordinates, eta, zb_cell, h, swe_1d_constants)
    # parameters for initial free surface profile setup
    bump_center_x = 5.0  # center of the bump
    bump_half_width_x = 1.0  # bump's half width
    
    h_small = 0.01

    #loop over cells
    @inbounds for i in 1:numOfCells
        if coordinates[i,1] < bump_center_x
            eta[i] = 1.0
        else
            eta[i] = 0.5 #0.5  #0.5
        end
    end

    #update water depth
    @inbounds for i in 1:numOfCells
        h[i] = eta[i] - zb_cell[i]
    end

    h[h.<0.0] .= h_small  #ensure positivity of water depth h

    #update the free surface elevation again in case h has been clipped
    @inbounds for i in 1:numOfCells
        eta[i] = h[i] + zb_cell[i]
    end

    #optionally plot the free surface for checking 
    #plot_free_surface_elevation(xCells, zb_cell, eta, h)
end


#zb_face = zeros(mesh.nFaces)      #zb at faces (points in 1D)
zb_cell = zeros(Float64, numOfCells)      #zb at cell centers 
#S0 = zeros(mesh.nCells)          #bed slope at cell centers 

setup_bed!(numOfCells, cell_centroids, zb_cell)

#make a copy of the bathymetry truth: zb_cell_truth
zb_cell_truth = deepcopy(zb_cell)

# setup initial conditions 
eta = zeros(numOfCells)          #free surface elevation at cells 
h = zeros(numOfCells)            #water depth at cells 
q_x = zeros(numOfCells)            #q_x=hu at cells 
q_y = zeros(numOfCells)            #q_y=hv at cells 
total_water_volume = Float64[]      #total volume of water in the domain 

setup_initial_eta!(mesh, eta, zb_cell, h, swe_1d_constants)

#setup boundary conditions
left_bcType = "inletQ"
#left_bcType = "wall"
left_bcValue = inletQ = 0.8
right_bcType = "exitH"
#right_bcType = "wall"
right_bcValue = exitH = 0.75

#fluxes on faces 
#fluxes = zeros(Number, 2, mesh.nFaces)
#temp flux
#flux = zeros(Number, 2)

#set up ODE parameters. In this problem, p=zb_cell
p = Float64.(zb_cell)
Q0 = hcat(h, q)   #initial condition: Q = [h q]

# time information 
tspan = (0.0, 1.0)
dt = 0.1
t = tspan[1]:dt:tspan[2]

dt_save = (tspan[2] - tspan[1])/10.0
t_save = tspan[1]:dt_save:tspan[2]

# define the ODE
ode_f = ODEFunction((dQdt, Q, p, t) -> swe_2D_rhs!(dQdt, Q, p, t,
        mesh, swe_1d_constants, left_bcType, right_bcType, left_bcValue, right_bcValue, ManningN), 
        jac_prototype=nothing)

prob = ODEProblem(ode_f, Q0, tspan, p)



println("All done!")
