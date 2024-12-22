#boundary conditions for 2D SWE

mutable struct BoundaryConditions2D
    # Numbers of each boundary type
    nInletQ_BCs::Int      #number of inlet-Q boundaries
    nExitH_BCs::Int       #number of exit-H boundaries
    nWall_BCs::Int        #number of wall boundaries
    nSymm_BCs::Int        #number of symmetry boundaries

    #indices of each boundary in the boundary list
    inletQ_BC_indices::Vector{Int}   #indices of inlet-Q boundaries in the boundary list
    exitH_BC_indices::Vector{Int}    #indices of exit-H boundaries in the boundary list
    wall_BC_indices::Vector{Int}     #indices of wall boundaries in the boundary list
    symm_BC_indices::Vector{Int}     #indices of symmetry (no slip) boundaries in the boundary list

     # Inlet Q boundary data
    inletQ_faceIDs::Vector{Vector{Int}}                     #face IDs of the inlet-q boundaries
    inletQ_ghostCellIDs::Vector{Vector{Int}}               #ghost cell IDs of the inlet-q boundaries
    inletQ_internalCellIDs::Vector{Vector{Int}}           #internal cell IDs of the inlet-q boundaries
    inletQ_faceCentroids::Vector{Matrix{Float64}}          #face centroids of the inlet-q boundaries
    inletQ_faceOutwardNormals::Vector{Matrix{Float64}}     #face outward normals of the inlet-q boundaries
    inletQ_TotalQ::Vector{Float64}                         #total discharge of the inlet-q boundaries
    inletQ_H::Vector{Vector{Float64}}                     #water depth of the inlet-q boundaries
    inletQ_A::Vector{Vector{Float64}}                     #cross-sectional area of the inlet-q boundaries
    inletQ_ManningN::Vector{Vector{Float64}}             #Manning's n of the inlet-q boundaries
    inletQ_Length::Vector{Vector{Float64}}               #length of the inlet-q boundaries
    inletQ_TotalA::Vector{Float64}                       #total cross-sectional area of the inlet-q boundaries
    inletQ_DryWet::Vector{Vector{Int}}                   #dry/wet flag of the inlet-q boundaries

     # Exit H boundary data
    exitH_faceIDs::Vector{Vector{Int}}                     #face IDs of the exit-h boundaries
    exitH_ghostCellIDs::Vector{Vector{Int}}                 #ghost cell IDs of the exit-h boundaries
    exitH_internalCellIDs::Vector{Vector{Int}}             #internal cell IDs of the exit-h boundaries
    exitH_faceCentroids::Vector{Matrix{Float64}}          #face centroids of the exit-h boundaries
    exitH_faceOutwardNormals::Vector{Matrix{Float64}}     #face outward normals of the exit-h boundaries
    exitH_WSE::Vector{Float64}                             #water surface elevation of the exit-h boundaries
    exitH_H::Vector{Vector{Float64}}                     #water depth of the exit-h boundaries
    exitH_A::Vector{Vector{Float64}}                     #cross-sectional area of the exit-h boundaries

     # Wall boundary data
    wall_faceIDs::Vector{Vector{Int}}                     #face IDs of the wall boundaries
    wall_ghostCellIDs::Vector{Vector{Int}}               #ghost cell IDs of the wall boundaries
    wall_internalCellIDs::Vector{Vector{Int}}           #internal cell IDs of the wall boundaries
    wall_faceCentroids::Vector{Matrix{Float64}}          #face centroids of the wall boundaries
    wall_outwardNormals::Vector{Matrix{Float64}}         #face outward normals of the wall boundaries
    wall_H::Vector{Vector{Float64}}                     #water depth of the wall boundaries
    wall_A::Vector{Vector{Float64}}                     #cross-sectional area of the wall boundaries

     # Symmetry boundary data
    symm_faceIDs::Vector{Vector{Int}}                     #face IDs of the symmetry boundaries
    symm_ghostCellIDs::Vector{Vector{Int}}               #ghost cell IDs of the symmetry boundaries
    symm_internalCellIDs::Vector{Vector{Int}}           #internal cell IDs of the symmetry boundaries
    symm_faceCentroids::Vector{Matrix{Float64}}          #face centroids of the symmetry boundaries
    symm_outwardNormals::Vector{Matrix{Float64}}         #face outward normals of the symmetry boundaries
    symm_H::Vector{Vector{Float64}}                     #water depth of the symmetry boundaries
    symm_A::Vector{Vector{Float64}}                     #cross-sectional area of the symmetry boundaries

    # Inner constructor with keyword arguments
    function BoundaryConditions2D(;
        nInletQ_BCs::Int=0,
        nExitH_BCs::Int=0,
        nWall_BCs::Int=0,
        nSymm_BCs::Int=0,
        inletQ_BC_indices::Vector{Int}=Int[],
        exitH_BC_indices::Vector{Int}=Int[],
        wall_BC_indices::Vector{Int}=Int[],
        symm_BC_indices::Vector{Int}=Int[],
        inletQ_faceIDs::Vector{Vector{Int}}=Vector{Int}[],
        inletQ_ghostCellIDs::Vector{Vector{Int}}=Vector{Int}[],
        inletQ_internalCellIDs::Vector{Vector{Int}}=Vector{Int}[],
        inletQ_faceCentroids::Vector{Matrix{Float64}}=Matrix{Float64}[],
        inletQ_faceOutwardNormals::Vector{Matrix{Float64}}=Matrix{Float64}[],
        inletQ_TotalQ::Vector{Float64}=Float64[],
        inletQ_H::Vector{Vector{Float64}}=Vector{Float64}[],
        inletQ_A::Vector{Vector{Float64}}=Vector{Float64}[],
        inletQ_ManningN::Vector{Vector{Float64}}=Vector{Float64}[],
        inletQ_Length::Vector{Vector{Float64}}=Vector{Float64}[],
        inletQ_TotalA::Vector{Float64}=Float64[],
        inletQ_DryWet::Vector{Vector{Int}}=Vector{Int}[],
        exitH_faceIDs::Vector{Vector{Int}}=Vector{Int}[],
        exitH_ghostCellIDs::Vector{Vector{Int}}=Vector{Int}[],
        exitH_internalCellIDs::Vector{Vector{Int}}=Vector{Int}[],
        exitH_faceCentroids::Vector{Matrix{Float64}}=Matrix{Float64}[],
        exitH_faceOutwardNormals::Vector{Matrix{Float64}}=Matrix{Float64}[],
        exitH_WSE::Vector{Float64}=Float64[],
        exitH_H::Vector{Vector{Float64}}=Vector{Float64}[],
        exitH_A::Vector{Vector{Float64}}=Vector{Float64}[],
        wall_faceIDs::Vector{Vector{Int}}=Vector{Int}[],
        wall_ghostCellIDs::Vector{Vector{Int}}=Vector{Int}[],
        wall_internalCellIDs::Vector{Vector{Int}}=Vector{Int}[],
        wall_faceCentroids::Vector{Matrix{Float64}}=Matrix{Float64}[],
        wall_outwardNormals::Vector{Matrix{Float64}}=Matrix{Float64}[],
        wall_H::Vector{Vector{Float64}}=Vector{Float64}[],
        wall_A::Vector{Vector{Float64}}=Vector{Float64}[],
        symm_faceIDs::Vector{Vector{Int}}=Vector{Int}[],
        symm_ghostCellIDs::Vector{Vector{Int}}=Vector{Int}[],
        symm_internalCellIDs::Vector{Vector{Int}}=Vector{Int}[],
        symm_faceCentroids::Vector{Matrix{Float64}}=Matrix{Float64}[],
        symm_outwardNormals::Vector{Matrix{Float64}}=Matrix{Float64}[],
        symm_H::Vector{Vector{Float64}}=Vector{Float64}[],
        symm_A::Vector{Vector{Float64}}=Vector{Float64}[]        
    )
        return new(nInletQ_BCs, nExitH_BCs, nWall_BCs, nSymm_BCs, inletQ_BC_indices, exitH_BC_indices, wall_BC_indices, symm_BC_indices, 
                   inletQ_faceIDs, inletQ_ghostCellIDs, inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, 
                   inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length, inletQ_TotalA, inletQ_DryWet, exitH_faceIDs, exitH_ghostCellIDs,
                   exitH_internalCellIDs, exitH_faceCentroids, exitH_faceOutwardNormals, exitH_WSE, exitH_H, exitH_A, 
                   wall_faceIDs, wall_ghostCellIDs, wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals, wall_H, wall_A, 
                   symm_faceIDs, symm_ghostCellIDs, symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals, symm_H, symm_A)
    end
end

 # create boundary conditions
 function initialize_boundary_conditions_2D(srh_all_Dict)
    
    # Get boundary counts
    #nInletQ = srh_all_Dict["nInletQ_BCs"]
    #nExitH = srh_all_Dict["nExitH_BCs"]
    #nWall = srh_all_Dict["nWall_BCs"]
    #nSymm = srh_all_Dict["nSymm_BCs"]

    #inletQ_BC_indices = Int[]
    #exitH_BC_indices = Int[]
    #wall_BC_indices = Int[]
    #symm_BC_indices = Int[]

    #create an empty BoundaryConditions object
    boundary_conditions = BoundaryConditions2D()

    #compute the indices of each boundary in the global boundary list
    compute_boundary_indices_2D!(boundary_conditions, srh_all_Dict)

    #preprocess all boundary conditions
    preprocess_all_boundaries_2D(boundary_conditions, srh_all_Dict)
    
    return boundary_conditions
 end

 #compute the indices of each boundary in the global boundary list
 #populate inletQ_BC_indices, exitH_BC_indices, wall_BC_indices, symm_BC_indices
 function compute_boundary_indices_2D!(boundary_conditions, srh_all_Dict)

    my_mesh_2D = srh_all_Dict["my_mesh_2D"]

    srhhydro_BC = srh_all_Dict["srhhydro_BC"]   #boundary conditions
    srhhydro_IQParams = srh_all_Dict["srhhydro_IQParams"]   #inlet-Q parameters
    srhhydro_EWSParamsC = srh_all_Dict["srhhydro_EWSParamsC"]   #exit-H parameters

    boundary_conditions.nInletQ_BCs = srh_all_Dict["nInletQ_BCs"]
    boundary_conditions.nExitH_BCs = srh_all_Dict["nExitH_BCs"]
    boundary_conditions.nWall_BCs = srh_all_Dict["nWall_BCs"]
    boundary_conditions.nSymm_BCs = srh_all_Dict["nSymm_BCs"]

    #loop over all boundaries to compute the indices of each boundary in the global boundary list
    for iBoundary in 1:my_mesh_2D.numOfBoundaries
        if haskey(srhhydro_BC, iBoundary)
            println("Key iBoundary exists in the Python dictionary srhhydro_BC.")
            println("The boundary: ", srhhydro_BC[iBoundary])
        else
            println("Key iBoundary does not exist in the Python dictionary srhhydro_BC.")
            readline()
            exit(-1)
        end

        boundaryType = srhhydro_BC[iBoundary]   #boundary type: wall, inletQ, exitH

        if lowercase(boundaryType) == "inlet-q"
            println("INLET-Q boundary condition is set for boundary ", iBoundary)

            push!(boundary_conditions.inletQ_BC_indices, iBoundary)

            #find the corresponding Q in IQParams
            if haskey(srhhydro_IQParams, iBoundary)
                println("Key IQParams exists in the Python dictionary srhhydro_IQParams.")
                println("Inlet specific discharge: ", srhhydro_IQParams[iBoundary])

                #update the boundary value: currenlty only support constant discharge
                #Q_value = parse(Float64, srhhydro_IQParams[iBoundary][1])
                #println("Inlet discharge: ", Q_value)

            else
                println("Key IQParams does not exist in the Python dictionary srhhydro_IQParams.")
                readline()
                exit(-1)
            end

        elseif lowercase(boundaryType) == "exit-h"
            println("EXIT-H boundary condition is set for boundary ", iBoundary)

            push!(boundary_conditions.exitH_BC_indices, iBoundary)

            #find the corresponding H in EWSParamsC
            if haskey(srhhydro_EWSParamsC, iBoundary)
                println("Key EWSParamsC exists in the Python dictionary srhhydro_EWSParamsC.")
                println("Exit water depth: ", srhhydro_EWSParamsC[iBoundary])

                #update the boundary value: currenlty only support constant water depth
                #H_value = parse(Float64, srhhydro_EWSParamsC[iBoundary][1])
                #println("Exit water depth: ", H_value)

            else
                println("Key EWSParamsC does not exist in the Python dictionary srhhydro_EWSParamsC.")
                readline()
                exit(-1)
            end

        elseif lowercase(boundaryType) == "wall"
            println("WALL boundary condition is set for boundary ", iBoundary)

            push!(boundary_conditions.wall_BC_indices, iBoundary)

        elseif lowercase(boundaryType) == "symm"
            println("SYMMETRY boundary condition is set for boundary ", iBoundary)

            push!(boundary_conditions.symm_BC_indices, iBoundary)
        end

    end

    #number of wall boundaries
    nWall_BCs = length(boundary_conditions.wall_BC_indices)
    srh_all_Dict["nWall_BCs"] = nWall_BCs  #update the number of wall boundaries in the dictionary

    #number of symmetry boundaries
    nSymm_BCs = length(boundary_conditions.symm_BC_indices)
    srh_all_Dict["nSymm_BCs"] = nSymm_BCs  #update the number of symmetry boundaries in the dictionary

    #println("inletQ_BC_indices: ", inletQ_BC_indices)
    #println("exitH_BC_indices: ", exitH_BC_indices) 
    #println("wall_BC_indices: ", wall_BC_indices)
    #println("symm_BC_indices: ", symm_BC_indices)

end


#preprocess all boundary conditions
function preprocess_all_boundaries_2D(boundary_conditions, srh_all_Dict)   

    my_mesh_2D = srh_all_Dict["my_mesh_2D"]

    #preprocess inlet-q boundaries
    preprocess_inlet_q_boundaries_2D(srh_all_Dict, boundary_conditions)

    #preprocess exit-h boundaries
    preprocess_exit_h_boundaries_2D(srh_all_Dict, boundary_conditions)

    #preprocess wall boundaries
    preprocess_wall_boundaries_2D(srh_all_Dict, boundary_conditions)

    #preprocess symmetry boundaries
    preprocess_symmetry_boundaries_2D(srh_all_Dict, boundary_conditions)

end

#preprocess inlet-q boundaries
function preprocess_inlet_q_boundaries_2D(srh_all_Dict, boundary_conditions)

    my_mesh_2D = srh_all_Dict["my_mesh_2D"]

    nInletQ_BCs = boundary_conditions.nInletQ_BCs
    srhhydro_IQParams = srh_all_Dict["srhhydro_IQParams"]
    inletQ_BC_indices = boundary_conditions.inletQ_BC_indices
    inletQ_faceIDs = boundary_conditions.inletQ_faceIDs
    inletQ_ghostCellIDs = boundary_conditions.inletQ_ghostCellIDs
    inletQ_internalCellIDs = boundary_conditions.inletQ_internalCellIDs

    #inletQ_faceCentroids = boundary_conditions.inletQ_faceCentroids
    #inletQ_faceOutwardNormals = boundary_conditions.inletQ_faceOutwardNormals
    #inletQ_TotalQ = boundary_conditions.inletQ_TotalQ
    #inletQ_H = boundary_conditions.inletQ_H
    #inletQ_A = boundary_conditions.inletQ_A
    #inletQ_ManningN = boundary_conditions.inletQ_ManningN
    #inletQ_Length = boundary_conditions.inletQ_Length
    #inletQ_TotalA = boundary_conditions.inletQ_TotalA
    #inletQ_DryWet = boundary_conditions.inletQ_DryWet

    face_normals = my_mesh_2D.face_normals   #face normals
    boundaryFaces_direction_Dict = my_mesh_2D.boundaryFaces_direction_Dict   #face directions of boundary faces

    #loop through all inlet-q boundaries
    for iInletQ in 1:nInletQ_BCs
        iBoundary = inletQ_BC_indices[iInletQ]

        println("Preprocessing INLET-Q boundary ", iInletQ, " with index in BC list ", iBoundary)

        #get the boundary face IDs of the current inlet-q boundary
        current_boundaryFaceIDs = my_mesh_2D.boundaryFaces_Dict[iBoundary]
        current_boundaryFace_directions = boundaryFaces_direction_Dict[iBoundary]
        current_ghostCellIDs = [my_mesh_2D.boundaryFaceID_to_ghostCellID_Dict[faceID] for faceID in current_boundaryFaceIDs]
        current_internalCellIDs = [my_mesh_2D.boundaryFaceID_to_internalCellID_Dict[faceID] for faceID in current_boundaryFaceIDs]

        boundary_conditions.inletQ_faceIDs[iInletQ] = current_boundaryFaceIDs   #the face ID of the current inlet-q boundary
        boundary_conditions.inletQ_ghostCellIDs[iInletQ] = current_ghostCellIDs   #the ghost cell ID of the current inlet-q boundary
        boundary_conditions.inletQ_internalCellIDs[iInletQ] = current_internalCellIDs   #the internal cell ID of the current inlet-q boundary    

        #number of bounary faces for the current inlet-q boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)

        boundary_conditions.inletQ_faceCentroids[iInletQ] = zeros(Float64, nBoundaryFaces, 3)   #face centroids of the current inlet-q boundary
        boundary_conditions.inletQ_faceOutwardNormals[iInletQ] = zeros(Float64, nBoundaryFaces, 2)   #face outward normals of the current inlet-q boundary
        boundary_conditions.inletQ_Length[iInletQ] = zeros(Float64, nBoundaryFaces)   #length for each face in the current inlet-q boundary

        #loop through all faces in the current inlet-q boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]

            #get the face centroid of the current inlet-q boundary
            faceCentroid = (my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[faceID][1], :] + my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[faceID][2], :]) / 2.0
            boundary_conditions.inletQ_faceCentroids[iInletQ][iFace, :] = faceCentroid

            #get the face direction of the current inlet-q boundary
            face_direction = current_boundaryFace_directions[iFace]

            #get the face outward normal of the current inlet-q boundary
            face_normal = face_normals[faceID]
            boundary_conditions.inletQ_faceOutwardNormals[iInletQ][iFace, :] = face_direction .* face_normal

            #get the length of the current face
            boundary_conditions.inletQ_Length[iInletQ][iFace] = my_mesh_2D.face_lengths[faceID]

            #println("Face ID: ", faceID, ", Ghost cell ID: ", ghostCellID, ", Internal cell ID: ", internalCellID, ", Face centroid: ", inletQ_faceCentroids[iInletQ][iFace, :], ", Face normal: ", inletQ_faceOutwardNormals[iInletQ][iFace, :], ", Face length: ", inletQ_Length[iInletQ][iFace])
        end

        #get the total discharge of the current inlet-q boundary
        boundary_conditions.inletQ_TotalQ[iInletQ] = parse(Float64, srhhydro_IQParams[iBoundary][1])

        #define some variables for the current inlet-q boundary
        boundary_conditions.inletQ_H[iInletQ] = zeros(Float64, nBoundaryFaces)   #inlet water depth for the current inlet-q boundary
        boundary_conditions.inletQ_A[iInletQ] = zeros(Float64, nBoundaryFaces)   #inlet cross-sectional area for the current inlet-q boundary
        boundary_conditions.inletQ_ManningN[iInletQ] = zeros(Float64, nBoundaryFaces)   #Manning's n for the current inlet-q boundary

        boundary_conditions.inletQ_TotalA[iInletQ] = 0.0   #total cross-sectional area for the current inlet-q boundary
        boundary_conditions.inletQ_DryWet[iInletQ] = zeros(Int, nBoundaryFaces)   #dry(=0)/wet(=1) flag for the current inlet-q boundary
    end

end

#preprocess exit-h boundaries
function preprocess_exit_h_boundaries_2D(srh_all_Dict, boundary_conditions)

    my_mesh_2D = srh_all_Dict["my_mesh_2D"]

    nExitH_BCs = boundary_conditions.nExitH_BCs
    srhhydro_EWSParamsC = srh_all_Dict["srhhydro_EWSParamsC"]
    exitH_BC_indices = boundary_conditions.exitH_BC_indices
    exitH_faceIDs = boundary_conditions.exitH_faceIDs
    exitH_ghostCellIDs = boundary_conditions.exitH_ghostCellIDs
    exitH_internalCellIDs = boundary_conditions.exitH_internalCellIDs

    #exitH_faceCentroids = boundary_conditions.exitH_faceCentroids
    #exitH_faceOutwardNormals = boundary_conditions.exitH_faceOutwardNormals
    #exitH_WSE = boundary_conditions.exitH_WSE
    #exitH_H = boundary_conditions.exitH_H
    #exitH_A = boundary_conditions.exitH_A

    face_normals = my_mesh_2D.face_normals   #face normals
    boundaryFaces_direction_Dict = my_mesh_2D.boundaryFaces_direction_Dict   #face directions of boundary faces

    for iExitH in 1:nExitH_BCs
        iBoundary = exitH_BC_indices[iExitH]

        println("Preprocessing EXIT-H boundary ", iExitH, " with index in BC list ", iBoundary)

        #get the boundary face IDs of the current exit-h boundary
        current_boundaryFaceIDs = my_mesh_2D.boundaryFaces_Dict[iBoundary]
        current_boundaryFace_directions = boundaryFaces_direction_Dict[iBoundary]
        current_ghostCellIDs = [my_mesh_2D.boundaryFaceID_to_ghostCellID_Dict[abs(faceID)] for faceID in current_boundaryFaceIDs]
        current_internalCellIDs = [my_mesh_2D.boundaryFaceID_to_internalCellID_Dict[abs(faceID)] for faceID in current_boundaryFaceIDs]

        boundary_conditions.exitH_faceIDs[iExitH] = current_boundaryFaceIDs   #the face ID of the current exit-h boundary
        boundary_conditions.exitH_ghostCellIDs[iExitH] = current_ghostCellIDs   #the ghost cell ID of the current exit-h boundary
        boundary_conditions.exitH_internalCellIDs[iExitH] = current_internalCellIDs   #the internal cell ID of the current exit-h boundary    

        #number of bounary faces for the current exit-h boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)

        boundary_conditions.exitH_faceCentroids[iExitH] = zeros(Float64, nBoundaryFaces, 3)   #face centroids of the current exit-h boundary
        boundary_conditions.exitH_faceOutwardNormals[iExitH] = zeros(Float64, nBoundaryFaces, 2)   #face outward normals of the current exit-h boundary

        #loop through all faces in the current exit-h boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]

            #get the face centroid of the current exit-h boundary
            faceCentroid = (my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][1], :] + my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][2], :]) / 2.0
            boundary_conditions.exitH_faceCentroids[iExitH][iFace, :] = faceCentroid

            #get the face direction of the current inlet-q boundary
            face_direction = current_boundaryFace_directions[iFace]

            #get the face outward normal of the current inlet-q boundary
            face_normal = face_normals[faceID]
            boundary_conditions.exitH_faceOutwardNormals[iExitH][iFace, :] = face_direction .* face_normal

            #println("Face ID: ", faceID, ", Ghost cell ID: ", ghostCellID, ", Internal cell ID: ", internalCellID, ", Face centroid: ", exitH_faceCentroids[iExitH][iFace, :], ", Face normal: ", exitH_faceOutwardNormals[iExitH][iFace, :])
        end

        #get the water surface elevation of the current exit-h boundary
        boundary_conditions.exitH_WSE[iExitH] = parse(Float64, srhhydro_EWSParamsC[iBoundary][1])

        #get the ghost cell ID of the current exit-h boundary
        boundary_conditions.exitH_H[iExitH] = zeros(Float64, nBoundaryFaces)   #water depth for the current exit-h boundary   
        boundary_conditions.exitH_A[iExitH] = zeros(Float64, nBoundaryFaces)   #cross-sectional area for the current exit-h boundary
    end

end

#preprocess wall boundaries
function preprocess_wall_boundaries_2D(srh_all_Dict, boundary_conditions)

    my_mesh_2D = srh_all_Dict["my_mesh_2D"]

    nWall_BCs = boundary_conditions.nWall_BCs
    wall_BC_indices = boundary_conditions.wall_BC_indices
    wall_faceIDs = boundary_conditions.wall_faceIDs
    wall_ghostCellIDs = boundary_conditions.wall_ghostCellIDs
    wall_internalCellIDs = boundary_conditions.wall_internalCellIDs

    #wall_faceCentroids = boundary_conditions.wall_faceCentroids
    #wall_outwardNormals = boundary_conditions.wall_outwardNormals
    #wall_H = boundary_conditions.wall_H
    #wall_A = boundary_conditions.wall_A

    face_normals = my_mesh_2D.face_normals   #face normals
    boundaryFaces_direction_Dict = my_mesh_2D.boundaryFaces_direction_Dict   #face directions of boundary faces

    #loop through all wall boundaries
    for iWall in 1:nWall_BCs
        iBoundary = wall_BC_indices[iWall]

        println("Preprocessing WALL boundary ", iWall, " with index in BC list ", iBoundary)

        #get the boundary face IDs of the current wall boundary
        current_boundaryFaceIDs = my_mesh_2D.boundaryFaces_Dict[iBoundary]
        current_boundaryFace_directions = boundaryFaces_direction_Dict[iBoundary]
        current_ghostCellIDs = [my_mesh_2D.boundaryFaceID_to_ghostCellID_Dict[abs(faceID)] for faceID in current_boundaryFaceIDs]
        current_internalCellIDs = [my_mesh_2D.boundaryFaceID_to_internalCellID_Dict[abs(faceID)] for faceID in current_boundaryFaceIDs]

        boundary_conditions.wall_faceIDs[iWall] = current_boundaryFaceIDs   #the face ID of the current wall boundary
        boundary_conditions.wall_ghostCellIDs[iWall] = current_ghostCellIDs   #the ghost cell ID of the current wall boundary
        boundary_conditions.wall_internalCellIDs[iWall] = current_internalCellIDs   #the internal cell ID of the current wall boundary    

        #number of bounary faces for the current wall boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)

        boundary_conditions.wall_faceCentroids[iWall] = zeros(Float64, nBoundaryFaces, 3)   #face centroids of the current wall boundary
        boundary_conditions.wall_outwardNormals[iWall] = zeros(Float64, nBoundaryFaces, 2)   #face outward normals of the current wall boundary


        #loop through all faces in the current wall boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]

            #get the face centroid of the current wall boundary
            faceCentroid = (my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][1], :] + my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][2], :]) / 2.0
            boundary_conditions.wall_faceCentroids[iWall][iFace, :] = faceCentroid

            #get the face direction of the current inlet-q boundary
            face_direction = current_boundaryFace_directions[iFace]

            #get the face outward normal of the current inlet-q boundary
            face_normal = face_normals[faceID]
            boundary_conditions.wall_outwardNormals[iWall][iFace, :] = face_direction .* face_normal

            #println("Face ID: ", faceID, ", Ghost cell ID: ", ghostCellID, ", Internal cell ID: ", internalCellID, ", Face centroid: ", wall_faceCentroids[iWall][iFace, :], ", Face normal: ", wall_outwardNormals[iWall][iFace, :])
        end

        #get the ghost cell ID of the current wall boundary
        boundary_conditions.wall_H[iWall] = zeros(Float64, nBoundaryFaces)   # water depth for the current wall boundary   
        boundary_conditions.wall_A[iWall] = zeros(Float64, nBoundaryFaces)   # cross-sectional area for the current wall boundary  
    end

end

#preprocess symmetry boundaries: called every time step to update the boundary condition
function preprocess_symmetry_boundaries_2D(srh_all_Dict, boundary_conditions)

    my_mesh_2D = srh_all_Dict["my_mesh_2D"]

    nSymm_BCs = boundary_conditions.nSymm_BCs
    symm_BC_indices = boundary_conditions.symm_BC_indices
    symm_faceIDs = boundary_conditions.symm_faceIDs
    symm_ghostCellIDs = boundary_conditions.symm_ghostCellIDs
    symm_internalCellIDs = boundary_conditions.symm_internalCellIDs

    #symm_faceCentroids = boundary_conditions.symm_faceCentroids
    #symm_outwardNormals = boundary_conditions.symm_outwardNormals
    #symm_H = boundary_conditions.symm_H
    #symm_A = boundary_conditions.symm_A

    face_normals = my_mesh_2D.face_normals   #face normals
    boundaryFaces_direction_Dict = my_mesh_2D.boundaryFaces_direction_Dict   #face directions of boundary faces

    #loop through all symmetry boundaries
    for iSymm in 1:nSymm_BCs
        iBoundary = symm_BC_indices[iSymm]

        println("Preprocessing SYMMETRY boundary ", iSymm, " with index in BC list ", iBoundary)

        #get the boundary face IDs of the current wall boundary
        current_boundaryFaceIDs = my_mesh_2D.boundaryFaces_Dict[iBoundary]
        current_boundaryFace_directions = boundaryFaces_direction_Dict[iBoundary]
        current_ghostCellIDs = [my_mesh_2D.boundaryFaceID_to_ghostCellID_Dict[abs(faceID)] for faceID in current_boundaryFaceIDs]
        current_internalCellIDs = [my_mesh_2D.boundaryFaceID_to_internalCellID_Dict[abs(faceID)] for faceID in current_boundaryFaceIDs]

        boundary_conditions.symm_faceIDs[iSymm] = current_boundaryFaceIDs   #the face ID of the current symmetry boundary
        boundary_conditions.symm_ghostCellIDs[iSymm] = current_ghostCellIDs   #the ghost cell ID of the current symmetry boundary
        boundary_conditions.symm_internalCellIDs[iSymm] = current_internalCellIDs   #the internal cell ID of the current symmetry boundary    

        #number of bounary faces for the current symmetry boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)

        boundary_conditions.symm_faceCentroids[iSymm] = zeros(Float64, nBoundaryFaces, 3)   #face centroids of the current symmetry boundary
        boundary_conditions.symm_outwardNormals[iSymm] = zeros(Float64, nBoundaryFaces, 2)   #face outward normals of the current symmetry boundary

        #loop through all faces in the current symmetry boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]

            #get the face centroid of the current symmetry boundary
            faceCentroid = (my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][1], :] + my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][2], :]) / 2.0
            boundary_conditions.symm_faceCentroids[iSymm][iFace, :] = faceCentroid

            #get the face direction of the current symmetry boundary
            face_direction = current_boundaryFace_directions[iFace]

            #get the face outward normal of the current symmetry boundary
            face_normal = face_normals[faceID]
            boundary_conditions.symm_outwardNormals[iSymm][iFace, :] = face_direction .* face_normal

            #println("Face ID: ", faceID, ", Ghost cell ID: ", ghostCellID, ", Internal cell ID: ", internalCellID, ", Face centroid: ", symm_faceCentroids[iSymm][iFace, :], ", Face normal: ", symm_outwardNormals[iSymm][iFace, :])
        end

        #get the ghost cell ID of the current symmetry boundary
        boundary_conditions.symm_H[iSymm] = zeros(Float64, nBoundaryFaces)   # water depth for the current symmetry boundary   
        boundary_conditions.symm_A[iSymm] = zeros(Float64, nBoundaryFaces)   # cross-sectional area for the current symmetry boundary  
    end

end
