#boundary conditions for 2D SWE (data within are immutable)

Base.@kwdef struct BoundaryConditions2D
    # Numbers of each boundary type
    nInletQ_BCs::Int = 0      #number of inlet-Q boundaries
    nExitH_BCs::Int = 0       #number of exit-H boundaries
    nWall_BCs::Int = 0        #number of wall boundaries
    nSymm_BCs::Int = 0        #number of symmetry boundaries

    #indices of each boundary in the boundary list
    inletQ_BC_indices::Vector{Int} = Int[]   #indices of inlet-Q boundaries in the boundary list
    exitH_BC_indices::Vector{Int} = Int[]    #indices of exit-H boundaries in the boundary list
    wall_BC_indices::Vector{Int} = Int[]     #indices of wall boundaries in the boundary list
    symm_BC_indices::Vector{Int} = Int[]     #indices of symmetry (no slip) boundaries in the boundary list

     # Inlet Q boundary data
    inletQ_faceIDs::Union{Vector{Vector{Int}}, Nothing} = nothing                     #face IDs of the inlet-q boundaries
    inletQ_ghostCellIDs::Union{Vector{Vector{Int}}, Nothing} = nothing               #ghost cell IDs of the inlet-q boundaries
    inletQ_internalCellIDs::Union{Vector{Vector{Int}}, Nothing} = nothing           #internal cell IDs of the inlet-q boundaries
    #inletQ_faceCentroids::Vector{Matrix{Float64}} = Matrix{Float64}[]          #face centroids of the inlet-q boundaries
    inletQ_faceOutwardNormals::Union{Vector{Matrix{Float64}}, Nothing} = nothing     #face outward normals of the inlet-q boundaries

    # Exit H boundary data
    exitH_faceIDs::Union{Vector{Vector{Int}}, Nothing} = nothing                     #face IDs of the exit-h boundaries
    exitH_ghostCellIDs::Union{Vector{Vector{Int}}, Nothing} = nothing                 #ghost cell IDs of the exit-h boundaries
    exitH_internalCellIDs::Union{Vector{Vector{Int}}, Nothing} = nothing             #internal cell IDs of the exit-h boundaries
    exitH_faceOutwardNormals::Union{Vector{Matrix{Float64}}, Nothing} = nothing     #face outward normals of the exit-h boundaries

    # Wall boundary data
    wall_faceIDs::Union{Vector{Vector{Int}}, Nothing} = nothing                     #face IDs of the wall boundaries
    wall_ghostCellIDs::Union{Vector{Vector{Int}}, Nothing} = nothing               #ghost cell IDs of the wall boundaries
    wall_internalCellIDs::Union{Vector{Vector{Int}}, Nothing} = nothing           #internal cell IDs of the wall boundaries
    wall_faceCentroids::Union{Vector{Matrix{Float64}}, Nothing} = nothing          #face centroids of the wall boundaries
    wall_outwardNormals::Union{Vector{Matrix{Float64}}, Nothing} = nothing         #face outward normals of the wall boundaries

     # Symmetry boundary data
    symm_faceIDs::Union{Vector{Vector{Int}}, Nothing} = nothing                     #face IDs of the symmetry boundaries
    symm_ghostCellIDs::Union{Vector{Vector{Int}}, Nothing} = nothing               #ghost cell IDs of the symmetry boundaries
    symm_internalCellIDs::Union{Vector{Vector{Int}}, Nothing} = nothing           #internal cell IDs of the symmetry boundaries
    symm_faceCentroids::Union{Vector{Matrix{Float64}}, Nothing} = nothing          #face centroids of the symmetry boundaries
    symm_outwardNormals::Union{Vector{Matrix{Float64}}, Nothing} = nothing         #face outward normals of the symmetry boundaries

end

 # create boundary conditions
 function initialize_boundary_conditions_2D(settings, srh_all_Dict, nodeCoordinates)

    inletQ_BC_indices = []
    exitH_BC_indices = []
    wall_BC_indices = []
    symm_BC_indices = []

    #compute the indices of each boundary in the global boundary list
    compute_boundary_indices_2D!(settings, inletQ_BC_indices, exitH_BC_indices, wall_BC_indices, symm_BC_indices, srh_all_Dict)

    # Get boundary counts
    nInletQ_BCs = srh_all_Dict["nInletQ_BCs"]
    nExitH_BCs = srh_all_Dict["nExitH_BCs"]
    nWall_BCs = srh_all_Dict["nWall_BCs"]
    nSymm_BCs = srh_all_Dict["nSymm_BCs"]
   
    #inlet_q arrays
    inletQ_faceIDs = nInletQ_BCs > 0 ? Array{Array{Int}}(undef, nInletQ_BCs) : [[0]]   #face IDs of the inlet-q boundaries
    inletQ_ghostCellIDs = nInletQ_BCs > 0 ? Array{Array{Int}}(undef, nInletQ_BCs) : [[0]]   #ghost cell IDs of the inlet-q boundaries
    inletQ_internalCellIDs = nInletQ_BCs > 0 ? Array{Array{Int}}(undef, nInletQ_BCs) : [[0]]   #internal cell IDs of the inlet-q boundaries

    inletQ_faceOutwardNormals = nInletQ_BCs > 0 ? Vector{Matrix{Float64}}(undef, nInletQ_BCs) : [[0.0 0.0; 0.0 0.0]]   #face outward normals of the inlet-q boundaries

    inletQ_TotalQ = nInletQ_BCs > 0 ? Array{Float64}(undef, nInletQ_BCs) : [0.0]   #total discharge for each inlet-q boundary
    inletQ_H = nInletQ_BCs > 0 ? Array{Array{Float64}}(undef, nInletQ_BCs) : [[0.0]]   #inlet water depth for each face in each inlet-q boundary
    inletQ_A = nInletQ_BCs > 0 ? Array{Array{Float64}}(undef, nInletQ_BCs) : [[0.0]]   #area for each face in each inlet-q boundary
    inletQ_ManningN = nInletQ_BCs > 0 ? Array{Array{Float64}}(undef, nInletQ_BCs) : [[0.0]]   #Manning's n for each face in each inlet-q boundary
    inletQ_Length = nInletQ_BCs > 0 ? Array{Array{Float64}}(undef, nInletQ_BCs) : [[0.0]]   #length for each face in each inlet-q boundary
    inletQ_TotalA = nInletQ_BCs > 0 ? Array{Float64}(undef, nInletQ_BCs) : [0.0]   #total cross-sectional area for each inlet-q boundary
    inletQ_DryWet = nInletQ_BCs > 0 ? Array{Array{Int}}(undef, nInletQ_BCs) : [[0]]   #dry(=0)/wet(=1) flag for each inlet-q boundary

    #exit_h arrays
    exitH_faceIDs = nExitH_BCs > 0 ? Array{Array{Int}}(undef, nExitH_BCs) : nothing   #face IDs of the exit-h boundaries
    exitH_ghostCellIDs = nExitH_BCs > 0 ? Array{Array{Int}}(undef, nExitH_BCs) : nothing   #ghost cell IDs of the exit-h boundaries
    exitH_internalCellIDs = nExitH_BCs > 0 ? Array{Array{Int}}(undef, nExitH_BCs) : nothing   #internal cell IDs of the exit-h boundaries

    exitH_faceOutwardNormals = nExitH_BCs > 0 ? Vector{Matrix{Float64}}(undef, nExitH_BCs) : nothing   #face outward normals of the exit-h boundaries

    exitH_WSE = nExitH_BCs > 0 ? Array{Float64}(undef, nExitH_BCs) : nothing   #WSE for each exit-h boundary
    exitH_H = nExitH_BCs > 0 ? Array{Array{Float64}}(undef, nExitH_BCs) : nothing   #inlet water depth for each exit-h boundary
    exitH_A = nExitH_BCs > 0 ? Array{Array{Float64}}(undef, nExitH_BCs) : nothing   #inlet cross-sectional area for each exit-h boundary

    #wall arrays
    wall_faceIDs = nWall_BCs > 0 ? Array{Array{Int}}(undef, nWall_BCs) : nothing   #face IDs of the wall boundaries
    wall_ghostCellIDs = nWall_BCs > 0 ? Array{Array{Int}}(undef, nWall_BCs) : nothing   #ghost cell IDs of the wall boundaries
    wall_internalCellIDs = nWall_BCs > 0 ? Array{Array{Int}}(undef, nWall_BCs) : nothing   #internal cell IDs of the wall boundaries

    wall_faceCentroids = nWall_BCs > 0 ? Vector{Matrix{Float64}}(undef, nWall_BCs) : nothing   #face centroids of the wall boundaries
    wall_outwardNormals = nWall_BCs > 0 ? Vector{Matrix{Float64}}(undef, nWall_BCs) : nothing   #face outward normals of the wall boundaries

    wall_H = nWall_BCs > 0 ? Array{Array{Float64}}(undef, nWall_BCs) : nothing   #water depth for each wall boundary
    wall_A = nWall_BCs > 0 ? Array{Array{Float64}}(undef, nWall_BCs) : nothing   #cross-sectional area for each wall boundary

    #symmetry arrays
    symm_faceIDs = nSymm_BCs > 0 ? Array{Array{Int}}(undef, nSymm_BCs) : nothing
    symm_ghostCellIDs = nSymm_BCs > 0 ? Array{Array{Int}}(undef, nSymm_BCs) : nothing
    symm_internalCellIDs = nSymm_BCs > 0 ? Array{Array{Int}}(undef, nSymm_BCs) : nothing

    symm_faceCentroids = nSymm_BCs > 0 ? Vector{Matrix{Float64}}(undef, nSymm_BCs) : nothing
    symm_outwardNormals = nSymm_BCs > 0 ? Vector{Matrix{Float64}}(undef, nSymm_BCs) : nothing

    symm_H = nSymm_BCs > 0 ? Array{Array{Float64}}(undef, nSymm_BCs) : nothing
    symm_A = nSymm_BCs > 0 ? Array{Array{Float64}}(undef, nSymm_BCs) : nothing
    
    #preprocess all boundary conditions
    preprocess_all_boundaries_2D(settings, srh_all_Dict, nInletQ_BCs, nExitH_BCs, nWall_BCs, nSymm_BCs, inletQ_BC_indices, exitH_BC_indices, wall_BC_indices, symm_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, inletQ_internalCellIDs, inletQ_faceOutwardNormals, inletQ_Length, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_TotalA, inletQ_DryWet, exitH_faceIDs, exitH_ghostCellIDs, exitH_internalCellIDs, exitH_faceOutwardNormals, exitH_WSE, exitH_H, exitH_A, wall_faceIDs, wall_ghostCellIDs, wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals, wall_H, wall_A, symm_faceIDs, symm_ghostCellIDs, symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals, symm_H, symm_A, nodeCoordinates)

    boundary_conditions = BoundaryConditions2D(
        nInletQ_BCs=nInletQ_BCs, 
        nExitH_BCs=nExitH_BCs, 
        nWall_BCs=nWall_BCs, 
        nSymm_BCs=nSymm_BCs, 
        inletQ_BC_indices=inletQ_BC_indices, 
        exitH_BC_indices=exitH_BC_indices, 
        wall_BC_indices=wall_BC_indices, 
        symm_BC_indices=symm_BC_indices, 
        inletQ_faceIDs=inletQ_faceIDs, 
        inletQ_ghostCellIDs=inletQ_ghostCellIDs, 
        inletQ_internalCellIDs=inletQ_internalCellIDs, 
        inletQ_faceOutwardNormals=inletQ_faceOutwardNormals, 
        exitH_faceIDs=exitH_faceIDs, 
        exitH_ghostCellIDs=exitH_ghostCellIDs, 
        exitH_internalCellIDs=exitH_internalCellIDs, 
        exitH_faceOutwardNormals=exitH_faceOutwardNormals, 
        wall_faceIDs=wall_faceIDs, 
        wall_ghostCellIDs=wall_ghostCellIDs, 
        wall_internalCellIDs=wall_internalCellIDs, 
        wall_faceCentroids=wall_faceCentroids, 
        wall_outwardNormals=wall_outwardNormals,
        symm_faceIDs=symm_faceIDs, 
        symm_ghostCellIDs=symm_ghostCellIDs, 
        symm_internalCellIDs=symm_internalCellIDs, 
        symm_faceCentroids=symm_faceCentroids, 
        symm_outwardNormals=symm_outwardNormals
        )

    return boundary_conditions, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length, inletQ_TotalA, inletQ_DryWet, exitH_WSE, exitH_H, exitH_A, wall_H, wall_A, symm_H, symm_A
 end

 #compute the indices of each boundary in the global boundary list
 #populate inletQ_BC_indices, exitH_BC_indices, wall_BC_indices, symm_BC_indices
 function compute_boundary_indices_2D!(settings, inletQ_BC_indices, exitH_BC_indices, wall_BC_indices, symm_BC_indices, srh_all_Dict)

    my_mesh_2D = srh_all_Dict["my_mesh_2D"]

    srhhydro_BC = srh_all_Dict["srhhydro_BC"]   #boundary conditions
    srhhydro_IQParams = srh_all_Dict["srhhydro_IQParams"]   #inlet-Q parameters
    srhhydro_EWSParamsC = srh_all_Dict["srhhydro_EWSParamsC"]   #exit-H parameters

    #loop over all boundaries to compute the indices of each boundary in the global boundary list
    for iBoundary in 1:my_mesh_2D.numOfBoundaries
        
        if haskey(srhhydro_BC, iBoundary)
            if settings.bVerbose
                println("Key iBoundary exists in the Python dictionary srhhydro_BC.")
                println("The boundary: ", srhhydro_BC[iBoundary])
            end
        else
            println("Key iBoundary does not exist in the Python dictionary srhhydro_BC.")
            readline()
            exit(-1)
        end

        boundaryType = srhhydro_BC[iBoundary]   #boundary type: wall, inletQ, exitH

        if lowercase(boundaryType) == "inlet-q"
            if settings.bVerbose
                println("INLET-Q boundary condition is set for boundary ", iBoundary)
            end

            push!(inletQ_BC_indices, iBoundary)

            #find the corresponding Q in IQParams
            if haskey(srhhydro_IQParams, iBoundary)
                if settings.bVerbose
                    println("Key IQParams exists in the Python dictionary srhhydro_IQParams.")
                    println("Inlet specific discharge: ", srhhydro_IQParams[iBoundary])
                end

                #update the boundary value: currenlty only support constant discharge
                #Q_value = parse(Float64, srhhydro_IQParams[iBoundary][1])
                #println("Inlet discharge: ", Q_value)

            else
                println("Key IQParams does not exist in the Python dictionary srhhydro_IQParams.")
                readline()
                exit(-1)
            end

        elseif lowercase(boundaryType) == "exit-h"
            if settings.bVerbose
                println("EXIT-H boundary condition is set for boundary ", iBoundary)
            end

            push!(exitH_BC_indices, iBoundary)

            #find the corresponding H in EWSParamsC
            if haskey(srhhydro_EWSParamsC, iBoundary)
                if settings.bVerbose
                    println("Key EWSParamsC exists in the Python dictionary srhhydro_EWSParamsC.")
                    println("Exit water depth: ", srhhydro_EWSParamsC[iBoundary])
                end

                #update the boundary value: currenlty only support constant water depth
                #H_value = parse(Float64, srhhydro_EWSParamsC[iBoundary][1])
                #println("Exit water depth: ", H_value)

            else
                println("Key EWSParamsC does not exist in the Python dictionary srhhydro_EWSParamsC.")
                readline()
                exit(-1)
            end

        elseif lowercase(boundaryType) == "wall"
            if settings.bVerbose
                println("WALL boundary condition is set for boundary ", iBoundary)
            end

            push!(wall_BC_indices, iBoundary)

        elseif lowercase(boundaryType) == "symm"
            if settings.bVerbose
                println("SYMMETRY boundary condition is set for boundary ", iBoundary)
            end

            push!(symm_BC_indices, iBoundary)
        end

    end

    #number of wall boundaries
    nWall_BCs = length(wall_BC_indices)
    srh_all_Dict["nWall_BCs"] = nWall_BCs  #update the number of wall boundaries in the dictionary
    #nWall_BCs = nWall_BCs

    #number of symmetry boundaries
    nSymm_BCs = length(symm_BC_indices)
    srh_all_Dict["nSymm_BCs"] = nSymm_BCs  #update the number of symmetry boundaries in the dictionary
    #nSymm_BCs = nSymm_BCs

    if settings.bVerbose
        println("inletQ_BC_indices: ", inletQ_BC_indices)
        println("exitH_BC_indices: ", exitH_BC_indices) 
        println("wall_BC_indices: ", wall_BC_indices)
        println("symm_BC_indices: ", symm_BC_indices)
    end

end


#preprocess all boundary conditions
function preprocess_all_boundaries_2D(settings, srh_all_Dict, nInletQ_BCs, nExitH_BCs, nWall_BCs, nSymm_BCs, inletQ_BC_indices, exitH_BC_indices, wall_BC_indices, symm_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, inletQ_internalCellIDs, inletQ_faceOutwardNormals, inletQ_Length, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_TotalA, inletQ_DryWet, exitH_faceIDs, exitH_ghostCellIDs, exitH_internalCellIDs, exitH_faceOutwardNormals, exitH_WSE, exitH_H, exitH_A, wall_faceIDs, wall_ghostCellIDs, wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals, wall_H, wall_A, symm_faceIDs, symm_ghostCellIDs, symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals, symm_H, symm_A, nodeCoordinates)   

    #preprocess inlet-q boundaries
    preprocess_inlet_q_boundaries_2D(settings, srh_all_Dict, nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, inletQ_internalCellIDs, inletQ_faceOutwardNormals, inletQ_Length, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_TotalA, inletQ_DryWet, nodeCoordinates)

    #preprocess exit-h boundaries
    preprocess_exit_h_boundaries_2D(settings, srh_all_Dict, nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, exitH_internalCellIDs,exitH_faceOutwardNormals, exitH_WSE, exitH_H, exitH_A, nodeCoordinates)

    #preprocess wall boundaries
    preprocess_wall_boundaries_2D(settings, srh_all_Dict, nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs, wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals, wall_H, wall_A, nodeCoordinates)

    #preprocess symmetry boundaries
    preprocess_symmetry_boundaries_2D(settings, srh_all_Dict, nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs, symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals, symm_H, symm_A, nodeCoordinates)

end

#preprocess inlet-q boundaries
function preprocess_inlet_q_boundaries_2D(settings, srh_all_Dict, nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, inletQ_internalCellIDs, inletQ_faceOutwardNormals, inletQ_Length, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_TotalA, inletQ_DryWet, nodeCoordinates)

    my_mesh_2D = srh_all_Dict["my_mesh_2D"]
    
    srhhydro_IQParams = srh_all_Dict["srhhydro_IQParams"]
 
    face_normals = my_mesh_2D.face_normals   #face normals
    boundaryFaces_direction_Dict = my_mesh_2D.boundaryFaces_direction_Dict   #face directions of boundary faces

    #loop through all inlet-q boundaries
    for iInletQ in 1:nInletQ_BCs
        iBoundary = inletQ_BC_indices[iInletQ]

        if settings.bVerbose
            println("Preprocessing INLET-Q boundary ", iInletQ, " with index in BC list ", iBoundary)
        end

        #get the boundary face IDs of the current inlet-q boundary
        current_boundaryFaceIDs = my_mesh_2D.boundaryFaces_Dict[iBoundary]
        current_boundaryFace_directions = boundaryFaces_direction_Dict[iBoundary]
        current_ghostCellIDs = [my_mesh_2D.boundaryFaceID_to_ghostCellID_Dict[faceID] for faceID in current_boundaryFaceIDs]
        current_internalCellIDs = [my_mesh_2D.boundaryFaceID_to_internalCellID_Dict[faceID] for faceID in current_boundaryFaceIDs]

        inletQ_faceIDs[iInletQ] = current_boundaryFaceIDs   #the face ID of the current inlet-q boundary
        inletQ_ghostCellIDs[iInletQ] = current_ghostCellIDs   #the ghost cell ID of the current inlet-q boundary
        inletQ_internalCellIDs[iInletQ] = current_internalCellIDs   #the internal cell ID of the current inlet-q boundary    

        #number of bounary faces for the current inlet-q boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)

        #inletQ_faceCentroids[iInletQ] = zeros(Float64, nBoundaryFaces, 3)   #face centroids of the current inlet-q boundary
        inletQ_faceOutwardNormals[iInletQ] = zeros(Float64, nBoundaryFaces, 2)   #face outward normals of the current inlet-q boundary
        inletQ_Length[iInletQ] = zeros(Float64, nBoundaryFaces)   #length for each face in the current inlet-q boundary

        #loop through all faces in the current inlet-q boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]

            #get the face centroid of the current inlet-q boundary
            #faceCentroid = (nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[faceID][1], :] + nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[faceID][2], :]) / 2.0
            #inletQ_faceCentroids[iInletQ][iFace, :] = faceCentroid

            #get the face direction of the current inlet-q boundary
            face_direction = current_boundaryFace_directions[iFace]

            #get the face outward normal of the current inlet-q boundary
            face_normal = face_normals[faceID]
            inletQ_faceOutwardNormals[iInletQ][iFace, :] = face_direction .* face_normal

            #get the length of the current face
            inletQ_Length[iInletQ][iFace] = my_mesh_2D.face_lengths[faceID]

            #println("Face ID: ", faceID, ", Ghost cell ID: ", ghostCellID, ", Internal cell ID: ", internalCellID, ", Face centroid: ", inletQ_faceCentroids[iInletQ][iFace, :], ", Face normal: ", inletQ_faceOutwardNormals[iInletQ][iFace, :], ", Face length: ", inletQ_Length[iInletQ][iFace])
        end

        #get the total discharge of the current inlet-q boundary
        inletQ_TotalQ[iInletQ] = parse(Float64, srhhydro_IQParams[iBoundary][1])

        #define some variables for the current inlet-q boundary
        inletQ_H[iInletQ] = zeros(Float64, nBoundaryFaces)   #inlet water depth for the current inlet-q boundary
        inletQ_A[iInletQ] = zeros(Float64, nBoundaryFaces)   #inlet cross-sectional area for the current inlet-q boundary
        inletQ_ManningN[iInletQ] = zeros(Float64, nBoundaryFaces)   #Manning's n for the current inlet-q boundary

        inletQ_TotalA[iInletQ] = 0.0   #total cross-sectional area for the current inlet-q boundary
        inletQ_DryWet[iInletQ] = zeros(Int64, nBoundaryFaces)   #dry(=0)/wet(=1) flag for the current inlet-q boundary
    end

end

#preprocess exit-h boundaries
function preprocess_exit_h_boundaries_2D(settings, srh_all_Dict, nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, exitH_internalCellIDs,  exitH_faceOutwardNormals, exitH_WSE, exitH_H, exitH_A, nodeCoordinates)

    my_mesh_2D = srh_all_Dict["my_mesh_2D"]

    srhhydro_EWSParamsC = srh_all_Dict["srhhydro_EWSParamsC"]

    face_normals = my_mesh_2D.face_normals   #face normals
    boundaryFaces_direction_Dict = my_mesh_2D.boundaryFaces_direction_Dict   #face directions of boundary faces

    for iExitH in 1:nExitH_BCs
        iBoundary = exitH_BC_indices[iExitH]

        if settings.bVerbose
            println("Preprocessing EXIT-H boundary ", iExitH, " with index in BC list ", iBoundary)
        end

        #get the boundary face IDs of the current exit-h boundary
        current_boundaryFaceIDs = my_mesh_2D.boundaryFaces_Dict[iBoundary]
        current_boundaryFace_directions = boundaryFaces_direction_Dict[iBoundary]
        current_ghostCellIDs = [my_mesh_2D.boundaryFaceID_to_ghostCellID_Dict[abs(faceID)] for faceID in current_boundaryFaceIDs]
        current_internalCellIDs = [my_mesh_2D.boundaryFaceID_to_internalCellID_Dict[abs(faceID)] for faceID in current_boundaryFaceIDs]

        exitH_faceIDs[iExitH] = current_boundaryFaceIDs   #the face ID of the current exit-h boundary
        exitH_ghostCellIDs[iExitH] = current_ghostCellIDs   #the ghost cell ID of the current exit-h boundary
        exitH_internalCellIDs[iExitH] = current_internalCellIDs   #the internal cell ID of the current exit-h boundary    

        #number of bounary faces for the current exit-h boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)

        #exitH_faceCentroids[iExitH] = zeros(Float64, nBoundaryFaces, 3)   #face centroids of the current exit-h boundary
        exitH_faceOutwardNormals[iExitH] = zeros(Float64, nBoundaryFaces, 2)   #face outward normals of the current exit-h boundary

        #loop through all faces in the current exit-h boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]

            #get the face centroid of the current exit-h boundary
            #faceCentroid = (nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][1], :] + nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][2], :]) / 2.0
            #exitH_faceCentroids[iExitH][iFace, :] = faceCentroid

            #get the face direction of the current inlet-q boundary
            face_direction = current_boundaryFace_directions[iFace]

            #get the face outward normal of the current inlet-q boundary
            face_normal = face_normals[faceID]
            exitH_faceOutwardNormals[iExitH][iFace, :] = face_direction .* face_normal

            #println("Face ID: ", faceID, ", Ghost cell ID: ", ghostCellID, ", Internal cell ID: ", internalCellID, ", Face centroid: ", exitH_faceCentroids[iExitH][iFace, :], ", Face normal: ", exitH_faceOutwardNormals[iExitH][iFace, :])
        end

        #get the water surface elevation of the current exit-h boundary
        exitH_WSE[iExitH] = parse(Float64, srhhydro_EWSParamsC[iBoundary][1])

        #get the ghost cell ID of the current exit-h boundary
        exitH_H[iExitH] = zeros(Float64, nBoundaryFaces)   #water depth for the current exit-h boundary   
        exitH_A[iExitH] = zeros(Float64, nBoundaryFaces)   #cross-sectional area for the current exit-h boundary
    end

end

#preprocess wall boundaries
function preprocess_wall_boundaries_2D(settings, srh_all_Dict, nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs, wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals, wall_H, wall_A, nodeCoordinates)

    my_mesh_2D = srh_all_Dict["my_mesh_2D"]
 
    face_normals = my_mesh_2D.face_normals   #face normals
    boundaryFaces_direction_Dict = my_mesh_2D.boundaryFaces_direction_Dict   #face directions of boundary faces

    #loop through all wall boundaries
    for iWall in 1:nWall_BCs
        iBoundary = wall_BC_indices[iWall]

        if settings.bVerbose
            println("Preprocessing WALL boundary ", iWall, " with index in BC list ", iBoundary)
        end

        #get the boundary face IDs of the current wall boundary
        current_boundaryFaceIDs = my_mesh_2D.boundaryFaces_Dict[iBoundary]
        current_boundaryFace_directions = boundaryFaces_direction_Dict[iBoundary]
        current_ghostCellIDs = [my_mesh_2D.boundaryFaceID_to_ghostCellID_Dict[abs(faceID)] for faceID in current_boundaryFaceIDs]
        current_internalCellIDs = [my_mesh_2D.boundaryFaceID_to_internalCellID_Dict[abs(faceID)] for faceID in current_boundaryFaceIDs]

        wall_faceIDs[iWall] = current_boundaryFaceIDs   #the face ID of the current wall boundary
        wall_ghostCellIDs[iWall] = current_ghostCellIDs   #the ghost cell ID of the current wall boundary
        wall_internalCellIDs[iWall] = current_internalCellIDs   #the internal cell ID of the current wall boundary    

        #number of bounary faces for the current wall boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)

        wall_faceCentroids[iWall] = zeros(Float64, nBoundaryFaces, 3)   #face centroids of the current wall boundary
        wall_outwardNormals[iWall] = zeros(Float64, nBoundaryFaces, 2)   #face outward normals of the current wall boundary


        #loop through all faces in the current wall boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]

            #get the face centroid of the current wall boundary
            faceCentroid = (nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][1], :] + nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][2], :]) / 2.0
            wall_faceCentroids[iWall][iFace, :] = faceCentroid

            #get the face direction of the current inlet-q boundary
            face_direction = current_boundaryFace_directions[iFace]

            #get the face outward normal of the current inlet-q boundary
            face_normal = face_normals[faceID]
            wall_outwardNormals[iWall][iFace, :] = face_direction .* face_normal

            #println("Face ID: ", faceID, ", Ghost cell ID: ", ghostCellID, ", Internal cell ID: ", internalCellID, ", Face centroid: ", wall_faceCentroids[iWall][iFace, :], ", Face normal: ", wall_outwardNormals[iWall][iFace, :])
        end

        #get the ghost cell ID of the current wall boundary
        wall_H[iWall] = zeros(Float64, nBoundaryFaces)   # water depth for the current wall boundary   
        wall_A[iWall] = zeros(Float64, nBoundaryFaces)   # cross-sectional area for the current wall boundary  
    end

end

#preprocess symmetry boundaries: called every time step to update the boundary condition
function preprocess_symmetry_boundaries_2D(settings, srh_all_Dict, nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs, symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals, symm_H, symm_A, nodeCoordinates)

    my_mesh_2D = srh_all_Dict["my_mesh_2D"]

    face_normals = my_mesh_2D.face_normals   #face normals
    boundaryFaces_direction_Dict = my_mesh_2D.boundaryFaces_direction_Dict   #face directions of boundary faces

    #loop through all symmetry boundaries
    for iSymm in 1:nSymm_BCs
        iBoundary = symm_BC_indices[iSymm]

        if settings.bVerbose
            println("Preprocessing SYMMETRY boundary ", iSymm, " with index in BC list ", iBoundary)
        end

        #get the boundary face IDs of the current wall boundary
        current_boundaryFaceIDs = my_mesh_2D.boundaryFaces_Dict[iBoundary]
        current_boundaryFace_directions = boundaryFaces_direction_Dict[iBoundary]
        current_ghostCellIDs = [my_mesh_2D.boundaryFaceID_to_ghostCellID_Dict[abs(faceID)] for faceID in current_boundaryFaceIDs]
        current_internalCellIDs = [my_mesh_2D.boundaryFaceID_to_internalCellID_Dict[abs(faceID)] for faceID in current_boundaryFaceIDs]

        symm_faceIDs[iSymm] = current_boundaryFaceIDs   #the face ID of the current symmetry boundary
        symm_ghostCellIDs[iSymm] = current_ghostCellIDs   #the ghost cell ID of the current symmetry boundary
        symm_internalCellIDs[iSymm] = current_internalCellIDs   #the internal cell ID of the current symmetry boundary    

        #number of bounary faces for the current symmetry boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)

        symm_faceCentroids[iSymm] = zeros(Float64, nBoundaryFaces, 3)   #face centroids of the current symmetry boundary
        symm_outwardNormals[iSymm] = zeros(Float64, nBoundaryFaces, 2)   #face outward normals of the current symmetry boundary

        #loop through all faces in the current symmetry boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]

            #get the face centroid of the current symmetry boundary
            faceCentroid = (nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][1], :] + nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][2], :]) / 2.0
            symm_faceCentroids[iSymm][iFace, :] = faceCentroid

            #get the face direction of the current symmetry boundary
            face_direction = current_boundaryFace_directions[iFace]

            #get the face outward normal of the current symmetry boundary
            face_normal = face_normals[faceID]
            symm_outwardNormals[iSymm][iFace, :] = face_direction .* face_normal

            #println("Face ID: ", faceID, ", Ghost cell ID: ", ghostCellID, ", Internal cell ID: ", internalCellID, ", Face centroid: ", symm_faceCentroids[iSymm][iFace, :], ", Face normal: ", symm_outwardNormals[iSymm][iFace, :])
        end

        #get the ghost cell ID of the current symmetry boundary
        symm_H[iSymm] = zeros(Float64, nBoundaryFaces)   # water depth for the current symmetry boundary   
        symm_A[iSymm] = zeros(Float64, nBoundaryFaces)   # cross-sectional area for the current symmetry boundary  
    end

end

#process all boundaries in just one function call to update Q_ghostCells: called every time step to update the boundary condition
#return a new Q_ghostCells
@noinline function process_all_boundaries_2d(
    settings::ControlSettings, 
    #Q_cells::Matrix{T}, 
    h::Vector{T},
    q_x::Vector{T},
    q_y::Vector{T},
    my_mesh_2D::mesh_2D, 
    boundary_conditions::BoundaryConditions2D, 
    ManningN_cells::Vector{T}, 
    zb_faces::Vector{T}, 
    swe_2D_constants::swe_2D_consts, 
    inletQ_Length::Vector{Array{T}}, 
    inletQ_TotalQ::Vector{T}, 
    exitH_WSE::Vector{T}) where {T}

    Zygote.ignore() do
        if settings.bVerbose
            #println("process_all_boundaries_2d")
            #println(" ")
        end
    end

    Zygote.ignore() do
        println("\nStarting process_all_boundaries_2d")
        println(" ")
        println("Input types:")
        println("h: ", typeof(h))
        println("q_x: ", typeof(q_x))
        println("q_y: ", typeof(q_y))
        println("my_mesh_2D: ", typeof(my_mesh_2D))
        println("boundary_conditions.inletQ_faceOutwardNormals: ", typeof(boundary_conditions.inletQ_faceOutwardNormals))
        println("ManningN_cells: ", typeof(ManningN_cells))
        println("zb_faces: ", typeof(zb_faces))
        println("inletQ_TotalQ: ", typeof(inletQ_TotalQ))
        println("inletQ_Length: ", typeof(inletQ_Length))
        println("exitH_WSE: ", typeof(exitH_WSE))
    end

    Zygote.ignore() do
        println("\nBoundary Conditions Types:")
        println("boundary_conditions: ", typeof(boundary_conditions))
        println("inletQ_BC_indices: ", typeof(boundary_conditions.inletQ_BC_indices))
        println("inletQ_faceOutwardNormals: ", typeof(boundary_conditions.inletQ_faceOutwardNormals))
        println("inletQ_faceIDs: ", typeof(boundary_conditions.inletQ_faceIDs))
        println("inletQ_ghostCellIDs: ", typeof(boundary_conditions.inletQ_ghostCellIDs))
        println("inletQ_internalCellIDs: ", typeof(boundary_conditions.inletQ_internalCellIDs))
        println("nInletQ_BCs: ", typeof(boundary_conditions.nInletQ_BCs))
    end

    #get views of the Q_cells array
    #h = @view Q_cells[:, 1]         
    #q_x = @view Q_cells[:, 2]
    #q_y = @view Q_cells[:, 3]

    #make copy of the Q_cells array
    #h = copy(Q_cells[:, 1])         
    #q_x = copy(Q_cells[:, 2])
    #q_y = copy(Q_cells[:, 3])

    h_ghost_local = zeros(eltype(h), my_mesh_2D.numOfAllBounaryFaces)
    q_x_ghost_local = zeros(eltype(q_x), my_mesh_2D.numOfAllBounaryFaces)
    q_y_ghost_local = zeros(eltype(q_y), my_mesh_2D.numOfAllBounaryFaces)

    # Before inlet-Q boundaries
    Zygote.ignore() do
        println("\nBefore inlet-Q boundaries")
        println("h_ghost_local: ", typeof(h_ghost_local))
        println("q_x_ghost_local: ", typeof(q_x_ghost_local))
        println("q_y_ghost_local: ", typeof(q_y_ghost_local))
    end

    #for inlet-q boundaries
    # Loop through all inlet-q boundaries
    for iInletQ in 1:boundary_conditions.nInletQ_BCs
        
        iBoundary = boundary_conditions.inletQ_BC_indices[iInletQ]

        if settings.bVerbose
            Zygote.ignore() do
                #println("Processing INLET-Q boundary ", iInletQ, " with index in BC list ", iBoundary)
            end
        end

        current_boundaryFaceIDs = boundary_conditions.inletQ_faceIDs[iInletQ]
        current_ghostCellIDs = boundary_conditions.inletQ_ghostCellIDs[iInletQ]
        current_internalCellIDs = boundary_conditions.inletQ_internalCellIDs[iInletQ]

        current_inletQ_Length = inletQ_Length[iInletQ]
        current_inletQ_TotalQ = inletQ_TotalQ[iInletQ]

        # Create new arrays for this boundary's updates
        h_new = zeros(eltype(h), length(current_ghostCellIDs))
        q_x_new = zeros(eltype(q_x), length(current_ghostCellIDs))
        q_y_new = zeros(eltype(q_y), length(current_ghostCellIDs))

        # Before total_A calculation
        Zygote.ignore() do
            println("\nBefore total_A calculation:")
            println("current_inletQ_Length type: ", typeof(current_inletQ_Length))
            println("h type: ", typeof(h))
            println("ManningN_cells type: ", typeof(ManningN_cells))
        end


        # First pass: compute total_A and dry/wet status
        # current_inletQ_DryWet: 1 for wet, 0 for dry
        #current_inletQ_DryWet = map(internalCellID -> h[internalCellID] > swe_2D_constants.h_small ? 1 : 0, current_internalCellIDs)  #does this break AD?
        #current_inletQ_DryWet = map(internalCellID -> 1, current_internalCellIDs)
        current_inletQ_DryWet = Int.(h[current_internalCellIDs] .> swe_2D_constants.h_small)

        # Compute total_A only for wet faces
        #total_A = sum(current_inletQ_Length[iFace]^(5.0 / 3.0) * h[current_internalCellIDs[iFace]] / ManningN_cells[current_internalCellIDs[iFace]]
                      #for iFace in 1:length(current_internalCellIDs) if current_inletQ_DryWet[iFace] == 1
        #)
        total_A = 100.0
       
        # After total_A calculation
        Zygote.ignore() do
            println("\nAfter total_A calculation:")
            println(" ")
            println("type of current_inletQ_DryWet: ", typeof(current_inletQ_DryWet))
            println("total_A type: ", typeof(total_A))
        end

        Zygote.ignore() do
            @assert total_A > 1e-10 "Total cross-sectional conveyance for inlet-q boundary $iInletQ is not positive: $total_A"
        end

        # Second pass: compute ghost cell values
        h_new = map(internalCellID -> h[internalCellID], current_internalCellIDs)

        face_normals = boundary_conditions.inletQ_faceOutwardNormals[iInletQ]
        velocity_normals = @. current_inletQ_TotalQ / total_A * 
                     current_inletQ_Length^(2.0 / 3.0) / 
                     ManningN_cells[current_internalCellIDs]

        Zygote.ignore() do
            println("face_normals: ", typeof(face_normals))
            println("face_normals: ", face_normals)
            println("face_normals size: ", size(face_normals))
            println("face_normals length: ", length(face_normals))

            println("velocity_normals: ", typeof(velocity_normals))
            println("velocity_normals: ", velocity_normals)
            println("velocity_normals size: ", size(velocity_normals))
            println("velocity_normals length: ", length(velocity_normals))

            println(" ")
        end

        # Calculate q components using list comprehensions
        q_x_new = [-h[internalCellID] * velocity_normals[iFace] * face_normals[iFace, 1] * 
        current_inletQ_DryWet[iFace] 
        for (iFace, internalCellID) in enumerate(current_internalCellIDs)]

        q_y_new = [-h[internalCellID] * velocity_normals[iFace] * face_normals[iFace, 2] * 
        current_inletQ_DryWet[iFace]
        for (iFace, internalCellID) in enumerate(current_internalCellIDs)]

    #     q_updates = [(
    #        let
    #            ManningN_face = ManningN_cells[internalCellID]
    #            face_normal = boundary_conditions.inletQ_faceOutwardNormals[iInletQ][iFace, :]
    #            velocity_normal = (current_inletQ_TotalQ / total_A *
    #                             current_inletQ_Length[iFace]^(2.0 / 3.0) / ManningN_face)
    #             if current_inletQ_DryWet[iFace] == 0
    #                (eltype(Q_cells)(0.0), eltype(Q_cells)(0.0))
    #            else
    #                (-h[internalCellID] * velocity_normal * face_normal[1],
    #                 -h[internalCellID] * velocity_normal * face_normal[2])
    #            end
    #        end
    #    ) for (iFace, internalCellID) in enumerate(current_internalCellIDs)]

    #     q_x_new = [q_updates[iFace][1] for iFace in 1:length(current_internalCellIDs)]
    #     q_y_new = [q_updates[iFace][2] for iFace in 1:length(current_internalCellIDs)]

        Zygote.ignore() do
            #@show typeof(current_inletQ_TotalQ)
            #@show typeof(total_A)
            #@show typeof(current_inletQ_Length)
            #@show typeof(ManningN_cells)

            #@show current_inletQ_TotalQ
            #@show total_A
            #@show current_inletQ_Length
            #@show ManningN_cells
        end
      
        Zygote.ignore() do
            @show typeof(h_new)
            @show typeof(q_x_new)
            @show typeof(q_y_new)
            @show size(q_x_new)
            @show size(q_y_new)
            @show size(h_new)
            @show h_new
            @show q_x_new
            @show q_y_new
        end

        # Before update_1d_array calls
        Zygote.ignore() do
            println(" ")
            println("\nBefore update_1d_array:")
            println("h_new type: ", typeof(h_new))
            println("current_ghostCellIDs type: ", typeof(current_ghostCellIDs))
        end

        # Update ghost arrays using update_1d_array
        h_ghost_local = update_1d_array(h_ghost_local, current_ghostCellIDs, h_new)
        q_x_ghost_local = update_1d_array(q_x_ghost_local, current_ghostCellIDs, q_x_new)
        q_y_ghost_local = update_1d_array(q_y_ghost_local, current_ghostCellIDs, q_y_new)
        
        #fake it for now
        #q_x_ghost_local = zeros(eltype(Q_cells), my_mesh_2D.numOfAllBounaryFaces)
        #q_y_ghost_local = zeros(eltype(Q_cells), my_mesh_2D.numOfAllBounaryFaces)
    end

    # After inlet-Q boundaries
    Zygote.ignore() do
        println("\nAfter inlet-Q boundaries")
        println(" ")
        println("h_ghost_local: ", typeof(h_ghost_local))
        println("q_x_ghost_local: ", typeof(q_x_ghost_local))
        println("q_y_ghost_local: ", typeof(q_y_ghost_local))
    end

    # Before exit-H boundaries
    Zygote.ignore() do
        println("\nBefore exit-H boundaries")
        println(" ")
    end

    #for exit-h boundaries
    #loop through all exit-h boundaries
    for iExitH in 1:boundary_conditions.nExitH_BCs
        iBoundary = boundary_conditions.exitH_BC_indices[iExitH]

        if settings.bVerbose
            Zygote.ignore() do
                #println("Processing EXIT-H boundary ", iExitH, " with index in BC list ", iBoundary)
            end
        end

        #iBoundary = exitH_BC_indices[iExitH]

        current_ghostCellIDs = boundary_conditions.exitH_ghostCellIDs[iExitH]  # ghost cell IDs for this boundary
        current_internalCellIDs = boundary_conditions.exitH_internalCellIDs[iExitH]  # internal cell IDs for this boundary
        current_faceIDs = boundary_conditions.exitH_faceIDs[iExitH]  # face IDs for this boundary

        #current_faceCentroids = exitH_faceCentroids[iExitH]  # face centroids for this boundary
        
        # Calculate new h_ghost values
        h_new = convert.(eltype(h), max.(swe_2D_constants.h_small,
            exitH_WSE[iExitH] .- zb_faces[current_faceIDs]))

        # Update arrays using update_1d_array
        h_ghost_local = update_1d_array(h_ghost_local, current_ghostCellIDs, h_new)
        q_x_ghost_local = update_1d_array(q_x_ghost_local, current_ghostCellIDs, q_x[current_internalCellIDs])
        q_y_ghost_local = update_1d_array(q_y_ghost_local, current_ghostCellIDs, q_y[current_internalCellIDs])
    end

    # After exit-H boundaries
    Zygote.ignore() do
        println("\nAfter exit-H boundaries")
        println(" ")
        println("h_ghost_local: ", typeof(h_ghost_local))
        println("q_x_ghost_local: ", typeof(q_x_ghost_local))
        println("q_y_ghost_local: ", typeof(q_y_ghost_local))

        println(" ")    
        println("\nBefore wall boundaries")
    end

    #for wall boundaries
    #loop through all wall boundaries
    for iWall in 1:boundary_conditions.nWall_BCs
        iBoundary = boundary_conditions.wall_BC_indices[iWall]

        if settings.bVerbose
            Zygote.ignore() do
                #println("Processing WALL boundary ", iWall, " with index in BC list ", iBoundary)
            end
        end

        # Get indices
        current_ghostCellIDs = boundary_conditions.wall_ghostCellIDs[iWall]
        current_internalCellIDs = boundary_conditions.wall_internalCellIDs[iWall]

        # updates (without in-place mutation)
        h_ghost_local = update_1d_array(h_ghost_local, current_ghostCellIDs, h[current_internalCellIDs])
        q_x_ghost_local = update_1d_array(q_x_ghost_local, current_ghostCellIDs, -q_x[current_internalCellIDs])
        q_y_ghost_local = update_1d_array(q_y_ghost_local, current_ghostCellIDs, -q_y[current_internalCellIDs])
    end

    # After wall boundaries
    Zygote.ignore() do
        println("\nAfter wall boundaries")
        println(" ")
        println("h_ghost_local: ", typeof(h_ghost_local))
        println("q_x_ghost_local: ", typeof(q_x_ghost_local))
        println("q_y_ghost_local: ", typeof(q_y_ghost_local))
    end

    #for symmetry boundaries
    #loop through all symmetry boundaries
    for iSymm in 1:boundary_conditions.nSymm_BCs
        iBoundary = boundary_conditions.symm_BC_indices[iSymm]

        if settings.bVerbose
            Zygote.ignore() do
                #println("Processing SYMMETRY boundary ", iSymm, " with index in BC list ", iBoundary)
            end
        end

        current_ghostCellIDs = boundary_conditions.symm_ghostCellIDs[iSymm]
        current_internalCellIDs = boundary_conditions.symm_internalCellIDs[iSymm]
        face_normals = boundary_conditions.symm_outwardNormals[iSymm]

        # Compute v·n (dot product of velocity and face normal)
        v_dot_n = q_x[current_internalCellIDs] .* face_normals[:, 1] .+
                  q_y[current_internalCellIDs] .* face_normals[:, 2]

        # Update depths using update_1d_array
        h_ghost_local = update_1d_array(h_ghost_local, current_ghostCellIDs, h[current_internalCellIDs])

        # Update velocities using update_1d_array
        q_x_new = q_x[current_internalCellIDs] .- 2.0 .* v_dot_n .* face_normals[:, 1]
        q_y_new = q_y[current_internalCellIDs] .- 2.0 .* v_dot_n .* face_normals[:, 2]

        q_x_ghost_local = update_1d_array(q_x_ghost_local, current_ghostCellIDs, q_x_new)
        q_y_ghost_local = update_1d_array(q_y_ghost_local, current_ghostCellIDs, q_y_new)
    end

    # After symmetry boundaries
    Zygote.ignore() do
        println("\nAfter symmetry boundaries")
        println(" ")
        println("h_ghost_local: ", typeof(h_ghost_local))
        println("q_x_ghost_local: ", typeof(q_x_ghost_local))
        println("q_y_ghost_local: ", typeof(q_y_ghost_local))
    end

    # Stack the updated ghost cell values and return new Q_ghostCells
    new_Q_ghostCells = hcat(h_ghost_local, q_x_ghost_local, q_y_ghost_local)

    Zygote.ignore() do
        #@show size(new_Q_ghostCells)
        #@show typeof(new_Q_ghostCells)
        #@show new_Q_ghostCells
    end

    Zygote.ignore() do
        println("\nReturning new_Q_ghostCells")
        println(" ")
        println("new_Q_ghostCells: ", typeof(new_Q_ghostCells))
        println(" ")
    end

    return new_Q_ghostCells
end

#update inletQ_TotalQ
function update_inletQ_TotalQ(inlet_discharges_current)
    #update inletQ_TotalQ based on the provided inlet_discharges_current

    Zygote.ignore() do
       println("update_inletQ_TotalQ")
    end

    
    #@show gradient(x -> sum([x[i] for i in 1:length(x)]), inlet_discharges_current)
    

    #inletQ_TotalQ =  [inlet_discharges_current[iInletQ] for iInletQ in 1:(length(inlet_discharges_current))]
    inletQ_TotalQ = copy(inlet_discharges_current)

    return inletQ_TotalQ

end


