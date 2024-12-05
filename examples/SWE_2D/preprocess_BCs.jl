#preprcess all boundary conditions, such as calculate the list of faces, internal cells, ghost cells, etc.
#These functions are supposed to be called only once before the time loop

using LinearAlgebra

function compute_boundary_indices!(my_mesh_2D, srh_all_Dict, inletQ_BC_indices, exitH_BC_indices, wall_BC_indices, symm_BC_indices)

    srhhydro_BC = srh_all_Dict["srhhydro_BC"]   #boundary conditions
    srhhydro_IQParams = srh_all_Dict["srhhydro_IQParams"]   #inlet-Q parameters
    srhhydro_EWSParamsC = srh_all_Dict["srhhydro_EWSParamsC"]   #exit-H parameters
    
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
            
            push!(inletQ_BC_indices, iBoundary)
            
            #find the corresponding Q in IQParams
            if haskey(srhhydro_IQParams, iBoundary)
                println("Key IQParams exists in the Python dictionary srhhydro_IQParams.")
                println("Inlet specific discharge: ", srhhydro_IQParams[iBoundary])
                
                #update the boundary value: currenlty only support constant discharge
                Q_value = parse(Float64, srhhydro_IQParams[iBoundary][1])
                println("Inlet discharge: ", Q_value)
                
            else
                println("Key IQParams does not exist in the Python dictionary srhhydro_IQParams.")
                readline()
                exit(-1)
            end
            
        elseif lowercase(boundaryType) == "exit-h"
            println("EXIT-H boundary condition is set for boundary ", iBoundary)
            
            push!(exitH_BC_indices, iBoundary)
            
            #find the corresponding H in EWSParamsC
            if haskey(srhhydro_EWSParamsC, iBoundary)
                println("Key EWSParamsC exists in the Python dictionary srhhydro_EWSParamsC.")
                println("Exit water depth: ", srhhydro_EWSParamsC[iBoundary])
                
                #update the boundary value: currenlty only support constant water depth
                H_value = parse(Float64, srhhydro_EWSParamsC[iBoundary][1])
                println("Exit water depth: ", H_value)
                
            else
                println("Key EWSParamsC does not exist in the Python dictionary srhhydro_EWSParamsC.")
                readline()
                exit(-1)
            end
            
        elseif lowercase(boundaryType) == "wall"
            println("WALL boundary condition is set for boundary ", iBoundary)
            
            push!(wall_BC_indices, iBoundary)
            
        elseif lowercase(boundaryType) == "symm"
            println("SYMMETRY boundary condition is set for boundary ", iBoundary)
            
            push!(symm_BC_indices, iBoundary)
        end
        
    end
    
    #number of wall boundaries
    nWall_BCs = length(wall_BC_indices)
    srh_all_Dict["nWall_BCs"] = nWall_BCs  #update the number of wall boundaries in the dictionary
    
    #number of symmetry boundaries
    nSymm_BCs = length(symm_BC_indices)
    srh_all_Dict["nSymm_BCs"] = nSymm_BCs  #update the number of symmetry boundaries in the dictionary
    
    println("inletQ_BC_indices: ", inletQ_BC_indices)
    println("exitH_BC_indices: ", exitH_BC_indices) 
    println("wall_BC_indices: ", wall_BC_indices)
    println("symm_BC_indices: ", symm_BC_indices)
    
end

#preprocess inlet-q boundaries
function preprocess_inlet_q_boundaries(my_mesh_2D, nInletQ_BCs, srhhydro_IQParams, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, 
    inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length, inletQ_TotalA, inletQ_DryWet)
        
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
        
        inletQ_faceIDs[iInletQ] = current_boundaryFaceIDs   #the face ID of the current inlet-q boundary
        inletQ_ghostCellIDs[iInletQ] = current_ghostCellIDs   #the ghost cell ID of the current inlet-q boundary
        inletQ_internalCellIDs[iInletQ] = current_internalCellIDs   #the internal cell ID of the current inlet-q boundary    
        
        #number of bounary faces for the current inlet-q boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)
        
        inletQ_faceCentroids[iInletQ] = zeros(Float64, nBoundaryFaces, 3)   #face centroids of the current inlet-q boundary
        inletQ_faceOutwardNormals[iInletQ] = zeros(Float64, nBoundaryFaces, 2)   #face outward normals of the current inlet-q boundary
        inletQ_Length[iInletQ] = zeros(Float64, nBoundaryFaces)   #length for each face in the current inlet-q boundary
        
        #loop through all faces in the current inlet-q boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]
            
            #get the face centroid of the current inlet-q boundary
            faceCentroid = (my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[faceID][1],:] + my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[faceID][2],:]) / 2.0
            inletQ_faceCentroids[iInletQ][iFace, :] = faceCentroid

            #get the face direction of the current inlet-q boundary
            face_direction = current_boundaryFace_directions[iFace]

            #get the face outward normal of the current inlet-q boundary
            face_normal = face_normals[faceID]
            inletQ_faceOutwardNormals[iInletQ][iFace, :] = face_direction .* face_normal  

            #get the length of the current face
            inletQ_Length[iInletQ][iFace] = my_mesh_2D.face_lengths[faceID]
            
            println("Face ID: ", faceID, ", Ghost cell ID: ", ghostCellID, ", Internal cell ID: ", internalCellID, ", Face centroid: ", inletQ_faceCentroids[iInletQ][iFace, :], ", Face normal: ", inletQ_faceOutwardNormals[iInletQ][iFace, :], ", Face length: ", inletQ_Length[iInletQ][iFace])
        end
        
        #get the total discharge of the current inlet-q boundary
        inletQ_TotalQ[iInletQ] = parse(Float64, srhhydro_IQParams[iBoundary][1])
        
        #get the ghost cell ID of the current inlet-q boundary
        inletQ_H[iInletQ] = zeros(Float64, nBoundaryFaces)   #inlet water depth for the current inlet-q boundary
        inletQ_A[iInletQ] = zeros(Float64, nBoundaryFaces)   #inlet cross-sectional area for the current inlet-q boundary
        inletQ_ManningN[iInletQ] = zeros(Float64, nBoundaryFaces)   #Manning's n for the current inlet-q boundary
        
        inletQ_TotalA[iInletQ] = 0.0   #total cross-sectional area for the current inlet-q boundary
        inletQ_DryWet[iInletQ] = zeros(Int, nBoundaryFaces)   #dry(=0)/wet(=1) flag for the current inlet-q boundary
    end
    
end

#preprocess exit-h boundaries
function preprocess_exit_h_boundaries(my_mesh_2D, nExitH_BCs, srhhydro_EWSParamsC, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, 
    exitH_internalCellIDs, exitH_faceCentroids, exitH_faceOutwardNormals, exitH_WSE, exitH_H, exitH_A)

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
        
        exitH_faceIDs[iExitH] = current_boundaryFaceIDs   #the face ID of the current exit-h boundary
        exitH_ghostCellIDs[iExitH] = current_ghostCellIDs   #the ghost cell ID of the current exit-h boundary
        exitH_internalCellIDs[iExitH] = current_internalCellIDs   #the internal cell ID of the current exit-h boundary    
        
        #number of bounary faces for the current exit-h boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)
        
        exitH_faceCentroids[iExitH] = zeros(Float64, nBoundaryFaces, 3)   #face centroids of the current exit-h boundary
        exitH_faceOutwardNormals[iExitH] = zeros(Float64, nBoundaryFaces, 2)   #face outward normals of the current exit-h boundary
        
        #loop through all faces in the current exit-h boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]
            
            #get the face centroid of the current exit-h boundary
            faceCentroid = (my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][1],:] + my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][2],:]) / 2.0
            exitH_faceCentroids[iExitH][iFace, :] = faceCentroid

            #get the face direction of the current inlet-q boundary
            face_direction = current_boundaryFace_directions[iFace]

            #get the face outward normal of the current inlet-q boundary
            face_normal = face_normals[faceID]
            exitH_faceOutwardNormals[iExitH][iFace, :] = face_direction .* face_normal  
            
            println("Face ID: ", faceID, ", Ghost cell ID: ", ghostCellID, ", Internal cell ID: ", internalCellID, ", Face centroid: ", exitH_faceCentroids[iExitH][iFace, :], ", Face normal: ", exitH_faceOutwardNormals[iExitH][iFace, :])
        end
        
        #get the water surface elevation of the current exit-h boundary
        exitH_WSE[iExitH] = parse(Float64, srhhydro_EWSParamsC[iBoundary][1])
        
        #get the ghost cell ID of the current exit-h boundary
        exitH_H[iExitH] = zeros(Float64, nBoundaryFaces)   #water depth for the current exit-h boundary   
        exitH_A[iExitH] = zeros(Float64, nBoundaryFaces)   #cross-sectional area for the current exit-h boundary
    end
    
end

#preprocess wall boundaries
function preprocess_wall_boundaries(my_mesh_2D, nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs, 
    wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals, wall_H, wall_A)

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
            faceCentroid = (my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][1],:] + my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][2],:]) / 2.0
            wall_faceCentroids[iWall][iFace, :] = faceCentroid

            #get the face direction of the current inlet-q boundary
            face_direction = current_boundaryFace_directions[iFace]

            #get the face outward normal of the current inlet-q boundary
            face_normal = face_normals[faceID]
            wall_outwardNormals[iWall][iFace, :] = face_direction .* face_normal  
            
            println("Face ID: ", faceID, ", Ghost cell ID: ", ghostCellID, ", Internal cell ID: ", internalCellID, ", Face centroid: ", wall_faceCentroids[iWall][iFace, :], ", Face normal: ", wall_outwardNormals[iWall][iFace, :])
        end
        
        #get the ghost cell ID of the current wall boundary
        wall_H[iWall] = zeros(Float64, nBoundaryFaces)   # water depth for the current wall boundary   
        wall_A[iWall] = zeros(Float64, nBoundaryFaces)   # cross-sectional area for the current wall boundary  
    end
    
end

#preprocess symmetry boundaries: called every time step to update the boundary condition
function preprocess_symmetry_boundaries(my_mesh_2D, nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs, 
    symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals, symm_H, symm_A)

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
            faceCentroid = (my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][1],:] + my_mesh_2D.nodeCoordinates[my_mesh_2D.faceNodes_r_Dict[abs(faceID)][2],:]) / 2.0
            symm_faceCentroids[iSymm][iFace, :] = faceCentroid

            #get the face direction of the current symmetry boundary
            face_direction = current_boundaryFace_directions[iFace]

            #get the face outward normal of the current symmetry boundary
            face_normal = face_normals[faceID]
            symm_outwardNormals[iSymm][iFace, :] = face_direction .* face_normal  
            
            println("Face ID: ", faceID, ", Ghost cell ID: ", ghostCellID, ", Internal cell ID: ", internalCellID, ", Face centroid: ", symm_faceCentroids[iSymm][iFace, :], ", Face normal: ", symm_outwardNormals[iSymm][iFace, :])
        end
        
        #get the ghost cell ID of the current symmetry boundary
        symm_H[iSymm] = zeros(Float64, nBoundaryFaces)   # water depth for the current symmetry boundary   
        symm_A[iSymm] = zeros(Float64, nBoundaryFaces)   # cross-sectional area for the current symmetry boundary  
    end

end

#process inlet-q boundaries: called every time step to update the boundary condition
function process_inlet_q_boundaries(nInletQ_BCs, inletQ_BC_indices, inletQ_faceIDs, inletQ_ghostCellIDs, 
    inletQ_internalCellIDs, inletQ_faceCentroids, inletQ_faceOutwardNormals, inletQ_TotalQ, inletQ_H, inletQ_A, inletQ_ManningN, inletQ_Length,
    inletQ_TotalA, inletQ_DryWet, Q_cells, Q_ghostCells, swe_2D_constants)

    h = @view Q_cells[:,1]
    q_x = @view Q_cells[:,2]
    q_y = @view Q_cells[:,3]

    h_ghost = @view Q_ghostCells[:,1]
    q_x_ghost = @view Q_ghostCells[:,2]
    q_y_ghost = @view Q_ghostCells[:,3]
        
    #loop through all inlet-q boundaries
    for iInletQ in 1:nInletQ_BCs
        iBoundary = inletQ_BC_indices[iInletQ]
        
        #println("Processing INLET-Q boundary ", iInletQ, " with index in BC list ", iBoundary)
        
        #get the boundary face IDs of the current inlet-q boundary
        
        current_boundaryFaceIDs = inletQ_faceIDs[iInletQ]          #the face ID of the current inlet-q boundary
        current_ghostCellIDs = inletQ_ghostCellIDs[iInletQ]        #the ghost cell ID of the current inlet-q boundary
        current_internalCellIDs = inletQ_internalCellIDs[iInletQ]  #the internal cell ID of the current inlet-q boundary    

        current_inletQ_H = inletQ_H[iInletQ]   #inlet water depth for each face in the current inlet-q boundary
        current_inletQ_A = inletQ_A[iInletQ]   #area for each face in the current inlet-q boundary
        current_inletQ_ManningN = inletQ_ManningN[iInletQ]   #Manning's n for each face in the current inlet-q boundary
        current_inletQ_Length = inletQ_Length[iInletQ]   #length for each face in the current inlet-q boundary
        current_inletQ_DryWet = inletQ_DryWet[iInletQ]   #dry(=0)/wet(=1) flag for each face in the current inlet-q boundary
        
        #number of bounary faces for the current inlet-q boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)
        
        current_inletQ_faceCentroids = inletQ_faceCentroids[iInletQ]  #face centroids of the current inlet-q boundary

        #get current total discharge for the current inlet-q boundary
        current_inletQ_TotalQ = inletQ_TotalQ[iInletQ] 

        inletQ_TotalA[iInletQ] = 0.0   #total cross-sectional area for the current inlet-q boundary
        
        #loop through all faces in the current inlet-q boundary and compute the total area (only account for the wet cells)
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]
            
            #get the face centroid of the current inlet-q boundary
            faceCentroid = current_inletQ_faceCentroids[iFace, :]

            current_inletQ_H[iFace] = h[internalCellID]   #get the water depth of the internal Cell
            current_inletQ_ManningN[iFace] = ManningN_cells[internalCellID]   #Manning's n for the current inlet-q boundary

            #get the water depth of the internal Cell
            if h[internalCellID] > swe_2D_constants.h_small
                current_inletQ_DryWet[iFace] = 1   #wet cell
                current_inletQ_A[iFace] = current_inletQ_Length[iFace]^(5/3) * h[internalCellID] / current_inletQ_ManningN[iFace]    #using conveyance method
                inletQ_TotalA[iInletQ] += current_inletQ_A[iFace]   #accumulate the total cross-sectional wet area
            else
                current_inletQ_DryWet[iFace] = 0   #dry cell
            end
        end

        if inletQ_TotalA[iInletQ] <= 1e-10
            println("Error: total cross-sectional conveyance for the current inlet-q boundary is not positive: ", iInletQ)
            readline()
            exit(-1)
        end
        
        #loop through all faces in the current inlet-q boundary and update the boundary condition
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]

            h_ghost[ghostCellID] = h[internalCellID]   #update the ghost cell water depth. Use the internal cell value no matter dry or wet

            if current_inletQ_DryWet[iFace] == 0   #dry cell
                q_x_ghost[ghostCellID] = 0.0   #update the ghost cell discharge
                q_y_ghost[ghostCellID] = 0.0   #update the ghost cell discharge
            
            else   #wet cell
                #face outward normal vector 
                face_normal = inletQ_faceOutwardNormals[iInletQ][iFace, :]
                
                #normal velocity
                velocity_normal = current_inletQ_TotalQ / inletQ_TotalA[iInletQ] * current_inletQ_Length[iFace]^(2/3) / current_inletQ_ManningN[iFace] 
               
                q_x_ghost[ghostCellID] = - h_ghost[ghostCellID] * velocity_normal * face_normal[1]   #update the ghost cell discharge
                q_y_ghost[ghostCellID] = - h_ghost[ghostCellID] * velocity_normal * face_normal[2]   #update the ghost cell discharge
            end
        end

    end
    
end

#process exit-h boundaries: called every time step to update the boundary condition
function process_exit_h_boundaries(nExitH_BCs, exitH_BC_indices, exitH_faceIDs, exitH_ghostCellIDs, 
    exitH_internalCellIDs, exitH_faceCentroids, exitH_WSE, Q_cells, Q_ghostCells, swe_2D_constants)

    h = @view Q_cells[:,1]
    q_x = @view Q_cells[:,2]
    q_y = @view Q_cells[:,3]

    h_ghost = @view Q_ghostCells[:,1]
    q_x_ghost = @view Q_ghostCells[:,2]
    q_y_ghost = @view Q_ghostCells[:,3]
    
    #loop through all exit-h boundaries
    for iExitH in 1:nExitH_BCs
        iBoundary = exitH_BC_indices[iExitH]
        
        #println("Processing EXIT-H boundary ", iExitH, " with index in BC list ", iBoundary)
        
        #get the boundary face IDs of the current exit-h boundary
        current_boundaryFaceIDs = exitH_faceIDs[iExitH]         #the face ID of the current exit-h boundary
        current_ghostCellIDs = exitH_ghostCellIDs[iExitH]       #the ghost cell ID of the current exit-h boundary
        current_internalCellIDs = exitH_internalCellIDs[iExitH] #the internal cell ID of the current exit-h boundary    
        
        #number of bounary faces for the current exit-h boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)
        
        #loop through all faces in the current exit-h boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]
            
            #get the face centroid of the current exit-h boundary
            faceCentroid = exitH_faceCentroids[iExitH][iFace, :]

            #update the ghost cell water depth. Use the exit water surface elevation
            h_ghost[ghostCellID] = max(swe_2D_constants.h_small, exitH_WSE[iExitH] - faceCentroid[3])  

            #update the ghost cell discharge                
            q_x_ghost[ghostCellID] = q_x[internalCellID]   #update the ghost cell discharge
            q_y_ghost[ghostCellID] = q_y[internalCellID]   #update the ghost cell discharge
        end

    end
    
end

#process wall boundaries: called every time step to update the boundary condition
function process_wall_boundaries(nWall_BCs, wall_BC_indices, wall_faceIDs, wall_ghostCellIDs, 
    wall_internalCellIDs, wall_faceCentroids, wall_outwardNormals, Q_cells, Q_ghostCells)

    h = @view Q_cells[:,1]
    q_x = @view Q_cells[:,2]
    q_y = @view Q_cells[:,3]

    h_ghost = @view Q_ghostCells[:,1]
    q_x_ghost = @view Q_ghostCells[:,2]
    q_y_ghost = @view Q_ghostCells[:,3]
    
    #loop through all wall boundaries
    for iWall in 1:nWall_BCs
        iBoundary = wall_BC_indices[iWall]
        
        #println("Prrocessing WALL boundary ", iWall, " with index in BC list ", iBoundary)
        
        #get the boundary face IDs of the current wall boundary
        current_boundaryFaceIDs = wall_faceIDs[iWall]       #the face ID of the current wall boundary
        current_ghostCellIDs = wall_ghostCellIDs[iWall]     #the ghost cell ID of the current wall boundary
        current_internalCellIDs = wall_internalCellIDs[iWall] #the internal cell ID of the current wall boundary    
        
        #number of bounary faces for the current wall boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)
        
        #loop through all faces in the current wall boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]
            
            #get the face centroid of the current wall boundary
            faceCentroid = wall_faceCentroids[iWall][iFace, :]

            #get the face outward normal of the current inlet-q boundary
            face_normal = wall_outwardNormals[iWall][iFace, :]
            
            #update the ghost cell water depth. Use the internal cell value
            h_ghost[ghostCellID] = h[internalCellID]   

            #update the ghost cell discharge                
            q_x_ghost[ghostCellID] = -q_x[internalCellID]   #update the ghost cell discharge
            q_y_ghost[ghostCellID] = -q_y[internalCellID]   #update the ghost cell discharge            
        end
       
    end
    
end

#process symmetry boundaries: called every time step to update the boundary condition
function process_symmetry_boundaries(nSymm_BCs, symm_BC_indices, symm_faceIDs, symm_ghostCellIDs, 
    symm_internalCellIDs, symm_faceCentroids, symm_outwardNormals, Q_cells, Q_ghostCells)

    h = @view Q_cells[:,1]
    q_x = @view Q_cells[:,2]
    q_y = @view Q_cells[:,3]

    h_ghost = @view Q_ghostCells[:,1]
    q_x_ghost = @view Q_ghostCells[:,2]
    q_y_ghost = @view Q_ghostCells[:,3]
    
    #loop through all symmetry boundaries
    for iSymm in 1:nSymm_BCs
        iBoundary = symm_BC_indices[iSymm]
        
        #println("Prrocessing SYMMETRY boundary ", iSymm, " with index in BC list ", iBoundary)
        
        #get the boundary face IDs of the current symmetry boundary
        current_boundaryFaceIDs = symm_faceIDs[iSymm]       #the face ID of the current symmetry boundary
        current_ghostCellIDs = symm_ghostCellIDs[iSymm]     #the ghost cell ID of the current symmetry boundary
        current_internalCellIDs = symm_internalCellIDs[iSymm] #the internal cell ID of the current symmetry boundary    
        
        #number of bounary faces for the current symmetry boundary
        nBoundaryFaces = length(current_boundaryFaceIDs)
        
        #loop through all faces in the current symmetry boundary
        for iFace in 1:nBoundaryFaces
            faceID = current_boundaryFaceIDs[iFace]
            ghostCellID = current_ghostCellIDs[iFace]
            internalCellID = current_internalCellIDs[iFace]
            
            #get the face centroid of the current symmetry boundary
            faceCentroid = symm_faceCentroids[iSymm][iFace, :]

            #get the face outward normal of the current symmetry boundary
            face_normal = symm_outwardNormals[iSymm][iFace, :]
            
            #update the ghost cell water depth. Use the internal cell value
            h_ghost[ghostCellID] = h[internalCellID]   

            #update the ghost cell discharge: symmetry about the face
            #Get internal cell velocity vector
            v_internal = [q_x[internalCellID], q_y[internalCellID]]
            
            #Compute reflection of velocity vector about face normal
            #v_reflected = v - 2(v·n)n where n is unit normal
            v_dot_n = dot(v_internal, face_normal)
            v_reflected = v_internal - 2 * v_dot_n * face_normal
            
            #Update ghost cell velocity components
            q_x_ghost[ghostCellID] = v_reflected[1]
            q_y_ghost[ghostCellID] = v_reflected[2]

        end
       
    end
    
end
