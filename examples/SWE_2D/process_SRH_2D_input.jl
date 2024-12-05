#process input files for SRH-2D case, e.g.. srhhydro, srhgeom, srhmat, srhmon files

function process_SRH_2D_input!(srhhydro_file_name::String)
    println("Python version: ", PyCall.pyversion)
    println("Python path: ", PyCall.pyprogramname)
    
    pyHMT2D = pyimport("pyHMT2D")
    pyHMT2D.gVerbose = true
    
    # Call a function from the custom package
    #file_path = joinpath(@__DIR__, "backwater.srhhydro" ) 
    #file_path = joinpath(@__DIR__, "simple.srhhydro" ) 
    file_path = joinpath(@__DIR__, srhhydro_file_name )
    my_srh_2d_data = pyHMT2D.SRH_2D.SRH_2D_Data(file_path)
    
    #ManningN_cell = my_srh_2d_data.ManningN_cell
    #ManningN_node = my_srh_2d_data.ManningN_node
    
    #println("ManningN_cell: ", ManningN_cell)
    #println("ManningN_node: ", ManningN_node)
    
    println("Result from Python function: ", my_srh_2d_data)
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
    
    srhhydro_Grid_fileName = srhhydro_obj.srhhydro_content["Grid"]
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
    
    #initialize 2D mesh
    my_mesh_2D = initialize_mesh_2D(srhgeom_obj, srhhydro_BC)
    
    # Access srhmat file data
    srhmat_obj = my_srh_2d_data.srhmat_obj
    srhmat_numOfMaterials = srhmat_obj.numOfMaterials
    srhmat_matNameList = srhmat_obj.matNameList
    srhmat_matZoneCells = srhmat_obj.matZoneCells

    #build material ID of each cell 
    matID_cells = zeros(Int, my_mesh_2D.numOfCells)

    for iCell in 1:my_mesh_2D.numOfCells
        bFound = false
        for (k, v) in srhmat_matZoneCells
            if iCell in v
                matID_cells[iCell] = k
                bFound = true
            end
        end

        if !bFound
            println("Cell ", iCell, " is not found in any material zone. Assign it to the default material zone 0.")
            matID_cells[iCell] = 0
        end
    end
    
    println("numOfMaterials: ", srhmat_numOfMaterials)
    println("matNameList: ", srhmat_matNameList)
    println("matZoneCells: ", srhmat_matZoneCells)
    println("matID_cells: ", matID_cells)

    #assemble all data into a dictionary and return
    srh_all_Dict = Dict()
    srh_all_Dict["srhhydro_obj"] = srhhydro_obj
    srh_all_Dict["srhhydro_ManningsN"] = srhhydro_ManningsN
    srh_all_Dict["srhhydro_BC"] = srhhydro_BC
    srh_all_Dict["srhhydro_MONITORINGLines"] = srhhydro_MONITORINGLines
    srh_all_Dict["srhhydro_Grid_fileName"] = srhhydro_Grid_fileName
    srh_all_Dict["srhhydro_RunType"] = srhhydro_RunType
    srh_all_Dict["srhhydro_TurbulenceModel"] = srhhydro_TurbulenceModel
    srh_all_Dict["srhhydro_OutputFormat"] = srhhydro_OutputFormat
    srh_all_Dict["srhhydro_OutputInterval"] = srhhydro_OutputInterval
    srh_all_Dict["srhhydro_Case"] = srhhydro_Case
    srh_all_Dict["srhhydro_ParabolicTurbulence"] = srhhydro_ParabolicTurbulence
    srh_all_Dict["srhhydro_Description"] = srhhydro_Description
    srh_all_Dict["srhhydro_OutputOption"] = srhhydro_OutputOption
    srh_all_Dict["srhhydro_HydroMat"] = srhhydro_HydroMat
    srh_all_Dict["srhhydro_InitCondOption"] = srhhydro_InitCondOption
    srh_all_Dict["srhhydro_SimTime"] = srhhydro_SimTime
    srh_all_Dict["srhhydro_SRHHYDRO"] = srhhydro_SRHHYDRO
    srh_all_Dict["srhhydro_ModelTemp"] = srhhydro_ModelTemp
    srh_all_Dict["srhhydro_UnsteadyOutput"] = srhhydro_UnsteadyOutput
    srh_all_Dict["srhhydro_IQParams"] = srhhydro_IQParams
    srh_all_Dict["srhhydro_EWSParamsC"] = srhhydro_EWSParamsC
    srh_all_Dict["nInletQ_BCs"] = nInletQ_BCs
    srh_all_Dict["nExitH_BCs"] = nExitH_BCs
    srh_all_Dict["nWall_BCs"] = 0  #nWall_BCs not defined yet. It is computed outside later. 
    srh_all_Dict["nSymm_BCs"] = 0  #nSymm_BCs not defined yet. It is computed outside later. 

    srh_all_Dict["srhgeom_obj"] = srhgeom_obj
    srh_all_Dict["srhmat_obj"] = srhmat_obj
    srh_all_Dict["my_mesh_2D"] = my_mesh_2D

    srh_all_Dict["srhmat_numOfMaterials"] = srhmat_numOfMaterials
    srh_all_Dict["srhmat_matNameList"] = srhmat_matNameList
    srh_all_Dict["srhmat_matZoneCells"] = srhmat_matZoneCells
    srh_all_Dict["matID_cells"] = matID_cells

    return srh_all_Dict  
end

#update swe_2d_constants based on the SRH-2D data
function update_swe_2D_constants!(swe_2D_constants, srh_all_Dict)

    srhhydro_SimTime = srh_all_Dict["srhhydro_SimTime"]  
    tStart = srhhydro_SimTime[1]
    tEnd = srhhydro_SimTime[2] * 3600 #convert hours to seconds
    dt = srhhydro_SimTime[3]

    #update the constants
    swe_2D_constants.tStart = tStart
    swe_2D_constants.tEnd = tEnd
    swe_2D_constants.dt = dt
    swe_2D_constants.t = tStart   #set the current time as the starting time 
end