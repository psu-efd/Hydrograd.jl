

"""
A struct to handle srhhydro file for SRH-2D
"""
mutable struct SRH_2D_SRHHydro
    srhhydro_filename::String
    srhhydro_content::Dict{String, Any}

    function SRH_2D_SRHHydro(settings, srhhydro_filename::String)
        # Constructor
        if !isfile(srhhydro_filename)
            throw(ErrorException("The SRHHYDRO file $srhhydro_filename does not exist. Exiting..."))
        end
        
        obj = new(srhhydro_filename, Dict{String, Any}())
        parse_srhhydro_file(settings, obj)  # Call parse immediately after creation
        return obj
    end
end

"""
Parse the SRHHydro file and build the dictionary self.srhhydro_content
"""
function parse_srhhydro_file(settings, hydro::SRH_2D_SRHHydro)
    res_all = Dict{String, Any}()
    res_ManningsN = Dict{Int, Float64}()
    res_BC = Dict{Int, String}()
    res_MONITORING = Dict{Int, String}()
    res_IQParams = Dict{Int, Vector{String}}()
    res_ISupCrParams = Dict{Int, String}()
    res_EWSParamsC = Dict{Int, Vector{String}}()
    res_EWSParamsRC = Dict{Int, Vector{String}}()
    res_EQParams = Dict{Int, String}()
    res_NDParams = Dict{Int, String}()

    # Read file line by line
    for line in eachline(hydro.srhhydro_filename)
        # Use Julia's built-in parsing for shell-like syntax
        parts = split(strip(line))
        
        if length(parts) ≤ 1
            continue
        end

        if parts[1] == "ManningsN"
            res_ManningsN[parse(Int, parts[2])] = parse(Float64, parts[3])
        elseif parts[1] == "BC"
            if parts[3] == "MONITORING"
                res_MONITORING[parse(Int, parts[2])] = parts[3]
            else
                res_BC[parse(Int, parts[2])] = parts[3]
            end
        elseif parts[1] == "IQParams"
            res_IQParams[parse(Int, parts[2])] = parts[3:5]
        elseif parts[1] == "ISupCrParams"
            res_ISupCrParams[parse(Int, parts[2])] = parts[3]
        elseif parts[1] == "EWSParamsC"
            res_EWSParamsC[parse(Int, parts[2])] = parts[3:5]
        elseif parts[1] == "EWSParamsRC"
            res_EWSParamsRC[parse(Int, parts[2])] = parts[3:5]
        elseif parts[1] == "EQParams"
            res_EQParams[parse(Int, parts[2])] = parts[3]
        elseif parts[1] == "NDParams"
            res_NDParams[parse(Int, parts[2])] = parts[3]
        elseif parts[1] == "OutputFormat"
            res_all[parts[1]] = parts[2:3]
        elseif parts[1] == "SimTime"
            res_all[parts[1]] = [parse(Float64, x) for x in parts[2:4]]
        elseif parts[1] == "ParabolicTurbulence"
            res_all[parts[1]] = parse(Float64, parts[2])
        elseif parts[1] == "OutputOption"
            res_all[parts[1]] = parse(Int, parts[2])
        elseif parts[1] == "OutputInterval"
            res_all[parts[1]] = parse(Float64, parts[2])
        else
            res_all[strip(parts[1])] = parts[2]
        end
    end

    # Add sub-dictionaries to main dictionary
    !isempty(res_ManningsN) && (res_all["ManningsN"] = res_ManningsN)
    !isempty(res_BC) && (res_all["BC"] = res_BC)
    !isempty(res_MONITORING) && (res_all["MONITORING"] = res_MONITORING)
    !isempty(res_IQParams) && (res_all["IQParams"] = res_IQParams)
    !isempty(res_ISupCrParams) && (res_all["ISupCrParams"] = res_ISupCrParams)
    !isempty(res_EWSParamsC) && (res_all["EWSParamsC"] = res_EWSParamsC)
    !isempty(res_EWSParamsRC) && (res_all["EWSParamsRC"] = res_EWSParamsRC)
    !isempty(res_EQParams) && (res_all["EQParams"] = res_EQParams)
    !isempty(res_NDParams) && (res_all["NDParams"] = res_NDParams)

    hydro.srhhydro_content = res_all
end

"""
Get the dictionary for Manning's n
"""
function srhhydro_get_manningN_dict(hydro::SRH_2D_SRHHydro)
    return hydro.srhhydro_content["ManningsN"]
end


"""
Get boundary condition dictionary
"""
function srhhydro_get_bc_dict(hydro::SRH_2D_SRHHydro)
    return hydro.srhhydro_content["BC"]
end

"""
Get InletQ dictionary
"""
function srhhydro_get_inletQ_dict(hydro::SRH_2D_SRHHydro)
    return get(hydro.srhhydro_content, "IQParams", Dict{Int, Vector{String}}())
end


"""
Get ExitH dictionary
"""
function srhhydro_get_exitH_dict(hydro::SRH_2D_SRHHydro)
    return get(hydro.srhhydro_content, "EWSParamsC", Dict{Int, Vector{String}}())
end


"""
Get case name
"""
function srhhydro_get_case_name(hydro::SRH_2D_SRHHydro)
    return hydro.srhhydro_content["Case"]
end

"""
Get HydroMat file name
"""
function srhhydro_get_hydroMat_file_name(hydro::SRH_2D_SRHHydro)
    return hydro.srhhydro_content["HydroMat"]
end

"""
Get the start and end time of simulation (in hours)
"""
function srhhydro_get_simulation_start_end_time(hydro::SRH_2D_SRHHydro)
    value = hydro.srhhydro_content["SimTime"]
    return value[1], value[3]  # start time and end time
end

"""
Get the time step size of simulation (in seconds)
"""
function srhhydro_get_simulation_time_step_size(hydro::SRH_2D_SRHHydro)
    return hydro.srhhydro_content["SimTime"][2]
end

"""
Get grid file name (should be like "xxx.srhgeom")
"""
function srhhydro_get_grid_file_name(hydro::SRH_2D_SRHHydro)
    return hydro.srhhydro_content["Grid"]
end

"""
Get material file name (should be like "xxx.srhmat")
"""
function srhhydro_get_mat_file_name(hydro::SRH_2D_SRHHydro)
    return hydro.srhhydro_content["HydroMat"]
end