#special postprocessing for the oneD_channel_with_bump case 


using JLD2 
using CSV
using DataFrames
using JSON3

# Postprocessing function for the checkpoint results
function postprocess_checkpoint_results(checkpoint_file)

    # Load checkpoint data
    checkpoint_results = load(checkpoint_file)
    
    ITER = checkpoint_results["ITER"]
    LOSS = checkpoint_results["LOSS"]
    PRED = checkpoint_results["PRED"]
    PARS = checkpoint_results["PARS"]
    
    #extract the Losses
    loss_total = []
    loss_pred = []
    loss_pred_WSE = []
    loss_pred_uv = []
    loss_bound = []
    loss_slope = []

    for curLoss in LOSS
        append!(loss_total, curLoss[1])
        append!(loss_pred, curLoss[2])
        append!(loss_pred_WSE, curLoss[3])
        append!(loss_pred_uv, curLoss[4])
        append!(loss_bound, curLoss[5])
        append!(loss_slope, curLoss[6])
    end 

    loss_history_df = DataFrame(iter_numbers = ITER, 
                                loss_total = loss_total, 
                                loss_pred = loss_pred, 
                                loss_pred_WSE = loss_pred_WSE, 
                                loss_pred_uv = loss_pred_uv, 
                                loss_bound = loss_bound, 
                                loss_slope = loss_slope)

    #save loss history to a file
    case_path = dirname(checkpoint_file)
    CSV.write(joinpath(case_path, "loss_history_checkpoint.csv"), loss_history_df)

    #save parameters to a file
    param_name = "zb"
    num_params = length(PARS[1])  # Get number of parameters from first vector
    param_names = ["$(param_name)_$i" for i in 1:num_params]

    pars_df = hcat(DataFrame(iter_numbers = ITER), DataFrame(hcat(PARS...)', param_names))

    # Write to CSV
    CSV.write(joinpath(case_path, "parameters_history_checkpoint.csv"), pars_df)

    #save the parameters at each checkpoint to file as a json file
    for i in eachindex(ITER)
        param_values = PARS[i]
        param_dict = Dict("zb_cells_param" => param_values)
        open(joinpath(case_path, "parameters_checkpoint_$(ITER[i]).json"), "w") do io
            JSON3.pretty(io, param_dict)
        end
    end
    
end

# Postprocessing function for the checkpoint results
checkpoint_file = "checkpoint_inversion_iter_350.jld2"
postprocess_checkpoint_results(checkpoint_file)

println("Checkpoint results postprocessed successfully.")