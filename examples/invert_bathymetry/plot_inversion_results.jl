#functions to make plots for the inversion 

using JLD2
using Plots
using LaTeXStrings
using Printf

plot_font = "sans-serif"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false,
        titlefont=font(16, plot_font),   #for title
        tickfont=font(12, plot_font),    #for tick label 
        guidefont=font(16, plot_font),   #for axis title label
        annotationfontsize=14,           #for annotate
        legendfont=font(10, plot_font))  #for legend


function  plot_inversion_results(save_path, mesh, zb_cell_truth)

    inversion_results = load(joinpath(save_path, "inversion_results.jld2"))
    
    LOSS = inversion_results["LOSS"]
    PRED = inversion_results["PRED"]
    PARS = inversion_results["PARS"]
    sol = inversion_results["sol"]

    plot_invesion_loss_time_history(LOSS)

    animate_inversion_process(LOSS, PRED, PARS, sol, mesh, zb_cell_truth)
    
end

function plot_invesion_loss_time_history(LOSS)

    iter_numbers = collect(range(1,size(LOSS)[1]))

    loss_total = []
    loss_pred = []
    loss_slope = []

    for curLoss in LOSS
        append!(loss_total, curLoss[1])
        append!(loss_pred, curLoss[2])
        append!(loss_slope, curLoss[3])
    end 

    #get rid of zeros in loss 
    @. loss_total = max(1e-12, loss_total)
    @. loss_pred = max(1e-12, loss_pred)
    @. loss_slope = max(1e-12, loss_slope)

    p1 = plot(iter_numbers, loss_total, c=:black, dpi=300, yscale=:log10, lablel="loss_total")
    plot!(iter_numbers, loss_pred, c=:blue, lablel="loss_pred")
    plot!(iter_numbers, loss_slope, c=:green, lablel="loss_slope")
                
    #xlims!(0, 30)
    #ylims!(-0.01, 2.01)
    xlabel!("Inversion iteration")
    ylabel!("Inversion loss")
    plot!(legend=:topright, fg_legend=:transparent, bg_legend=:transparent)

    savefig(p1, joinpath(save_path, "inversion_loss_history.png"))
end 

function animate_inversion_process(LOSS, PRED, PARS, sol, mesh, zb_cell_truth)
    anim = Animation()

    let animCounter = 1

        for (curLoss, curPred, curPars) in zip(LOSS, PRED, PARS)

            animCounter = animCounter + 1

            if (animCounter % 10) == 0

                #curPars is zb_cell
                if !isa(curPars, Vector{Float64})  #NLopt returns an optimization object, not an arrary
                    curPars = curPars.u
                end
                #interploate_zb_from_face_to_cell_and_compute_S0!(mesh, curPars, zb_cell, S0)
                zb_cell = curPars

                plt = plot(mesh.xCells, curPred[:, 1, end] + zb_cell, c=:aqua, label="free surface (predicted)", dpi=300)
                scatter!(mesh.xCells, Array(sol)[:,1,end]+zb_cell_truth, mc=:blue, ms=2, ma=0.5, label="free surface (truth)")
                plot!(mesh.xCells, zb_cell, c=:red, label="bed (inverted)")
                plot!(mesh.xCells, zb_cell_truth, c=:black, label="bed (truth)")

                xlims!(-0.1, 10.1)
                ylims!(-0.2, 2.01)
                xlabel!("x (m)")
                ylabel!("z (m)")

                #plot!(legend=:bottomleft, fg_legend=:transparent, bg_legend=:transparent)
                plot!(legend=(0.1, 0.4), fg_legend=:transparent, bg_legend=:transparent)

                annotate!(5.0, 1.9,
                          text(
                            "Inversion iteration = $(animCounter)"
                            )
                         )

                lossString = @sprintf("Loss = %.3e", curLoss[1])
                annotate!(5.0, 1.7,
                          text(
                            #"Loss = $(round(curLoss, digits = 2))"
                            #"Loss = $(curLoss)"
                            lossString
                            )
                         )

                frame(anim, plt)

            end

        end
    end

    #gif(anim, joinpath(save_path, "free_surface_1d.gif"), fps=15)
    gif(anim, joinpath(save_path, "inversion_process.mp4"), fps=1)
end