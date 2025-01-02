#functions to make plots for the inversion 

using JLD2
using Plots
using LaTeXStrings
using Printf

plot_font = "Times New Roman"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false,
        titlefont=font(16, plot_font),   #for title
        tickfont=font(12, plot_font),    #for tick label 
        guidefont=font(16, plot_font),   #for axis title label
        annotationfontsize=14,           #for annotate
        legendfont=font(10, plot_font))  #for legend


function  plot_inversion_results(save_path, mesh, zb_cell)

    inversion_results = load(joinpath(save_path, "inversion_results.jld2"))
    
    LOSS = inversion_results["LOSS"]
    PRED = inversion_results["PRED"]
    PARS = inversion_results["PARS"]
    sol = inversion_results["sol"]

    #plot_invesion_loss_time_history(LOSS)

    animate_inversion_process(LOSS, PRED, PARS, sol, mesh, zb_cell)
    
end

function plot_invesion_loss_time_history(LOSS)

    iter_numbers = collect(range(1,size(LOSS)[1]))

    p1 = plot(iter_numbers, LOSS, c=:black, dpi=300, yscale=:log10)
                
    xlims!(0, 30)
    #ylims!(-0.01, 2.01)
    xlabel!("Inversion iteration")
    ylabel!("Inversion loss")
    #plot!(legend=:bottomleft, fg_legend=:transparent, bg_legend=:transparent)

    savefig(p1, joinpath(save_path, "inversion_loss_history.png"))
end 

function animate_inversion_process(LOSS, PRED, PARS, sol, mesh, zb_cell)
    anim = Animation()

    let animCounter = 1

        for (curLoss, curPred, curPars) in zip(LOSS, PRED, PARS)

            animCounter = animCounter + 1

            if (animCounter % 1) == 0

                plt = plot(mesh.xCells, curPred[:, 1, end] + zb_cell, c=:aqua, label="free surface (predicted)", dpi=300)
                scatter!(mesh.xCells, Array(sol)[:,1,end]+zb_cell, mc=:blue, ms=2, ma=0.5, label="free surface (truth)")
                plot!(mesh.xCells, zb_cell, c=:black, label="bed")
                xlims!(-0.1, 10.1)
                ylims!(-0.01, 2.01)
                xlabel!("x (m)")
                ylabel!("z (m)")
                plot!(legend=:bottomleft, fg_legend=:transparent, bg_legend=:transparent)

                annotate!(5.0, 1.9,
                          text(
                            "Inversion iteration = $(animCounter)"
                            )
                         )

                lossString = @sprintf("Loss = %.3e", curLoss)
                annotate!(5.0, 1.7,
                          text(
                            #"Loss = $(round(curLoss, digits = 2))"
                            #"Loss = $(curLoss)"
                            lossString
                            )
                         )

                parsString = @sprintf("inletQ, exitH, n = %.5f, %.5f, %.5f", curPars[1], curPars[2], curPars[3])
                annotate!(5.0, 1.5,
                         #latexstring(
                         text(
                           #"inletQ, exitH, n = $(round(curPars[1], digits = 5)), $(round(curPars[2], digits = 5)), $(round(curPars[3], digits = 5))"
                           parsString
                           )
                        )

                frame(anim, plt)

            end

        end
    end

    #gif(anim, joinpath(save_path, "free_surface_1d.gif"), fps=15)
    gif(anim, joinpath(save_path, "inversion_process.mp4"), fps=1)
end