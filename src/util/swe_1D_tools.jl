#some tools for 1D SWE

using Printf
using DifferentialEquations
using CSV 
using DataFrames
using Plots

#calculate the total water volume in the domain
function swe_1D_calc_total_water_volume(h, dx)
    total_water_volume = sum(h) * dx
    return total_water_volume
end

#save results 
function swe_1D_save_results(sol, total_water_volume, mesh, zb_cell, save_path)

    dx=mesh.dx 
    nCells=mesh.nCells
    xCells=mesh.xCells

    #solution at the end of the simulation
    sol_final = sol.u[end]

    #calculate total volume of water 
    for (Q, t) in tuples(sol)
        push!(total_water_volume, swe_1D_calc_total_water_volume(Q[:, 1], dx))
    end

    #save results to files the same directory of this jl file
    open(joinpath(save_path, "results.csv"), "w") do fo
        println(fo, "x, h, q, eta, zb")
        for i in 1:nCells
            @printf(fo, "%.6f, %.6f, %.6f, %.6f, %.6f\n", xCells[i], sol_final[i, 1], sol_final[i, 2], sol_final[i, 1] + zb_cell[i], zb_cell[i])
        end
    end

    open(joinpath(save_path, "total_water_volume.csv"), "w") do fo
        println(fo, "total_water_volume")
        for volume in total_water_volume
            println(fo, volume)
        end
    end
end


#make plots of the results
function swe_1D_make_plots(save_path)
    data_dataframe = CSV.read(joinpath(save_path, "results.csv"), DataFrame)
    data = Matrix(data_dataframe)  #x, h, q, eta, zb

    p1 = plot(data[:, 1], data[:, 4], c=:aqua, label="free surface", dpi=300)
    plot!(data[:, 1], data[:, 5], c=:black, label="bed")
    
    xlims!(-0.1, 10.1)
    ylims!(-0.01, 1.51)
    xlabel!("x (m)")
    ylabel!("z (m)")
    plot!(legend=:bottomleft, fg_legend=:transparent, bg_legend=:transparent)

    data_dataframe = CSV.read(joinpath(save_path, "total_water_volume.csv"), DataFrame)
    data = Matrix(data_dataframe)  #total_water_volume

    p2 = plot(data[:, 1])

    display(p1)
    #display(p2)

    savefig(p1, joinpath(save_path, "simulation_results_tEnd.png"))
end


#make animation
function swe_1D_make_animation(sol, mesh, zb_cell, save_path)

    anim = Animation()

    let animCounter = 1

        for (Q, t) in tuples(sol)

            animCounter = animCounter + 1

            if (animCounter % 10) == 0

                plt = plot(mesh.xCells, Q[:, 1] + zb_cell, c=:aqua, label="free surface", dpi=300)
                plot!(mesh.xCells, zb_cell, c=:black, label="bed")
                xlims!(-0.1, 10.1)
                ylims!(-0.01, 1.51)

                xlabel!("x (m)")
                ylabel!("z (m)")
                plot!(legend=:bottomleft, fg_legend=:transparent, bg_legend=:transparent)

                annotate!(5.0, 1.35,
                          text(
                            "Time = $(t) s"
                            )
                         )

                frame(anim, plt)

            end

        end
    end

    #gif(anim, joinpath(save_path, "free_surface_1d.gif"), fps=15)
    gif(anim, joinpath(save_path, "free_surface_1d.mp4"), fps=15)
end

