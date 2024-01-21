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
function swe_1D_save_results(sol, total_water_volume, mesh, zb_cell)

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
    open(joinpath(dirname(@__FILE__), "results.csv"), "w") do fo
        println(fo, "x, h, q, eta, zb")
        for i in 1:nCells
            @printf(fo, "%.6f, %.6f, %.6f, %.6f, %.6f\n", xCells[i], sol_final[i, 1], sol_final[i, 2], sol_final[i, 1] + zb_cell[i], zb_cell[i])
        end
    end

    open(joinpath(dirname(@__FILE__), "total_water_volume.csv"), "w") do fo
        println(fo, "total_water_volume")
        for volume in total_water_volume
            println(fo, volume)
        end
    end
end


#make plots of the results
function swe_1D_make_plots()
    data_dataframe = CSV.read(joinpath(dirname(@__FILE__), "results.csv"), DataFrame)
    data = Matrix(data_dataframe)  #x, h, q, eta, zb

    p1 = plot(data[:, 1], data[:, 4], label="free surface")
    plot!(data[:, 1], data[:, 5], label="bed")
    plot!(legend=:outerbottom, fg_legend=:transparent)

    data_dataframe = CSV.read(joinpath(dirname(@__FILE__), "total_water_volume.csv"), DataFrame)
    data = Matrix(data_dataframe)  #total_water_volume

    p2 = plot(data[:, 1])

    display(p1)
    #display(p2)
end


#make animation
function swe_1D_make_animation(sol, mesh, fields)

    anim = Animation()

    let animCounter = 1

        for (Q, t) in tuples(sol)

            animCounter = animCounter + 1

            if (animCounter % 10) == 0

                plt = plot(mesh.xCells, Q[:, 1] + fields.zb_cell, label="free surface")
                plot!(mesh.xCells, fields.zb_cell, label="bed")
                xlims!(-0.1, 10.1)
                ylims!(-0.1, 2.0)
                plot!(legend=:bottomleft, fg_legend=:transparent)
                frame(anim, plt)

            end

        end
    end

    #gif(anim, joinpath(dirname(@__FILE__), "free_surface_1d.gif"), fps=15)
    gif(anim, joinpath(dirname(@__FILE__), "free_surface_1d.mp4"), fps=15)
end

