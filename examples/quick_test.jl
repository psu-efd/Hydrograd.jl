using Zygote
using Hydrograd

h_ghost_local = rand(Float64, 10)
current_ghostCellIDs = [1, 2, 3, 4, 5]
current_internalCellIDs = [6, 7, 8, 9, 10]
h = ones(Float64, 10)

h_ghost_local = update_1d_array(h_ghost_local, current_ghostCellIDs, h[current_internalCellIDs])