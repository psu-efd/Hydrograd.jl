#define variables 

#control file
control_file = nothing

#settings
settings = nothing

#define a swe_2D_constants object with some values from the control file
swe_2D_constants = nothing

#read data from SRH-2D hydro, geom, and material files; it aslo create the 2D mesh.
srh_all_Dict = nothing

#Get the 2D mesh 
my_mesh_2D = nothing

#mesh related variables which might be mutable (should not be in struct mesh_2D)
#nodeCoordinates: Float64 2D array [numOfNodes, 3]
nodeCoordinates = nothing

zb_cells =nothing
zb_ghostCells = nothing 
zb_faces = nothing
S0 = nothing

zb_cells_truth = nothing

#get the true Manning's n and inlet discharges
srhhydro_ManningsN_Dict = nothing

ManningN_values_truth = nothing

#get the true inlet discharges (could be nothing if no inletQ_BCs)
srhhydro_inletQ_Dict = nothing

inlet_discharges_truth = nothing

#define the parameters array (nothing for forward simulation)
zb_cells_param = nothing
ManningN_list_param = nothing
inlet_discharges_param = nothing

#preprocess: create a array of model parameters for 2D shallow water equations
#params_array: the 1D array of all parameters (zb_cells_param, ManningN_list_param, inlet_discharges_param)
#active_range: the range of active parameters 
#param_ranges: the range of each parameter
params_array =nothing
active_range = nothing
param_ranges = nothing

#Initial setup of Manning's n using the SRH-2D data (if performing inversion on Manning's n, ManningN_cells will be updated later in the inversion process)
ManningN_cells = nothing
ManningN_ghostCells = nothing

# setup initial conditions 
wse = nothing          #free surface elevation at cells 
h = nothing            #water depth at cells 
q_x = nothing          #q_x=hu at cells 
q_y = nothing          #q_y=hv at cells 

wse_ghostCells = nothing          #free surface elevation at ghost cells 
h_ghostCells = nothing            #water depth at ghost cells 
q_x_ghostCells = nothing          #q_x=hu at ghost cells 
q_y_ghostCells = nothing          #q_y=hv at ghost cells 

total_water_volume = []   #total volume of water in the domain 

#setup initial condition for wse, h, q_x, q_y
initial_condition_options = nothing
initial_condition_constant_values = nothing
initial_condition_values_from_file = nothing


#create and preprocess boundary conditions: boundary_conditions only contains the static information of the boundaries.
boundary_conditions = nothing
inletQ_TotalQ = nothing
inletQ_H = nothing
inletQ_A = nothing
inletQ_ManningN = nothing
inletQ_Length = nothing
inletQ_TotalA = nothing
inletQ_DryWet = nothing
exitH_WSE = nothing
exitH_H = nothing
exitH_A = nothing
wall_H = nothing
wall_A = nothing
symm_H = nothing
symm_A = nothing

#set up initial condition for ODE solver
Q0 = nothing

Q_ghost = nothing   #ghost cell values of Q

# Create the ODEFunction with the typed function
ode_f = nothing

jac_sparsity = nothing

# time information (the same for forward simulation, inversion, and sensitivity analysis)
tspan = nothing
dt = nothing
t = nothing

dt_save = nothing
t_save = nothing