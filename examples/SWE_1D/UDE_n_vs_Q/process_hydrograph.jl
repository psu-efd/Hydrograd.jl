using CSV, DataFrames

function load_hydrograph(save_path)
    # Reading a CSV file
    df = CSV.read(joinpath(save_path, "hydrograph2.csv"), DataFrame)
    #println(data)

    # Convert the dataframe into an array
    # using the Matrix function
    array_data = Matrix(df)
    #println(array_data)

    #make a linear regression
    hydrograph_t = array_data[:,1]
    hydrograph_Q = array_data[:,2]

    p1=plot(hydrograph_t, hydrograph_Q)
    display(p1)
end