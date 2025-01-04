using Lux, Random, Zygote, Optimisers, ComponentArrays, Plots, Statistics, NNlib

# Generate data
x_train = reshape(collect(range(0, 1, 100)), 1, :)  # Shape: (1, 100)
y_train = reshape(0.01 .+ 0.03 .* sin.(2π * x_train) .+ 0.005 .* randn(size(x_train)), 1, :)

@show size(x_train)
@show size(y_train)

# Define model
n_low, n_upper = -0.02, 0.04
model = Chain(
    Dense(1 => 16, relu),
    Dense(16 => 16, relu),
    Dense(16 => 1)
)

# Forward pass
function predict(model, x, ps, st)
    y, st = model(x, ps, st)
    scaled_output = n_low .+ (n_upper - n_low) .* sigmoid.(y)
    return scaled_output, st
end

# Loss function
function loss(model, ps, st, x, y)
    pred, _ = predict(model, x, ps, st)
    return mean((y .- pred).^2)
end

# Training loop
function train_model!(model, x_train, y_train; epochs=1000, lr=0.01)
    opt = Optimisers.Adam(lr)
    rng = Random.default_rng()
    ps, st = Lux.setup(rng, model)
    opt_state = Optimisers.setup(opt, ps)

    for epoch in 1:epochs
        grads = Zygote.gradient(ps -> loss(model, ps, st, x_train, y_train), ps)[1]
        opt_state, ps = Optimisers.update(opt_state, ps, grads)
        
        if epoch % 100 == 0
            println("Epoch $epoch, Loss: ", loss(model, ps, st, x_train, y_train))
        end
    end
    return ps, st
end

# Train
ps, st = train_model!(model, x_train, y_train)

# Predict using the trained model
x_test = [0.5]  # Example test input: water depth = 2.0
n_predicted, _ = predict(model, x_test, ps, st)
println("Predicted Manning's n for h = $x_test: $n_predicted")

# Plot predictions vs truth data
function plot_predictions(model, x_train, y_train)
    y_pred, _ = predict(model, x_train, ps, st)
    scatter(vec(x_train), vec(y_train), label="Truth", 
           xlabel="Water Depth (h)", ylabel="Manning's n", 
           title="Prediction vs Truth", legend=:bottomright)
    plot!(vec(x_train), vec(y_pred), label="Prediction", color=:red)
end

plot_predictions(model, x_train, y_train)
