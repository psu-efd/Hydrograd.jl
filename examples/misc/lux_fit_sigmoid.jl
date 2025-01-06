using Lux
using Optim
using Random
using Zygote
using Plots

# Define the sigmoid function
my_sigmoid(x) = 1 / (1 + exp(-x))

# Generate training data
rng = Random.default_rng()
x_train = range(-10.0, 10.0, length=100)  # 100 points in the range [-10, 10]
y_train = my_sigmoid.(x_train)              # True values using the sigmoid function

# Define the neural network model
model = Lux.Chain(
    Lux.Dense(1, 16, tanh),  # First hidden layer: 1 input -> 16 units with tanh activation
    Lux.Dense(16, 16, tanh), # Second hidden layer: 16 units -> 16 units with tanh activation
    Lux.Dense(16, 1)         # Output layer: 16 units -> 1 output
)

# Initialize model parameters and state
ps, st = Lux.setup(rng, model)

# Define the loss function (mean squared error)
function loss(ps, st, x, y)
    y_pred, st_new = model(ps, st, x)
    return mean((y .- y_pred).^2), st_new
end

# Flatten parameters for optimization
flat_params, re = Flux.destructure(model)  # Flatten model parameters
flat_loss = function(p)                    # Loss function compatible with Optim
    ps_new = re(p)                         # Restructure parameters
    l, _ = loss(ps_new, st, x_train', y_train')  # Evaluate loss
    return l
end

# Train the model using Optim.jl
result = Optim.optimize(flat_loss, flat_params, Optim.BFGS(); maxiters=500)

# Restructure parameters after optimization
trained_ps = re(result.minimizer)

# Make predictions
y_pred, _ = model(trained_ps, st, x_train')

# Plot results
plot(
    x_train, y_train, label="True Sigmoid", lw=2, linestyle=:dash,
    xlabel="x", ylabel="f(x)", title="Approximating Sigmoid with Neural Network"
)
plot!(x_train, y_pred[:], label="NN Approximation", lw=2, color=:red)
