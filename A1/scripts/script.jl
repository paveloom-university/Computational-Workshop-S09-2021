# This script solves a Fredholm equation of the first type

# Print about the packages
println('\n', " "^4, "> Loading the packages...")

using LinearAlgebra # Norm
using Optim # Optimization
using Plots # Plotting
using QuadGK # Gauss–Kronrod integration

# Use the GR backend for plots
gr()

# Change the default font for plots
default(fontfamily = "Computer Modern", dpi=300, legend=:outerright)

# Print about the computing
println(" "^4, "> Computing the solution (Test)...")

# Define the kernel (Test)
K(x, s) = exp(s * x)

# Define the right part of the equation (Test)
u(x) = (exp(x + 1) - 1) / (x + 1)

# Set the integration intervals
a = 0
b = 1
c = 0
d = 1

# Set the number of nodes
n = 100

# Set the initial value of the regularization parameter
α = 0.001

# Calculate the step
h = (b - a) / n

# Calculate the nodes for the s argument
s = [ a + 0.5 * h + (i - 1) * h for i in 1:n ]

# Calculate the nodes for the t argument
t = copy(s)

# Compute the g vector
g = [ quadgk(x -> K(x, s[i]) * u(x), c, d; rtol=1e-8)[1] for i in 1:n ]

# Prepare a matrix for the computation
Bα = Matrix{Float64}(undef, n, n)

# Prepare a vector for the solution
z = Vector{Float64}(undef, n)

# Prepare a range of nodes for the residual calculation
xr = range(c, d; length=1000)

# Compute the residual of the solution with the specified regularization parameter
function residual(θ::Vector{Float64})::Float64
    # Unpack the parameters
    α = θ[1]

    # Compute the Bα matrix
    for i in 1:n, j in 1:n
        αc = if (i == j == 1) || (i == j == n)
            α * (1 + 1 / h^2)
        elseif i == j
            α * (1 + 2 / h^2)
        elseif (i == j + 1) || (i + 1 == j)
            -α / h^2
        else
            0
        end
        Bα[i, j] = quadgk(x -> K(x, s[i]) * K(x, t[j]), c, d; rtol=1e-8)[1] * h + αc
    end

    # Compute the solution
    z .= Symmetric(Bα) \ g

    # Calculate the residual
    r = norm([ sum(K.(x, s) .* z .* h) for x in xr ] .- u.(xr))

    return r
end

# Optimize over the regularization parameter
res = Optim.optimize(
    residual,
    [α,],
    LBFGS(),
    Optim.Options(
        show_trace = false,
        extended_trace = true,
        store_trace = true,
    );
    inplace = false,
)

# Save the trace and the results
open("optimization.trace", "w") do io
    println(io, res.trace)
    println(
        io,
        " * Parameters:\n",
        " "^4, "α = $(res.minimum)\n",
    )
    show(io, res)
end

# Print the residual
println(
    '\n',
    " "^6, "Regularization parameter: ", res.minimizer[1], '\n',
    " "^6, "Residual: ", res.minimum, '\n',
)

# Print about the plotting
println(" "^4, "> Plotting the solution (Test)...")

# Plot the solution
p = plot(s, z; label="Приближенное решение", xlabel="s", ylabel="z(s)", legend=:outerright);

# Add the solution to the plot
plot!(p, x -> exp(x), a, b; label="Точное решение");

# Save the figure
savefig("$(@__DIR__)/../plots/test.pdf")

# Print about the figure
println('\n', " "^6, "* The figure `test.pdf` is saved. *", '\n')

# Print about the computing
println(" "^4, "> Computing the solution (Assignment)...")

# Define the kernel
K(x, s) = 1 / (1 + x + s)

# Define the right part of the equation
u(x) = 1 / sqrt(2 - x) * (log((2 + x) / (1 + x)) +
       2 * log((sqrt(3) + sqrt(2 - x)) / (2 + sqrt(2 - x))))

# Compute the g vector
g = [ quadgk(x -> K(x, s[i]) * u(x), c, d; rtol=1e-8)[1] for i in 1:n ]

# Optimize over the regularization parameter
res = Optim.optimize(
    residual,
    [α,],
    LBFGS(),
    Optim.Options(
        show_trace = false,
        extended_trace = true,
        store_trace = true,
    );
    inplace = false,
)

# Save the trace and the results
open("optimization.trace", "w") do io
    println(io, res.trace)
    println(
        io,
        " * Parameters:\n",
        " "^4, "α = $(res.minimum)\n",
    )
    show(io, res)
end

# Print the residual
println(
    '\n',
    " "^6, "Regularization parameter: ", res.minimizer[1], '\n',
    " "^6, "Residual: ", res.minimum, '\n',
)

# Print about the plotting
println(" "^4, "> Plotting the solution (Assignment)...")

# Plot the solution
p = plot(s, z; label="", xlabel="s", ylabel="z(s)");

# Save the figure
savefig("$(@__DIR__)/../plots/assignment.pdf")

# Print about the figure
println('\n', " "^6, "* The figure `assignment.pdf` is saved. *", '\n')
