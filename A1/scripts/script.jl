# This script solves a Fredholm equation of the first type

println('\n', " "^4, "> Loading the packages...")

using LinearAlgebra # Norm
using Optim # Optimization
using Plots # Plotting
using Printf # Formatted printing
using QuadGK # Gauss–Kronrod integration

# Use the GR backend for plots
gr()

# Change the default font for plots
default(fontfamily="Computer Modern", dpi=300, legend=:outerright)

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

# Calculate the step
h = (b - a) / n

# Calculate the nodes for the s argument
s = [ a + 0.5 * h + (i - 1) * h for i in 1:n ]

# Set the initial value of the regularization parameter
α = 0.001

# Prepare a range of nodes for the residual calculation
xr = range(c, d; length=1000)

# Compute the residual of the solution with the specified number of nodes and the regularization parameter
function calculate(n::Int, α::Float64)::Tuple{Vector{Float64},Float64}
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

    # Compute the Bα matrix
    for i in 1:n, j in 1:n
        Bα[i, j] = quadgk(x -> K(x, s[i]) * K(x, t[j]), c, d; rtol=1e-8)[1] * h + (i == j ? α : 0)
    end

    # Compute the solution
    z = Symmetric(Bα) \ g

    # Calculate the residual
    r = norm([ sum(K.(x, s) .* z .* h) for x in xr ] .- u.(xr))
    
    return z, r
end

# Compute the residual of the solution with the specified regularization parameter (optimization only)
function residual(θ::Vector{Float64})::Float64
    return calculate(100, θ[1])[2]
end

# Optimize over the regularization parameter
res = Optim.optimize(
    residual,
    [α,],
    LBFGS(),
    Optim.Options(
        show_trace=false,
        extended_trace=true,
        store_trace=true,
    );
    inplace=false,
)

# Save the trace and the results
open("test_optimization.trace", "w") do io
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

println(" "^4, "> Plotting the solution (Test)...")

# Recalculate the solution
z = calculate(n, res.minimizer[1])

# Plot the solution
p = plot(s, z; label="Приближенное решение", xlabel="s", ylabel="z(s)");

# Add the solution to the plot
plot!(p, x -> exp(x), a, b; label="Точное решение");

# Save the figure
savefig(p, "$(@__DIR__)/../plots/test.pdf")

println('\n', " "^6, "* The figure `test.pdf` is saved. *", '\n')

println(" "^4, "> Plotting the residuals against the number of nodes (Test)...")

# Define the `n`s to plot against
ns = [ i for i in 10:300 ]

# Calculate the residuals for these `n`s
rs = @. getindex(calculate(ns, res.minimizer[1]), 2)

# Make a plot of the residual against `n`
p = plot(ns, rs; label="", xlabel="n", ylabel="r");

# Save the figure
savefig(p, "$(@__DIR__)/../plots/test_residual.pdf")

println('\n', " "^6, "* The figure `test_residual.pdf` is saved. *", '\n')

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
        show_trace=false,
        extended_trace=true,
        store_trace=true,
    );
    inplace=false,
)

# Save the trace and the results
open("assignment_optimization.trace", "w") do io
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

println(" "^4, "> Plotting the solution (Assignment)...")

# Recalculate the solution
z = calculate(n, res.minimizer[1])

# Plot the solution
p = plot(s, z; label="", xlabel="s", ylabel="z(s)");

# Save the figure
savefig(p, "$(@__DIR__)/../plots/assignment.pdf")

println('\n', " "^6, "* The figure `assignment.pdf` is saved. *", '\n')

println(" "^4, "> Plotting the residuals against the number of nodes (Assignment)...")

# Define the `n`s to plot against
ns = [ i for i in 10:300 ]

# Calculate the residuals for these `n`s
rs = @. getindex(calculate(ns, res.minimizer[1]), 2)

# Make a plot of the residual against `n`
p = plot(ns, rs; label="", xlabel="n", ylabel="r");

# Save the figure
savefig(p, "$(@__DIR__)/../plots/assignment_residual.pdf")

println('\n', " "^6, "* The figure `assignment_residual.pdf` is saved. *", '\n')

println(" "^4, "> Calculating the solutions for tables (Assignment)...")

# Define the header
header = "α = 1e-4" * " "^7 * "1e-5" * " "^11 * "1e-6" * " "^11 * "1e-7" *
         " "^11 * "1e-8" * " "^11 * "1e-9" * " "^11 * "1e-10" * " "^10 * "1e-11" *
         " "^10 * "1e-12" * " "^10 * "1e-13" * " "^10 * "1e-14" * " "^10 * "1e-15"

# Calculate the solutions for different numbers of nodes and regularization parameters
for n in [100, 200, 300]
    step = Int(n * 0.05)
    data = Matrix{Float64}(undef, 20, 12)
    format = Printf.Format("%.8e "^11 * "%.8e\n")
    open("$(@__DIR__)/../tables/$(n).dat", "w") do io
        println(io, header)
        for (i, α) in enumerate([1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15])
            sol, _ = calculate(n, α)
            data[:, i] = sol[1:step:end]
        end
        for i in 1:20
            Printf.format(io, format, data[i, :]...)
        end
    end
end

println('\n', " "^6, "* Tables `100.dat`, `200.dat`, and `300.dat` are saved. *", '\n')
