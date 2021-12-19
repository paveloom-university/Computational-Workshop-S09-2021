# This script solves a parabolic partial differential equation

# println('\n', " "^4, "> Loading the packages...")

# using Plots # Plotting

# # Use the GR backend for plots
# gr()

# # Change the default font for plots
# default(fontfamily = "Computer Modern", dpi = 300, legend = :outerright)

# Define the paths
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
TABLES = joinpath(ROOT_DIR, "tables")

# Prepare the output directories
mkpath(TABLES)

# Define the equation's functions
p(x) = x + 3
b(_, _) = 0
c(x, _) = -x
α₁(t) = 0
α₂(t) = -1
β₁(t) = 1
β₂(t) = 0

# Define the maximum time
T = 0.1

# Define the solution's functions
u(x, t) = x + 3t
f(x, t) = x^2 + 3 * x * t + 2
φ(x) = x
α(t) = 1
β(t) = 1 + 3t

"Check if a solution will be stable on the specified grid"
function is_stable(N, M; implicit = false)::Bool
    A = maximum(p, 0:0.0001:1)
    h = 1 / N
    τ = T / M
    ν = τ / h^2
    return A * ν ≤ 0.5
end

"Using provided N, find M, which gives a stable solution"
function find_stable(N)::Int
    M = 5
    while !is_stable(N, M)
        M *= 2
    end
    return M
end

# Define sets for numbers of nodes
N = [5, 10, 20]
M = find_stable.(N)

"Calculate the value of the differential operator L"
function L(um, x, t, h, i, k)::Float64
    return p(x[i] + h / 2) * (um[k, i+1] - um[k, i]) / h^2 -
           p(x[i] - h / 2) * (um[k, i] - um[k, i-1]) / h^2 +
           b(x[i], t[k]) * (um[k, i+1] - um[k, i-1]) / (2h) +
           c(x[i], t[k]) * um[k, i]
end

"Solve the problem for the specified grid, return the solution matrix"
function solve(N, M, h, τ, x, t)::Matrix{Float64}
    # Prepare a matrix for the solution
    # (number of columns is the number of nodes for x's)
    um = Matrix{Float64}(undef, M + 1, N + 1)
    # Compute the first layer
    um[1, :] .= φ(x)
    # Compute other layers
    for k = 2:M+1
        # Compute for i in 2:N
        for i = 2:N
            um[k, i] = um[k-1, i] +
                       τ * (L(um, x, t, h, i, k - 1) +
                            f(x[i], t[k-1]))
        end
        # Compute the rest
        um[k, 1] = (α(t[k]) + α₂(t[k]) * (4 * um[k, 2] - um[k, 3]) / (2h)) /
                   (α₁(t[k]) + 3 * α₂(t[k]) / (2h))
        um[k, N+1] = (β(t[k]) + β₂(t[k]) * (4 * um[k, N] - um[k, N-1]) / (2h)) /
                     (β₁(t[k]) + 3 * β₂(t[k]) / (2h))
    end
    return um
end

"""
Compare the approximate solution to the exact one,
find the biggest absolute difference between the nodes
"""
function max_difference(um, x, t)::Float64
    M1, N1 = size(um)
    # Find the maximum
    Δu = 0
    for i = 1:N1, k = 1:M1
        _Δu = abs(u(x[i], t[k]) - um[k, i])
        if _Δu > Δu
            Δu = _Δu
        end
    end
    return Δu
end

"Create the TeX tables and write them to disk"
function tables(um, N, M)
    open(joinpath(TABLES, "$N, $M.tex"), "w") do io
        # Create and write the table with the solution grid
        for k in Int.([1:M/5:M+1]...)
            s = "& "
            for i in Int.([1:N/5:N+1]...)
                s = "$(s)$(sprint(show, round(um[k, i]; digits=2))) & "
            end
            s = "$(s)\\\\"
            println(io, s)
        end
    end
    # Find the biggest difference between the nodes
    # max_difference(um, x, t)
end

# For each grid
for l in eachindex(N)
    # Compute the steps
    h = 1 / N[l]
    τ = T / M[l]
    # Compute the nodes
    x = [j * h for j = 0:N[l]]
    t = [j * τ for j = 0:M[l]]
    # Solve the problem
    um = solve(N[l], M[l], h, τ, x, t)
    # Create and write the tables
    tables(um, N[l], M[l])
end
