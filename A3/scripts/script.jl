# This script solves an elliptic partial differential equation

using Printf
using LinearAlgebra

# Define the paths
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
TABLES = joinpath(ROOT_DIR, "tables")

# Prepare the output directories
mkpath(TABLES)

# Define the problem
l_x = 1
l_y = 1
p(x, y) = 1
q(x, y) = 1
u(x, y) = x^3 * y + y^2 * x
f(x, y) = -(6 * x * y + 2 * x)

# Define the grids
N = [5, 10, 20]
M = copy(N)

# Set the expected precision of the results
ε = 1e-4
digits = 4

"Compute the value of the differential operator L"
function L(U, x, y, h_x, h_y, i, j)::Float64
    return Λ₁(U, x, y, h_x, i, j) + Λ₂(U, x, y, h_y, i, j)
end

"Compute the value of the differential operator Λ₁"
function Λ₁(U, x, y, h_x, i, j)::Float64
    return p(x[i] + h_x / 2, y[j]) * (U[i+1, j] - U[i, j]) / h_x^2 -
           p(x[i] - h_x / 2, y[j]) * (U[i, j] - U[i-1, j]) / h_x^2
end

"Compute the value of the differential operator Λ₂"
function Λ₂(U, x, y, h_y, i, j)::Float64
    return q(x[i], y[j] + h_y / 2) * (U[i, j+1] - U[i, j]) / h_y^2 -
           q(x[i], y[j] - h_y / 2) * (U[i, j] - U[i, j-1]) / h_y^2
end

"Compute ξ, which is needed to compute an optimal number of iterations"
function compute_xi(l_x, l_y, h_x, h_y)::Float64
    rx = 0:0.001:l_x
    ry = 0:0.001:l_y
    pairs = [(x, y) for x in rx, y in ry]
    c₁ = minimum(args -> p(args...), pairs)
    c₂ = maximum(args -> p(args...), pairs)
    d₁ = minimum(args -> q(args...), pairs)
    d₂ = maximum(args -> q(args...), pairs)
    k₁ = 4 / h_x^2
    k₂ = 4 / h_y^2
    arg₁ = (π * h_x) / (2 * l_x)
    arg₂ = (π * h_y) / (2 * l_y)
    δ = c₁ * k₁ * sin(arg₁)^2 +
        d₁ * k₂ * sin(arg₂)^2
    Δ = c₂ * k₁ * cos(arg₁)^2 +
        d₂ * k₂ * cos(arg₂)^2
    return δ / Δ
end

"Initialize a new matrix and compute the boundary values"
function boundary_values(N, M, x, y)::Matrix{Float64}
    # Prepare a matrix for the solution
    U = zeros(N + 1, M + 1)
    # Compute the boundary values
    U[:, 1] .= u.(x, 0)
    U[:, M+1] .= u.(x, l_y)
    U[1, 2:M] .= u.(0, y[2:M])
    U[N+1, 2:M] .= u.(l_x, y[2:M])
    return U
end

"Iterate using the simple iteration method"
function simple_iteration(Uₖ₋₁, ξ, ε, N, M, h_x, h_y, x, y)::Tuple{Matrix{Float64},Int}
    # Prepare a new matrix
    Uₖ = copy(Uₖ₋₁)
    # Compute the optimal number of iterations
    m = ceil(Int, log(1 / ε) / (2ξ))
    # Iterate enough times to achieve desired precision
    for _ in 1:m
        # Compute values at the inner nodes
        for i = 2:N, j = 2:M
            k₁ = p(x[i] - h_x / 2, y[j]) / h_x^2
            k₂ = p(x[i] + h_x / 2, y[j]) / h_x^2
            k₃ = q(x[i], y[j] - h_y / 2) / h_y^2
            k₄ = q(x[i], y[j] + h_y / 2) / h_y^2
            Uₖ[i, j] = (k₁ * Uₖ₋₁[i-1, j] +
                        k₂ * Uₖ₋₁[i+1, j] +
                        k₃ * Uₖ₋₁[i, j-1] +
                        k₄ * Uₖ₋₁[i, j+1] +
                        f(x[i], y[j])) /
                       (k₁ + k₂ + k₃ + k₄)
        end
        # Reassign the current iteration as the previous one
        Uₖ₋₁ = Uₖ
    end
    return Uₖ, m
end

"Iterate using the Seidel's method"
function seidel(Uₖ₋₁, ξ, ε, N, M, h_x, h_y, x, y)::Tuple{Matrix{Float64},Int}
    # Prepare a new matrix
    Uₖ = copy(Uₖ₋₁)
    # Compute the optimal number of iterations
    m = ceil(Int, log(1 / ε) / (4ξ))
    # Iterate enough times to achieve desired precision
    for _ in 1:m
        # Compute values at the inner nodes
        for i = N:-1:2, j = M:-1:2
            k₁ = p(x[i] - h_x / 2, y[j]) / h_x^2
            k₂ = p(x[i] + h_x / 2, y[j]) / h_x^2
            k₃ = q(x[i], y[j] - h_y / 2) / h_y^2
            k₄ = q(x[i], y[j] + h_y / 2) / h_y^2
            Uₖ[i, j] = (k₁ * Uₖ[i-1, j] +
                        k₂ * Uₖ₋₁[i+1, j] +
                        k₃ * Uₖ[i, j-1] +
                        k₄ * Uₖ₋₁[i, j+1] +
                        f(x[i], y[j])) /
                       (k₁ + k₂ + k₃ + k₄)
        end
        # Reassign the current iteration as the previous one
        Uₖ₋₁ = Uₖ
    end
    return Uₖ, m
end

"Iterate using the upper relaxation method"
function upper_relaxation(Uₖ₋₁, ξ, ε, N, M, h_x, h_y, x, y)::Tuple{Matrix{Float64},Int}
    # Prepare a new matrix
    Uₖ = copy(Uₖ₋₁)
    # Compute the optimal number of iterations
    m = ceil(Int, log(1 / ε) / sqrt(ξ))
    # Define ω in (0,2), which affects the speed of convergence
    ω = 1.0
    # Iterate enough times to achieve desired precision
    for _ in 1:m
        # Compute values at the inner nodes
        for i = N:-1:2, j = M:-1:2
            k₁ = p(x[i] + h_x / 2, y[j]) / h_x^2
            k₂ = p(x[i] - h_x / 2, y[j]) / h_x^2
            k₃ = q(x[i], y[j] + h_y / 2) / h_y^2
            k₄ = q(x[i], y[j] - h_y / 2) / h_y^2
            Uₖ[i, j] = Uₖ₋₁[i, j] +
                       ω * (k₁ * (Uₖ₋₁[i+1, j] - Uₖ₋₁[i, j]) -
                            k₂ * (Uₖ₋₁[i, j] - Uₖ[i-1, j]) +
                            k₃ * (Uₖ₋₁[i, j+1] - Uₖ₋₁[i, j]) -
                            k₄ * (Uₖ₋₁[i, j] - Uₖ[i, j-1]) +
                            f(x[i], y[j])) /
                       (k₁ + k₂ + k₃ + k₄)
        end
        # Reassign the current iteration as the previous one
        Uₖ₋₁ = Uₖ
    end
    return Uₖ, m
end

"Iterate, using the alternately triangular iterative method"
function triangular(Uₖ₋₁, ε, N, M, h_x, h_y, x, y)::Tuple{Matrix{Float64},Int}
    # Prepare a new matrix
    Uₖ = copy(Uₖ₋₁)
    # Compute the coefficients
    pairs = [(x, y) for x in 0:0.001:l_x, y in 0:0.001:l_y]
    c₁ = minimum(args -> p(args...), pairs)
    c₂ = maximum(args -> p(args...), pairs)
    d₁ = minimum(args -> q(args...), pairs)
    d₂ = maximum(args -> q(args...), pairs)
    k₁ = 4 / h_x^2
    k₂ = 4 / h_y^2
    δ = c₁ * k₁ * sin((π * h_x) / (2 * l_x))^2 +
        d₁ * k₂ * sin((π * h_y) / (2 * l_y))^2
    Δ = c₂ * k₁ + d₂ * k₂
    ω = 2 / sqrt(δ * Δ)
    η = δ / Δ
    γ₁ = δ / (2 + 2 * sqrt(η))
    γ₂ = δ / (4 * sqrt(η))
    ξ = γ₁ / γ₂
    κ₁ = ω / h_x^2
    κ₂ = ω / h_y^2
    τ = 2 / (γ₁ + γ₂)
    # Compute the optimal number of iterations
    m = ceil(Int, log(1 / ε) / log((1 + ξ) / (1 - ξ)))
    # Iterate enough times to achieve desired precision
    for _ in 1:m
        # Prepare intermediate matrices
        Ũ = zeros(size(Uₖ₋₁)...)
        U̅ = zeros(size(Uₖ₋₁)...)
        # Compute values at the inner nodes of the first intermediate matrix
        for i = 2:N, j = 2:M
            k₁ = κ₁ * p(x[i] - h_x / 2, y[j])
            k₂ = κ₂ * q(x[i], y[j] - h_y / 2)
            Ũ[i, j] = (k₁ * Ũ[i-1, j] + k₂ * Ũ[i, j-1] +
                        L(Uₖ₋₁, x, y, h_x, h_y, i, j) + f(x[i], y[j])) /
                       (1 + k₁ + k₂)
        end
        # Compute values at the inner nodes of the second intermediate matrix
        for i = N:-1:2, j = M:-1:2
            k₁ = κ₁ * p(x[i] + h_x / 2, y[j])
            k₂ = κ₂ * q(x[i], y[j] + h_y / 2)
            U̅[i, j] = (k₁ * U̅[i+1, j] + k₂ * U̅[i, j+1] + Ũ[i, j]) /
                       (1 + k₁ + k₂)
        end
        # Add up the matrices
        Uₖ[2:N, 2:M] .= Uₖ₋₁[2:N, 2:M] .+ τ .* U̅[2:N, 2:M]
        # Reassign the current iteration as the previous one
        Uₖ₋₁ = Uₖ
    end
    return Uₖ, m
end

"Iterate using the variable directions method"
function variable_directions(Uₖ₋₁, ε, N, M, h_x, h_y, x, y)::Tuple{Matrix{Float64},Int}
    # Prepare a new matrix and an intermediate matrix
    Uₖ = copy(Uₖ₋₁)
    U̅ = copy(Uₖ₋₁)
    # Compute the optimal number of iterations
    m = ceil(Int, N / (2 * π) * log(1 / ε))
    # Compute the coefficients
    pairs = [(x, y) for x in 0:0.001:l_x, y in 0:0.001:l_y]
    c₁ = minimum(args -> p(args...), pairs)
    c₂ = maximum(args -> p(args...), pairs)
    d₁ = minimum(args -> q(args...), pairs)
    d₂ = maximum(args -> q(args...), pairs)
    k₁ = 4 / h_x^2
    k₂ = 4 / h_y^2
    arg₁ = (π * h_x) / (2 * l_x)
    arg₂ = (π * h_y) / (2 * l_y)
    δ₁ = c₁ * k₁ * sin(arg₁)^2
    δ₂ = d₁ * k₂ * sin(arg₂)^2
    Δ₁ = c₂ * k₁ * cos(arg₁)^2
    Δ₂ = d₂ * k₁ * cos(arg₂)^2
    δ = min(δ₁, δ₂)
    Δ = max(Δ₁, Δ₂)
    τ = 2 / sqrt(δ * Δ)
    # Iterate enough times to achieve desired precision
    for _ in 1:m
        # Compute values at the inner nodes of the intermediate matrix
        for j = 2:M
            # Compute the linear system's matrix
            Ũ = Tridiagonal(
                # A's
                [repeat([τ / (2 * h_x^2)], N - 1); 0],
                # B's
                [1; repeat([-τ / h_x^2 + 1], N - 1); 1],
                # C's
                [0; repeat([τ / (2 * h_x^2)], N - 1)],
            )
            # Compute the linear system's right-hand side vector
            g = [
                Uₖ₋₁[1, j]
                [-Uₖ₋₁[i, j] -
                 τ / 2 * (Λ₂(Uₖ₋₁, x, y, h_y, i, j) +
                          f(x[i], y[j]))
                 for i in 2:N
                ]
                Uₖ₋₁[N+1, j]]
            # Compute the solution of the linear system
            U̅[:, j] .= Ũ \ g
        end
        # Compute values at the inner nodes of the input matrix
        for i = 2:N
            # Compute the linear system's matrix
            Ũ = Tridiagonal(
                # A's
                [repeat([τ / (2 * h_y^2)], N - 1); 0],
                # B's
                [1; repeat([-τ / h_y^2 + 1], N - 1); 1],
                # C's
                [0; repeat([τ / (2 * h_y^2)], N - 1)],
            )
            # Compute the linear system's right-hand side vector
            g = [
                Uₖ₋₁[i, 1]
                [-U̅[i, j] -
                 τ / 2 * (Λ₁(U̅, x, y, h_x, i, j) +
                          f(x[i], y[j]))
                 for j in 2:M
                ]
                Uₖ₋₁[i, M+1]]
            # Compute the solution of the linear system
            Uₖ[i, :] .= Ũ \ g
        end
        # Reassign the current iteration as the previous one
        Uₖ₋₁ = Uₖ
    end
    return Uₖ, m
end

"Create the TeX tables and write them to disk"
function tables(dir, U, N, M; digits = 8)
    # Prepare the output directory
    OUTPUT_DIR = joinpath(TABLES, dir)
    mkpath(OUTPUT_DIR)
    # Prepare the ranges
    ri = Int.([1:N/5:N+1]...)
    rj = Int.([1:M/5:M+1]...)
    # Create and write the table with the solution grid
    open(joinpath(OUTPUT_DIR, "$N, $M.tex"), "w") do io
        for i in ri
            s = ""
            for j in rj
                number = round(U[i, j]; digits)
                s = "$(s)& \$ $(Printf.format(Printf.Format("%.0$(digits)f"), number == 0 ? 0 : number)) \$ "
            end
            s = "$(s)\\\\"
            println(io, s)
        end
    end
end

"""
Compare the approximate solution to the exact one,
find the biggest absolute difference between the nodes
"""
function max_difference(U, x, y)::Float64
    M1, N1 = size(U)
    # Find the maximum
    Δu = 0
    for i = 1:N1, j = 1:M1
        _Δu = abs(u(x[i], y[j]) - U[i, j])
        if _Δu > Δu
            Δu = _Δu
        end
    end
    return Δu
end

# For each grid
for l in eachindex(N)
    # Compute the steps
    h_x = l_x / N[l]
    h_y = l_y / M[l]
    # Print the numbers of nodes
    println(
        '\n',
        " "^5, "N: ", N[l], '\n',
        " "^5, "M: ", M[l]
    )
    # Compute the nodes
    x = collect(0:h_x:1.0)
    y = collect(0:h_y:1.0)
    # Compute ξ
    ξ = compute_xi(l_x, l_y, h_x, h_y)
    # Create a new matrix and compute the boundary values
    U = boundary_values(N[l], M[l], x, y)
    # Iterate, using the simple iteration method
    U, m = simple_iteration(U, ξ, ε, N[l], M[l], h_x, h_y, x, y)
    # Create and write the TeX tables
    tables("simple_iteration", U, N[l], M[l]; digits)
    # Print the info
    println(
        '\n',
        " "^5, "> Simple iteration:", '\n',
        " "^5, "  m = ", m, '\n',
        " "^5, "  Δu = ", max_difference(U, x, y)
    )
    # Iterate, using the Seidel's method
    U, m = seidel(U, ξ, ε, N[l], M[l], h_x, h_y, x, y)
    # Create and write the TeX tables
    tables("seidel", U, N[l], M[l]; digits)
    # Print the info
    println(
        '\n',
        " "^5, "> Seidel:", '\n',
        " "^5, "  m = ", m, '\n',
        " "^5, "  Δu = ", max_difference(U, x, y)
    )
    # Iterate, using the upper relaxation method
    U, m = upper_relaxation(U, ξ, ε, N[l], M[l], h_x, h_y, x, y)
    # Create and write the TeX tables
    tables("upper_relaxation", U, N[l], M[l]; digits)
    # Print the info
    println(
        '\n',
        " "^5, "> Upper relaxation:", '\n',
        " "^5, "  m = ", m, '\n',
        " "^5, "  Δu = ", max_difference(U, x, y)
    )
    # Iterate, using the alternately triangular iterative method
    U, m = triangular(U, ε, N[l], M[l], h_x, h_y, x, y)
    # Create and write the TeX tables
    tables("triangular", U, N[l], M[l]; digits)
    # Print the info
    println(
        '\n',
        " "^5, "> Alternately triangular iterative method:", '\n',
        " "^5, "  m = ", m, '\n',
        " "^5, "  Δu = ", max_difference(U, x, y)
    )
    # Iterate, using the variable directions method
    U, m = variable_directions(U, ε, N[l], M[l], h_x, h_y, x, y)
    # Create and write the TeX tables
    tables("variable_directions", U, N[l], M[l]; digits)
    # Print the info
    println(
        '\n',
        " "^5, "> Variable directions method:", '\n',
        " "^5, "  m = ", m, '\n',
        " "^5, "  Δu = ", max_difference(U, x, y)
    )
end

println()
