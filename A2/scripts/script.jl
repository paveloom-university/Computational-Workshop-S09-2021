# This script solves a parabolic partial differential equation

using LinearAlgebra

# Define the paths
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
TABLES = joinpath(ROOT_DIR, "tables")

# Prepare the output directories
mkpath(TABLES)

# Define the maximum time
T = 0.1

# Define the differential equation's functions
p(x) = x + 3
b(_, _) = 0
c(x, _) = -x
α₁(t) = 0
α₂(t) = -1
β₁(t) = 1
β₂(t) = 0

# Define the first problem's functions
u(x, t) = x + 3t
f(x, t) = x^2 + 3 * x * t + 2
φ(x) = x
α(t) = 1
β(t) = 1 + 3t

"Check if a solution will be stable on the specified grid"
function is_stable(N, M)::Bool
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

"Calculate the value of the differential operator L"
function L(um, x, t, h, i, k)::Float64
    return p(x[i] + h / 2) * (um[k, i+1] - um[k, i]) / h^2 -
           p(x[i] - h / 2) * (um[k, i] - um[k, i-1]) / h^2 +
           b(x[i], t[k]) * (um[k, i+1] - um[k, i-1]) / (2h) +
           c(x[i], t[k]) * um[k, i]
end

"Solve the problem for the specified grid, return the solution matrix"
function solve(N, M, h, τ, x, t; σ = 0)::Matrix{Float64}
    # Prepare a matrix for the solution
    # (number of columns is the number of nodes for x's)
    um = Matrix{Float64}(undef, M + 1, N + 1)
    # Compute the first layer
    um[1, :] .= φ.(x)
    if σ == 0
        # Use the explicit scheme for other layers
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
    else
        # Use the implicit scheme for other layers
        for k = 2:M+1
            # Compute the linear system's matrix
            tm = Tridiagonal(
                # A's
                [
                    [σ * (p(x[i] - h / 2) / h^2 - b(x[i], t[k]) / (2h)) for i = 2:N]
                    -β₂(t[k]) / h
                ],
                # B's
                [
                    α₁(t[k]) + α₂(t[k]) / h
                    -[1 / τ + σ * ((p(x[i] + h / 2) + p(x[i] - h / 2)) / h^2 - c(x[i], t[i])) for i = 2:N]
                    β₁(t[k]) + β₂(t[k]) / h
                ],
                # C's
                [
                    -α₂(t[k]) / h
                    [σ * (p(x[i] + h / 2) / h^2 + b(x[i], t[k]) / (2h)) for i = 2:N]
                ],
            )
            # Compute the linear system's right-hand side vector
            g = [
                α(t[k])
                [-um[k-1, i] / τ - (1 - σ) * L(um, x, t, h, i, k - 1) - f(x[i], t[k] - (1 - σ) * τ) for i = 2:N]
                β(t[k])
            ]
            # Compute the solution of the linear system
            um[k, :] = tm \ g
        end
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
function tables(dir, digits, um, N, M)
    # Prepare the output directories
    OUTPUT_DIR = joinpath(TABLES, dir)
    mkpath(OUTPUT_DIR)
    # Create and write the table with the solution grid
    open(joinpath(OUTPUT_DIR, "$N, $M.tex"), "w") do io
        for k in Int.([1:M/5:M+1]...)
            s = ""
            for i in Int.([1:N/5:N+1]...)
                s = "$(s)& \$ $(sprint(show, round(um[k, i]; digits))) \$ "
            end
            s = "$(s)\\\\"
            println(io, s)
        end
    end
end

"""
Solve the specified problem for different grids,
create the TeX tables and write them to disk
"""
function solve_all(dir, digits)
    # Print the output directory
    println('\n', " "^5, "> Output for u(x, t) = ", dir)
    # Define grids
    N = [5, 10, 20]
    Mₑ = find_stable.(N)
    Mᵢ = repeat([Int(1 / (0.1 / 100))], length(N))
    # For each grid
    for l in eachindex(N)
        # Compute the steps
        h = 1 / N[l]
        τₑ = T / Mₑ[l]
        τᵢ = T / Mᵢ[l]
        # Compute the nodes
        x = [j * h for j = 0:N[l]]
        tₑ = [j * τₑ for j = 0:Mₑ[l]]
        tᵢ = [j * τᵢ for j = 0:Mᵢ[l]]
        # Solve the problem with the explicit scheme
        um = solve(N[l], Mₑ[l], h, τₑ, x, tₑ)
        # Print info about the solution
        println(
            '\n',
            " "^5, "σ: 0", '\n',
            " "^5, "N: ", N[l], '\n',
            " "^5, "M: ", Mₑ[l], '\n',
            " "^5, "Δu: ", max_difference(um, x, tₑ)
        )
        # Create and write the tables
        tables(joinpath(dir, "σ = 0"), digits, um, N[l], Mₑ[l])
        # Solve the problem with the implicit scheme (σ = 0.5)
        um = solve(N[l], Mᵢ[l], h, τᵢ, x, tᵢ; σ = 0.5)
        # Print info about the solution
        println(
            '\n',
            " "^5, "σ: 0.5", '\n',
            " "^5, "N: ", N[l], '\n',
            " "^5, "M: ", Mᵢ[l], '\n',
            " "^5, "Δu: ", max_difference(um, x, tᵢ)
        )
        # Create and write the tables
        tables(joinpath(dir, "σ = 0.5"), digits, um, N[l], Mᵢ[l])
        # Solve the problem with the implicit scheme (σ = 1)
        um = solve(N[l], Mᵢ[l], h, τᵢ, x, tᵢ; σ = 1.0)
        # Print info about the solution
        println(
            '\n',
            " "^5, "σ: 1.0", '\n',
            " "^5, "N: ", N[l], '\n',
            " "^5, "M: ", Mᵢ[l], '\n',
            " "^5, "Δu: ", max_difference(um, x, tᵢ)
        )
        # Create and write the tables
        tables(joinpath(dir, "σ = 1.0"), digits, um, N[l], Mᵢ[l])
    end
end

# Solve the first problem
solve_all("x + 3t", 2)

# Define the second problem's functions
u(x, t) = x^3 + t^3
f(x, t) = x^4 + x * t^3 - 9 * x^2 + 3 * t^2 - 18 * x
φ(x) = x^3
α(t) = 0
β(t) = 1 + t^3

# Solve the second problem
solve_all("x^3 + t^3", 6)

println()
