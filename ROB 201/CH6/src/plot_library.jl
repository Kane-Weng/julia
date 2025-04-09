#===============MESHGRID HELPER FUNCTION=======================#
function meshgrid(x::AbstractVector, y::AbstractVector)
    X = repeat(x, 1, length(y))
    Y = repeat(y', length(x), 1)
    return X, Y 
end

#===============NOT SUGGESTED TO RUN, ALMOST CRASHED====================#
function plot_circle_hyperbola(r1, r2)
    # Define the functions for f1 and f2
    f1(x1, x2) = x1^2 + x2^2 - r1^2 # Circle
    f2(x1, x2) = x1^2 - x2^2 - r2^2 # Hyperbola with two branches

    # Define the grid for plotting
    x1_vals = LinRange(-2, 2, 5000)
    x2_vals = LinRange(-2, 2, 5000)

    # Create meshgrid-like structure using broadcasting
    X1, X2 = meshgrid(x1_vals, x2_vals)

    # Compute the values of f1 and f2 over the grid
    F1 = f1.(X1, X2)
    F2 = f2.(X1, X2)

    # Plot contour lines for both f1 and f2, ensure levels=[0]
    plt1 = contour(x1_vals, x2_vals, F1, levels=[0], linewidth=2, label="f1(x1, x2) = 0", color=:blue, colorbar=false, aspect_ratio=1)
    contour!(x1_vals, x2_vals, F2, levels=[0], linewidth=2, label="f2(x1, x2) = 0", color=:red, colorbar=false)

    # Mark the root found by NLsolve
    scatter!([root[1]], [root[2]], marker=:star, markersize=10, color=:green, label="Root")

    # Add labels and title
    xlabel!("x1")
    ylabel!("x2")
    title!("Solution of System of Equations")

    # Display the plot
    display(plt1)
end