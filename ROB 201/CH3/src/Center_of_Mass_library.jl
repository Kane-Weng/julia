#======================== KEY FUNCTION: COMPUTE COM ===============================#
function computeCoM(f,g, xmin, xmax)
    # set up functions
    areaIntegrand(x) = f(x) - g(x)
    MxIntegrand(x) = x* areaIntegrand(x)
    MyIntegrand(x) = 0.5*(f(x)^2 - g(x)^2)

    #quadgk
    A, = quadgk(areaIntegrand, xmin, xmax)
    Mx, = quadgk(MxIntegrand, xmin, xmax)
    My, = quadgk(MyIntegrand, xmin, xmax)

    xc = Mx / A    
    yc = My / A

    # return values    
    return (xc, yc)
end

function plot_custom_robot_link(f, g, L, joint_offset, includeCoordinates=false, joint_markersize=20)
    # f = upper shape of the link
    # g = lower shape
    # L = overall length of the link
    # joint_offset is distance from the left and right ends of the link for placing the revolute joints
    # Origin of the link is placed at the leftmost end of the body
    # includeCoordinates = boolean variable to plot coordinate axes or not
    # joint_markersize = self-evident. Good value is 20

    # Define the revolute joint locations
    joint1 = (joint_offset, 0.0)
    joint2 = (L - joint_offset, 0.0)

    # Compute x_vals and y_vals
    x_range = 0:(L / 1000.0):L # high resolution
    x_vals = collect(x_range) 
    y_vals_upper = f.(x_vals) 
    y_vals_lower = g.(x_vals)

    ymax = maximum(y_vals_upper)
    ymin = minimum(y_vals_lower)

    xlim_min = -joint_offset
    xlim_max = L + joint_offset

    ylim_min = minimum([1.25 * ymin, -2 * joint_offset])
    ylim_max = maximum([1.25 * ymax,  2 * joint_offset])

    # Create the plot
    p1 = plot(
        legend = false,
        aspect_ratio = :equal,
        ylims = (ylim_min, ylim_max),
        xlims = (xlim_min, xlim_max)
    )

    # Draw the custom robot link shape
    plot!(
        [x_vals; reverse(x_vals)],
        [y_vals_upper; reverse(y_vals_lower)],
        color = :blue,
        linewidth = 5,
        seriestype = :shape,
        linecolor = :blue,
        fillalpha = 1,
        fillcolor = :lightgray,
        label = nothing
    )

    # Draw the revolute joints
    scatter!(
        [joint1[1], joint2[1]],
        [joint1[2], joint2[2]],
        color = :red,
        markersize = joint_markersize,
        shape = :circle,
        label = nothing
    )

    if includeCoordinates
        # Draw the coordinate system centered at the left end of object (longer and thicker axes)
        arrow_length = 0.9 * ylim_max
        quiver!([0.0], [0.0], quiver = ([arrow_length], [0]), color = :black, linewidth = 4)
        quiver!([0.0], [0.0], quiver = ([0], [arrow_length]), color = :black, linewidth = 4)

        # Place a dot at the origin for enhanced visibility
        scatter!([0, 0], [0, 0], color = :gray, markersize = joint_markersize / 2, shape = :circle, label = nothing)

        # Add x and y labels on the axes
        annotate!((0.9 * arrow_length, -arrow_length / 6), L"x") 
        annotate!((-arrow_length / 6, 0.9 * arrow_length), text(L"y", 20)) 
    end

    # Return the plot and some useful values
    return (
        plot = p1,
        xmin = x_vals[1],
        xmax = x_vals[end],
        ymin = ymin,
        ymax = ymax,
        joint1 = joint1,
        joint2 = joint2
    )
end

function center_of_mass_symbol(center::Tuple{Real, Real}, radius::Real, thickness::Real)
    # Create the outer circle
    plot!(
        θ -> center[1] + radius * cos(θ),
        θ -> center[2] + radius * sin(θ),
        0, 2π,
        linewidth = thickness,
        color = :black,
        seriestype = :path,
        fillalpha = 0,
        label = nothing
    )

    # Create the plus sign
    plot!(
        [center[1], center[1]],
        [center[2] - radius, center[2] + radius],
        linewidth = thickness,
        color = :black,
        label = nothing
    )
    plot!(
        [center[1] - radius, center[1] + radius],
        [center[2], center[2]],
        linewidth = thickness,
        color = :black,
        label = nothing
    )

    # Fill in quadrant 1 (upper right)
    θ_vals = range(0, π/2, length = 100)
    x_vals = [center[1]; center[1] .+ radius .* cos.(θ_vals); center[1]]
    y_vals = [center[2]; center[2] .+ radius .* sin.(θ_vals); center[2]]
    plot!(x_vals, y_vals, fillcolor = :black, seriestype = :shape, linecolor = :transparent, label = nothing)

    # Fill in quadrant 3 (lower left)
    θ_vals = range(π, 3π/2, length = 100)
    x_vals = [center[1]; center[1] .+ radius .* cos.(θ_vals); center[1]]
    y_vals = [center[2]; center[2] .+ radius .* sin.(θ_vals); center[2]]
    plot!(x_vals, y_vals, fillcolor = :black, seriestype = :shape, linecolor = :transparent, label = nothing)

    # Fill in quadrant 2 (upper left)
    θ_vals = range(π/2, π, length = 100)
    x_vals = [center[1]; center[1] .+ radius .* cos.(θ_vals); center[1]]
    y_vals = [center[2]; center[2] .+ radius .* sin.(θ_vals); center[2]]
    plot!(x_vals, y_vals, fillcolor = :white, seriestype = :shape, linecolor = :transparent, label = nothing)

    # Fill in quadrant 4 (lower right)
    θ_vals = range(3π/2, 2π, length = 100)
    x_vals = [center[1]; center[1] .+ radius .* cos.(θ_vals); center[1]]
    y_vals = [center[2]; center[2] .+ radius .* sin.(θ_vals); center[2]]
    plot!(x_vals, y_vals, fillcolor = :white, seriestype = :shape, linecolor = :transparent, label = nothing)
end
