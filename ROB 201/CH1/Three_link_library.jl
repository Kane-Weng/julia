function plot_pendulum()
    # constants for the pendulums
    lengths = [2,1.5,1]
    angles = [20,20,20]
    ball_radius = 0.1
    line_extension = 0.5
    bar_thickness = 0.05

    # convert angles to radians
    angles_rad = [a * (pi/180) for a in angles]

    # initialize starting point
    x,y = 0, 0

    # plot
    p1 = plot(; xlim=(-0.5,4), ylim=(-0.5,3), aspect_ratio=:equal, legend=false)

    x_prev, y_prev = x,y
    for i in 1:3
        # calculate end point of the pendulum
        x_end = x + lengths[i] * cos(sum(angles_rad[1:i]))
        y_end = y + lengths[i] * sin(sum(angles_rad[1:i]))

        # plot pendulum bar
        plot!([x,x_end], [y,y_end], color=:blue, linewidth=bar_thickness*100, line=:solid)

        # plot the ball at the end
        scatter!([x_end], [y_end], color=:red, markersize=ball_radius*100)

        # plot the dotted line extension
        if i<3
            x_dotted = x_end + line_extension * cos(sum(angles_rad[1:i]))
            y_dotted = y_end + line_extension * sin(sum(angles_rad[1:i]))
            plot!([x_end,x_dotted],[y_end,y_dotted], color=:black, linewidth=1, line=:solid)
        end

        # Annotate angles with arcs 
        arc_radius = 0.5
        if i==1
            t = range(0, angles_rad[i], length=100)
            x_arc = arc_radius .* cos.(t)
            y_arc = arc_radius .* sin.(t)
            plot!(x_arc, y_arc, color=:red, linewidth=1.5)
            annotate!(0.8, 0.13, text(L"\theta_1", 12))
            annotate!(x_end, y_end+.3, text(L"p_1", 15))
        elseif i==2
            t = range(angles_rad[1], sum(angles_rad[1:2]), length=100)
            x_arc = x .+ arc_radius .* cos.(t)
            y_arc = y .+ arc_radius .* sin.(t)
            plot!(x_arc, y_arc, color=:red, linewidth=1.5)
            annotate!(x+0.6, y+0.3, text(L"\theta_2", 12))
            annotate!(x_end-.1, y_end+.3, text(L"p_2", 15))
        else
            t = range(sum(angles_rad[1:2]), sum(angles_rad[1:3]), length=100)
            x_arc = x .+ arc_radius .* cos.(t)
            y_arc = y .+ arc_radius .* sin.(t)
            plot!(x_arc, y_arc, color=:red, linewidth=1.5)
            annotate!(x+0.5, y+0.5, text(L"\theta_3", 12))
            annotate!(x_end, y_end+.3, text(L"p_3", 15))
        end

        x_prev, y_prev = x,y
        x,y = x_end, y_end
    end

    # add base line
    x_base_line = 2*line_extension*cos(0)
    y_base_line = 2*line_extension*sin(0)
    plot!([0,x_base_line], [0,y_base_line], color=:black, linewidth=2, line=:arrow, arrow=:closed)
    annotate!(x_base_line + 0.1, -0.1, text("x", 12))

    x_base_line = 2*line_extension*cos(pi/2)
    y_base_line = 2*line_extension*sin(pi/2)
    plot!([0,x_base_line], [0,y_base_line], color=:black, linewidth=2, line=:arrow, arrow=:closed)
    annotate!(-0.2, y_base_line+0.05, text("y", 12))

    # add encoders
    annotate!(-0.3, 2.6, text("Encoders installed here", 14, :left, fontweight=:bold))

    # arrow to connection
    plot!([0.5,0.1], [2.4,0.2], line=:arrow, arrow=:closed, color=:green, linewidth=1.5)

    x_conn1 = lengths[1] * cos(angles_rad[1])
    y_conn1 = lengths[1] * sin(angles_rad[1])
    plot!([0.5,x_conn1-0.15], [2.4,y_conn1+0.2], line=:arrow, arrow=:closed, color=:green, linewidth=1.5)

    x_conn2 = x_conn1 + lengths[2] * cos(sum(angles_rad[1:2]))
    y_conn2 = y_conn1 + lengths[2] * sin(sum(angles_rad[1:2]))
    plot!([0.5,x_conn2-0.15], [2.4,y_conn2+0.1], line=:arrow, arrow=:closed, color=:green, linewidth=1.5)

    display(p1)
    #png(p1, "robotKinematicChainWithEncoderAnnotation")
end

function modelParameters()
    # modelParameters
    L1, L2, L3 = [2,1.5,1]
    return (L1=L1,L2=L2,L3=L3)
end

function linkPositions(th1,th2,th3)
    params = modelParameters()
    p0 = [0;0]
    p1 = p0 + [params.L1 * cos(th1); params.L1 * sin(th1)]
    p2 = p1 + [params.L2 * cos(th1+th2); params.L2 * sin(th1+th2)]
    p3 = p2 + [params.L3 * cos(th1+th2+th3); params.L3 * sin(th1+th2+th3)]
    return (p0=p0, p1=p1, p2=p2, p3=p3)
end

function plot_points(positions; line_thickness=5, ball_size=10)
    call_count = 0
    p0 = positions.p0
    p1 = positions.p1
    p2 = positions.p2
    p3 = positions.p3

    line_colors = [:blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :black]
    # Increment the call count and determine the color for this call
    call_count += 1
    current_color = line_colors[mod(call_count, length(line_colors)) + 1]
    
    # Plot the lines
    plot!([p0[1], p1[1], p2[1], p3[1]], [p0[2], p1[2], p2[2], p3[2]], 
        linewidth=line_thickness, color=current_color, label="Call $(call_count)")
    
    # Plot the balls
    scatter!([p0[1]], [p0[2]], color=:black, markersize=ball_size, label=nothing)
    scatter!([p1[1], p2[1], p3[1]], [p1[2], p2[2], p3[2]], color=:red, 
        markersize=ball_size, label=nothing)
    
    display(current())
end