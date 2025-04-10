# Define Linear Drag Model
struct Params
    m::Float64
    kappa::Float64
    g::Float64
end

function parallel_model(option)
    # Initialize parameters
    m = 78e3
    g = 9.81
    v0 = 120
    
    if (option=="linear")
        kappa = (0.9*v0*m)/3000; @show kappa
    elseif (option=="nonlinear")
        kappa = (m/3000)*log(10); @show kappa
    end

    params = Params(m, kappa, g)

    # Define the ODE function
    function f(v, params, t)
        if (option=="linear")
            dv = -(params.kappa / params.m) * v 
        elseif (option=="nonlinear")
            dv = -(params.kappa / params.m) * v.^2
        end
        return dv
    end

    # Set initial condition as a vector
    v0 = [v0] 

    # Set time span
    tspan = (0.0, 100.0)

    # Create the ODE problem
    prob = ODEProblem{false}(f, v0, tspan, params)

    # Solve the ODE problem
    sol = solve(prob, Tsit5())
    if (option=="linear")
        Tf = -(m/kappa)*log(0.1)
        return Tf, sol
    elseif (option=="nonlinear")
        Tf = 9*m/(kappa*v0[1])
        return Tf, sol
    end 
end

function vertical_model(option)
    # Initialize the Params structure with specific values
    m = 75 # kg
    g = 9.81 # m/s^2
    vTerminal = -2 # m/s
    if (option=="linear")
        kappa = -m*g/vTerminal; @show kappa
    elseif (option=="nonlinear")
        kappa = -m*g/(vTerminal^2); @show kappa
    end 
    
    params = Params(m, kappa, g)

    # define the ODE dv/dt = f(v,params,t)
    function vertical_linear(v, params, t)
        if (option=="linear")
            dv = -params.g .-(params.kappa / params.m) * v
        elseif (option=="nonlinear")
            dv = -params.g .-(params.kappa / params.m) * v.^2
        end 
        return dv
    end

    # Set the initial condition as a vector
    v0 = -25 # m/s speed when chute is opened
    v0= [v0]

    # Set the time interval
    T = (0.0, 1.5) # 

    # Setup the ODE problem with out-of-place function
    problem = ODEProblem{false}(vertical_linear, v0, T, params)

    # solve the ODE problem using the Runge Kutta Tsitouras 5/4 Integrator
    sol = solve(problem, Tsit5());
    return sol
end