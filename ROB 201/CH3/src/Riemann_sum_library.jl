function riemannLowerUpperSums(f,a,b,N=100, nRefine=4, aTol=1e-6)
    if abs(a - b) < aTol
        println("abs(a-b) < $aTol. Estimated integral is zero.")
        return 0.0
    end
    if nRefine < 2
        nRefine = 2
        println("nRefine was set to its minimum value of $nRefine")
    end
    N = max(ceil(Int64, N), 2)
    Deltax = (b-a)/N
    lowerSum = 0
    upperSum = 0
    xi = a
    for i = 1:N
        # x_{i+1}
        xip1 = xi + Deltax
        
        # add in extra points so that we can find an estimate of the max and min
        xRefine = range(xi, xip1, length=ceil(Int64, nRefine))
        fRefined = f.(xRefine)
        fMax = maximum(fRefined)
        fMin = minimum(fRefined)
        
        # Now that we have max and min for f, we can compute upper and lower sums
        upperSum = upperSum + fMax * Deltax
        lowerSum = lowerSum + fMin * Deltax
        
        # Update x_i
        xi = xip1 
    end
    # @show b - xi # for debugging
    estIntegral = (upperSum + lowerSum) / 2.0
    pmError = (upperSum - lowerSum) / 2

    return (lowerSum = lowerSum, upperSum = upperSum, estIntegral = estIntegral, pmError = pmError)
end

function trapezoidalRule(f, a, b, N, aTol = 1e-6)
    if abs(a - b) < aTol
        println("abs(a-b) < $aTol. Estimated integral is zero.")
        return 0.0
    end
    
    # Make sure N is an integer greater or equal to 2
    N = max(ceil(Int64, N), 2)
    
    # Define \Delta x
    Deltax = (b - a) / N
    
    # Initialize Sum
    estIntegralTrapezoidal = 0.0
    
    # Initialize x_i
    xi = a
    for i = 1:N
        # Define x_{i+1}
        xip1 = xi + Deltax
        
        # Update integral estimate using trapezoidal rule
        estIntegralTrapezoidal = estIntegralTrapezoidal + Deltax * (f(xi) + f(xip1)) / 2.0
        
        # Update xi
        xi = xip1
    end
    
    return estIntegralTrapezoidal
end

#================================ KEY FUNCTION: SIMPSON RULE ======================================#
function simpsonRuleWithErrorBounds(f, a, b, N = 100, aTol = 1e-6)
    # Uses a centered quadratic
    if abs(a - b) < aTol
        println("abs(a-b) < $aTol. Estimated integral is zero.")
        return 0.0
    end
    
    # Make sure N is an integer greater or equal to 2
    N = max(ceil(Int64, N), 2)
    
    # Define \Delta x
    Deltax = (b - a) / N
    
    # Initialize Sum
    estIntegralSimpson = 0.0
    
    # Initialize total error
    pmError = 0.0
    
    # Initialize x_i
    xi = a
    for i = 1:N
        xip1 = xi + Deltax
        xc = (xi + xip1) / 2
        a1 = f(xi)
        a2 = f(xip1)
        ac = f(xc)
        
        # Simpson's rule estimate
        estIntegralSimpson = estIntegralSimpson + Deltax * (a1 + a2) / 6.0 + 2.0 * Deltax * ac / 3.0
        
        # Refine the intervals to estimate error
        nRefine = 50
        xRefined = range(xi, xip1, length = nRefine)
        fRefined = f.(xRefined)
        
        # Coefficients of the centered quadratic
        alpha = (2 / Deltax^2) * (a1 + a2 - 2 * ac)
        beta = (1 / Deltax) * (a2 - a1)
        gamma = ac
        
        # Define the centered quadratic
        g(x) = alpha * (x - xc)^2 + beta * (x - xc) + gamma
        gRefined = g.(xRefined)
        
        # Compute error
        errorLocal = fRefined - gRefined
        errorChoice = "avg"
        
        if errorChoice == "max"
            maxError = maximum(abs.(fRefined - gRefined))
            pmError = pmError + Deltax * maxError
        elseif errorChoice == "avg"
            avgError = sum(abs.(fRefined - gRefined)) / length(errorLocal)
            pmError = pmError + Deltax * avgError
        else
            meanError = sum(errorLocal) / length(errorLocal)
            pmError = pmError + Deltax * meanError
        end
        
        # Update xi
        xi = xip1
    end
    
    return (estIntegralSimpson = estIntegralSimpson, pmError = pmError)
end