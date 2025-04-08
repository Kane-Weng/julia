#=================BERNOULLI'S DEFINITION OF EULER'S CONSTANT=======================#
function bernoulli_e(n)
    n=floor(Int,n)
    if n<2
        n=2
    end
    e_LowerApprox = (1.0+1.0/n)^n 
    e_UpperApprox = 1.0
    for k = 1:n
        e_UpperApprox = e_UpperApprox + 1/factorial(big(k))
    end
    e_UpperApprox = e_UpperApprox + 1/(n*factorial(big(n)))
    e_Approx = Float64((e_LowerApprox+e_UpperApprox)/2)
    e_ErrorBound = Float64((e_LowerApprox-e_UpperApprox)/2)
    println("$e_LowerApprox < e < $e_UpperApprox for n = $n")
    return (est=e_Approx, low=e_LowerApprox, up=e_UpperApprox, error=e_ErrorBound)
end

F = bernoulli_e(52*7*24)
