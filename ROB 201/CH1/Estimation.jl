#===================ESTIMATE PI AS ARCHIMEDES========================#
function archimedes_pi(n,radius=0.5)
    # circle with radius 0.5 => diamter = 1
    diameter = 2*radius
    # force n to be an integer
    n = floor(Int, n)
    if n<3
        n=3
    end
    theta = 360/n

    # using trig functions based on degrees (radians are defined through pi)
    circumInnerApprox = (2*radius*sind(theta/2)) * n
    circumOuterApprox = (2*radius*tand(theta/2)) * n

    # estimate pi
    piLowerBound = circumInnerApprox/diameter
    piUpperBound = circumOuterApprox/diameter
    println("$piLowerBound < pi < $piUpperBound for n = $n")
    piApprox = (piLowerBound+piUpperBound)/2
    piErrorBound = (piLowerBound-piUpperBound)/2
    return (piEst=piApprox, piErr=piErrorBound, piL=piLowerBound, piU=piUpperBound)
end

F = archimedes_pi(1000)
#println(F.piEst)   # or F[1]
println(" ")

#====================APPROXIMATING SQRT(2) USING BISECTION ALGORITHM===========================#
include("Bisection_Algorithm.jl")

f(x) = x^2-2
F_2 = bisection(f,1,2,1e-10)
#@show F_2
println("The root is approximately $(F_2.est)")
println("The square of the root is $(F_2.est^2)")
println(" ")

F_3 = bisection2(f,1,2,50)
#@show F_3
println("The root is approximately $(F_3.est)")
println("The square of the root is $(F_3.est^2)")
println(" ")