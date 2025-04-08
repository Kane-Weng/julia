include("Bisection_Algorithm.jl")
using IJulia

f(x) = atan(sin(x^3+4)-2.0^(-1-x^2))
y = -0.5
h(x) = f(x)-y  # find the root of h(x)

xmin = 1
xmax = 1.5

F = bisection(h, xmin, xmax, 1e-5)
#@show F
latex_string = "The value of \$f^{-1}(y = 0.5) = $(F.est)\$"
println(latex_string)