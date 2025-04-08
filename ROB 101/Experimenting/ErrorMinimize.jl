# Solving no solution
using LinearAlgebra

Phi = [1 1;3 1;7 1]
Y = [2;4;5]

alphaStar = inv(Phi'*Phi) * Phi' *Y
@show alphaStar
@show Phi * alphaStar

error = norm(Phi * alphaStar-Y)