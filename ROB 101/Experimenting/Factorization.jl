using LinearAlgebra
A = [5 3 7 4; 0 3 8 4;5 2 0 2;1 2 3 2]
b = [151/8;16;18;551/40]
F = lu(A)
@show F.L
@show F.U

function forwardsub(L, b)
    n = length(b)
    x = Vector{Float64}(undef, n); 
    x[1] = b[1]/L[1,1] 
    for i = 2:n 
        x[i] = (b[i]- (L[i,1:i-1])'*x[1:i-1] )/L[i,i] 
    end
    return x
end
@show y = forwardsub(F.L,b)

function backwardsub(U, b)

    if minimum(abs.(diag(U))) < 1e-6
        println(" U is nearly singular. ")
        return x = NaN
    end
    
    n = length(b)
    x = Vector{Float64}(undef, n)
    x[n] = b[n] / U[n,n]
    for i = n-1:-1:1
        x[i] = (b[i] - (U[i,(i+1):n])' * x[(i+1):n]) / U[i,i]
    end
    
    return x    
end
@show x = backwardsub(F.U,y)
@show A*x-b

#QR

function gram_schmidt(U)
    nRowsU, nColsU = size(U)
    V = Array{Float64,2}(undef, nRowsU, 0)
    for k = 1:nColsU 
      uk = U[:, k]
      vk = copy(uk)
        for i = 1:size(V,2)
          vi = V[:,i]
          vk = vk - ( dot(uk, vi)/dot(vi,vi) )*vi 
        end 
        V = [V vk]
    end
    return V 
  end
#A2 = [2 1 -3;4 2 -6;1 -1 1]
#b2 = [0;0;0]

#display(nullspace(A2))

#Q = gram_schmidt(A2)
#print("Q:","\n")
#display(Q)


#=print("det(A) is:",det(A2),"\n")

Q = gram_schmidt(A2')
print("Q:","\n")
display(Q)

R = Q'*A2'
print("R:","\n")
display(R)

beta = forwardsub(R',b2)
print("beta:","\n")
display(beta)

xStar = Q*beta 
print("minimum norm solution:","\n")
display(xStar)

#A3 = [1 -3 1;-1 2 -5; 5 -13 13]
#@show Q = gram_schmidt(A3)
#@show R = Q'*A3 =#