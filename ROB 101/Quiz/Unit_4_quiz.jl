using LinearAlgebra

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

function gram_schmidt_nr(U)
  nRowsU, nColsU = size(U)
  V = Array{Float64,2}(undef, nRowsU, 0)
  for k = 1:nColsU 
    uk = U[:, k]
    vk = copy(uk)
      for i = 1:size(V,2)
        vi = V[:,i]
        vk = vk - ( dot(uk, vi)/dot(vi,vi) )*vi 
      end 
      V = [V vk/norm(vk)]
  end
  return V 
end

function forwardsub(L, b)
    n = length(b)
    x = Vector{Float64}(undef, n); 
    x[1] = b[1]/L[1,1] 
    for i = 2:n 
        x[i] = (b[i]- (L[i,1:i-1])'*x[1:i-1] )/L[i,i] 
    end
    return x
end

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

#=
A = [2 1 -3;4 2 -6;1 -1 1]
b = [0; 0; 0]

#display(nullspace(A))
#display(gram_schmidt(A))  # This function gram_schmidt() produce orthogonal columns only

display(nullspace(A'))

A2 = [2 1 -3;1 -1 1]
b2 = [0;0]

display(nullspace(A2'))

Q = gram_schmidt_nr(A2') # This function gram_schmidt_nr() produce orthonormal columns
print("Q:","\n")
display(Q)

R = Q'*A2'
print("R:","\n")
display(R)

beta = forwardsub(R',b2)
print("beta:","\n")
display(beta) =#

A = [1 -3 1;-1 2 -5;5 -13 13]
b = [4; 3; 8]
Q = gram_schmidt(A)
#display(Q)

A2 = [1 -3;-1 2;5 -13]
Q = gram_schmidt_nr(A2)
#display(Q)
R = Q'*A2
#display(R)

x = backwardsub(R,Q'*b)
display(x)

error = b - A2*x
display(error)
display(norm(error)) 
