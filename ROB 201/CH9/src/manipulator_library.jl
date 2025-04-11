#============================ SIMPLIFYING EXPRESSION =========================#
function simple(expr)
    # If the expression is a matrix, handle each element separately
    if expr isa Matrix
        rows, cols = size(expr)
        new_mat = Matrix{Any}(undef, rows, cols)
        for i in 1:rows
            for j in 1:cols
                new_mat[i, j] = simple(expr[i, j])
            end
        end
        return new_mat
    end

    # Convert Symbolics.jl expression to a string
    expr_str = string(expr)

    # Replace implicit multiplication (like "0.5(") with explicit multiplication ("0.5*(")
    # expr_str = replace(expr_str, r"([\d\.]+)\(" => s"\1*(")
    expr_str = replace(expr_str, r"([\d\.]+)([a-zA-Z\(])" => s"\1*\2")
    
    # Remove unwanted artifacts (e.g., "Num")
    expr_str = replace(expr_str, "Num[" => "[")
    expr_str = replace(expr_str, "Symbol (" => "")
    expr_str = replace(expr_str, "Integer (" => "")
    # Replace double divide
    expr_str = replace(expr_str, "//" => "/")
 
    # Convert the modified string to SymPy.jl expression
    sympy_expr = SymPy.sympify(expr_str)
 
    # Use SymPy to simplify the expression
    simplified_expr = SymPy.simplify(sympy_expr)
 
    # Convert simplified expression back to string
    simplified_str = string(simplified_expr)
 
    # Remove 1.0* from the simplified expression string
    final = replace(simplified_str, r"(?<!\d)1\.0\*(?![\d\.])" => "")
 
    return final
end

#============================ PROF'S MOD =========================#
function dyn_mod_2LinkManipulator(q, dq)
    # DYN_MOD_2LINKMANIPULATOR
    # 2024-07-22 14:44:12
    #
    # Author: Grizzle
    #
    # Model NOTATION: D(q)ddq + C(q,dq)*dq + G(q) = B*tau 
    # The Robot Equations: From Lagrange's Equations of Motion
    #
    g, L1, L2, m1, m2 = modelParameters()
    #
    # Variable names for the model
    q1, q2 = q 
    dq1, dq2 = dq
    #
    D = zeros(2, 2)
      D[1, 1] = m1*(L1^2) + 0.5m2*((L1*cos(q1) + L2*cos(q1 + q2))*(2L1*cos(q1) +
                 2L2*cos(q1 + q2)) + (-L1*sin(q1) - L2*sin(q1 + q2))*(-2L1*
                sin(q1) - 2L2*sin(q1 + q2)))
      D[1, 2] = 0.5m2*(2L2*(L1*cos(q1) + L2*cos(q1 + q2))*cos(q1 + q2) - 2L2*(-
                L1*sin(q1) - L2*sin(q1 + q2))*sin(q1 + q2))
      D[2, 1] = 0.5m2*(2L2*(L1*cos(q1) + L2*cos(q1 + q2))*cos(q1 + q2) - 2L2*(-
                L1*sin(q1) - L2*sin(q1 + q2))*sin(q1 + q2))
      D[2, 2] = 0.5m2*(L2^2)*(2(cos(q1 + q2)^2) + 2(sin(q1 + q2)^2))
    #
    C = zeros(2, 2)
      C[1, 1] = 0.25dq2*m2*(-2L2*(L1*cos(q1) + L2*cos(q1 + q2))*sin(q1 + q2) -
                 2L2*(-L1*sin(q1) - L2*sin(q1 + q2))*cos(q1 + q2) - L2*(2L1*
                cos(q1) + 2L2*cos(q1 + q2))*sin(q1 + q2) - L2*(-2L1*sin(q1) -
                 2L2*sin(q1 + q2))*cos(q1 + q2))
      C[1, 2] = 0.25dq1*m2*(-2L2*(L1*cos(q1) + L2*cos(q1 + q2))*sin(q1 + q2) -
                 2L2*(-L1*sin(q1) - L2*sin(q1 + q2))*cos(q1 + q2) - L2*(2L1*
                cos(q1) + 2L2*cos(q1 + q2))*sin(q1 + q2) - L2*(-2L1*sin(q1) -
                 2L2*sin(q1 + q2))*cos(q1 + q2)) + 0.5dq2*m2*(-2L2*(L1*cos(q1) +
                 L2*cos(q1 + q2))*sin(q1 + q2) - 2L2*(-L1*sin(q1) - L2*sin(q1 +
                 q2))*cos(q1 + q2))
      C[2, 1] = dq1*(0.5m2*(-2L2*(L1*cos(q1) + L2*cos(q1 + q2))*sin(q1 + q2) -
                 2L2*(-L1*cos(q1) - L2*cos(q1 + q2))*sin(q1 + q2)) - 0.25m2*(-
                2L2*(L1*cos(q1) + L2*cos(q1 + q2))*sin(q1 + q2) - 2L2*(-L1*
                sin(q1) - L2*sin(q1 + q2))*cos(q1 + q2) - L2*(2L1*cos(q1) + 2L2*
                cos(q1 + q2))*sin(q1 + q2) - L2*(-2L1*sin(q1) - 2L2*sin(q1 +
                 q2))*cos(q1 + q2)))
    #
    G = zeros(2)
      G[1] = g*m2*(L1*cos(q1) + L2*cos(q1 + q2)) + L1*g*m1*cos(q1)
      G[2] = L2*g*m2*cos(q1 + q2)
    #
    B = zeros(2, 2)
      B[1, 1] = 1
      B[2, 2] = 1
    #
    JacG = zeros(2, 2)
      JacG[1, 1] = g*m2*(-L1*sin(q1) - L2*sin(q1 + q2)) - L1*g*m1*sin(q1)
      JacG[1, 2] = -L2*g*m2*sin(q1 + q2)
      JacG[2, 1] = -L2*g*m2*sin(q1 + q2)
      JacG[2, 2] = -L2*g*m2*sin(q1 + q2)
    #
    return (D=D, C=C, G=G, B=B, JacG=JacG)
end
    
function dyn_mod_3LinkManipulator(q, dq)
    # DYN_MOD_3LINKMANIPULATOR
    # 2025-04-10 00:55:46
    #
    # Author: Grizzle
    #
    # Model NOTATION: D(q)ddq + C(q,dq)*dq + G(q) = B*tau 
    # The Robot Equations: From Lagrange's Equations of Motion
    #
    g, L1, L2, L3, m1, m2, m3 = modelParameters()
    #
    # Variable names for the model
    q1, q2, q3 = q 
    dq1, dq2, dq3 = dq
    #
    D = zeros(3, 3)
        D[1, 1] = L1^2*m1 + m2*(L1^2 + 2*L1*L2*cos(q2) + L2^2) + m3*(L1^2 + 2*L1*
                L2*cos(q2) + 2*L1*L3*cos(q2 + q3) + L2^2 + 2*L2*L3*cos(q3) +
                    L3^2)
        D[1, 2] = L2*m2*(L1*cos(q2) + L2) + m3*(L1*L2*cos(q2) + L1*L3*cos(q2 +
                    q3) + L2^2 + 2*L2*L3*cos(q3) + L3^2)
        D[1, 3] = L3*m3*(L1*cos(q2 + q3) + L2*cos(q3) + L3)
        D[2, 1] = L2*m2*(L1*cos(q2) + L2) + m3*(L1*L2*cos(q2) + L1*L3*cos(q2 +
                    q3) + L2^2 + 2*L2*L3*cos(q3) + L3^2)
        D[2, 2] = L2^2*m2 + m3*(L2^2 + 2*L2*L3*cos(q3) + L3^2)
        D[2, 3] = L3*m3*(L2*cos(q3) + L3)
        D[3, 1] = L3*m3*(L1*cos(q2 + q3) + L2*cos(q3) + L3)
        D[3, 2] = L3*m3*(L2*cos(q3) + L3)
        D[3, 3] = L3^2*m3
    #
    C = zeros(3, 3)
        C[1, 1] = -L1*L2*dq2*m2*sin(q2) - L1*L2*dq2*m3*sin(q2) - L1*L3*dq2*m3*
                sin(q2 + q3) - L1*L3*dq3*m3*sin(q2 + q3) - L2*L3*dq3*m3*sin(q3)
        C[1, 2] = -L1*L2*dq1*m2*sin(q2) - L1*L2*dq1*m3*sin(q2) - L1*L2*dq2*m2*
                sin(q2) - L1*L2*dq2*m3*sin(q2) - L1*L3*dq1*m3*sin(q2 + q3) - L1*
                L3*dq2*m3*sin(q2 + q3) - L1*L3*dq3*m3*sin(q2 + q3) - L2*L3*dq3*
                m3*sin(q3)
        C[1, 3] = -L3*m3*(L1*sin(q2 + q3) + L2*sin(q3))*(dq1 + dq2 + dq3)
        C[2, 1] = L1*L2*dq1*m2*sin(q2) + L1*L2*dq1*m3*sin(q2) + L1*L3*dq1*m3*
                sin(q2 + q3) - L2*L3*dq3*m3*sin(q3)
        C[2, 2] = -L2*L3*dq3*m3*sin(q3)
        C[2, 3] = -L2*L3*m3*(dq1 + dq2 + dq3)*sin(q3)
        C[3, 1] = L3*m3*(L1*dq1*sin(q2 + q3) + L2*dq1*sin(q3) + L2*dq2*sin(q3))
        C[3, 2] = L2*L3*m3*(dq1 + dq2)*sin(q3)
        C[3, 3] = 0.0
    #
    G = zeros(3)
        G[1] = g*m2*(L1*cos(q1) + L2*cos(q1 + q2)) + g*m3*(L1*cos(q1) + L2*
                cos(q1 + q2) + L3*cos(q1 + q2 + q3)) + L1*g*m1*cos(q1)
        G[2] = g*m3*(L2*cos(q1 + q2) + L3*cos(q1 + q2 + q3)) + L2*g*m2*cos(q1 +
                    q2)
        G[3] = L3*g*m3*cos(q1 + q2 + q3)
    #
    B = zeros(3, 3)
        B[1, 1] = 1
        B[2, 2] = 1
        B[3, 3] = 1
    #
    JacG = zeros(3, 3)
        JacG[1, 1] = g*m2*(-L1*sin(q1) - L2*sin(q1 + q2)) + g*m3*(-L1*sin(q1) - L2*
                sin(q1 + q2) - L3*sin(q1 + q2 + q3)) - L1*g*m1*sin(q1)
        JacG[1, 2] = g*m3*(-L2*sin(q1 + q2) - L3*sin(q1 + q2 + q3)) - L2*g*m2*sin(q1 +
                    q2)
        JacG[1, 3] = -L3*g*m3*sin(q1 + q2 + q3)
        JacG[2, 1] = g*m3*(-L2*sin(q1 + q2) - L3*sin(q1 + q2 + q3)) - L2*g*m2*sin(q1 +
                    q2)
        JacG[2, 2] = g*m3*(-L2*sin(q1 + q2) - L3*sin(q1 + q2 + q3)) - L2*g*m2*sin(q1 +
                    q2)
        JacG[2, 3] = -L3*g*m3*sin(q1 + q2 + q3)
        JacG[3, 1] = -L3*g*m3*sin(q1 + q2 + q3)
        JacG[3, 2] = -L3*g*m3*sin(q1 + q2 + q3)
        JacG[3, 3] = -L3*g*m3*sin(q1 + q2 + q3)
    #
    return (D=D, C=C, G=G, B=B, JacG=JacG)
end

function dyn_mod_2LinkManipulatorWithSprings(q, dq)
  # DYN_MOD_2LINKMANIPULATORWITHSPRINGS
  # 2024-07-22 14:52:28
  #
  # Author: Grizzle
  #
  # Model NOTATION: D(q)ddq + C(q,dq)*dq + G(q) = B*tau 
  # The Robot Equations: From Lagrange's Equations of Motion
  #
   
  #
  # Variable names for the model
  q1, q2 = q 
  dq1, dq2 = dq
  #
  D = zeros(2, 2)
    D[1, 1] = 3.0*cos(q2) + 4.75
    D[1, 2] = 1.5*cos(q2) + 0.75
    D[2, 1] = 1.5*cos(q2) + 0.75
    D[2, 2] = 0.750000000000000
  #
  C = zeros(2, 2)
    C[1, 1] = -1.5*dq2*sin(q2)
    C[1, 2] = -1.5*(dq1 + dq2)*sin(q2)
    C[2, 1] = 1.5*dq1*sin(q2)
    C[2, 2] = 0.0
  #
  G = zeros(2)
    G[1] = 14.715cos(q1 + q2) + 24.0q1 + 39.24cos(q1)
    G[2] = 14.715cos(q1 + q2) + 24.0q2
  #
  B = zeros(2, 2)
    B[1, 1] = 1
    B[2, 2] = 1
  #
  JacG = zeros(2, 2)
    JacG[1, 1] = 24.0 - 14.715sin(q1 + q2) - 39.24sin(q1)
    JacG[1, 2] = -14.715sin(q1 + q2)
    JacG[2, 1] = -14.715sin(q1 + q2)
    JacG[2, 2] = 24.0 - 14.715sin(q1 + q2)
  #
  return (D=D, C=C, G=G, B=B, JacG=JacG)
end

function cleanUp(A, tol=1e-10)
  B = copy(A)
  for i in eachindex(B)
      if isa(B[i], Complex)
          # Clean up the real part
          real_part = abs(real(B[i])) < tol ? 0.0 : real(B[i])
          # Clean up the imaginary part
          imag_part = abs(imag(B[i])) < tol ? 0.0 : imag(B[i])
          # Reconstruct the complex number
          B[i] = real_part + imag_part * im
      else
          # Original cleanup for non-complex numbers
          if abs(B[i]) < tol
              B[i] = 0.0
          end
      end
  end
  return B
end