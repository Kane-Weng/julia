function modelParameters()
    L1, L2, L3 = [2, 1.5, 1] # Link lengths
    return (L1=L1, L2=L2, L3=L3)
end

function linkPositions(th1, th2, th3)
    params = modelParameters()
    p0 = [0, 0] # Base position
    p1 = p0 + [params.L1 * cos(th1); params.L1 * sin(th1)] # Position of link 1
    p2 = p1 + [params.L2 * cos(th1 + th2); params.L2 * sin(th1 + th2)] # Position of link 2
    p3 = p2 + [params.L3 * cos(th1 + th2 + th3); params.L3 * sin(th1 + th2 + th3)] # Position of link 3
    return (p0=p0, p1=p1, p2=p2, p3=p3)
end

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

function dyn_mod_3LinkWalker(q, dq)
    # DYN_MOD_3LINKWALKER
    # 2024-07-22 14:48:18
    #
    # Author: Grizzle
    #
    # Model NOTATION: D(q)ddq + C(q,dq)*dq + G(q) = B*tau 
    # The Robot Equations: From Lagrange's Equations of Motion
    #
    g, r, L, m, Mh, Mt = modelParameters3LinkWalker()
    #
    # Variable names for the model
    th1, th2, th3 = q 
    dth1, dth2, dth3 = dq
    #
    D = zeros(3, 3)
      D[1, 1] = r^2*(Mh + Mt + 5*m/4)
      D[1, 2] = -0.5*m*r^2*cos(th1 - th2)
      D[1, 3] = L*Mt*r*cos(th1 - th3)
      D[2, 1] = -0.5*m*r^2*cos(th1 - th2)
      D[2, 2] = 0.25*m*r^2
      D[2, 3] = 0
      D[3, 1] = L*Mt*r*cos(th1 - th3)
      D[3, 2] = 0
      D[3, 3] = L^2*Mt
    #
    C = zeros(3, 3)
      C[1, 1] = 0.0
      C[1, 2] = -0.5*dth2*m*r^2*sin(th1 - th2)
      C[1, 3] = L*Mt*dth3*r*sin(th1 - th3)
      C[2, 1] = 0.5*dth1*m*r^2*sin(th1 - th2)
      C[2, 2] = 0.0
      C[2, 3] = 0.0
      C[3, 1] = -L*Mt*dth1*r*sin(th1 - th3)
      C[3, 2] = 0.0
      C[3, 3] = 0.0
    #
    G = zeros(3)
      G[1] = -Mh*g*r*sin(th1) - Mt*g*r*sin(th1) - (3/2)*g*m*r*sin(th1)
      G[2] = (1/2)*g*m*r*sin(th2)
      G[3] = -L*Mt*g*sin(th3)
    #
    B = zeros(3, 2)
      B[1, 1] = -1
      B[2, 2] = -1
      B[3, 1] = 1
      B[3, 2] = 1
    #
    JacG = zeros(3, 3)
      JacG[1, 1] = -Mh*g*r*cos(th1) - Mt*g*r*cos(th1) - (3/2)*g*m*r*cos(th1)
      JacG[2, 2] = (1/2)*g*m*r*cos(th2)
      JacG[3, 3] = -L*Mt*g*cos(th3)
    #
    return (D=D, C=C, G=G, B=B, JacG=JacG)
end