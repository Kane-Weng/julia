function bisection(f,a,b,tol)
    # Both a and b are above/below x-axis
    if f(a) * f(b) > 0
        error("f(a) and f(b) must have opposite signs")
        return Na
    end

    kmax = 1e5
    k=1
    c=NaN
    while ((b-a)/2>tol || abs(f(c))>tol) && k<kmax
        c = (a+b)/2
        if f(a)*f(c) <0
            b=c
        else
            a=c
        end
        k = k+1
    end
    root_ErrorBound = (b-a)/2
    return (est=(a+b)/2, low=a, up=b, error=root_ErrorBound)
end

function bisection2(f,a,b,n)
    # Both a and b are above/below x-axis
    if f(a) * f(b) > 0
        error("f(a) and f(b) must have opposite signs")
        return NaN
    end
    c=NaN
    for k=1:n
        c = (a+b)/2
        if f(a)*f(c) <0
            b=c
        else
            a=c
        end
        k = k+1
    end
    root_ErrorBound = (b-a)/2
    return (est=(a+b)/2, low=a, up=b, error=root_ErrorBound)
end