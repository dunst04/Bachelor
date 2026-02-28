
function bisection(f, a, b, xtol = 1e-8, ftol = 1e-8, maxiter = 100)

    @assert f(a) * f(b) < 0 "f(a) and f(b) must have opposite signs"


    for i in 1:maxiter

        c = (a + b) / 2
        fc = f(c)

        if abs(b-a)<xtol || abs(fc)<ftol
            return (c, fc, i)
        end

        if f(a) * fc < 0
            b = c
        else
            a = c
        end

    end
    c = (a + b) / 2
    println("maximum number of iterations has been achieved")
    return (c, f(c), maxiter)
end

g(x) = x^2 - 2
x̂ = bisection(g, 0.0, 2.0)[1]
abs(x̂ - sqrt(2))