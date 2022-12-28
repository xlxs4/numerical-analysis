##
using Plots
##

##
const f(x) = ℯ^sin(x)^3 + x^6 - 2x^4 - x^3 - 1
const f′(x) = 3sin(x)^2 * cos(x) * log(ℯ)ℯ^sin(x)^3 + 6x^5 - 8x^3 - 3x^2
const tol = 1e-5
##

##
plot(f, -2, 2)
##

##
function bisection(f, interval, tol)
    a, b = interval
    fₐ = f(a)
    fₐ * f(b) < 0 || throw(DomainError("The intermediate value theorem has to hold for the given interval"))
    iter = 0
    while (b - a) / 2 > tol
        iter += 1
        m = (a + b) / 2
        f(m) == 0 && return m
        fₐ * f(m) < 0 ? b = m : a = m
    end
    (a + b) / 2, iter
end

function bisection_error_bound(interval, tol)
    a, b = interval
    floor(log((b - a) / tol) / log(2))
end
##

##
bisection(f, (-1.3, -1), tol)
bisection(f, (1.5, 1.55), tol)
##
