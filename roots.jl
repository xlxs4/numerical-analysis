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

##
function newton_raphson(f, f′, x₀, tol, max_iter=1000)
    x = x₀
    iter = 0
    for i ∈ 1:max_iter
        iter += 1
        x′ = x - f(x) / f′(x)
        abs(x′ - x) < tol && return x′, iter
        x = x′
    end
    throw(DomainError("Failed to converge after $max_iter iterations."))
end
##

##
newton_raphson(f, f′, -1.3, tol)
newton_raphson(f, f′, 0.1, tol)
newton_raphson(f, f′, 1.5, tol)
##

##
function secant(f, x₀, x₁, tol, max_iter=1000)
    iter = 0
    for i ∈ 1:max_iter
        iter += 1
        x₂ = x₁ - f(x₁) * (x₁ - x₀) / (f(x₁) - f(x₀))
        abs(x₂ - x₁) < tol && return x₂, iter
        x₀, x₁ = x₁, x₂
    end
    throw(DomainError("Failed to converge after $max_iter iterations."))
end
##

##
secant(f, -1.5, -1.3, tol)
secant(f, -0.1, 0.1, tol)
secant(f, 1.4, 1.5, tol)
##
