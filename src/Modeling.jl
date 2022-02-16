using LinearAlgebra

struct StateSpace
    A::Matrix{Float64}
    B::Vector{Float64}
    C::Matrix{Float64}
    D::Float64

    function StateSpace(a::Matrix, b::Vector, c::Matrix, d::Float64 = 0.)
	return new(a, b, c, d)
    end

    # Creates state space from float values
    function StateSpace(a::Float64, b::Float64, c::Float64, d::Float64 = 0.)
	return new(fill(a, (1,1)), [b], fill(c, (1, 1)), d)
    end
end

function msd(;m::Float64 = 1., b::Float64 = 1., k::Float64 = 1.)
    A = [
	 0 1;
	 (-k/m) (-b/m)
	 ]
    return StateSpace(A, [1., 0.], [1. 0.])
end

# Maybe add low- and high-pass filters in here so that I can easily define them. 

function low_pass(wc::Float64)
    return StateSpace(-wc, 1., wc)
end

function high_pass(wc::Float64)
    return StateSpace(-wc, 1., -wc, 1.)
end

mutable struct DiscreteSS
    Ad::Matrix{Float64}
    Bd::Vector{Float64}
    C::Matrix{Float64}
    D::Float64

    function DiscreteSS(ss::StateSpace, dt::Float64)
	this = new()
	this.Ad = exp(ss.A * dt)
	if det(ss.A) == 0
	    this.Bd = dt / 2 * (this.Ad + I) * ss.B
	else
	    this.Bd = inv(ss.A) * (this.Ad - I) * ss.B
	end
	this.C = ss.C
	this.D = ss.D
	return this
    end
end

function get_y(ss::StateSpace, x::Vector{Float64}, u::Float64 = 0.)
    return (ss.C * x)[1, 1] + ss.D * u
end

function get_y(dss::DiscreteSS, x::Vector{Float64}, u::Float64 = 0.)
    return (dss.C * x)[1, 1] + dss.D * u
end

function next_x(dss::DiscreteSS, x::Vector{Float64}, u::Float64 = 0.)
    return dss.Ad * x + dss.Bd * u
end

function rk4(f::Function, x::Vector{Float64}, dt::Float64, params, u::Float64 = 0.)
    k1 = f(x, params, u)
    k2 = f(x + dt * k1 / 2, params, u)
    k3 = f(x + dt * k2 / 2, params, u)
    k4 = f(x + dt * k3, params, u)
    return x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
end

function rk4(ss::StateSpace, x::Vector{Float64}, dt::Float64, u::Float64 = 0.)
    k1 = ss.A * x + ss.B * u
    k2 = ss.A * (x + dt * k1 / 2) + ss.B * u
    k3 = ss.A * (x + dt * k2 / 2) + ss.B * u
    k4 = ss.A * (x + dt * k3) + ss.B * u
    return x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
end

function eulers(f::Function, x::Vector{Float64}, dt::Float64, params, u::Float64 = 0.)
    return f(x, params, u) * dt + x
end

function eulers(ss::StateSpace, x::Vector{Float64}, dt::Float64, u::Float64 = 0.)
    return (ss.A * x + ss.B * u) * dt + x
end
