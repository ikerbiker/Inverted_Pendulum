using LinearAlgebra
using Polynomials
include("Modeling.jl")

mutable struct Controller
    C::Matrix{Float64}
    a::Matrix{Float64}
    A::Matrix{Float64}
    α::Matrix{Float64}
    K::Matrix{Float64}
    kr::Float64

    function get_C(ss::StateSpace)
	C = zeros(size(ss.A)...)
	for i in 1:length(ss.B)
	    C[:,i] = ss.A ^ (i - 1) * ss.B
	end
	if det(C) == 0
	    throw("NOT CONTROLLABLE")
	end
	return C
    end

    function get_a(A::Matrix{Float64})
	return Matrix(reverse(coeffs(fromroots(A)))[2:end]')
    end

    function get_A(a::Matrix{Float64})
	A = Matrix(1.0I, length(a), length(a))
	for i in 1:(length(a) - 1)
	    for j in (i + 1):length(a)
		A[i, j] = a[j - i]
	    end
	end
	return A
    end

    function Controller(ss::StateSpace)
	this = new()
	this.C = get_C(ss)
	this.a = get_a(ss.A)
	this.A = get_A(this.a)
	return this
    end

    function Controller(ss::StateSpace, poles::Vector)
	this = Controller(ss)
	set_poles!(this, ss, poles)
	return this
    end
end

function set_poles!(c::Controller, ss::StateSpace, poles::Vector{Float64})
    c.α = Matrix(reverse(coeffs(fromroots(poles)))[2:end]')
    c.K = (c.α - c.a) * inv(c.A) * inv(c.C)
    c.kr = -1 / (ss.C * inv(ss.A - ss.B * c.K) * ss.B)[1]
end

function set_poles!(c::Controller, ss::StateSpace, poles::Vector{ComplexF64})
    Δ = coeffs(fromroots(poles))
    if all(isreal, Δ)
	c.α = Matrix(reverse(Δ)[2:end]')
	c.K = (c.α - c.a) * inv(c.A) * inv(c.C)
	c.kr = -1 / (ss.C * inv(ss.A - ss.B * c.K) * ss.B)[1]
    else
	throw("COMPLEX POLES MUST BE CONJUGATE PAIRS")
    end
end

function get_u(c::Controller, x::Vector{Float64}, r::Float64)
    return c.kr * r - (c.K * x)[1]
end
