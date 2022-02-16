include("src/Modeling.jl")
include("src/Controls.jl")
using Plots
ENV["GKSwstype"]="nul" # Call this if using WSL


mutable struct Params
    mc::Float64
    mp::Float64
    g::Float64
    l::Float64
    function Params()
        return new(10., 1., 10., 2.)
    end
end

# Differential equation used in model of inverted pendulum
function xdot(x::Vector{Float64}, params::Params, force::Float64)
    return [
        x[2],
        (-params.mp * params.g * sin(x[3]) * cos(x[3]) + params.mp * params.l * sin(x[3]) * x[4] ^ 2 + force) / 
        (params.mc + params.mp - params.mp * sin(x[3])^2),
        x[4],
        (-params.g * (params.mc + params.mp) * sin(x[3]) + params.mp * params.l * sin(x[3]) * cos(x[3]) * x[4] ^ 2 + force * cos(x[3])) / 
        (-params.l * (params.mc + params.mp) + params.mp * params.l * cos(x[3]) ^ 2)
    ]
end

params = Params()

# Linearized State-Space model for making the controller
A = [
    0 1 0 0;
    0 0 (params.mp * params.g / params.mc) 0;
    0 0 0 1;
    0 0 (params.g * (params.mc + params.mp) / params.mc / params.l) 0
]

B = [0; 1 / params.mc; 0; -1 / params.mc / params.l]

C = [1 0 0 0]

ss = StateSpace(A, B, C)
controller = Controller(ss, [-1, -1.1, -1.2, -1.3])

# Setting up time constraints
dt = 0.01
t = 0:dt:20
x = [[0., 0., 0., 0.]] #Initial conditions of pendulum

# Iterating through all the time steps
for i in t[2:end]
    u = get_u(controller, last(x), 0.0)
    push!(x, rk4(xdot, last(x), dt, params, u))
    if i == 1.0
        last(x)[4] = -2.0 # Instantaneous change in angular velocity, simulating flicking the ball
    end
end

# Making GIF

anim = @animate for i in 1:length(t)
    plot([x[i][1], x[i][1] + params.l * sin(x[i][3])], 
        [0, params.l * cos(x[i][3])], 
        leg = false, 
        xlims = (-10, 1), 
        ylims = (-1, 2.5), 
        aspect_ratio = 1)
    scatter!([x[i][1]], [0], m = (:rect, 12))
    scatter!([x[i][1] + params.l * sin(x[i][3])], [params.l * cos(x[i][3])], m = (:circle, 8))
end every 10;
gif(anim, "gifs/Inverted_pendulum.gif");
