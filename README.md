# Inverted_Pendulum
Controlling an inverted Pendulum

Uses Plots.jl package for plotting. If you are using WSL to generate the GIFs, you will need to set the graphical output to `NULL`. See the example below

```julia
using Pkg
Pkg.add("Plots")
ENV["GKSwstype"]="nul" # Call this if using WSL (2)
gr()
```

You can run the code by excecuting `julia Inverted_Pendulum.jl`

The `Modeling.jl` and `Controls.jl` files are for providing various modeling techniques such as the Runge-Kutta 4th order method for solving initial value problems, State-Space models for linear systems, and providing full-state feedback control for given State-Space systems

Here is an example output GIF

![model_gif](https://raw.githubusercontent.com/ikerbiker/Inverted_Pendulum/main/gifs/Inverted_pendulum.gif)
