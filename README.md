# Inverted_Pendulum
Controlling an inverted Pendulum

Uses Plots.jl package for plotting, along with custom packages I have created entitled Modeling.jl and Controls.jl. If you are using WSL to generate the GIFs, you will need to set the graphical output to `NULL`. See the example below

```julia
using Pkg
Pkg.add("Plots")
ENV["GKSwstype"]="nul" # Call this if using WSL (2)
gr()
```

You can run the code by excecuting `julia --project=Deps Inverted_Pendulum.jl`, although you may need to instantiate the REPL and download the dependencies by using `julia --project=Deps -e 'using Pkg; Pkg.instantiate()'`

The `Modeling.jl` and `Controls.jl` packages are for providing various modeling techniques such as the Runge-Kutta 4th order method for solving initial value problems, State-Space models for linear systems, and providing full-state feedback control for given State-Space systems

Here is an example output GIF

![model_gif](https://raw.githubusercontent.com/smith-isaac/Inverted_Pendulum/main/gifs/Inverted_pendulum.gif)
