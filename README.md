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

Here is an example output GIF

![](./gifs/Inverted_Pendulum.gif)
