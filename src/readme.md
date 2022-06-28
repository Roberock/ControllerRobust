## Functions for data analysis and to run the solver

```
# run monte carlo analysis on the system model g
MonteCarloController 
```

```
https://github.com/Roberock/ControllerRobust/edit/master/src/EstimateReliabilityScores.jl
```

```
g_controller

    # INPUTS 
    # design - the controller gains design = [coeffs numerator ceffs denominator]
    # θ - system parameters (uncertain parameters) 
    # tvec - discretized time
    # do_linear - [true,false]
    
    # OUTPUTS: the reliabiliy performance funtions of the controller
    # g[1] : hurwitz stability
    # g[2] : control effort 
    # g[3] : settling time 
    
```
