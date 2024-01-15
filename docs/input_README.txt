nx:                   Number of grid points in x (the interval is [0, 2pi]). Positive integer.
ny:                   Number of grid points in y (the interval is [0, 2pi]). Positive integer.
nu_1:                 value of nu_1. Positivie real number.
nu_2:                 value of nu_2. Positivie real number.
dt:                   time step. Positive real number.
Tfinal:               final time of simulation. Positive real number.
cutoff_time:          time at which to cut off the plotting of the energy from below. It is made to remove transient regime. Nonnegative real number.
anim_solution:       whether to plot an animation of the solution (takes considerably more time) (1 = true, 0 = false)
plot_energy:          whether to plot the energy of the solution (1 = true, 0 = false)
anim_freq:    whether to plot the Fourier coefficients of the solution (1 = true, 0 = false)
average_solution:     whether to average the solution over space (1 = true, 0 = false)
freq_anim_sol:       frequency of plotting the animation (every freq_anim_sol-th time step). The higher the value, the less smooth the animation. Positive integer.
freq_plot_e:          frequency of plotting the energy (every freq_plot_e-th time step). The higher the value, the less smooth the energy plot. Positive integer.
freq_anim_freq:       frequency of plotting the Fourier coefficients (every freq_anim_freq-th time step). The higher the value, the less smooth the Fourier coefficients plot. Positive integer.
method:               method to use for solving the PDE. 1 = IMEX-Euler, 2 = IMEX-RK4.