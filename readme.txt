demo.m  
demo how to run our simpler gradient methods VGD-VHL, ScalGD-VHL , in 1-D signal case. 
   

solverPgd_fh.m  
PGD-VHL algorithm, 1D case, fast computing via Hankel structure

solverVgd_fh.m  
VGD-VHL algorithm, 1D case, fast computing via Hankel structure


SolverScaledgd_fh.m 
ScalGD-VHL algorithm, 1D case, fast computing via Hankel structure

generateSignals_bdft_withsep.m
generate the true signal and measurements.

--------------------------------------------------------------------------------
./2Dcase
demo2d.m
 run a fxied 2D case via our algorithms, and plot the delay-doppler location via 2D MUSIC

real2Dcase.mat
ture signal matrix; a fixed case.

solverPgd2d.m
PGD-VHL algorithm, 2D case

solverVgd2d.m
VGD-VHL algorithm, 2D case

solverScaledgd2d.m
ScalGD-VHL algorithm, 2D case

music_2d.m
MUSIC for 2D  signal super-resolution

getSignals_ofdm.m
generate 2-D  signal