# Simpler GDs for Blind Super-Resolution with Lower Iteration Complexity 

This is the code for "Simpler Gradient Methods for Blind Super-Resolution with Lower Iteration Complexity" by Jinsheng Li, Wei Cui, Xu Zhang, 
in IEEE Transactions on Signal Processing, DOI: 10.1109/TSP.2024.3470071, [arxiv](https://arxiv.org/abs/2403.09031).
## Abstract
We study the problem of blind super-resolution, which can be formulated as a low-rank matrix recovery problem via vectorized Hankel lift (VHL). The previous gradient descent
 method based on VHL named PGD-VHL relies on additional
 regularization such as the projection and balancing penalty,
 exhibiting a suboptimal iteration complexity. In this paper, we
 propose a simpler unconstrained optimization problem without
 the above two types of regularization and develop two new
 and provable gradient methods named VGD-VHL and ScalGD
VHL. A novel and sharp analysis is provided for the theoretical
 guarantees of our algorithms, which demonstrates that our
 methods offer lower iteration complexity than PGD-VHL. In
 addition, ScalGD-VHL has the lowest iteration complexity while
 being independent of the condition number. Furthermore, our
 novel analysis reveals that the blind super-resolution problem is
 less incoherence-demanding, thereby eliminating the necessity for
 incoherent projections to achieve linear convergence. Empirical
 results illustrate that our methods exhibit superior computational
 efficiency while achieving comparable recovery performance to
 prior arts.
## Algorithms
<div align=left> <img src=fig/algs.JPG width="500" height="500" align="center" />  </div>  

## Experiments   
Some experiments in our paper.  
### Phase transitions  
![](fig/pst.png)
###  Runtime comparisons
<div align=center> <img src=fig/Time_comparisons.png width="500" height="250" align="center" />  </div>
 
###  Robust recovery  
 <div align=center> <img src=fig/robust_recovery.png width="500" height="250" align="center" /> </div>


## Citation
If you find this code useful for your research, please consider citing:
```bibtex
@ARTICLE{10474161,
  author={Li, Jinsheng and Cui, Wei and Zhang, Xu},
  journal={IEEE Transactions on Signal Processing}, 
  title={Simpler Gradient Methods for Blind Super-Resolution with Lower Iteration Complexity}, 
  year={2024},
  volume={72},
  number={},
  pages={},
  keywords={—Blind super-resolution, vanilla gradient descent,
 scaled gradient descent, low-rank matrix factorization},
  doi={10.1109/TSP.2024.3470071}}
```
## Code descriptions
### main, 1D case 
demo.m   
demo how to run our simpler gradient methods VGD-VHL, ScalGD-VHL , in 1-D signal case.  
   

solverPgd_fh.m   
PGD-VHL algorithm, 1D case, fast computing via Hankel structure 

solverVgd_fh.m   
VGD-VHL algorithm, 1D case, fast computing via Hankel structure 


solverScaledgd_fh.m  
ScalGD-VHL algorithm, 1D case, fast computing via Hankel structure 

generateSignals_bdft_withsep.m  
generate the true signal and measurements.  

--------------------------------------------------------------------------------
### ./2Dcase  

demo2d.m 
 run a fixed 2D case via our algorithms, and plot the delay-doppler location via 2D MUSIC 

real2Dcase.mat  
true signal matrix; a fixed case. 

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
