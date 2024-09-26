# Symmetric Hankel projected gradient descent (SHGD)

This is the code for "Projected Gradient Descent for Spectral Compressed Sensing via Symmetric Hankel Factorization" by Jinsheng Li, Wei Cui, Xu Zhang, 
in IEEE Transactions on Signal Processing, doi: 10.1109/TSP.2024.3378004, [arxiv](https://arxiv.org/abs/2403.09031).
## Abstract
Current spectral compressed sensing methods via Hankel matrix completion employ symmetric factorization to demonstrate the low-rank property of the Hankel matrix. However, previous non-convex gradient methods only utilize asymmetric factorization to achieve spectral compressed sensing. In this paper, we propose a novel nonconvex projected gradient descent method for spectral compressed sensing via symmetric factorization named Symmetric Hankel Projected Gradient Descent (SHGD), which updates only one matrix and avoids a balancing regularization term. SHGD reduces about half of the computation and storage costs compared to the prior gradient method based on asymmetric factorization. Besides, the symmetric factorization employed in our work is completely novel to the prior low-rank factorization model, introducing a new factorization ambiguity under complex orthogonal transformation. Novel distance metrics are designed for our factorization method and a linear convergence guarantee to the desired signal is established with $O (r^2\log(n))$ observations. Numerical simulations demonstrate the superior performance of the proposed SHGD method in phase transitions and computation efficiency compared to state-of-the-art methods.
## Algorithms
<div align=left> <img src=fig/Alg.JPG width="500" height="300" align="center" />  </div>  

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
  title={Projected Gradient Descent for Spectral Compressed Sensing via Symmetric Hankel Factorization}, 
  year={2024},
  volume={72},
  number={},
  pages={1590-1606},
  keywords={Symmetric matrices;Compressed sensing;Matrix decomposition;Sparse matrices;Costs;Gradient methods;Convergence;Spectral compressed sensing;Hankel matrix completion;symmetric matrix factorization},
  doi={10.1109/TSP.2024.3378004}}
```
## Code descriptions
