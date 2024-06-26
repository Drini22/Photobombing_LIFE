# Photobombing_LIFE

Contains C++ code (.hpp and .cpp files) used to find the SNR of Q and target confusion value $\mathcal{D}$ from the SNR of P.

Where, 
- SNR = Signal to noise ratio
- P={p_1,p_2} is a configuration of 2 point sources located in the field of view of LIFE. The SNR of P is calculated using LIFEsim. 
- Q is a configuration of 1 point source that minimises the L2 loss L(P,Q). To calculate this quantity one needs to know the SNR of P (given by LIFEsim).

The main steps of the algorithm :
1. Initial point sources -> Start with $Q_1 = \{p_1\}$
2. Linear regression -> The luminosity of $Q$ is updated by calculating integrals. In practice one uses the Simpson formula or similar to approximate integrals to compute $$l_Q = \frac{<\mu_p, \mu_q>}{<\mu_q, \mu_q>} $$
  if $l_Q < 0$ then we set the position of $Q$ to the opposite side and its luminosity to $|l_Q|$.
3. Newton-Raphson method -> The position of $Q$ is updated by computing the gradient $\nabla_Q L$ and the Hessian $H_Q L$ of the loss with respect to $Q$. In practice the analytical expressions are generated and optimised for c++ using sympy library in python \cite{10.7717/peerj-cs.103}.
    $$(\delta_{x,q,i+1} ,\delta_{y,q,i+1}) = (\delta_{x,q,i} ,\delta_{y,q,i}) - H_Q L^{-1}\nabla_Q L $$ 
4. Repeat 2 and 3 alternatively until good enough convergence. -> Since the Newton-Raphson method converges exponentially fast, in practice only 5 steps are enough for most cases.
5. Repeat 1 but start with $Q_2=p_2$. Keep the best of $Q_1,Q_2$ -> We choose between $Q_1$ and $Q_2$ by setting $$Q^{\star} =  \mbox{argmin}_{Q \in \{Q_1,Q_2\}} \mathcal{L}(Q,P)$$
6. Return the value of target confusion map and the SNR of $Q^\star*
    return : $$\mathcal{\mathcal{D}}(P,\sigma) = \min_{P_{s} \in \mathcal{P}(P)-\{P\}} \mathcal{L} (P_{s}, Q^{\star},\sigma) \notag - \mathcal{L}(P,Q^\star,\sigma)$$
    if $\mathcal{D} > 0$, then by definition there is contamination.
