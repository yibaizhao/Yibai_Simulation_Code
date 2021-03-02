# Yibai_Simulation_Code

This repo shows Yibai's simulation code based on R.

The goal of this simulation is to test the robustness of our proposed model under various settings.

## Simulation I
In this setting, time to event $T$ follows Weibull PH model and the probability of time to event at birth($\tau_0$) is a function of covariate $\mathbf z$.

## Simulation II
In this scenario, we keep everything the same as Simulation I but relax the proportional hazards assumption in the Weibull regression model, so the baseline hazard depends on covariate $\mathbf{z}$. For strata $j$, our mixture model can be expressed as
\begin{equation*}
	P(T\leq t; \mathbf{z} = j, \mathbf\theta) = \pi(\mathbf{z};\mathbf\theta_1) + (1 - \pi(\mathbf{z};\mathbf\theta_1))[1 - S(t;\mathbf\theta_2 | \mathbf{c}=0, \mathbf{z} = j)]
\end{equation*},
where $S(t;\mathbf\theta_2 | \mathbf{c}=0, \mathbf{z} = j) = \exp(-\lambda_j\exp(\mathbf{x}^T\mathbf\beta)t^{\gamma_j})$.