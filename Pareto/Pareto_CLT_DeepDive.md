# CLT and the Pareto Distribution: Deep Dive

**Classical CLT:**  
If X₁, X₂, … are i.i.d. with mean μ and finite variance σ², then:  
Zₙ = √n (X̄ₙ – μ)/σ → 𝒩(0,1) in distribution.

**Pareto heavy tails:**  
Pareto(α, xₘ) has tail P(X > x) ∼ (xₘ/x)^α.  
- Mean finite iff α > 1, μ = αxₘ/(α–1).  
- Variance finite iff α > 2, σ² = αxₘ² / ((α–1)² (α–2)).

**Failure mode:**  
For 1 < α ≤ 2, the mean exists but variance is infinite → classical CLT fails.  
The generalized CLT applies: properly scaled sums converge to an α-stable law (non-Gaussian).  
Studentizing by sample sd does not restore normality: t-statistics remain heavy-tailed.

**Practical consequences:**  
- t-intervals for the mean under-cover severely when α ≤ 2.  
- Robust estimators (trimmed means, Winsorized means, Huber M-estimators) regain asymptotic normality.

**Takeaways for demonstrations:**  
1. α = 2.5 (finite variance): LLN + CLT work as expected → sample mean normalized looks Gaussian.  
2. α = 1.5 (finite mean, infinite variance): LLN holds but CLT fails → histograms remain heavy-tailed, Q–Q vs Normal bends, coverage of 95% CIs is poor.  
3. Trimmed means dramatically stabilize inference in heavy tails.

This contrast makes Pareto the textbook demonstration of where LLN still works but CLT fails spectacularly.
