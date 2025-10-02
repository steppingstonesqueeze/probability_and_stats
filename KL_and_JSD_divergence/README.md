
# Bayesian Coin Tossing: Likelihood Ratio, Posterior Odds, and KL

This Streamlit app lets you:
- Set **prior probability** that a coin is **fair** vs **biased**.
- Choose the **biased coin parameter** θ (the probability of heads under the biased hypothesis).
- Add heads/tails manually or simulate random tosses from a **true p**.
- Watch how the **likelihood ratio (LR)** and **posterior odds** evolve after each toss.
- Inspect **KL divergences** to connect expected log-LR with information gain.

## Quickstart

```bash
# 1) Create a fresh environment (optional but recommended)
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

# 2) Install requirements
pip install -r requirements.txt

# 3) Run the app
streamlit run app.py
```

Then open the local URL Streamlit prints (usually http://localhost:8501).

## How it works

Two **simple hypotheses**:
- H₀ (Fair): p = 0.5
- H₁ (Biased): p = θ (chosen in the sidebar)

Given H heads and T tails after N = H + T tosses, the **likelihood ratio** is:
LR = [θ^H (1-θ)^T] / [0.5^N].

**Posterior odds (Fair vs Biased)** = **Prior odds** × LR⁻¹ (or equivalently, using LR as L(H₁)/L(H₀) we track posterior odds of Fair vs Biased via its reciprocal).  
In the app, we directly compute posterior **P(Fair | data)** from odds.

**KL Divergence** for Bernoulli distributions connects to **expected log-LR**:
- If the true coin has parameter p*, then E[log LR in favor of that true model] = KL(Bern(p*) || Bern(other)).
- We show both **model-based KLs** (assuming H₀ or H₁ is true) and a **plug-in KL** using the empirical p̂ = H/N.
