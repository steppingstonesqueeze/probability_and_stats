
import streamlit as st
import math
import numpy as np
import matplotlib.pyplot as plt

st.set_page_config(page_title="Bayesian Coin: Odds, LR, and KL", layout="wide")

# ---------- Utilities ----------
def safe_log(x):
    if x <= 0:
        return float("-inf")
    return math.log(x)

def kl_bernoulli(p, q):
    # KL(p || q) for Bernoulli; handle edge cases
    eps = 1e-12
    p = min(max(p, eps), 1 - eps)
    q = min(max(q, eps), 1 - eps)
    return p * math.log(p / q) + (1 - p) * math.log((1 - p) / (1 - q))

def likelihood_binomial(h, t, p):
    # Binomial likelihood up to binomial coefficient (common to both models); 
    # for LR this cancels out, so we can use p^h (1-p)^t directly.
    return (p ** h) * ((1 - p) ** t)

def odds_from_prob(p):
    eps = 1e-12
    p = min(max(p, eps), 1 - eps)
    return p / (1 - p)

def prob_from_odds(o):
    if o == float("inf"):
        return 1.0
    if o <= 0:
        return 0.0
    return o / (1 + o)

def cumulative_arrays(seq):
    # seq is list of 0/1 where 1=heads
    heads = np.cumsum(seq)
    tails = np.arange(1, len(seq) + 1) - heads
    return heads, tails

# ---------- Session State ----------
if "seq" not in st.session_state:
    st.session_state.seq = []  # 1 for H, 0 for T
if "true_p" not in st.session_state:
    st.session_state.true_p = 0.5

# ---------- Sidebar Controls ----------
st.sidebar.header("Model & Prior")
p_fair_prior = st.sidebar.slider("Prior P(Fair)", min_value=0.0, max_value=1.0, value=0.5, step=0.01)
p_biased_prior = 1 - p_fair_prior
theta_alt = st.sidebar.slider("Biased coin parameter θ (H-prob under H₁)", 0.01, 0.99, 0.7, 0.01)

st.sidebar.header("Data Generation")
true_p = st.sidebar.slider("True coin p (for simulation)", 0.0, 1.0, st.session_state.true_p, 0.01)
st.session_state.true_p = true_p
n_to_add = st.sidebar.number_input("Add N random tosses", min_value=1, max_value=1000, value=10, step=1)
seed_opt = st.sidebar.text_input("Random seed (optional)", value="")

col_buttons = st.sidebar.columns(3)
with col_buttons[0]:
    if st.button("Add Head"):
        st.session_state.seq.append(1)
with col_buttons[1]:
    if st.button("Add Tail"):
        st.session_state.seq.append(0)
with col_buttons[2]:
    if st.button("Reset"):
        st.session_state.seq = []

if st.sidebar.button(f"Simulate +{n_to_add}"):
    if seed_opt.strip() != "":
        try:
            np.random.seed(int(seed_opt))
        except:
            pass
    draws = np.random.binomial(1, true_p, size=n_to_add).tolist()
    st.session_state.seq.extend(draws)

# ---------- Main Layout ----------
st.title("Bayesian Coin Tossing: Likelihood Ratio, Posterior Odds, and KL Divergence")

st.markdown("""
This app compares two **simple hypotheses** about a coin:
- **H₀ (Fair):** p = 0.5
- **H₁ (Biased):** p = θ (chosen in sidebar)

You choose a **prior** over the two hypotheses, observe tosses, and watch how the **likelihood ratio (LR)** updates the **posterior odds** in real time.  
We also report **KL divergences** between Bernoulli distributions to connect expected log-likelihood ratios with information gain.
""")

seq = st.session_state.seq
N = len(seq)
h_cum, t_cum = cumulative_arrays(seq) if N > 0 else (np.array([]), np.array([]))

# ---------- Summary Cards ----------
c1, c2, c3, c4 = st.columns(4)
with c1:
    st.metric("Tosses (N)", N)
with c2:
    heads = int(h_cum[-1]) if N > 0 else 0
    st.metric("Heads (H)", heads)
with c3:
    tails = int(t_cum[-1]) if N > 0 else 0
    st.metric("Tails (T)", tails)
with c4:
    f_hat = heads / N if N > 0 else 0.0
    st.metric("Empirical p̂ (Heads)", f_hat)

# ---------- Compute LR / Odds evolution ----------
# Prior odds: P(F)/P(B)
prior_odds = odds_from_prob(p_fair_prior)

def lr_at(h, t, theta):
    # LR = L(H1)/L(H0) = [theta^h (1-theta)^t] / [0.5^h (0.5)^t] = [theta^h (1-theta)^t] / [0.5^N]
    if h + t == 0:
        return 1.0
    num = likelihood_binomial(h, t, theta)
    den = likelihood_binomial(h, t, 0.5)
    if den == 0:
        return float("inf")
    return num / den

lrs = []
post_odds_list = []
post_fair_list = []
log_lr_list = []

if N > 0:
    for i in range(N):
        h_i = int(h_cum[i])
        t_i = int(t_cum[i])
        lr_i = lr_at(h_i, t_i, theta_alt)
        lrs.append(lr_i)
        log_lr_list.append(safe_log(lr_i))
        post_odds_i = prior_odds * lr_i  # Posterior odds(F/B) = Prior odds * LR
        post_odds_list.append(post_odds_i)
        post_fair_list.append(prob_from_odds(post_odds_i))
else:
    lrs = []
    post_odds_list = []
    post_fair_list = []
    log_lr_list = []

# Final posterior probabilities
final_post_odds = post_odds_list[-1] if len(post_odds_list) else prior_odds
post_fair = prob_from_odds(final_post_odds)
post_biased = 1 - post_fair

st.subheader("Posterior after N tosses")
colA, colB, colC = st.columns(3)
with colA:
    st.write(f"**Prior P(Fair)** = {p_fair_prior:0.4f}")
    st.write(f"**Prior P(Biased)** = {p_biased_prior:0.4f}")
with colB:
    st.write(f"**Posterior P(Fair | data)** = {post_fair:0.6f}")
with colC:
    st.write(f"**Posterior P(Biased | data)** = {post_biased:0.6f}")

# ---------- Plots ----------
if N > 0:
    # Likelihood ratio (cumulative)
    st.subheader("Evolution: Likelihood Ratio, Log-LR, and Posterior P(Fair)")
    fig1, ax1 = plt.subplots()
    ax1.plot(range(1, N+1), lrs)
    ax1.set_xlabel("Toss index")
    ax1.set_ylabel("Likelihood Ratio  L(H₁)/L(H₀)")
    ax1.set_title("Cumulative Likelihood Ratio")
    st.pyplot(fig1)

    fig2, ax2 = plt.subplots()
    ax2.plot(range(1, N+1), log_lr_list)
    ax2.set_xlabel("Toss index")
    ax2.set_ylabel("log LR")
    ax2.set_title("Cumulative log-Likelihood Ratio")
    st.pyplot(fig2)

    fig3, ax3 = plt.subplots()
    ax3.plot(range(1, N+1), post_fair_list)
    ax3.set_xlabel("Toss index")
    ax3.set_ylabel("Posterior P(Fair)")
    ax3.set_ylim(0, 1)
    ax3.set_title("Posterior Probability of Fair Coin")
    st.pyplot(fig3)

# ---------- KL Divergences ----------
st.subheader("KL Divergence Connections")

col1, col2 = st.columns(2)

with col1:
    st.markdown("**Plug-in (Empirical) KLs using p̂**")
    if N == 0:
        st.info("Add some tosses to compute plug-in KL.")
    else:
        kl_hat_vs_fair = kl_bernoulli(f_hat, 0.5)
        kl_hat_vs_alt = kl_bernoulli(f_hat, theta_alt)
        st.write(f"KL( Bern(p̂) || Bern(0.5) ) = **{kl_hat_vs_fair:0.6f}**")
        st.write(f"KL( Bern(p̂) || Bern(θ)) with θ={theta_alt:0.2f} = **{kl_hat_vs_alt:0.6f}**")
        st.caption("Interpretation: With p̂ as a stand-in for the unknown truth, these KLs approximate expected log-likelihood ratios against each hypothesis.")

with col2:
    st.markdown("**Model-based KLs (if one hypothesis were the truth)**")
    kl_alt_vs_fair = kl_bernoulli(theta_alt, 0.5)
    kl_fair_vs_alt = kl_bernoulli(0.5, theta_alt)
    st.write(f"KL( Bern(θ) || Bern(0.5) ) with θ={theta_alt:0.2f} = **{kl_alt_vs_fair:0.6f}**")
    st.write(f"KL( Bern(0.5) || Bern(θ) ) with θ={theta_alt:0.2f} = **{kl_fair_vs_alt:0.6f}**")
    st.caption("Interpretation: If H₁ (θ) were true, the expected log-LR favoring H₁ over fair is KL(θ || 0.5); vice versa for the other direction.")

st.markdown("""---""")
st.markdown("""
**Notes**
- **Posterior odds** = **Prior odds** × **Likelihood ratio** (Bayes’ rule for two simple hypotheses).  
- **Expected log-LR under the true model** equals the **KL divergence** to the false model.  
- The plug-in KL with p̂ is a data-driven approximation when the truth is unknown.
""")

st.markdown("Made with ❤️ for hands-on Bayesian intuition.")
