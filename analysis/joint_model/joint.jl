# Joint generation-interval model using all traced pairs.
#
#   incubation : every infectee's onset - infection window
#   chain GIs  : infector & infectee infection windows (direct GI)
#   onset GIs  : infector onset + infectee infection -> GI = (E_j - O_i) + incubation,
#                with the infector's incubation integrated out and the incubation
#                distribution shared with the incubation stream.
#
# All likelihoods integrate the *pdf* by quadrature, which keeps everything
# autodiff-friendly (the gamma CDF's shape derivative is unsupported).
# Forward (phase) correction and multiple-infector marginalisation: later layers.

using CensoredDistributions, Distributions, Turing, Random, ForwardDiff
include(joinpath(@__DIR__, "data.jl"))

d = load_pairs()

I = d.incubation                                   # incubation stream
inc_pw = I.inf_hi .- I.inf_lo
inc_left = I.onset_lo .- I.inf_lo
inc_sw = I.onset_hi .- I.onset_lo

C = d.chain_gi                                      # direct (chain) GI stream
gi_pw = C.pwin_hi .- C.pwin_lo
gi_left = max.(C.swin_lo .- C.pwin_lo, 0.0)
gi_sw = C.swin_hi .- C.swin_lo

N = d.onset_gi                                      # onset GI stream: X = E_j - O_i
x_lo = N.swin_lo .- N.or_onset_hi
x_hi = N.swin_hi .- N.or_onset_lo

# forward correction: condition on the INFECTOR's timing and right-truncate at the
# tracing cutoff. The infector window is censored (primary_censored), not a midpoint.
T_obs = 31.0                                        # cutoff (max observed onset)

# The forward correction evaluates a convolution CDF per infector. With integer-day
# windows the (window, cutoff) inputs repeat heavily, so we evaluate each *unique*
# combination once and weight by its multiplicity. This is exact (the per-individual
# log-terms are identical within a group) and collapses 262 CDF calls to ~46, which
# dominates the per-iteration cost. Chain: cutoff = T_obs - infector_infection_lo.
ch_pairs = collect(zip(C.pwin_hi .- C.pwin_lo, T_obs .- C.pwin_lo))
ch_u = unique(ch_pairs); ch_w = first.(ch_u); ch_cut = last.(ch_u)
ch_n = Float64[count(==(p), ch_pairs) for p in ch_u]
# Onset: cutoff = (T_obs - infector_onset_lo) + mu_i (mu_i added inside the model).
on_pairs = collect(zip(N.or_onset_hi .- N.or_onset_lo, T_obs .- N.or_onset_lo))
on_u = unique(on_pairs); on_w = first.(on_u); on_base = last.(on_u)
on_n = Float64[count(==(p), on_pairs) for p in on_u]

# log P(delay in [a, a+w)) with the primary event uniform over [0, pw]
function censln(dist, a, w, pw)
    pc = primary_censored(dist, Uniform(0.0, pw))
    log(max(cdf(pc, a + w) - cdf(pc, a), 1e-12))
end

# P(delay <= z) via interval_censored (AD-safe CDF; raw gamma CDF is not)
cdfle(dist, z) = z <= 1e-3 ? 1e-6 : pdf(interval_censored(dist, [1e-3, z]), 1e-3)


# onset pair: integrate the infector incubation out of GI = (E_j - O_i) + incubation.
# The inner P(GI in [a, b)) is an exact CDF difference via interval_censored, which
# differentiates under a ChainRules backend (the bare gamma CDF does not).
function onsetln(gi, incub, x_lo, x_hi)
    n = 24
    h = 18.0 / n
    s = 0.0
    for k in 0:(n - 1)
        inc = (k + 0.5) * h
        a = max(x_lo + inc, 1e-3)
        b = x_hi + inc
        b <= a && continue
        s += pdf(incub, inc) * pdf(interval_censored(gi, [a, b]), a)
    end
    log(max(s * h, 1e-12))
end

@model function joint(inc_pw, inc_left, inc_sw, gi_pw, gi_left, gi_sw, x_lo, x_hi)
    mu_g ~ LogNormal(log(4.0), 0.5);  k_g ~ LogNormal(log(2.0), 0.5)   # generation interval
    mu_i ~ LogNormal(log(4.5), 0.4);  k_i ~ LogNormal(log(3.0), 0.4)   # incubation period
    gi    = Gamma(k_g, mu_g / k_g)
    incub = Gamma(k_i, mu_i / k_i)

    for j in eachindex(inc_left)
        Turing.@addlogprob! censln(incub, inc_left[j], inc_sw[j], inc_pw[j])
    end
    for j in eachindex(gi_left)
        Turing.@addlogprob! censln(gi, gi_left[j], gi_sw[j], gi_pw[j])
    end
    for j in eachindex(x_lo)
        Turing.@addlogprob! onsetln(gi, incub, x_lo[j], x_hi[j])
    end

    # forward right-truncation, conditioned on the infector. Infectee onset <= T means
    # infectee infection <= T - incubation, so GI <= (T - mu_i) - E_infector. For chain
    # pairs E_infector is the infection midpoint; for onset pairs it is onset - mu_i, so
    # the GI cutoff is T - infector_onset.
    conv = convolve_distributions(gi, incub)           # GI + infectee incubation, built once
    for j in eachindex(ch_w)                           # chain: infector window censored
        pc = primary_censored(conv, Uniform(0.0, ch_w[j]))
        Turing.@addlogprob! -ch_n[j] * log(max(cdf(pc, ch_cut[j]), 1e-6))
    end
    for j in eachindex(on_w)                           # onset: onset window censored; infector
        pc = primary_censored(conv, Uniform(0.0, on_w[j]))     # incubation ~ mu_i (a difference
        Turing.@addlogprob! -on_n[j] * log(max(cdf(pc, on_base[j] + mu_i), 1e-6))  # convolution)
    end
end

Random.seed!(1)
model = joint(inc_pw, inc_left, inc_sw, gi_pw, gi_left, gi_sw, x_lo, x_hi)
init = Turing.DynamicPPL.InitFromParams((mu_g = 4.0, k_g = 2.0, mu_i = 4.5, k_i = 3.0))

# live progress: print a one-line summary every 50 iterations and flush, so the
# output file is readable while sampling runs in the background.
const N_ITER = 600
function progress_cb(rng, model, sampler, transition, state, i; kwargs...)
    (i == 1 || i % 50 == 0 || i == N_ITER) || return
    try
        nt = NamedTuple(transition.params)
        lp = transition.stats.logjoint
        println(rpad("iter $i/$N_ITER", 16),
                "mu_g=", round(nt.mu_g, digits = 2), " k_g=", round(nt.k_g, digits = 2),
                "  mu_i=", round(nt.mu_i, digits = 2), " k_i=", round(nt.k_i, digits = 2),
                "  lp=", round(lp, digits = 1))
    catch e
        println("iter $i/$N_ITER  [progress: ", typeof(transition), "]")
    end
    flush(stdout)
end
chn = sample(model, NUTS(; adtype = AutoForwardDiff()), N_ITER;
             initial_params = init, progress = false, callback = progress_cb)

# Persist the posterior draws so the Quarto report can read them without re-running
# the ~25 min sampler. Each gamma is parameterised by its mean (mu) and shape (k),
# so mu_g/mu_i are directly the mean generation interval / incubation period.
draws = DataFrame(mu_g = vec(chn[:mu_g]), k_g = vec(chn[:k_g]),
                  mu_i = vec(chn[:mu_i]), k_i = vec(chn[:k_i]))
CSV.write(joinpath(@__DIR__, "joint_gi_draws.csv"), draws)

for (nm, sym) in (("generation interval", :mu_g), ("incubation period", :mu_i))
    q = quantile(vec(chn[sym]), [0.5, 0.025, 0.975])
    println(rpad(nm, 22), "mean = ", round(q[1], digits = 2),
            " d  (", round(q[2], digits = 2), ", ", round(q[3], digits = 2), ")")
end
println("\nforward-corrected, all ", d.n_chain, " chain + ", d.n_onset, " onset pairs + ",
        d.n_infectees, " incubation obs  (backward joint GI was 2.98; chain-only 3.1)")
