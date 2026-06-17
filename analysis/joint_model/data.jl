# Load and shape the contact-tracing pairs for the joint generation-interval model.
#
# Each traced individual has an infection *window* (the reported contact window)
# and a symptom-onset day. From these we build three data streams:
#   incubation : infectee onset - infectee infection      (every traced infectee)
#   chain GIs  : infectee infection - infector infection   (infector also traced)
#   onset GIs  : infectee infection - infector onset        (infector onset only)
# Days are integers, so each observed time is itself a one-day interval.

using CSV, DataFrames

const DAY = 1.0  # daily censoring width

function load_pairs(path = joinpath(@__DIR__, "ct_pairs.csv"))
    df = CSV.read(path, DataFrame; missingstring = "NA")

    # one row per traced infectee -> incubation period observations
    incub = unique(df, :ee_id)
    incubation = (
        inf_lo = Float64.(incub.ee_inf_lo),           # infection window lower
        inf_hi = Float64.(incub.ee_inf_hi) .+ DAY,    # upper (+ a day, daily resolution)
        onset_lo = Float64.(incub.ee_onset),          # onset day
        onset_hi = Float64.(incub.ee_onset) .+ DAY,
    )

    # chain-of-three pairs: both infection windows known -> direct GI
    chain = df[.!ismissing.(df.or_inf_lo), :]
    chain_gi = (
        pwin_lo = Float64.(chain.or_inf_lo),                       # infector infection window
        pwin_hi = Float64.(chain.or_inf_hi) .+ DAY,
        swin_lo = Float64.(chain.ee_inf_lo),                       # infectee infection window
        swin_hi = Float64.(chain.ee_inf_hi) .+ DAY,
        ee_id = chain.ee_id,
    )

    # non-chain pairs with a known infector onset -> GI minus infector incubation
    nonchain = df[ismissing.(df.or_inf_lo) .& .!ismissing.(df.or_onset), :]
    onset_gi = (
        or_onset_lo = Float64.(nonchain.or_onset),                # infector onset day
        or_onset_hi = Float64.(nonchain.or_onset) .+ DAY,
        swin_lo = Float64.(nonchain.ee_inf_lo),                   # infectee infection window
        swin_hi = Float64.(nonchain.ee_inf_hi) .+ DAY,
        ee_id = nonchain.ee_id,
        n_cand = nonchain.n_cand,                                 # candidate infectors per infectee
    )

    return (; incubation, chain_gi, onset_gi,
            n_infectees = nrow(incub),
            n_chain = nrow(chain),
            n_onset = nrow(nonchain))
end

if abspath(PROGRAM_FILE) == @__FILE__
    d = load_pairs()
    println("infectees (incubation): ", d.n_infectees)
    println("chain GI pairs:         ", d.n_chain)
    println("onset GI pairs:         ", d.n_onset)
end
