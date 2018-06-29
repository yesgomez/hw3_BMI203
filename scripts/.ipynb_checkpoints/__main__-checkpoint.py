from .sw import optimize_gap_penalties, import_pairs, temporary

negpairs_a, negpairs_b, pospairs_a, pospairs_b = import_pairs()
a, b = temporary()

# Run optimization first
results = optimize_gap_penalties(pospairs_a, pospairs_b)
results = optimize_gap_penalties(a, b)

# Then move on to...