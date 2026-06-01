### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 11111111-5555-1111-1111-111111111111
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
	Pkg.instantiate()
	Pkg.add("Revise")
    using Revise, FTAN
end

# ╔═╡ 22222222-5555-2222-2222-222222222222
md"# Source-State Matrix Workflow"

# ╔═╡ 33333333-5555-3333-3333-333333333333
begin
    dt = 0.25
    distance = 130.0
    t = collect((1:2048) .* dt)
    s1 = @. exp(-0.5 * ((t - distance / 3.3) / 8.0)^2) * cos(2π * 0.14 * t)
    s2 = @. exp(-0.5 * ((t - distance / 2.9) / 8.0)^2) * cos(2π * 0.18 * t)
    Wc = hcat(s1, s2, 0.7 .* s1 .+ 0.3 .* s2)
    Wac = hcat(0.95 .* s1, 0.8 .* s2, 0.6 .* s1 .+ 0.4 .* s2)
    periods = collect(range(4.0, 16.0, length=24))
    bank = MFTFilterBank(dt, size(Wc, 1), periods;
        velocity_range=(2.0, 5.0),
        bandwidth_factor=0.30,
        storage_mode=:full,
        N_initial=2 * size(Wc, 2))
    batch = analyze_causal_acausal_branches(Wc, Wac, distance, bank;
        state_labels=["state 1", "state 2", "mix"],
        compute_phase=true)
end

# ╔═╡ 44444444-5555-4444-4444-444444444444
plot_branch_correlation(batch; title="Synthetic source-state branch correlation")

# ╔═╡ Cell order:
# ╠═11111111-5555-1111-1111-111111111111
# ╟─22222222-5555-2222-2222-222222222222
# ╠═33333333-5555-3333-3333-333333333333
# ╠═44444444-5555-4444-4444-444444444444
