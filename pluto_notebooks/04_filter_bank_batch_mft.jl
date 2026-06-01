### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 11111111-4444-1111-1111-111111111111
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
	Pkg.instantiate()
	Pkg.add("Revise")
    using Revise, FTAN
end

# ╔═╡ 22222222-4444-2222-2222-222222222222
md"# Filter-Bank Batch MFT"

# ╔═╡ 33333333-4444-3333-3333-333333333333
begin
    dt = 0.25
    distance = 140.0
    t = collect((1:2048) .* dt)
    base = @. exp(-0.5 * ((t - distance / 3.1) / 8.0)^2) * cos(2π * 0.14 * t)
    W = hcat(base, 0.8 .* base, base .+ 0.1 .* sin.(2π .* 0.05 .* t))
    periods = collect(range(4.0, 18.0, length=28))
    bank = MFTFilterBank(dt, size(W, 1), periods;
        velocity_range=(2.0, 5.0),
        bandwidth_factor=0.30,
        storage_mode=:picks_only,
        N_initial=size(W, 2))
    results = perform_mft_analysis_batch!(bank, W, distance; compute_phase=true)
end

# ╔═╡ 44444444-4444-4444-4444-444444444444
plot_uc_consistency_comparison(results;
    labels=["state 1", "state 2", "state 3"],
    pair_label="batch synthetic",
    title="Batch MFT U-c consistency")

# ╔═╡ Cell order:
# ╠═11111111-4444-1111-1111-111111111111
# ╟─22222222-4444-2222-2222-222222222222
# ╠═33333333-4444-3333-3333-333333333333
# ╠═44444444-4444-4444-4444-444444444444
