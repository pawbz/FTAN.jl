### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 11111111-2222-1111-1111-111111111111
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
	Pkg.instantiate()
	Pkg.add("Revise")
    using Revise, FTAN
end

# ╔═╡ 22222222-2222-2222-2222-222222222222
md"# Causal/Acausal U-c Consistency"

# ╔═╡ 33333333-2222-3333-3333-333333333333
begin
    dt = 0.25
    distance = 150.0
    t = collect((1:2048) .* dt)
    causal = @. exp(-0.5 * ((t - distance / 3.4) / 8.0)^2) * cos(2π * 0.16 * t)
    acausal = @. exp(-0.5 * ((t - distance / 3.25) / 8.5)^2) * cos(2π * 0.16 * t + 0.2)
    periods = collect(range(4.0, 16.0, length=26))
    result = analyze_causal_acausal_branches(causal, acausal, dt, distance, periods;
        velocity_range=(2.0, 5.0),
        bandwidth_factor=0.30,
        compute_phase=true)
end

# ╔═╡ 44444444-2222-4444-4444-444444444444
begin
    rows = vcat(
        FTAN._uc_consistency_rows_from_mft_result(result.causal_result;
            pair_label="synthetic", label="causal", branch="causal"),
        FTAN._uc_consistency_rows_from_mft_result(result.acausal_result;
            pair_label="synthetic", label="acausal", branch="acausal"))
    plot_uc_consistency_comparison(rows; pair_label="synthetic",
        title="Causal/acausal U-c consistency")
end

# ╔═╡ Cell order:
# ╠═11111111-2222-1111-1111-111111111111
# ╟─22222222-2222-2222-2222-222222222222
# ╠═33333333-2222-3333-3333-333333333333
# ╠═44444444-2222-4444-4444-444444444444
