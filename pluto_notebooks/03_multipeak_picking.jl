### A Pluto.jl notebook ###
# v0.20.27

using Markdown
using InteractiveUtils

# ╔═╡ 11111111-3333-1111-1111-111111111111
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
	Pkg.instantiate()
	Pkg.add("Revise")
    using Revise, FTAN
end

# ╔═╡ 22222222-3333-2222-2222-222222222222
md"# Multipeak Picking"

# ╔═╡ 33333333-3333-3333-3333-333333333333
begin
    dt = 0.25
    distance = 120.0
    t = collect((1:2048) .* dt)
    data = @. exp(-0.5 * ((t - distance / 3.4) / 6.0)^2) * cos(2π * 0.18 * t)
    data .+= @. 0.8 * exp(-0.5 * ((t - distance / 2.5) / 7.0)^2) * cos(2π * 0.18 * t)
    periods = collect(range(4.0, 14.0, length=24))
    res = perform_mft_analysis(data, dt, distance, periods;
        velocity_range=(2.0, 5.0),
        bandwidth_factor=0.30,
        compute_phase=true,
        storage_mode=:full)
    modes = extract_all_modes(res; distance=distance, max_modes=4)
end

# ╔═╡ 44444444-3333-4444-4444-444444444444
plot_multipeak_dispersion(modes; title="Multiple MFT envelope peaks")

# ╔═╡ Cell order:
# ╠═11111111-3333-1111-1111-111111111111
# ╟─22222222-3333-2222-2222-222222222222
# ╠═33333333-3333-3333-3333-333333333333
# ╠═44444444-3333-4444-4444-444444444444
