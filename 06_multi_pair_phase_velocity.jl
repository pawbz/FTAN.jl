### A Pluto.jl notebook ###
# v0.20.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 66660001-0000-0000-0000-000000000001
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
    Pkg.add("Revise")
    Pkg.add("PlutoUI")
    Pkg.add("CSV")
    Pkg.add("DataFrames")
    using Revise, FTAN
    using PlutoUI
    using CSV
    using DataFrames
    using Dates
    using PlutoPlotly
end

# ╔═╡ 66660002-0000-0000-0000-000000000002
md"""
# Multi-Pair Phase Velocity

This notebook picks independent per-pair group and phase velocity dispersion
curves from many station-pair noise cross-correlations in one batched MFT run.

Workflow:
1. Load station coordinates and noise cross-correlation waveforms for many
   station pairs.
2. For each pair, pick whichever branch (causal or acausal) the upstream
   pipeline already selected, giving one waveform per pair.
3. Compute per-pair great-circle distance.
4. Run batched `perform_mft_analysis_batch!` calls across all pairs.
5. Export selected pairs' group velocities and independent `ph_to_vel` phase
   velocities to pDSurfTomo.
"""

# ╔═╡ 66660003-0000-0000-0000-000000000003
md"## Station coordinates"

# ╔═╡ 66660004-0000-0000-0000-000000000004
begin
    stations_csv_path = "/mnt/NAS/Sanket_DRDO/Minneapolis_2011_2013/data/stationlists/Stations_XI_2011_13_SN.csv"
    stations_df = CSV.read(stations_csv_path, DataFrame)
    station_coords = Dict{String,Tuple{Float64,Float64}}(
        String(row["Station Code"]) => (Float64(row.Latitude), Float64(row.Longitude))
        for row in eachrow(stations_df)
    )
    length(station_coords)
end

# ╔═╡ 66660005-0000-0000-0000-000000000005
md"""
## Waveforms

Each CSV holds one station pair's noise cross-correlation function (NCF) as
`time_lag_s` from -500 to 500 s. Lags ≥ 0 are the **causal** branch
(`global_average_c`); lags ≤ 0 are the **acausal** branch (`global_average_ac`).
The `selected_branch` column records which half the upstream pipeline judged
usable for that pair — we use exactly that one branch per pair.
"""

# ╔═╡ 35b2d963-1341-4590-bc16-1eb59a37a43f
folder = "/mnt/sanket1/SN_waveforms/months-2012-03_to_2012-08"

# ╔═╡ 6666000c-0000-0000-0000-00000000000c
md"## Analysis parameters"

# ╔═╡ 6666000d-0000-0000-0000-00000000000d
md"""
**Period range (s):** $(@bind period_min Slider(4.0:1.0:15.0, default=10.0, show_value=true)) to $(@bind period_max Slider(15.0:1.0:40.0, default=30.0, show_value=true))

**Number of periods:** $(@bind n_periods Slider(8:1:256, default=128, show_value=true))

**Gaussian bandwidth (% of centre frequency):** $(@bind bandwidth_percent Slider(5.0:1.0:60.0, default=10.0, show_value=true))

**Velocity search window (km/s):** $(@bind vmin Slider(2.0:0.1:3.5, default=2.5, show_value=true)) to $(@bind vmax Slider(3.5:0.1:6.0, default=4.5, show_value=true))
"""

# ╔═╡ 66660006-0000-0000-0000-000000000006
function read_sn_waveform_csvs(folder)
    isdir(folder) || error("Folder does not exist: $(folder)")

    files = sort(filter(f -> endswith(lowercase(f), ".csv"), readdir(folder; join=true)))
    isempty(files) && return DataFrame()

    dfs = DataFrame[]
    for path in files
        df = CSV.read(path, DataFrame)
        df.source_file .= basename(path)
        push!(dfs, df)
    end

    vcat(dfs...; cols=:union)
end

# ╔═╡ 66660007-0000-0000-0000-000000000007
waveforms_df = read_sn_waveform_csvs(folder)

# ╔═╡ 6666002b-0000-0000-0000-00000000002b
begin
    waveform_value_columns = [
        String(name) for name in names(waveforms_df)
        if name ∉ (:time_lag_s, :pair, :source_file, :selected_branch, :selected_state, :selected_kind, :best_selected_state) &&
           eltype(skipmissing(waveforms_df[!, name])) <: Number
    ]
    waveform_column_options = vcat(
        ["selected_branch_average" => "Selected branch average (global_average_c/ac)"],
        [col => col for col in waveform_value_columns],
    )
end

# ╔═╡ 6666002c-0000-0000-0000-00000000002c
md"""
**Waveform column:** $(@bind selected_waveform_column Select(waveform_column_options, default="global_average_mean"))
"""

# ╔═╡ 66660008-0000-0000-0000-000000000008
md"""
## One waveform per pair

The source folder can contain more than one file for the same `pair` (e.g. a
re-run of the upstream pipeline). When that happens we keep only the most
recently written file per pair (by `source_file`'s position in the sorted file
list) and warn about which pairs were affected, so a stale duplicate never
silently doubles a pair's row count.

For the surviving rows per pair, extract the single branch (`selected_branch`)
already chosen for that pair: causal samples are lags `1.0:1.0:500.0`; acausal
samples are lags `-1.0:-1.0:-500.0` (time read backward from zero lag). The
waveform dropdown controls which CSV amplitude column is read. The default
branch-aware option uses `global_average_c` for causal and `global_average_ac`
for acausal; choosing an exact column uses that column for whichever lag branch
is selected. Either way we get a length-500 vector starting at `dt = 1.0 s`,
matching the convention that sample 1 occurs at time `dt`, not `0`.
"""

# ╔═╡ 66660009-0000-0000-0000-000000000009
begin
    grouped_raw = groupby(waveforms_df, :pair)
    pair_labels = String.(sort(collect(first.(keys(grouped_raw)))))
    npairs = length(pair_labels)
    nt = 500

    W = zeros(Float64, nt, npairs)
    duplicate_pairs = String[]
    for (j, lbl) in enumerate(pair_labels)
        g_all = grouped_raw[(pair=lbl,)]
        source_files = sort(unique(g_all.source_file))
        if length(source_files) > 1
            push!(duplicate_pairs, lbl)
        end
        keep_file = last(source_files)  # most recent by sorted filename/listing order
        g = sort(g_all[g_all.source_file .== keep_file, :], :time_lag_s)

        branch = first(g.selected_branch)
        waveform_column = if selected_waveform_column == "selected_branch_average"
            branch == "causal" ? :global_average_c : :global_average_ac
        else
            Symbol(selected_waveform_column)
        end

        if branch == "causal"
            rows = sort(g[g.time_lag_s .>= 1.0, :], :time_lag_s)
            W[:, j] .= Float64.(rows[!, waveform_column])
        else
            rows = sort(g[g.time_lag_s .<= -1.0, :], :time_lag_s, rev=true)
            W[:, j] .= Float64.(rows[!, waveform_column])
        end
    end
    isempty(duplicate_pairs) ||
        @warn "Multiple source files found for $(length(duplicate_pairs)) pair(s); kept only the most recent file per pair" duplicate_pairs
    (npairs, size(W))
end

# ╔═╡ 6666000a-0000-0000-0000-00000000000a
md"""
## Per-pair distance

`haversine_distance_km` is the package's exported great-circle distance helper,
used here for the analysis distance fed into `perform_mft_analysis_batch!`.
"""

# ╔═╡ 6666000b-0000-0000-0000-00000000000b
begin
    parse_pair_label(lbl) = (s = split(lbl, "-"); (String(s[1]), String(s[2])))

    pair_distances_km = Float64[]
    for lbl in pair_labels
        sta, stb = parse_pair_label(lbl)
        lat1, lon1 = station_coords[sta]
        lat2, lon2 = station_coords[stb]
        push!(pair_distances_km, haversine_distance_km(lat1, lon1, lat2, lon2))
    end
    length(pair_distances_km)
end

# ╔═╡ 6666000e-0000-0000-0000-00000000000e
periods_analysis = collect(range(period_min, period_max, length=n_periods))

# ╔═╡ 6666000f-0000-0000-0000-00000000000f
md"## Batched MFT analysis — one call for all pairs"

# ╔═╡ b44ea373-3ad3-40f7-94f6-bead3f5b1eab
W

# ╔═╡ 66660010-0000-0000-0000-000000000010
begin
    bank = MFTFilterBank(1.0, size(W, 1), periods_analysis;
        velocity_range=(vmin, vmax),
        bandwidth_factor=bandwidth_percent / 100.0,
        storage_mode=:picks_only,
        N_initial=size(W, 2))
    results = perform_mft_analysis_batch!(bank, W, pair_distances_km; compute_phase=true, wavelength_ref_velocity=3.2, wavelength_fraction=2.5)
end

# ╔═╡ 66660011-0000-0000-0000-000000000011
md"### Sanity check: group-velocity dispersion across all pairs"

# ╔═╡ 66660012-0000-0000-0000-000000000012
plot_dispersion_curve(results; names=pair_labels)

# ╔═╡ 66660026-0000-0000-0000-000000000026
md"### Independent per-pair phase velocity from ph_to_vel"

# ╔═╡ 66660027-0000-0000-0000-000000000027
independent_phtovel_results = perform_mft_analysis_batch!(bank, W, pair_distances_km;
    compute_phase=true,
    use_phtovel=true,
    phase_anchor_velocity=3.3,
    phase_velocity_range=(vmin, vmax), wavelength_ref_velocity=3.2, wavelength_fraction=2.5)

# ╔═╡ 66660013-0000-0000-0000-000000000013
md"## Per-pair MFT results"

# ╔═╡ 66660014-0000-0000-0000-000000000014
begin
    mft_results_by_pair = Dict(pair_labels[i] => results[i] for i in eachindex(pair_labels))
    length(mft_results_by_pair)
end

# ╔═╡ 66660028-0000-0000-0000-000000000028
begin
    independent_phtovel_results_by_pair = Dict(pair_labels[i] => independent_phtovel_results[i] for i in eachindex(pair_labels))
    independent_phtovel_phase_velocities = Dict(lbl => independent_phtovel_results_by_pair[lbl].phase_velocities for lbl in pair_labels)
    length(independent_phtovel_phase_velocities)
end

# ╔═╡ 66660029-0000-0000-0000-000000000029
let
    dist_min = minimum(independent_phtovel_results_by_pair[lbl].distance for lbl in pair_labels)
    dist_max = maximum(independent_phtovel_results_by_pair[lbl].distance for lbl in pair_labels)

    traces = [scatter()]
    colorbar_shown = false
    for lbl in pair_labels
        res = independent_phtovel_results_by_pair[lbl]
        v = independent_phtovel_phase_velocities[lbl]
        valid = findall(isfinite, v)
        isempty(valid) && continue
        push!(traces, scatter(
            x=res.periods[valid], y=v[valid],
            mode="lines+markers",
            marker=attr(
                size=6,
                color=fill(res.distance, length(valid)),
                colorscale="RdBu",
                cmin=dist_min, cmax=dist_max,
                showscale=!colorbar_shown,
                colorbar=attr(title="Distance (km)"),
            ),
            line=attr(color="rgba(0,0,0,0.0)", width=1.5),
            name=lbl,
            hovertemplate="$(lbl)<br>Period: %{x:.2f} s<br>c: %{y:.3f} km/s<br>d=$(round(res.distance; digits=1)) km<extra></extra>",
            showlegend=false,
        ))
        colorbar_shown = true
    end
    layout = Layout(
        title="Per-Pair Phase Velocity (independent ph_to_vel, colored by distance)",
        xaxis=attr(title="Period (s)"),
        yaxis=attr(title="Phase Velocity (km/s)", range=[vmin, vmax]),
        width=950, height=600,
        plot_bgcolor="white", paper_bgcolor="white",
        showlegend=false,
    )
    WideCell(PlutoPlotly.plot(traces, layout))
end

# ╔═╡ 6666001f-0000-0000-0000-00000000001f
md"""
## Export to pDSurfTomo

Write flat 6-column files (`period lat1 lon1 lat2 lon2 velocity`, whitespace
separated, no header) that `pDSurfTomo_v1.jl` reads directly in its Section 1
("Flat dispersion file"). Group and phase velocities are written separately,
and the pair picker below controls which station pairs are exported.
"""

# ╔═╡ 66660020-0000-0000-0000-000000000020
md"""
**Pairs to export**

$(@bind dsurftomo_selected_pairs MultiCheckBox(pair_labels; default=pair_labels, select_all=true))
"""

# ╔═╡ 66660021-0000-0000-0000-000000000021
begin
    dsurftomo_pairs_to_export = dsurftomo_selected_pairs isa AbstractVector ? String.(dsurftomo_selected_pairs) : pair_labels
    dsurftomo_pair_set = Set(dsurftomo_pairs_to_export)

    dsurftomo_results_by_pair = Dict(lbl => mft_results_by_pair[lbl]
        for lbl in pair_labels if lbl in dsurftomo_pair_set)

    dsurftomo_group_velocities = Dict(lbl => mft_results_by_pair[lbl].group_velocities
        for lbl in keys(dsurftomo_results_by_pair))

    dsurftomo_phase_velocities = Dict(lbl => independent_phtovel_phase_velocities[lbl]
        for lbl in keys(dsurftomo_results_by_pair))

    dsurftomo_group_rows = pdsurftomo_dispersion_rows(dsurftomo_group_velocities, dsurftomo_results_by_pair, station_coords)
    dsurftomo_phase_rows = pdsurftomo_dispersion_rows(dsurftomo_phase_velocities, dsurftomo_results_by_pair, station_coords)

    dsurftomo_output_dir = folder
    dsurftomo_waveform_suffix = replace(String(selected_waveform_column), r"[^A-Za-z0-9_.-]+" => "_")
    dsurftomo_datetime_suffix = Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HHMMSS")
    dsurftomo_group_output_path = joinpath(dsurftomo_output_dir, "group_velocity_$(dsurftomo_waveform_suffix)_$(dsurftomo_datetime_suffix)_for_dsurftomo.txt")
    dsurftomo_phase_output_path = joinpath(dsurftomo_output_dir, "phase_velocity_independent_phtovel_$(dsurftomo_waveform_suffix)_$(dsurftomo_datetime_suffix)_for_dsurftomo.txt")

    write_pdsurftomo_dispersion(dsurftomo_group_output_path, dsurftomo_group_rows)
    write_pdsurftomo_dispersion(dsurftomo_phase_output_path, dsurftomo_phase_rows)
end

# ╔═╡ 6666002a-0000-0000-0000-00000000002a
let
    row_key(r) = (r.lat1, r.lon1, r.lat2, r.lon2)
    pair_key(lbl) = begin
        sta, stb = parse_pair_label(lbl)
        lat1, lon1 = station_coords[sta]
        lat2, lon2 = station_coords[stb]
        (lat1, lon1, lat2, lon2)
    end

    dist_min = minimum(res.distance for res in values(dsurftomo_results_by_pair))
    dist_max = maximum(res.distance for res in values(dsurftomo_results_by_pair))

    rows_by_kind = [
        ("Group", dsurftomo_group_rows, "circle"),
        ("Phase (independent ph_to_vel)", dsurftomo_phase_rows, "diamond"),
    ]

    traces = [scatter()]
    colorbar_shown = false
    for lbl in pair_labels
        haskey(dsurftomo_results_by_pair, lbl) || continue
        res = dsurftomo_results_by_pair[lbl]
        branch_x = Float64[]
        branch_y = Float64[]
        branch_text = String[]
        has_branch_point = false
        for ib in axes(res.phase_velocity_branches, 2)
            branch_number = res.phase_branch_numbers[ib]
            for ip in eachindex(res.periods)
                c = res.phase_velocity_branches[ip, ib]
                isfinite(c) && vmin <= c <= vmax || continue
                push!(branch_x, res.periods[ip])
                push!(branch_y, c)
                push!(branch_text, "N=$(branch_number)")
                has_branch_point = true
            end
            push!(branch_x, NaN)
            push!(branch_y, NaN)
            push!(branch_text, "")
        end
        has_branch_point || continue
        push!(traces, scatter(
            x=branch_x,
            y=branch_y,
            mode="markers",
            name="$(lbl) branches",
            legendgroup=lbl,
            showlegend=false,
            marker=attr(
                size=4,
                symbol="x",
                color=fill(res.distance, length(branch_y)),
                colorscale="RdBu",
                cmin=dist_min, cmax=dist_max,
                showscale=!colorbar_shown,
                colorbar=attr(title="Distance (km)"),
                opacity=0.28,
            ),
            hovertext=branch_text,
            hovertemplate="$(lbl)<br>candidate branch %{hovertext}<br>Period: %{x:.2f} s<br>c: %{y:.3f} km/s<br>d=$(round(res.distance; digits=1)) km<extra></extra>",
        ))
        colorbar_shown = true
    end
    for (ikind, (kind, rows, symbol)) in enumerate(rows_by_kind)
        for lbl in pair_labels
            haskey(dsurftomo_results_by_pair, lbl) || continue
            res = dsurftomo_results_by_pair[lbl]
            key = pair_key(lbl)
            pair_rows = sort(filter(r -> row_key(r) == key, rows), by=r -> r.period)
            isempty(pair_rows) && continue
            push!(traces, scatter(
                x=[r.period for r in pair_rows],
                y=[r.velocity for r in pair_rows],
                mode="lines+markers",
                name=lbl,
                legendgroup=lbl,
                showlegend=ikind == 1,
                marker=attr(
                    size=7,
                    symbol=symbol,
                    color=fill(res.distance, length(pair_rows)),
                    colorscale="RdBu",
                    cmin=dist_min, cmax=dist_max,
                    showscale=!colorbar_shown,
                    colorbar=attr(title="Distance (km)"),
                ),
                line=attr(color="rgba(0,0,0,0.0)", width=1.5),
                hovertemplate="$(lbl)<br>$(kind)<br>Period: %{x:.2f} s<br>Velocity: %{y:.3f} km/s<br>d=$(round(res.distance; digits=1)) km<extra></extra>",
            ))
            colorbar_shown = true
        end
    end

    layout = Layout(
        title="Exported pDSurfTomo Velocities ($(length(dsurftomo_results_by_pair)) selected pairs)",
        xaxis=attr(title="Period (s)"),
        yaxis=attr(title="Velocity (km/s)", range=[vmin, vmax]),
        width=950, height=600,
        plot_bgcolor="white", paper_bgcolor="white",
    )
    WideCell(PlutoPlotly.plot(traces, layout))
end

# ╔═╡ 66660022-0000-0000-0000-000000000022
md"""
**Wrote group:** $(length(dsurftomo_group_rows)) rows · $(length(unique(r.period for r in dsurftomo_group_rows))) periods ·
$(length(dsurftomo_results_by_pair)) selected pairs
to `$(dsurftomo_group_output_path)`

**Wrote phase (independent ph_to_vel):** $(length(dsurftomo_phase_rows)) rows · $(length(unique(r.period for r in dsurftomo_phase_rows))) periods ·
$(length(dsurftomo_results_by_pair)) selected pairs
to `$(dsurftomo_phase_output_path)`
"""

# ╔═╡ Cell order:
# ╠═66660001-0000-0000-0000-000000000001
# ╟─66660002-0000-0000-0000-000000000002
# ╟─66660003-0000-0000-0000-000000000003
# ╠═66660004-0000-0000-0000-000000000004
# ╟─66660005-0000-0000-0000-000000000005
# ╠═35b2d963-1341-4590-bc16-1eb59a37a43f
# ╟─6666002c-0000-0000-0000-00000000002c
# ╟─6666000c-0000-0000-0000-00000000000c
# ╟─6666000d-0000-0000-0000-00000000000d
# ╠═66660006-0000-0000-0000-000000000006
# ╠═66660007-0000-0000-0000-000000000007
# ╠═6666002b-0000-0000-0000-00000000002b
# ╟─66660008-0000-0000-0000-000000000008
# ╠═66660009-0000-0000-0000-000000000009
# ╟─6666000a-0000-0000-0000-00000000000a
# ╠═6666000b-0000-0000-0000-00000000000b
# ╠═6666000e-0000-0000-0000-00000000000e
# ╟─6666000f-0000-0000-0000-00000000000f
# ╠═b44ea373-3ad3-40f7-94f6-bead3f5b1eab
# ╠═66660010-0000-0000-0000-000000000010
# ╟─66660011-0000-0000-0000-000000000011
# ╠═66660012-0000-0000-0000-000000000012
# ╟─66660026-0000-0000-0000-000000000026
# ╠═66660027-0000-0000-0000-000000000027
# ╟─66660013-0000-0000-0000-000000000013
# ╠═66660014-0000-0000-0000-000000000014
# ╠═66660028-0000-0000-0000-000000000028
# ╟─66660029-0000-0000-0000-000000000029
# ╟─6666001f-0000-0000-0000-00000000001f
# ╠═66660020-0000-0000-0000-000000000020
# ╠═66660021-0000-0000-0000-000000000021
# ╠═6666002a-0000-0000-0000-00000000002a
# ╟─66660022-0000-0000-0000-000000000022
