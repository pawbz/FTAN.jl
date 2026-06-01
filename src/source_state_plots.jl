function plot_source_state_waveforms(item; dt, velocity_range, period_min, period_max)
    isnothing(item) && return PlutoPlotly.plot(PlutoPlotly.scatter(x=[0], y=[0], text=["No run selected"]))
    cluster_avg_ac = _column_normalise(item.acausal)
    cluster_avg_c = _column_normalise(item.causal)
    nth = size(cluster_avg_ac, 1)
    t_neg = [-(nth - i + 1) * dt for i in 1:nth]
    t_pos = [i * dt for i in 1:nth]
    t_full = [t_neg; t_pos]

    global_avg_ac = _vector_normalise(item.global_avg_ac)
    global_avg_c = _vector_normalise(item.global_avg_c)
    global_full = [reverse(global_avg_ac); global_avg_c]
    global_ncc = _ncc(global_avg_ac, global_avg_c)

    combo_labels_local = item.combo_labels
    ncomb = size(cluster_avg_ac, 2)
    traces = AbstractTrace[]
    colors = begin
        nc = max(ncomb, 1)
        cs = ColorSchemes.rainbow
        [Colors.hex(get(cs, (i - 1) / max(1, nc - 1))) for i in 1:nc]
    end

    total_ac = sum(item.counts_ac)
    total_c = sum(item.counts_c)
    amp_peak = maximum(abs.(vcat(vec(cluster_avg_ac), vec(cluster_avg_c), global_full)))
    vertical_spacing = amp_peak * 2.5 + 1e-3

    for combo_idx in 1:ncomb
        c = colors[mod1(combo_idx, length(colors))]
        a = cluster_avg_ac[:, combo_idx]
        b = cluster_avg_c[:, combo_idx]
        full_k = [reverse(a); b]
        ncc = _ncc(a, b)
        pct_ac = 100 * item.counts_ac[combo_idx] / max(total_ac, 1)
        pct_c = 100 * item.counts_c[combo_idx] / max(total_c, 1)
        state_label = combo_idx <= length(combo_labels_local) ? combo_labels_local[combo_idx] : string(combo_idx)
        legend_label = "State $(state_label) (ac: $(round(pct_ac; digits=1))%, c: $(round(pct_c; digits=1))%, corr=$(round(ncc; digits=3)))"
        offset = (combo_idx - 1) * vertical_spacing
        push!(traces, PlutoPlotly.scatter(x=t_full, y=global_full .+ offset, mode="lines",
            name=combo_idx == 1 ? "Global mean (corr=$(round(global_ncc; digits=3)))" : "Global mean",
            showlegend=combo_idx == 1,
            line=attr(color="rgba(0,0,0,0.18)", width=3)))
        push!(traces, PlutoPlotly.scatter(x=t_full, y=full_k .+ offset, mode="lines",
            name=legend_label, line=attr(color=c, width=2)))
    end

    shapes = if isnothing(item.distance)
        []
    else
        vmin, vmax = velocity_range
        t_fast = item.distance / vmax
        t_slow = item.distance / vmin
        [attr(type="line", x0=t, x1=t, y0=0, y1=1, yref="paper",
              line=attr(color="rgba(0,0,0,0.25)", width=1, dash="dash"))
         for t in (-t_slow, -t_fast, t_fast, t_slow)]
    end

    distance_label = isnothing(item.distance) ? "distance unavailable" : "$(round(Int, item.distance))km"
    title = "Source State Average Waveforms ($(item.pair_label) seed=$(item.seed) $(distance_label) $(_compact_number(period_min))-$(_compact_number(period_max))s)"
    return _plotly_plot(traces, PlutoPlotly.Layout(
        title=attr(text=title, font=attr(size=18, family="Computer Modern, serif")),
        height=500 * max(1, cld(ncomb, 5)),
        width=900,
        xaxis=attr(title="Lag (s)", zeroline=true, zerolinecolor="rgba(0,0,0,0.3)"),
        yaxis=attr(title="Amplitude"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=attr(x=0.5, xanchor="center", y=-0.2, orientation="h",
            font=attr(size=12, family="Computer Modern, serif")),
        shapes=shapes,
    ))
end

# ╔═╡ a1000006-0000-0000-0000-000000000001
function plot_cluster_histogram(counts_ac, counts_c; title="Cluster Usage", labels=nothing)
    K = length(counts_ac)
    total_ac = max(sum(counts_ac), 1)
    total_c  = max(sum(counts_c),  1)
    pct_ac = 100.0 .* Float64.(counts_ac) ./ total_ac
    pct_c  = 100.0 .* Float64.(counts_c)  ./ total_c
    xlabels = isnothing(labels) ? string.(1:K) : string.(labels)
    traces = [
        PlutoPlotly.bar(x=xlabels, y=pct_ac, name="Acausal",
            marker=attr(color="rgba(31,119,180,0.7)")),
        PlutoPlotly.bar(x=xlabels, y=pct_c,  name="Causal",
            marker=attr(color="rgba(214,39,40,0.7)")),
    ]
    layout = Layout(
        title=attr(text=title, font=attr(size=18)),
        barmode="group", height=400, width=700,
        xaxis=attr(title="Source state"),
        yaxis=attr(title="Usage (%)"),
        plot_bgcolor="white", paper_bgcolor="white",
    )
    return _plotly_plot(traces, layout)
end

# ╔═╡ a1000007-0000-0000-0000-000000000001
function plot_state_ncc_heatmap(acausal::AbstractMatrix, causal::AbstractMatrix;
        labels=nothing, title="State-State Normalised Correlation")
    n = size(acausal, 2)
    xlabels = isnothing(labels) ? string.(1:n) : string.(labels)

    function ncc_matrix(A)
        C = Matrix{Float32}(undef, n, n)
        cols = [begin v = vec(Float64.(A[:, i])); v .- mean(v) end for i in 1:n]
        norms = [norm(c) + 1e-8 for c in cols]
        for i in 1:n, j in 1:n
            C[i, j] = dot(cols[i], cols[j]) / (norms[i] * norms[j])
        end
        return C
    end

    C_ac = ncc_matrix(acausal)
    C_c  = ncc_matrix(causal)
    trace_ac = PlutoPlotly.heatmap(
        z=C_ac, x=xlabels, y=xlabels,
        colorscale="RdBu", zmid=0, zmin=-1, zmax=1,
        colorbar=attr(title="Corr", len=0.9, x=0.46),
        xaxis="x1", yaxis="y1",
    )
    trace_c = PlutoPlotly.heatmap(
        z=C_c, x=xlabels, y=xlabels,
        colorscale="RdBu", zmid=0, zmin=-1, zmax=1,
        colorbar=attr(title="Corr", len=0.9, x=1.01),
        xaxis="x2", yaxis="y2",
    )
    sz = max(350, n * 40)
    layout = Layout(
        title=attr(text=title, font=attr(size=16)),
        grid=attr(rows=1, columns=2, pattern="independent"),
        annotations=[
            attr(text="Acausal", x=0.22, xref="paper", y=1.05, yref="paper",
                 showarrow=false, font=attr(size=14)),
            attr(text="Causal",  x=0.78, xref="paper", y=1.05, yref="paper",
                 showarrow=false, font=attr(size=14)),
        ],
        xaxis=attr(title="State", tickangle=-45),
        yaxis=attr(title="State"),
        xaxis2=attr(title="State", tickangle=-45),
        yaxis2=attr(title="State"),
        width=900, height=sz + 80,
        plot_bgcolor="white", paper_bgcolor="white",
        margin=attr(t=80, b=80, l=80, r=80),
    )
    return _plotly_plot([trace_ac, trace_c], layout)
end
