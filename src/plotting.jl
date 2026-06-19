_plotly_traces(traces::AbstractVector) = AbstractTrace[trace for trace in traces]
_plotly_plot(traces::AbstractVector, layout; kwargs...) =
    PlutoPlotly.plot(_plotly_traces(traces), layout; kwargs...)

const DEFAULT_VELOCITY_RANGE = (1.0, 7.0)

_velocity_range_or_default(velocity_range) = isnothing(velocity_range) ?
    [Float64(DEFAULT_VELOCITY_RANGE[1]), Float64(DEFAULT_VELOCITY_RANGE[2])] :
    [Float64(velocity_range[1]), Float64(velocity_range[2])]

_compact_label(x; maxlen::Int=72) = begin
    s = String(x)
    lastindex(s) <= maxlen ? s : string(first(s, maxlen - 1), "...")
end

function _join_compact(items; max_items::Int=4)
    vals = String.(collect(items))
    isempty(vals) && return ""
    head = first(vals, min(max_items, length(vals)))
    suffix = length(vals) > max_items ? " + $(length(vals) - max_items) more" : ""
    join(head, ", ") * suffix
end

function _title_with_subtitle(title::AbstractString, subtitle=nothing;
                              font_family::AbstractString="Arial, sans-serif",
                              font_size::Integer=16)
    clean_title = _compact_label(title; maxlen=95)
    isnothing(subtitle) || isempty(String(subtitle)) ?
        attr(text=clean_title, font=attr(size=font_size + 3, family=font_family)) :
        attr(
            text="$(clean_title)<br><sup>$(_compact_label(subtitle; maxlen=140))</sup>",
            font=attr(size=font_size + 3, family=font_family)
        )
end

_bottom_legend(font_size::Integer; y::Real=-0.22) = attr(
    orientation="h",
    x=0.5,
    xanchor="center",
    y=y,
    yanchor="top",
    font=attr(size=max(font_size - 2, 9)),
    bgcolor="rgba(255,255,255,0.92)",
    bordercolor="rgba(0,0,0,0.18)",
    borderwidth=1,
    tracegroupgap=8
)

_publication_margin(; l=80, r=60, t=92, b=150) = attr(l=l, r=r, t=t, b=b)

function plot_envelopes(res;
                  colorscale="Reds",
                  width=800,
                  height=500,
                  font_family="Arial, sans-serif",
                  font_size=20,
                  title="Multiple Filter Technique Analysis",
                  show_picks=true,
                  show_branches=true,
                  pick_color="white",
                  pick_width=3)

    mft_has_full_storage(res) ||
        return PlutoPlotly.plot(scatter(x=[0.0], y=[0.0],
            text=["Envelope storage unavailable for picks-only MFT result"]))

	vsticks = range(1, stop=8, step=0.5)
	tticks = res.distance ./ vsticks
	tmin = res.distance ./ 2.0
	tmax = res.distance ./ 6.0
	
    layout = Layout(
        title = _title_with_subtitle(title, "velocity=$(vel_min)-$(vel_max) km/s; n=$(length(results))";
            font_family=font_family, font_size=font_size),
        yaxis = attr(
            title = "Velocity (km/s)",
			range = (tmin, tmax),
            showgrid = true,
            gridcolor = "rgba(128,128,128,0.2)",
            zeroline = false,
            titlefont = attr(size = font_size, family = font_family),
            tickfont = attr(size = font_size - 2, family = font_family),
			tickmode = "array",
        	tickvals = tticks,
        	ticktext = vsticks
        ),
        xaxis = attr(
            title = "Period (s)",
            type = "linear",
            showgrid = true,
            gridcolor = "rgba(128,128,128,0.2)",
            zeroline = false,
            titlefont = attr(size = font_size, family = font_family),
            tickfont = attr(size = font_size - 2, family = font_family)
        ),
        coloraxis = attr(
            colorbar = attr(
                title = attr(
                    text = "Envelope<br>Amplitude",
                    font = attr(size = font_size - 1, family = font_family)
                ),
                tickfont = attr(size = font_size - 2, family = font_family),
                len = 0.85,
                thickness = 20
            ),
            colorscale = colorscale,
            cmin = 0,  # Start from 0 for envelope
            cmax = quantile(vec(res.envelopes), 1.0)  # Clip to a percentile
        ),
        width = width,
        height = height,
        plot_bgcolor = "white",
        paper_bgcolor = "white",
        margin = attr(l=80, r=120, t=80, b=80),
        showlegend = show_picks
    )
    E = res.envelopes ./ maximum(res.envelopes, dims=1)
    # Create heatmap trace
    heatmap_trace = contour(
        y = res.time,
        x = res.periods,
        z = E,
        coloraxis = "coloraxis",
        hovertemplate = "Time: %{x:.2f} s<br>Period: %{y:.2f} s<br>Amplitude: %{z:.3f}<extra></extra>"
    )
    
    traces = [heatmap_trace]

    # Add overtone / multiple-branch picks if requested
    if show_branches && !isempty(res.all_peaks)
        branch_colors = ["cyan", "deepskyblue", "lime", "yellow"]
        nbranches = maximum(length.(res.all_peaks); init=0)
        for branch in 1:nbranches
            x_branch = Float64[]
            y_branch = Float64[]
            for (ip, peaks) in enumerate(res.all_peaks)
                if length(peaks) >= branch
                    push!(x_branch, res.periods[ip])
                    push!(y_branch, peaks[branch][1])
                end
            end
            if !isempty(x_branch)
                push!(traces, scatter(
                    x = x_branch,
                    y = y_branch,
                    mode = "lines+markers",
                    line = attr(color = branch_colors[mod1(branch, length(branch_colors))], width = 1.5, dash = branch == 1 ? "solid" : "dot"),
                    marker = attr(size = branch == 1 ? 6 : 5, color = branch_colors[mod1(branch, length(branch_colors))], symbol = branch == 1 ? "circle-open" : "x"),
                    name = branch == 1 ? "Peak branch 1 (max)" : "Peak branch $(branch)",
                    hovertemplate = "Period: %{x:.2f} s<br>Arrival: %{y:.2f} s<br>Branch: $(branch)<extra></extra>",
                    showlegend = true
                ))
            end
        end
    end
    
    # Add arrival time picks if requested
    if show_picks
        # Filter out NaN values
        valid_indices = findall(!isnan, res.arrival_times)
        if !isempty(valid_indices)
            # Color the picks by the per-period quality factor and add a colorbar
            pick_colors = res.quality_factors[valid_indices]
            # Use grayscale: black -> white (reversed Greys so that higher quality is lighter)
            pick_trace = scatter(
                x = res.periods[valid_indices],
                y = res.arrival_times[valid_indices],
                mode = "lines+markers",
                line = attr(color = pick_color, width = pick_width),
                marker = attr(
                    color = pick_colors,
                    colorscale = "Greys",
                    size = 8,
                    symbol = "circle",
                    line = attr(color = "black", width = 1),
                    showscale = false,
                    colorbar = attr(
                        title = attr(text = "Quality Factor", font = attr(size = font_size - 2)),
                        tickfont = attr(size = font_size - 4),
                        len = 0.7,
                        thickness = 15
                    )
                ),
                name = "Group Velocity Picks",
                hovertemplate = "Period: %{x:.2f} s<br>Arrival: %{y:.2f} s<br>Velocity: " .* 
                               string.(round.(res.group_velocities[valid_indices], digits=2)) .* 
                               " km/s<br>Quality: %{marker.color:.2f}<extra></extra>",
                showlegend = true
            )
            push!(traces, pick_trace)
        end
    end
    
    return _plotly_plot(traces, layout)
end

# ╔═╡ df761225-7a00-4a7c-b56f-611243a75e97
function plot_dispersion_curve(results::AbstractVector{MFTResult};
						 names = [string("Result ", i) for i in 1:length(results)],
                         colorscale="Viridis",
                         width=800,
                         height=500,
                         font_family="Arial, sans-serif",
                         font_size=20,
                         title="Group Velocity Dispersion (multiple)",
                         color_by_quality=true,
                         velocity_range=nothing,  # (vmin, vmax) or nothing for default 1-7 km/s
                         show_phase=false)  # optionally overlay phase velocity

    # Allow passing a single result as a convenience
    if length(results) == 1
        return plot_dispersion_curve(results[1]; colorscale=colorscale, width=width, height=height,
                                     font_family=font_family, font_size=font_size, title=title,
                                     color_by_quality=color_by_quality, velocity_range=velocity_range,
                                     show_phase=show_phase)
    end

    # Collect valid group velocities across all results to determine autoscaling
    all_group_vals = Float64[]
    for r in results
        append!(all_group_vals, filter(!isnan, r.group_velocities))
    end

    if isempty(all_group_vals)
        @warn "No valid group velocity measurements found across results"
        return PlutoPlotly.plot(scatter(x=[0], y=[0], text=["No valid data"]))
    end

    # Determine velocity range
    vel_min, vel_max = _velocity_range_or_default(velocity_range)

    layout = Layout(
        title = _title_with_subtitle(title, "velocity=$(vel_min)-$(vel_max) km/s; n=$(length(results))";
            font_family=font_family, font_size=font_size),
        xaxis = attr(
            title = "Period (s)",
            type = "linear",
            showgrid = true,
            gridcolor = "rgba(128,128,128,0.2)",
            zeroline = false,
            titlefont = attr(size = font_size, family = font_family),
            tickfont = attr(size = font_size - 2, family = font_family)
        ),
        yaxis = attr(
            title = "Group Velocity (km/s)",
            showgrid = true,
            gridcolor = "rgba(128,128,128,0.2)",
            zeroline = false,
            range = [vel_min, vel_max],
            titlefont = attr(size = font_size, family = font_family),
            tickfont = attr(size = font_size - 2, family = font_family)
        ),
        width = width,
        height = height,
        plot_bgcolor = "white",
        paper_bgcolor = "white",
        margin = _publication_margin(l=80, r=60, t=95, b=135),
        showlegend = true,
        legend = _bottom_legend(font_size; y=-0.24)
    )

    traces = [scatter()]

    # line color palette to cycle through for different results
    palette = ["steelblue", "firebrick", "darkgreen", "purple", "orange", "brown", "teal", "magenta"]

    # Plot each result using a distinct line color (no markers)
    for (i, r) in enumerate(results)
        valid_idx = findall(!isnan, r.group_velocities)
        if isempty(valid_idx)
            continue
        end

        periods_valid = r.periods[valid_idx]
        group_vel_valid = r.group_velocities[valid_idx]
        color = palette[mod1(fld(i+1, 2), length(palette))]

        # alternate between solid and dotted line styles while keeping the same color
        dash_style = isodd(i) ? "solid" : "dash"
        group_trace = scatter(
            x = periods_valid,
            y = group_vel_valid,
            mode = "lines",
            line = attr(color = color, width = 2, dash = dash_style),
            name = names[i],
            hovertemplate = "Period: %{x:.2f} s<br>Velocity: %{y:.3f} km/s<extra></extra>"
        )

        push!(traces, group_trace)
    end

    return _plotly_plot(traces, layout)
end

# ╔═╡ ee08ecac-5010-4a64-820e-43ee80823939
function plot_dispersion_curve(res::MFTResult;
                         colorscale="Viridis",
                         width=800,
                         height=500,
                         font_family="Arial, sans-serif",
                         font_size=20,
                         title="Group Velocity Dispersion",
                         color_by_quality=true,
                         velocity_range=nothing,  # (vmin, vmax) or nothing for default 1-7 km/s
                         show_phase=false)  # optionally overlay phase velocity
    
    # Filter valid measurements (non-NaN velocities)
    valid_idx = findall(!isnan, res.group_velocities)
    
    if isempty(valid_idx)
        @warn "No valid group velocity measurements found"
        return PlutoPlotly.plot(scatter(x=[0], y=[0], text=["No valid data"]))
    end
    
    periods_valid = res.periods[valid_idx]
    group_vel_valid = res.group_velocities[valid_idx]
    quality_valid = res.quality_factors[valid_idx]
    
    # Determine velocity range
    vel_min, vel_max = _velocity_range_or_default(velocity_range)
    
    layout = Layout(
        title = _title_with_subtitle(title, "velocity=$(vel_min)-$(vel_max) km/s; periods=$(length(valid_idx))";
            font_family=font_family, font_size=font_size),
        xaxis = attr(
            title = "Period (s)",
            type = "linear",
            showgrid = true,
            gridcolor = "rgba(128,128,128,0.2)",
            zeroline = false,
            titlefont = attr(size = font_size, family = font_family),
            tickfont = attr(size = font_size - 2, family = font_family)
        ),
        yaxis = attr(
            title = "Group Velocity (km/s)",
            showgrid = true,
            gridcolor = "rgba(128,128,128,0.2)",
            zeroline = false,
            range = [vel_min, vel_max],
            titlefont = attr(size = font_size, family = font_family),
            tickfont = attr(size = font_size - 2, family = font_family)
        ),
        width = width,
        height = height,
        plot_bgcolor = "white",
        paper_bgcolor = "white",
        margin = _publication_margin(l=80, r=60, t=95, b=135),
        showlegend = true,
        legend = _bottom_legend(font_size; y=-0.24)
    )
    
    traces = [scatter()]
    
    if color_by_quality
        # Normalize quality factors for color mapping
        q_norm = (quality_valid .- minimum(quality_valid)) ./ 
                 (maximum(quality_valid) - minimum(quality_valid) .+ 1e-8)
        
        group_trace = scatter(
            x = periods_valid,
            y = group_vel_valid,
            mode = "markers+lines",
            marker = attr(
                size = 10,
                color = quality_valid,
                colorscale = colorscale,
                showscale = true,
                colorbar = attr(
                    title = "Quality<br>Factor",
                    titlefont = attr(size = font_size - 2),
                    tickfont = attr(size = font_size - 4),
                    len = 0.7,
                    thickness = 15
                ),
                line = attr(color = "white", width = 1)
            ),
            line = attr(color = "gray", width = 1, dash = "dot"),
            name = "Group Velocity",
            hovertemplate = "Period: %{x:.2f} s<br>Velocity: %{y:.3f} km/s<br>Quality: %{marker.color:.1f}<extra></extra>"
        )
    else
        group_trace = scatter(
            x = periods_valid,
            y = group_vel_valid,
            mode = "markers+lines",
            marker = attr(
                size = 10,
                color = "steelblue",
                line = attr(color = "white", width = 1)
            ),
            line = attr(color = "steelblue", width = 2),
            name = "Group Velocity",
            hovertemplate = "Period: %{x:.2f} s<br>Velocity: %{y:.3f} km/s<extra></extra>"
        )
    end
    
    push!(traces, group_trace)

    # --- Phase velocity overlay (sacmft96 dophvel implementation) ---
    if show_phase
        phv_idx = findall(i -> !isnan(res.phase_velocities[i]), 1:length(res.periods))
        if !isempty(phv_idx)
            phase_trace = scatter(
                x = res.periods[phv_idx],
                y = res.phase_velocities[phv_idx],
                mode = "markers+lines",
                marker = attr(
                    size = 8, color = "firebrick", symbol = "diamond",
                    line = attr(color = "white", width = 1)
                ),
                line = attr(color = "firebrick", width = 2, dash = "dash"),
                name = "Phase Velocity",
                hovertemplate = "Period: %{x:.2f} s<br>Phase Vel: %{y:.3f} km/s<extra></extra>"
            )
            push!(traces, phase_trace)
        end
    end

    return _plotly_plot(traces, layout)
end

# ╔═╡ 3bab4eab-92da-4539-b3df-d5c637d1ce77
# Plot all individual CCs with the stack overlaid
function plot_crosscorr_bundle(result; title="Cross-Correlations")
    traces = [scatter()]
    
    # Individual CCs (semi-transparent)
    for (i, cc) in enumerate(result.cc_all)
        push!(traces, scatter(
            x = result.lags,
            y = cc,
            mode = "lines",
            line = attr(color = "lightgray", width = 0.5),
            showlegend = false,
            hoverinfo = "skip"
        ))
    end
    
    # Stacked average (bold)
    push!(traces, scatter(
        x = result.lags,
        y = result.cc_average,
        mode = "lines",
        line = attr(color = "red", width = 3),
        name = "Stack (n=$(result.n_events))"
    ))
    
    layout = Layout(
        title = title,
        xaxis = attr(title = "Time Lag (s)", zeroline = true),
        yaxis = attr(title = "Normalized Amplitude", zeroline = true),
        width = 900,
        height = 500,
        plot_bgcolor = "white"
    )
    
    return _plotly_plot(traces, layout)
end

function compute_velocity_error_metrics(res::MFTResult,
                                        periods_true::Vector{Float64},
                                        cg_true::Vector{Float64},
                                        cp_true::Vector{Float64})
    group_errors = Float64[]
    phase_errors = Float64[]

    for i in eachindex(res.periods)
        idx_ref = argmin(abs.(periods_true .- res.periods[i]))

        gv = res.group_velocities[i]
        if !isnan(gv)
            gref = cg_true[idx_ref]
            push!(group_errors, 100.0 * abs(gv - gref) / max(abs(gref), 1e-8))
        end

        pv = res.phase_velocities[i]
        if !isnan(pv)
            pref = cp_true[idx_ref]
            push!(phase_errors, 100.0 * abs(pv - pref) / max(abs(pref), 1e-8))
        end
    end

    group_mean = isempty(group_errors) ? NaN : mean(group_errors)
    group_max = isempty(group_errors) ? NaN : maximum(group_errors)
    phase_mean = isempty(phase_errors) ? NaN : mean(phase_errors)
    phase_max = isempty(phase_errors) ? NaN : maximum(phase_errors)

    return (; group_mean, group_max, phase_mean, phase_max,
            n_group=length(group_errors), n_phase=length(phase_errors))
end

function plot_phase_velocities(res::MFTResult;
                         width=900,
                         height=550,
                         font_family="Arial, sans-serif",
                         font_size=20,
                         title="Phase Velocity Branches",
                         velocity_range=nothing,
                         show_group=true)

    valid_phase = .!isnan.(res.phase_velocity_branches)
    if !any(valid_phase)
        @warn "No valid phase velocity branch measurements found"
        return PlutoPlotly.plot(scatter(x=[0], y=[0], text=["No valid phase branch data"]))
    end

    all_phase_vals = vec(res.phase_velocity_branches[valid_phase])
    all_vel_vals = copy(all_phase_vals)
    if show_group
        append!(all_vel_vals, filter(!isnan, res.group_velocities))
    end

    if isempty(all_vel_vals)
        @warn "No valid velocity measurements found"
        return PlutoPlotly.plot(scatter(x=[0], y=[0], text=["No valid velocity data"]))
    end

    vel_min, vel_max = _velocity_range_or_default(velocity_range)

    layout = Layout(
        title = _title_with_subtitle(title, "velocity=$(vel_min)-$(vel_max) km/s";
            font_family=font_family, font_size=font_size),
        xaxis = attr(
            title = "Period (s)",
            type = "linear",
            showgrid = true,
            gridcolor = "rgba(128,128,128,0.2)",
            zeroline = false,
            titlefont = attr(size = font_size, family = font_family),
            tickfont = attr(size = font_size - 2, family = font_family)
        ),
        yaxis = attr(
            title = "Velocity (km/s)",
            showgrid = true,
            gridcolor = "rgba(128,128,128,0.2)",
            zeroline = false,
            range = [vel_min, vel_max],
            titlefont = attr(size = font_size, family = font_family),
            tickfont = attr(size = font_size - 2, family = font_family)
        ),
        width = width,
        height = height,
        plot_bgcolor = "white",
        paper_bgcolor = "white",
        margin = _publication_margin(l=80, r=60, t=95, b=140),
        showlegend = true,
        legend = _bottom_legend(font_size; y=-0.24)
    )

    traces = [scatter()]
    branch_palette = ["#2c7bb6", "#00a6ca", "#00ccbc", "#d7191c", "#fdae61", "#abd9e9", "#7b3294"]

    for (ib, branch) in enumerate(res.phase_branch_numbers)
        valid_idx = findall(i -> !isnan(res.phase_velocity_branches[i, ib]), 1:length(res.periods))
        isempty(valid_idx) && continue

        is_canonical = any(i -> res.selected_phase_branches[i] == branch &&
                                isfinite(res.phase_velocities[i]) &&
                                res.phase_velocities[i] == res.phase_velocity_branches[i, ib],
                           valid_idx)
        push!(traces, scatter(
            x = res.periods[valid_idx],
            y = res.phase_velocity_branches[valid_idx, ib],
            mode = "lines+markers",
            marker = attr(
                size = is_canonical ? 8 : 6,
                color = branch_palette[mod1(ib, length(branch_palette))],
                symbol = is_canonical ? "diamond" : "circle-open",
                line = attr(color = "white", width = 1)
            ),
            line = attr(
                color = branch_palette[mod1(ib, length(branch_palette))],
                width = is_canonical ? 3 : 1.5,
                dash = is_canonical ? "solid" : (branch < 0 ? "dot" : "dash")
            ),
            name = branch == 0 ? "Phase N=0" : "Phase N=$(branch > 0 ? "+" : "")$(branch)",
            hovertemplate = "Period: %{x:.2f} s<br>Phase Vel: %{y:.3f} km/s<br>Branch: $(branch)<extra></extra>"
        ))
    end

    if show_group
        valid_group = findall(!isnan, res.group_velocities)
        if !isempty(valid_group)
            push!(traces, scatter(
                x = res.periods[valid_group],
                y = res.group_velocities[valid_group],
                mode = "lines+markers",
                marker = attr(size = 7, color = "black", symbol = "x"),
                line = attr(color = "black", width = 2, dash = "dot"),
                name = "Group velocity",
                hovertemplate = "Period: %{x:.2f} s<br>Group Vel: %{y:.3f} km/s<extra></extra>"
            ))
        end
    end

    return _plotly_plot(traces, layout)
end

function plot_multipeak_dispersion(modes::Vector{MultimodalDispersion}; 
                                  title="Multimodal Group Velocity Dispersion",
                                  colorscale="Viridis",
                                  width=900,
                                  height=600,
                                  font_family="Arial, sans-serif",
                                  font_size=16,
                                  show_amplitudes=true,
                                  velocity_range=nothing)
    
    if isempty(modes)
        @warn "No modes provided for plotting"
        return PlutoPlotly.plot(scatter(x=[0], y=[0], text=["No modes"]))
    end
    
    # Collect all valid velocity measurements to determine axis range
    all_vels = Float64[]
    for mode in modes
        append!(all_vels, filter(!isnan, mode.group_velocities))
    end
    
    if isempty(all_vels)
        @warn "No valid velocity measurements found in modes"
        return PlutoPlotly.plot(scatter(x=[0], y=[0], text=["No valid data"]))
    end
    
    vel_min, vel_max = _velocity_range_or_default(velocity_range)
    
    # Generate color palette for modes
    n_modes = length(modes)
    colors = _sample_scheme_colors(colorscale, n_modes)
    
    # Build traces for each mode
    traces = []
    
    for (mode_idx, mode) in enumerate(modes)
        # Filter valid measurements
        valid_idx = findall(!isnan, mode.group_velocities)
        
        if isempty(valid_idx)
            continue
        end
        
        periods_valid = mode.periods[valid_idx]
        vels_valid = mode.group_velocities[valid_idx]
        amps_valid = mode.peak_amplitudes[valid_idx]
        
        # Normalize amplitudes for marker sizing
        amp_min = minimum(amps_valid)
        amp_max = maximum(amps_valid)
        marker_sizes = if amp_max > amp_min
            8 .+ 12 .* (amps_valid .- amp_min) ./ (amp_max - amp_min)
        else
            fill(12.0, length(amps_valid))
        end
        
        # Create trace for this mode
        color = colors[mode_idx]
        trace = scatter(
            x=periods_valid,
            y=vels_valid,
            mode="markers+lines",
            marker=attr(
                size=marker_sizes,
                color=color,
                opacity=0.8,
                line=attr(color="white", width=1)
            ),
            line=attr(
                color=color,
                width=2,
                dash="solid"
            ),
            name="Mode $(mode.mode_index) (n=$(length(valid_idx)))",
            hovertemplate="<b>Mode $(mode.mode_index)</b><br>" *
                          "Period: %{x:.2f} s<br>" *
                          "Group velocity: %{y:.3f} km/s<br>" *
                          "Amplitude: %{customdata:.2f}<extra></extra>",
            customdata=amps_valid
        )
        
        push!(traces, trace)
    end
    
    if isempty(traces)
        @warn "No valid traces to plot"
        return PlutoPlotly.plot(scatter(x=[0], y=[0], text=["No valid traces"]))
    end
    
    # Create layout
    layout = Layout(
        title=_title_with_subtitle(title, "velocity=$(vel_min)-$(vel_max) km/s; modes=$(length(modes))";
            font_family=font_family, font_size=font_size),
        xaxis=attr(
            title="Period (s)", type="linear",
            showgrid=true,
            gridcolor="rgba(128,128,128,0.2)",
            zeroline=false,
            titlefont=attr(size=font_size, family=font_family),
            tickfont=attr(size=font_size-2, family=font_family)
        ),
        yaxis=attr(
            title="Group Velocity (km/s)",
            showgrid=true,
            gridcolor="rgba(128,128,128,0.2)",
            zeroline=false,
            range=[vel_min, vel_max],
            titlefont=attr(size=font_size, family=font_family),
            tickfont=attr(size=font_size-2, family=font_family)
        ),
        width=width,
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=_publication_margin(l=80, r=60, t=95, b=145),
        showlegend=true,
        legend=_bottom_legend(font_size; y=-0.24)
    )
    
    return _plotly_plot(traces, layout)
end


function _uc_consistency_rows_from_mft_result(res::MFTResult;
        pair_label::AbstractString="MFT result",
        label::AbstractString="MFT result",
        branch="")
    u_pred = any(isfinite, res.u_predicted_from_phase) ?
        res.u_predicted_from_phase : compute_group_velocity_from_phase(res)
    rows = NamedTuple[]
    for ip in eachindex(res.periods)
        period = Float64(res.periods[ip])
        u_meas = Float64(res.group_velocities[ip])
        u_hat = Float64(u_pred[ip])
        c = Float64(res.phase_velocities[ip])
        quality = Float64(res.quality_factors[ip])
        isfinite(period) && period > 0.0 || continue
        isfinite(u_meas) && u_meas > 0.0 || continue
        isfinite(u_hat) && u_hat > 0.0 || continue
        isfinite(c) && c > 0.0 || continue
        relerr = abs(u_meas - u_hat) / max(abs(u_hat), eps(Float64))
        selected_phase_branch = ip <= length(res.selected_phase_branches) ?
            Int(res.selected_phase_branches[ip]) : 0
        phase_suspect = ip <= length(res.phase_suspect) ? Bool(res.phase_suspect[ip]) : false
        push!(rows, (; pair_label=String(pair_label), label=String(label), branch,
            period,
            group_velocity=u_meas,
            phase_velocity=c,
            predicted_group_velocity=u_hat,
            relative_error=relerr,
            relative_agreement=1.0 / (1.0 + relerr),
            quality,
            selected_phase_branch,
            phase_suspect))
    end
    sort(rows, by=r -> (String(r.label), String(r.branch), r.period))
end

# ╔═╡ b7000002-0000-0000-0000-000000000001
function plot_uc_consistency_comparison(rows;
        pair_label=nothing,
        overlay_rows=NamedTuple[],
        codebook_rows=overlay_rows,
        velocity_range=nothing,
        title::AbstractString="U-c Consistency",
        subtitle=nothing,
        error_cap::Float64=0.5,
        width::Int=1100,
        height::Int=680)
    base_rows = isnothing(pair_label) ? collect(rows) :
        [r for r in rows if hasproperty(r, :pair_label) && String(r.pair_label) == String(pair_label)]
    overlay_all = isnothing(pair_label) ? collect(codebook_rows) :
        [r for r in codebook_rows if hasproperty(r, :pair_label) && String(r.pair_label) == String(pair_label)]
    isempty(base_rows) && isempty(overlay_all) &&
        return PlutoPlotly.plot(scatter(x=[0.0], y=[0.0], text=["No finite U-c rows"]))

    branches = ["causal", "acausal"]
    traces = AbstractTrace[]
    colors = Dict("causal" => "#1f77b4", "acausal" => "#d62728")
    symbols = Dict("causal" => "circle", "acausal" => "diamond")
    error_scales = Dict("causal" => "Blues", "acausal" => "Reds")

    for branch in branches
        b = sort([r for r in base_rows
            if hasproperty(r, :branch) && String(r.branch) == branch], by=r -> r.period)
        isempty(b) && continue
        selected_branch_text(r) = hasproperty(r, :selected_phase_branch) ?
            "<br>phase branch=$(r.selected_phase_branch)" : ""
        hover = ["$(branch)<br>T=$(round(r.period; digits=3)) s<br>U=$(round(r.group_velocity; digits=3)) km/s<br>c=$(round(r.phase_velocity; digits=3)) km/s<br>U_pred=$(round(r.predicted_group_velocity; digits=3)) km/s<br>rel err=$(round(r.relative_error; digits=3))$(selected_branch_text(r))<br>quality=$(round(r.quality; digits=2))" for r in b]
        push!(traces, scatter(
            x=[r.period for r in b], y=[r.group_velocity for r in b],
            yaxis="y", mode="markers+lines", name="$(branch) U", legendgroup=branch,
            marker=attr(size=8, color=[min(Float64(r.relative_error), error_cap) for r in b],
                cmin=0.0, cmax=error_cap, colorscale=error_scales[branch], showscale=false,
                symbol=symbols[branch], line=attr(color="black", width=0.8)),
            line=attr(color=colors[branch], width=1.5), text=hover, hoverinfo="text"))
        push!(traces, scatter(
            x=[r.period for r in b], y=[r.predicted_group_velocity for r in b],
            yaxis="y", mode="lines+markers", name="$(branch) U from c(T)", legendgroup=branch,
            marker=attr(size=5, color=colors[branch], symbol="circle-open"),
            line=attr(color=colors[branch], width=1.8, dash="dash"),
            hovertemplate="$(branch) U_pred<br>T=%{x:.3f} s<br>U_pred=%{y:.3f} km/s<extra></extra>"))
        push!(traces, scatter(
            x=[r.period for r in b], y=[r.phase_velocity for r in b],
            yaxis="y", mode="lines+markers", name="$(branch) c", legendgroup=branch,
            marker=attr(size=4.5, color=colors[branch], symbol="square-open"),
            line=attr(color=colors[branch], width=1.3, dash="dot"),
            hovertemplate="$(branch) phase velocity<br>T=%{x:.3f} s<br>c=%{y:.3f} km/s<extra></extra>"))
        push!(traces, scatter(
            x=[r.period for r in b], y=[min(Float64(r.relative_error), error_cap) for r in b],
            xaxis="x2", yaxis="y2", mode="lines+markers", name="$(branch) relative error",
            legendgroup=branch, showlegend=false,
            marker=attr(size=7, color=[min(Float64(r.relative_error), error_cap) for r in b],
                cmin=0.0, cmax=error_cap, colorscale=error_scales[branch], showscale=false,
                symbol=symbols[branch], line=attr(color="black", width=0.5)),
            line=attr(color=colors[branch], width=1.4),
            text=["$(branch)<br>T=$(round(r.period; digits=3)) s<br>relative error=$(round(r.relative_error; digits=3))$(r.relative_error > error_cap ? "<br>plotted at $(error_cap) cap" : "")" for r in b],
            hoverinfo="text"))
    end

    unlabeled_base = sort([r for r in base_rows if !hasproperty(r, :branch) || isempty(String(r.branch))],
        by=r -> (hasproperty(r, :label) ? String(r.label) : "", r.period))
    append!(overlay_all, unlabeled_base)
    overlay_labels = unique(String(hasproperty(r, :label) ? r.label : "overlay") for r in overlay_all)
    overlay_palette = ["#111111", "#7b3294", "#008837", "#e66101", "#5e3c99", "#4d4d4d"]
    overlay_groups = NamedTuple[]
    for label in overlay_labels
        label_rows = [r for r in overlay_all if String(hasproperty(r, :label) ? r.label : "overlay") == label]
        if any(r -> hasproperty(r, :branch) && !isempty(String(r.branch)), label_rows)
            for branch in branches
                rb = sort([r for r in label_rows
                    if hasproperty(r, :branch) && String(r.branch) == branch], by=r -> r.period)
                isempty(rb) || push!(overlay_groups, (; label, branch, rows=rb))
            end
        else
            push!(overlay_groups, (; label, branch="", rows=sort(label_rows, by=r -> r.period)))
        end
    end
    for (igroup, group) in enumerate(overlay_groups)
        rs = group.rows
        isempty(rs) && continue
        branch = String(group.branch)
        branch_suffix = isempty(branch) ? "" : " $(branch)"
        line_color = overlay_palette[mod1(igroup, length(overlay_palette))]
        marker_color = branch == "acausal" ? "#7b3294" : "#b8860b"
        label = String(group.label)
        legend_prefix = isempty(branch) ? label : "$(label) $(branch)"
        hover = ["$(label)$(isempty(branch) ? "" : "<br>$(branch)")<br>T=$(round(r.period; digits=3)) s<br>U=$(round(r.group_velocity; digits=3)) km/s<br>c=$(round(r.phase_velocity; digits=3)) km/s<br>U_pred=$(round(r.predicted_group_velocity; digits=3)) km/s<br>rel err=$(round(r.relative_error; digits=3))<br>quality=$(round(r.quality; digits=2))" for r in rs]
        push!(traces, scatter(
            x=[r.period for r in rs], y=[r.group_velocity for r in rs],
            yaxis="y", mode="markers+lines", name="$(legend_prefix) U",
            legendgroup="overlay $(igroup)",
            marker=attr(size=8.5, color=marker_color, symbol="star",
                line=attr(color=line_color, width=0.8)),
            line=attr(color=line_color, width=2.0), text=hover, hoverinfo="text"))
        push!(traces, scatter(
            x=[r.period for r in rs], y=[r.predicted_group_velocity for r in rs],
            yaxis="y", mode="lines+markers", name="$(legend_prefix) U(c)",
            legendgroup="overlay $(igroup)",
            marker=attr(size=5, color=line_color, symbol="star-open"),
            line=attr(color=line_color, width=1.8, dash="dash"),
            hovertemplate="overlay$(branch_suffix) U_pred<br>T=%{x:.3f} s<br>U_pred=%{y:.3f} km/s<extra></extra>"))
        push!(traces, scatter(
            x=[r.period for r in rs], y=[r.phase_velocity for r in rs],
            yaxis="y", mode="lines+markers", name="$(legend_prefix) c",
            legendgroup="overlay $(igroup)",
            marker=attr(size=4.8, color=line_color, symbol="square-open"),
            line=attr(color=line_color, width=1.3, dash="dot"),
            hovertemplate="overlay$(branch_suffix) phase velocity<br>T=%{x:.3f} s<br>c=%{y:.3f} km/s<extra></extra>"))
        push!(traces, scatter(
            x=[r.period for r in rs], y=[min(Float64(r.relative_error), error_cap) for r in rs],
            xaxis="x2", yaxis="y2", mode="lines+markers",
            name="overlay$(branch_suffix) relative error", legendgroup="overlay $(igroup)",
            showlegend=false,
            marker=attr(size=7.5, color=marker_color, symbol="star",
                line=attr(color=line_color, width=0.6)),
            line=attr(color=line_color, width=1.5),
            text=["$(label)$(isempty(branch) ? "" : "<br>$(branch)")<br>T=$(round(r.period; digits=3)) s<br>relative error=$(round(r.relative_error; digits=3))$(r.relative_error > error_cap ? "<br>plotted at $(error_cap) cap" : "")" for r in rs],
            hoverinfo="text"))
    end

    y_range = _velocity_range_or_default(velocity_range)
    subtitle_parts = String[]
    isnothing(subtitle) || isempty(String(subtitle)) || push!(subtitle_parts, String(subtitle))
    isnothing(pair_label) || push!(subtitle_parts, "pair=$(pair_label)")
    isempty(overlay_labels) || push!(subtitle_parts, "overlays=$(_join_compact(overlay_labels; max_items=3))")
    push!(subtitle_parts, "rows=$(length(base_rows) + length(overlay_all))")
    push!(subtitle_parts, "error cap=$(error_cap)")
    push!(subtitle_parts, "velocity=$(y_range[1])-$(y_range[2]) km/s")
    plot_subtitle = join(subtitle_parts, " | ")
    _plotly_plot(traces, Layout(
        title=_title_with_subtitle(title, plot_subtitle; font_family="Arial, sans-serif", font_size=15),
        xaxis=attr(type="log", domain=[0.0, 0.86], anchor="y",
            showticklabels=false, showgrid=true, gridcolor="rgba(0,0,0,0.10)", zeroline=false),
        xaxis2=attr(title="Period (s)", type="log", domain=[0.0, 0.86], anchor="y2",
            matches="x", showgrid=true, gridcolor="rgba(0,0,0,0.10)", zeroline=false),
        yaxis=attr(title="Velocity (km/s)", range=y_range, domain=[0.33, 1.0],
            showgrid=true, gridcolor="rgba(0,0,0,0.10)", zeroline=false),
        yaxis2=attr(title="relative |U - U(c)| / U(c) (capped at $(error_cap))",
            range=[0.0, error_cap], domain=[0.0, 0.23],
            showgrid=true, gridcolor="rgba(0,0,0,0.10)", zeroline=false),
        width=width, height=height,
        plot_bgcolor="white", paper_bgcolor="white",
        margin=_publication_margin(l=78, r=60, t=105, b=175),
        legend=_bottom_legend(13; y=-0.24)))
end

# ╔═╡ b7000003-0000-0000-0000-000000000001
function plot_uc_consistency_comparison(res::MFTResult;
        label::AbstractString="MFT result",
        pair_label::AbstractString=label,
        kwargs...)
    rows = _uc_consistency_rows_from_mft_result(res; pair_label=pair_label, label=label)
    plot_uc_consistency_comparison(rows; pair_label=pair_label, kwargs...)
end

# ╔═╡ b7000004-0000-0000-0000-000000000001
function plot_uc_consistency_comparison(results::AbstractVector{<:MFTResult};
        labels=nothing,
        pair_label::AbstractString="MFT results",
        kwargs...)
    labels0 = isnothing(labels) ? ["MFT result $(i)" for i in eachindex(results)] : String.(labels)
    rows = NamedTuple[]
    for (i, res) in enumerate(results)
        label = i <= length(labels0) ? labels0[i] : "MFT result $(i)"
        append!(rows, _uc_consistency_rows_from_mft_result(res;
            pair_label=pair_label, label=label))
    end
    plot_uc_consistency_comparison(rows; pair_label=pair_label, kwargs...)
end

function plot_group_phase_consistency(res::MFTResult;
                                      width=900,
                                      height=550,
                                      font_family="Arial, sans-serif",
                                      font_size=18,
                                      title="Group and Phase Consistency",
                                      velocity_range=nothing,
                                      show_suspect=true)
    u_pred = any(isfinite, res.u_predicted_from_phase) ?
        res.u_predicted_from_phase : compute_group_velocity_from_phase(res)
    valid_group = findall(i -> isfinite(res.group_velocities[i]), eachindex(res.periods))
    valid_pred = findall(i -> isfinite(u_pred[i]), eachindex(res.periods))
    isempty(valid_group) && isempty(valid_pred) &&
        return PlutoPlotly.plot(scatter(x=[0.0], y=[0.0], text=["No velocity data"]))

    vals = Float64[]
    append!(vals, res.group_velocities[valid_group])
    append!(vals, u_pred[valid_pred])
    vel_min, vel_max = _velocity_range_or_default(velocity_range)

    traces = [scatter()]
    if !isempty(valid_group)
        push!(traces, scatter(
            x=res.periods[valid_group], y=res.group_velocities[valid_group],
            mode="lines+markers",
            marker=attr(size=8, color="black"),
            line=attr(color="black", width=2),
            name="Measured U",
            hovertemplate="Period: %{x:.2f} s<br>U: %{y:.3f} km/s<extra></extra>"
        ))
    end
    if !isempty(valid_pred)
        push!(traces, scatter(
            x=res.periods[valid_pred], y=u_pred[valid_pred],
            mode="lines+markers",
            marker=attr(size=8, color="firebrick", symbol="diamond"),
            line=attr(color="firebrick", width=2, dash="dash"),
            name="U from c(T)",
            hovertemplate="Period: %{x:.2f} s<br>U_pred: %{y:.3f} km/s<extra></extra>"
        ))
    end
    if show_suspect
        suspect_idx = findall(i -> res.phase_suspect[i] && isfinite(u_pred[i]), eachindex(res.periods))
        if !isempty(suspect_idx)
            push!(traces, scatter(
                x=res.periods[suspect_idx], y=u_pred[suspect_idx],
                mode="markers",
                marker=attr(size=13, color="orange", symbol="x"),
                name="Phase suspect",
                hovertemplate="Period: %{x:.2f} s<br>Suspect U_pred: %{y:.3f} km/s<extra></extra>"
            ))
        end
    end

    layout = Layout(
        title=_title_with_subtitle(title, "velocity=$(vel_min)-$(vel_max) km/s; periods=$(length(unique(vcat(valid_group, valid_pred))))";
            font_family=font_family, font_size=font_size),
        xaxis=attr(title="Period (s)", type="linear", showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        yaxis=attr(title="Group Velocity (km/s)", range=[vel_min, vel_max], showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        width=width, height=height,
        plot_bgcolor="white", paper_bgcolor="white",
        margin=_publication_margin(l=80, r=60, t=95, b=135),
        showlegend=true,
        legend=_bottom_legend(font_size; y=-0.24)
    )
    return _plotly_plot(traces, layout)
end
