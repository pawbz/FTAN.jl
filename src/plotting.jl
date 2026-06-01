_plotly_traces(traces::AbstractVector) = AbstractTrace[trace for trace in traces]
_plotly_plot(traces::AbstractVector, layout; kwargs...) =
    PlutoPlotly.plot(_plotly_traces(traces), layout; kwargs...)

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
        title = attr(
            text = title,
            font = attr(size = font_size + 2, family = font_family)
        ),
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
                         velocity_range=nothing,  # (vmin, vmax) or nothing for auto
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
    if isnothing(velocity_range)
        vel_min = 0.9 * minimum(all_group_vals)
        vel_max = 1.1 * maximum(all_group_vals)
    else
        vel_min, vel_max = velocity_range
    end

    layout = Layout(
        title = attr(
            text = title,
            font = attr(size = font_size + 2, family = font_family)
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
        margin = attr(l=80, r=120, t=80, b=80),
        showlegend = true
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
        color = palette[fld(i+1, 2)]

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
                         velocity_range=nothing,  # (vmin, vmax) or nothing for auto
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
    if isnothing(velocity_range)
        vel_min = 0.9 * minimum(group_vel_valid)
        vel_max = 1.1 * maximum(group_vel_valid)
    else
        vel_min, vel_max = velocity_range
    end
    
    layout = Layout(
        title = attr(
            text = title,
            font = attr(size = font_size + 2, family = font_family)
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
        margin = attr(l=80, r=120, t=80, b=80),
        showlegend = true
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

    if isnothing(velocity_range)
        vel_min = 0.9 * minimum(all_vel_vals)
        vel_max = 1.1 * maximum(all_vel_vals)
    else
        vel_min, vel_max = velocity_range
    end

    layout = Layout(
        title = attr(
            text = title,
            font = attr(size = font_size + 2, family = font_family)
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
        margin = attr(l=80, r=120, t=80, b=80),
        showlegend = true
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
                                  show_amplitudes=true)
    
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
    
    vel_min = 0.9 * minimum(all_vels)
    vel_max = 1.1 * maximum(all_vels)
    
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
        title=attr(
            text=title,
            font=attr(size=font_size+2, family=font_family)
        ),
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
        margin=attr(l=80, r=120, t=80, b=80),
        showlegend=true,
        legend=attr(
            x=1.02,
            y=1.0,
            font=attr(size=font_size-2),
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="gray",
            borderwidth=1
        )
    )
    
    return _plotly_plot(traces, layout)
end

function plot_branch_comparison(result::BranchAnalysisResult; 
                               title="Causal vs Acausal Dispersion Comparison",
                               colorscale="Viridis",
                               width=1200,
                               height=600,
                               font_family="Arial, sans-serif",
                               font_size=14)
    
    # Collect all valid velocities for axis scaling
    all_vels = Float64[]
    for modes in [result.causal_modes, result.acausal_modes]
        for mode in modes
            append!(all_vels, filter(!isnan, mode.group_velocities))
        end
    end
    
    if isempty(all_vels)
        @warn "No valid velocity measurements in branch results"
        return PlutoPlotly.plot(scatter(x=[0], y=[0], text=["No data"]))
    end
    
    vel_min = 0.9 * minimum(all_vels)
    vel_max = 1.1 * maximum(all_vels)
    
    # Generate colors for modes
    n_max_modes = max(length(result.causal_modes), length(result.acausal_modes))
    colors = _sample_scheme_colors(colorscale, max(n_max_modes, 2))
    
    # Build causal branch traces
    traces_causal = [scatter()]
    for (mode_idx, mode) in enumerate(result.causal_modes)
        valid_idx = findall(!isnan, mode.group_velocities)
        if !isempty(valid_idx)
            periods_valid = mode.periods[valid_idx]
            vels_valid = mode.group_velocities[valid_idx]
            amps_valid = mode.peak_amplitudes[valid_idx]
            
            amp_min = minimum(amps_valid)
            amp_max = maximum(amps_valid)
            marker_sizes = (amp_max > amp_min) ? 
                8 .+ 8 .* (amps_valid .- amp_min) ./ (amp_max - amp_min) : 
                fill(10.0, length(amps_valid))
            
            trace = scatter(
                x=periods_valid, y=vels_valid,
                xaxis="x1", yaxis="y1",
                mode="markers+lines",
                marker=attr(size=marker_sizes, color=colors[mode_idx], opacity=0.7, 
                           line=attr(color="white", width=1)),
                line=attr(color=colors[mode_idx], width=2),
                name="C-Mode $(mode.mode_index)",
                legendgroup="mode_$(mode_idx)",
                hovertemplate="<b>Causal Mode $(mode.mode_index)</b><br>" *
                             "Period: %{x:.2f} s<br>Velocity: %{y:.3f} km/s<extra></extra>",
                showlegend=true
            )
            push!(traces_causal, trace)
        end
    end
    
    # Build acausal branch traces
    traces_acausal = [scatter()]
    for (mode_idx, mode) in enumerate(result.acausal_modes)
        valid_idx = findall(!isnan, mode.group_velocities)
        if !isempty(valid_idx)
            periods_valid = mode.periods[valid_idx]
            vels_valid = mode.group_velocities[valid_idx]
            amps_valid = mode.peak_amplitudes[valid_idx]
            
            amp_min = minimum(amps_valid)
            amp_max = maximum(amps_valid)
            marker_sizes = (amp_max > amp_min) ? 
                8 .+ 8 .* (amps_valid .- amp_min) ./ (amp_max - amp_min) : 
                fill(10.0, length(amps_valid))
            
            trace = scatter(
                x=periods_valid, y=vels_valid,
                xaxis="x2", yaxis="y1",
                mode="markers+lines",
                marker=attr(size=marker_sizes, color=colors[mode_idx], opacity=0.7, 
                           line=attr(color="white", width=1)),
                line=attr(color=colors[mode_idx], width=2, dash="dot"),
                name="A-Mode $(mode.mode_index)",
                legendgroup="mode_$(mode_idx)",
                hovertemplate="<b>Acausal Mode $(mode.mode_index)</b><br>" *
                             "Period: %{x:.2f} s<br>Velocity: %{y:.3f} km/s<extra></extra>",
                showlegend=true
            )
            push!(traces_acausal, trace)
        end
    end
    
    all_traces = vcat(traces_causal, traces_acausal)
    
    layout = Layout(
        title=attr(text=title, font=attr(size=font_size+2, family=font_family)),
        xaxis1=attr(title="Period (s)", type="linear", domain=[0, 0.45], showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        xaxis2=attr(title="Period (s)", type="linear", domain=[0.55, 1], showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        yaxis1=attr(title="Group Velocity (km/s)", range=[vel_min, vel_max], 
                   showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        width=width, height=height,
        plot_bgcolor="white", paper_bgcolor="white",
        margin=attr(l=80, r=80, t=80, b=80),
        showlegend=true,
        legend=attr(x=1.02, y=1.0, font=attr(size=font_size-2))
    )
    
    return _plotly_plot(all_traces, layout)
end

# ╔═╡ 4c915e4e-3cdc-11f1-9235-33f20a13a658
"""
    plot_branch_correlation(result::BranchAnalysisResult;
                           title="Causal-Acausal Envelope Correlation",
                           width=900, height=500,
                           font_family="Arial, sans-serif",
                           font_size=14,
                           correlation_threshold=0.85)

Plot the zero-lag correlation between causal and acausal branches across periods.

**Interpretation:**
- Correlation > threshold (default 0.85): Consistent, reliable arrivals on both branches
- Correlation < threshold: Possible mode/branch inconsistency; manually check
- NaN: Unable to compute (e.g., zero-amplitude envelope)

# Arguments
- `result` : BranchAnalysisResult from analyze_causal_acausal_branches
- `title` : Plot title
- `width`, `height` : Plot dimensions
- `correlation_threshold` : Horizontal threshold line to highlight acceptable correlations (default: 0.85)

# Returns
- PlutoPlotly plot object
"""
function plot_branch_correlation(result::BranchAnalysisResult;
                                title="Causal-Acausal Correlation",
                                width=900,
                                height=500,
                                font_family="Arial, sans-serif",
                                font_size=14,
                                correlation_threshold=0.85,
                                wavelength_ref_velocity::Union{Nothing,Real}=nothing,
                                wavelength_fraction::Union{Nothing,Real}=nothing)
    
    # Filter out NaN values for plotting
    wavelength_idx = Set(_wavelength_valid_indices(result.periods, result.distance;
                                                    wavelength_ref_velocity=wavelength_ref_velocity,
                                                    wavelength_fraction=wavelength_fraction))
    valid_idx = findall(ip -> !isnan(result.branch_correlation[ip]) && ip in wavelength_idx,
                        eachindex(result.branch_correlation))
    
    if isempty(valid_idx)
        @warn "No valid correlation values found"
        return PlutoPlotly.plot(scatter(x=[0], y=[0], text=["No correlation data"]))
    end
    
    periods_valid = result.periods[valid_idx]
    corr_valid = result.branch_correlation[valid_idx]
    
    # Color code by threshold
    colors = [c >= correlation_threshold ? "green" : "orange" for c in corr_valid]
    
    traces = [
        scatter(
            x=periods_valid,
            y=corr_valid,
            mode="markers+lines",
            marker=attr(size=10, color=colors, line=attr(color="darkgray", width=1)),
            line=attr(color="lightgray", width=2),
            name="Correlation",
            hovertemplate="Period: %{x:.2f} s<br>Correlation: %{y:.3f}<extra></extra>"
        ),
        scatter(
            x=[minimum(periods_valid), maximum(periods_valid)],
            y=[correlation_threshold, correlation_threshold],
            mode="lines",
            line=attr(color="red", width=2, dash="dash"),
            name="Threshold ($correlation_threshold)",
            hoverinfo="skip"
        )
    ]
    
    layout = Layout(
        title=attr(text=title, font=attr(size=font_size+2, family=font_family)),
        xaxis=attr(title="Period (s)", type="linear", range=[minimum(periods_valid), maximum(periods_valid)], showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        yaxis=attr(title="Correlation Coefficient", range=[-0.1, 1.05], 
                  showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        width=width, height=height,
        plot_bgcolor="white", paper_bgcolor="white",
        margin=attr(l=80, r=80, t=80, b=80),
        showlegend=true,
        legend=attr(x=0.02, y=0.98, font=attr(size=font_size-2),
                   bgcolor="rgba(255,255,255,0.8)", borderwidth=1)
    )
    
    return _plotly_plot(traces, layout)
end

function plot_phase_regression_fit(fit::PhaseVelocityPeriodFit;
        width::Int=950,
        height::Int=560,
        title::String="Robust phase-velocity regression")
    isempty(fit.observations) &&
        return PlutoPlotly.plot(scatter(x=[0.0], y=[0.0], text=["No phase observations"]))
    xs = [o.distance for o in fit.observations]
    status = [fit.inlier[i] ? "inlier" : fit.azimuth_outlier[i] ? "azimuth outlier" : "outlier"
              for i in eachindex(fit.observations)]
    colors = [fit.inlier[i] ? "royalblue" : fit.azimuth_outlier[i] ? "darkorange" : "crimson"
              for i in eachindex(fit.observations)]
    line_x = collect(range(minimum(xs), maximum(xs), length=100))
    traces = AbstractTrace[
        scatter(
            x=xs,
            y=fit.unwrapped_phase,
            mode="markers",
            name="observations",
            marker=attr(size=10, color=colors),
            text=[let az_label = isfinite(o.azimuth_deg) ? string(round(o.azimuth_deg; digits=1)) : "NA"
                    "$(o.pair_label)<br>N=$(fit.cycle_numbers[i])<br>$(status[i])<br>res=$(round(fit.residuals[i]; digits=3)) rad<br>sigma_phi=$(round(o.sigma_phi; digits=3))<br>q=$(round(o.quality; digits=2))<br>az=$(az_label)"
                  end for (i, o) in enumerate(fit.observations)],
            hovertemplate="%{text}<br>d=%{x:.1f} km<br>Φ=%{y:.3f} rad<extra></extra>",
        ),
        scatter(
            x=line_x,
            y=fit.slope .* line_x,
            mode="lines",
            name="WLS fit c=$(isfinite(fit.phase_velocity) ? round(fit.phase_velocity; digits=3) : NaN) km/s",
            line=attr(color="black", width=2.5),
        ),
    ]
    return _plotly_plot(traces, Layout(
        title=title,
        xaxis=attr(title="Distance d (km)"),
        yaxis=attr(title="Unwrapped phase Φ (rad)"),
        width=width,
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
    ))
end

function plot_branch_correlation(result::BranchBatchAnalysisResult;
                                 title="Causal-Acausal Correlation Across Source States",
                                 colorscale="Viridis",
                                 width=1000,
                                 height=550,
                                 font_family="Arial, sans-serif",
                                 font_size=14,
                                 correlation_threshold=0.85,
                                 reference_results::AbstractVector{<:BranchAnalysisResult}=BranchAnalysisResult[],
                                 reference_labels::AbstractVector{<:AbstractString}=String[],
                                 wavelength_ref_velocity::Union{Nothing,Real}=nothing,
                                 wavelength_fraction::Union{Nothing,Real}=nothing)
    nstates = size(result.branch_correlation, 2)
    nstates == 0 && return PlutoPlotly.plot(scatter(x=[0], y=[0], text=["No source states"]))

    colors = _sample_scheme_colors(colorscale, max(nstates, 2))
    traces = [scatter()]
    visible_periods = Float64[]

    for i in 1:nstates
        st = result.state_results[i]
        corr_i = result.branch_correlation[:, i]
        wavelength_idx = Set(_wavelength_valid_indices(result.periods, st.distance;
                                                        wavelength_ref_velocity=wavelength_ref_velocity,
                                                        wavelength_fraction=wavelength_fraction))
        valid_idx = findall(ip -> !isnan(corr_i[ip]) && ip in wavelength_idx, eachindex(corr_i))
        isempty(valid_idx) && continue
        append!(visible_periods, result.periods[valid_idx])

        push!(traces, scatter(
            x = result.periods[valid_idx],
            y = corr_i[valid_idx],
            mode = "lines+markers",
            marker = attr(size = 7, color = colors[i], line = attr(color = "white", width = 0.8)),
            line = attr(color = colors[i], width = 1.8),
            name = result.state_labels[i],
            hovertemplate = "State: $(result.state_labels[i])<br>Period: %{x:.2f} s<br>Correlation: %{y:.3f}<extra></extra>"
        ))
    end

    for (iref, ref) in enumerate(reference_results)
        corr_ref = ref.branch_correlation
        wavelength_idx = Set(_wavelength_valid_indices(ref.periods, ref.distance;
                                                        wavelength_ref_velocity=wavelength_ref_velocity,
                                                        wavelength_fraction=wavelength_fraction))
        valid_idx = findall(ip -> !isnan(corr_ref[ip]) && ip in wavelength_idx, eachindex(corr_ref))
        isempty(valid_idx) && continue
        append!(visible_periods, ref.periods[valid_idx])
        ref_label = iref <= length(reference_labels) ? String(reference_labels[iref]) : "Global average"
        push!(traces, scatter(
            x = ref.periods[valid_idx],
            y = corr_ref[valid_idx],
            mode = "lines+markers",
            marker = attr(size = 8, color = "#222222", symbol = "circle-open", line = attr(color = "#222222", width = 1.5)),
            line = attr(color = "#222222", width = 3.0),
            name = ref_label,
            hovertemplate = "$(ref_label)<br>Period: %{x:.2f} s<br>Correlation: %{y:.3f}<extra></extra>"
        ))
    end

    isempty(visible_periods) && return PlutoPlotly.plot(scatter(x=[0], y=[0], text=["No correlation data"]))

    # Threshold guide
    push!(traces, scatter(
        x = [minimum(visible_periods), maximum(visible_periods)],
        y = [correlation_threshold, correlation_threshold],
        mode = "lines",
        line = attr(color = "red", width = 2, dash = "dash"),
        name = "Threshold ($(correlation_threshold))",
        hoverinfo = "skip"
    ))

    layout = Layout(
        title = attr(text = title, font = attr(size = font_size + 2, family = font_family)),
        xaxis = attr(title = "Period (s)",
            type = "linear", showgrid = true, gridcolor = "rgba(128,128,128,0.2)"),
        yaxis = attr(title = "Correlation Coefficient", range = [-0.1, 1.05], showgrid = true, gridcolor = "rgba(128,128,128,0.2)"),
        width = width,
        height = height,
        plot_bgcolor = "white",
        paper_bgcolor = "white",
        margin = attr(l = 80, r = 80, t = 80, b = 80),
        showlegend = true,
        legend = attr(x = 1.02, y = 1.0, font = attr(size = font_size - 2), bgcolor = "rgba(255,255,255,0.8)", borderwidth = 1)
    )

    return _plotly_plot(traces, layout)
end

function plot_all_highcorr_groupvelocity_picks(result::BranchAnalysisResult;
                                               correlation_threshold::Float64=0.85,
                                               pair_and_average::Bool=false,
                                               velocity_tolerance_fraction::Float64=0.10,
                                               width::Int=1000,
                                               height::Int=600,
                                               font_family::String="Arial, sans-serif",
                                               font_size::Int=14,
                                               title::String="All High-Correlation Group-Velocity Picks",
                                               wavelength_ref_velocity::Union{Nothing,Real}=nothing,
                                               wavelength_fraction::Union{Nothing,Real}=nothing)
    wavelength_idx = Set(_wavelength_valid_indices(result.periods, result.distance;
                                                    wavelength_ref_velocity=wavelength_ref_velocity,
                                                    wavelength_fraction=wavelength_fraction))
    hi_idx = [ip for ip in high_correlation_indices(result; correlation_threshold=correlation_threshold)
             if ip in wavelength_idx]
    if isempty(hi_idx)
        @warn "No periods satisfy the correlation threshold $(correlation_threshold)"
        return PlutoPlotly.plot(scatter(x=[0.0], y=[0.0], text=["No high-correlation periods"]))
    end

    traces = [scatter()]
    all_vels = Float64[]

    if pair_and_average
        matched_by_period = Dict{Int,Vector{Float64}}()
        max_matched = 0

        for ip in hi_idx
            causal_peaks = result.causal_result.all_peaks[ip]
            acausal_peaks = result.acausal_result.all_peaks[ip]
            vavg = _matched_peak_average_velocities(causal_peaks, acausal_peaks, result.distance;
                                                    velocity_tolerance_fraction=velocity_tolerance_fraction)
            isempty(vavg) && continue
            matched_by_period[ip] = vavg
            max_matched = max(max_matched, length(vavg))
            append!(all_vels, vavg)
        end

        for pidx in 1:max_matched
            xp = Float64[]
            yp = Float64[]
            cp = Float64[]

            for ip in hi_idx
                if haskey(matched_by_period, ip)
                    vals = matched_by_period[ip]
                    if pidx <= length(vals)
                        push!(xp, result.periods[ip])
                        push!(yp, vals[pidx])
                        push!(cp, result.branch_correlation[ip])
                    end
                end
            end

            isempty(xp) && continue
            push!(traces, scatter(
                x=xp,
                y=yp,
                mode="markers+lines",
                marker=attr(size=8, color="#2ca02c", symbol="circle", line=attr(color="white", width=0.8)),
                line=attr(color="#2ca02c", width=1.8),
                name="Avg matched peak $(pidx)",
                hovertemplate="Matched Avg<br>Period: %{x:.2f} s<br>v_g: %{y:.3f} km/s<br>Corr: %{customdata:.3f}<extra></extra>",
                customdata=cp
            ))
        end
    else
        branches = [
            ("Causal", result.causal_result, "#1f77b4", "circle"),
            ("Acausal", result.acausal_result, "#d62728", "diamond"),
        ]

        for (branch_name, mft_res, color, marker_symbol) in branches
            max_peaks = maximum(length.(mft_res.all_peaks); init=0)
            max_peaks == 0 && continue

            for pidx in 1:max_peaks
                xp = Float64[]
                yp = Float64[]
                cp = Float64[]

                for ip in hi_idx
                    peaks = mft_res.all_peaks[ip]
                    if pidx <= length(peaks)
                        tpk, _ = peaks[pidx]
                        if isfinite(tpk) && tpk > 0.0
                            vpk = result.distance / tpk
                            if isfinite(vpk) && vpk > 0.0
                                push!(xp, mft_res.periods[ip])
                                push!(yp, vpk)
                                push!(cp, result.branch_correlation[ip])
                                push!(all_vels, vpk)
                            end
                        end
                    end
                end

                isempty(xp) && continue

                push!(traces, scatter(
                    x=xp,
                    y=yp,
                    mode="markers+lines",
                    marker=attr(size=8, color=color, symbol=marker_symbol, line=attr(color="white", width=0.8)),
                    line=attr(color=color, width=1.8),
                    name="$(branch_name) peak $(pidx)",
                    hovertemplate="$(branch_name)<br>Period: %{x:.2f} s<br>v_g: %{y:.3f} km/s<br>Corr: %{customdata:.3f}<extra></extra>",
                    customdata=cp
                ))
            end
        end
    end

    if length(traces) == 1 || isempty(all_vels)
        @warn "No valid peak-derived group-velocity picks found in high-correlation periods"
        return PlutoPlotly.plot(scatter(x=[0.0], y=[0.0], text=["No valid high-correlation picks"]))
    end

    y_min = 0.9 * minimum(all_vels)
    y_max = 1.1 * maximum(all_vels)
    n_hi = length(hi_idx)

    layout = Layout(
        title=attr(text="$(title) (N=$(n_hi), threshold=$(correlation_threshold))", font=attr(size=font_size + 2, family=font_family)),
        xaxis=attr(title="Period (s)", type="linear", range=[minimum(result.periods[hi_idx]), maximum(result.periods[hi_idx])], showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        yaxis=attr(title="Group Velocity (km/s)", range=[y_min, y_max], showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        width=width,
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=attr(l=80, r=80, t=80, b=80),
        showlegend=true,
        legend=attr(x=1.02, y=1.0, font=attr(size=font_size - 2), bgcolor="rgba(255,255,255,0.8)", borderwidth=1)
    )

    return _plotly_plot(traces, layout)
end

# ╔═╡ 58d6ecb6-2c42-4a53-9cbc-b72ebf48356c
"""
    plot_all_highcorr_groupvelocity_picks(result::BranchBatchAnalysisResult;
                                          correlation_threshold=0.85,
                                          colorscale="Viridis",
                                          width=1200,
                                          height=700,
                                          font_family="Arial, sans-serif",
                                          font_size=13,
                                          title="All High-Correlation Group-Velocity Picks Across Source States")

Plot all envelope-peak group-velocity picks for every source state, restricted to
periods where that state's branch correlation is above `correlation_threshold`.
"""
function plot_all_highcorr_groupvelocity_picks(result::BranchBatchAnalysisResult;
                                               correlation_threshold::Float64=0.85,
                                               pair_and_average::Bool=false,
                                               velocity_tolerance_fraction::Float64=0.10,
                                               reference_results::AbstractVector{<:BranchAnalysisResult}=BranchAnalysisResult[],
                                               reference_labels::AbstractVector{<:AbstractString}=String[],
                                               colorscale::String="Viridis",
                                               width::Int=1200,
                                               height::Int=700,
                                               font_family::String="Arial, sans-serif",
                                               font_size::Int=13,
                                               title::String="All High-Correlation Group-Velocity Picks Across Source States",
                                               wavelength_ref_velocity::Union{Nothing,Real}=nothing,
                                               wavelength_fraction::Union{Nothing,Real}=nothing)
    nstates = length(result.state_results)
    nstates == 0 && return PlutoPlotly.plot(scatter(x=[0.0], y=[0.0], text=["No source states"]))

    state_colors = _sample_scheme_colors(colorscale, max(nstates, 2))
    traces = [scatter()]
    all_vels = Float64[]
    total_hi = 0
    visible_periods = Float64[]

    for i in 1:nstates
        st = result.state_results[i]
        label = result.state_labels[i]
        wavelength_idx = Set(_wavelength_valid_indices(st.periods, st.distance;
                                                        wavelength_ref_velocity=wavelength_ref_velocity,
                                                        wavelength_fraction=wavelength_fraction))
        hi_idx = [ip for ip in high_correlation_indices(st; correlation_threshold=correlation_threshold)
                 if ip in wavelength_idx]
        total_hi += length(hi_idx)
        isempty(hi_idx) && continue
        append!(visible_periods, st.periods[hi_idx])

        if pair_and_average
            matched_by_period = Dict{Int,Vector{Float64}}()
            max_matched = 0

            for ip in hi_idx
                causal_peaks = st.causal_result.all_peaks[ip]
                acausal_peaks = st.acausal_result.all_peaks[ip]
                vavg = _matched_peak_average_velocities(causal_peaks, acausal_peaks, st.distance;
                                                        velocity_tolerance_fraction=velocity_tolerance_fraction)
                isempty(vavg) && continue
                matched_by_period[ip] = vavg
                max_matched = max(max_matched, length(vavg))
            end

            for pidx in 1:max_matched
                xp = Float64[]
                yp = Float64[]
                cp = Float64[]

                for ip in hi_idx
                    if haskey(matched_by_period, ip)
                        vals = matched_by_period[ip]
                        if pidx <= length(vals)
                            push!(xp, st.periods[ip])
                            push!(yp, vals[pidx])
                            push!(cp, st.branch_correlation[ip])
                            push!(all_vels, vals[pidx])
                        end
                    end
                end

                isempty(xp) && continue

                push!(traces, scatter(
                    x=xp,
                    y=yp,
                    mode="markers",
                    marker=attr(
                        size=7,
                        color=state_colors[i],
                        opacity=0.82,
                        symbol="circle",
                        line=attr(color="white", width=0.7)
                    ),
                    name="$(label) | Avg matched p$(pidx)",
                    hovertemplate="State: $(label)<br>Type: Avg matched<br>Peak: $(pidx)<br>Period: %{x:.2f} s<br>v_g: %{y:.3f} km/s<br>Corr: %{customdata:.3f}<extra></extra>",
                    customdata=cp
                ))
            end

            continue
        end

        branch_defs = [
            ("Causal", st.causal_result, "circle"),
            ("Acausal", st.acausal_result, "diamond"),
        ]

        for (branch_name, mft_res, marker_symbol) in branch_defs
            max_peaks = maximum(length.(mft_res.all_peaks); init=0)
            max_peaks == 0 && continue

            for pidx in 1:max_peaks
                xp = Float64[]
                yp = Float64[]
                cp = Float64[]

                for ip in hi_idx
                    peaks = mft_res.all_peaks[ip]
                    if pidx <= length(peaks)
                        tpk, _ = peaks[pidx]
                        if isfinite(tpk) && tpk > 0.0
                            vpk = st.distance / tpk
                            if isfinite(vpk) && vpk > 0.0
                                push!(xp, mft_res.periods[ip])
                                push!(yp, vpk)
                                push!(cp, st.branch_correlation[ip])
                                push!(all_vels, vpk)
                            end
                        end
                    end
                end

                isempty(xp) && continue

                push!(traces, scatter(
                    x=xp,
                    y=yp,
                    mode="markers",
                    marker=attr(
                        size=7,
                        color=state_colors[i],
                        opacity=0.78,
                        symbol=marker_symbol,
                        line=attr(color="white", width=0.7)
                    ),
                    name="$(label) | $(branch_name) p$(pidx)",
                    hovertemplate="State: $(label)<br>Branch: $(branch_name)<br>Peak: $(pidx)<br>Period: %{x:.2f} s<br>v_g: %{y:.3f} km/s<br>Corr: %{customdata:.3f}<extra></extra>",
                    customdata=cp
                ))
            end
        end
    end

    reference_styles = [
        ("Causal", "#08306b", "circle-open"),
        ("Acausal", "#7f0000", "diamond-open"),
    ]

    for (iref, ref) in enumerate(reference_results)
        ref_label = iref <= length(reference_labels) ? String(reference_labels[iref]) : "Global avg"
        ref_idx = _wavelength_valid_indices(ref.periods, ref.distance;
                                             wavelength_ref_velocity=wavelength_ref_velocity,
                                             wavelength_fraction=wavelength_fraction)
        isempty(ref_idx) && continue
        append!(visible_periods, ref.periods[ref_idx])
        for (branch_name, mft_res, color, marker_symbol) in [
                (reference_styles[1][1], ref.causal_result, reference_styles[1][2], reference_styles[1][3]),
                (reference_styles[2][1], ref.acausal_result, reference_styles[2][2], reference_styles[2][3])]
            max_peaks = maximum(length.(mft_res.all_peaks); init=0)
            max_peaks == 0 && continue
            for pidx in 1:max_peaks
                xp = Float64[]
                yp = Float64[]
                cp = Float64[]
                for ip in ref_idx
                    peaks = mft_res.all_peaks[ip]
                    if pidx <= length(peaks)
                        tpk, _ = peaks[pidx]
                        if isfinite(tpk) && tpk > 0.0
                            vpk = ref.distance / tpk
                            if isfinite(vpk) && vpk > 0.0
                                push!(xp, mft_res.periods[ip])
                                push!(yp, vpk)
                                push!(cp, ref.branch_correlation[ip])
                                push!(all_vels, vpk)
                            end
                        end
                    end
                end
                isempty(xp) && continue
                push!(traces, scatter(
                    x=xp,
                    y=yp,
                    mode="lines+markers",
                    marker=attr(
                        size=9,
                        color=color,
                        opacity=1.0,
                        symbol=marker_symbol,
                        line=attr(color=color, width=1.4)
                    ),
                    line=attr(color=color, width=2.6, dash="dash"),
                    name="$(ref_label) $(lowercase(branch_name)) peak $(pidx)",
                    hovertemplate="$(ref_label)<br>Branch: $(branch_name)<br>Peak: $(pidx)<br>Period: %{x:.2f} s<br>v_g: %{y:.3f} km/s<br>Corr: %{customdata:.3f}<extra></extra>",
                    customdata=cp
                ))
            end
        end
    end

    if length(traces) == 1 || isempty(all_vels)
        @warn "No valid high-correlation picks across source states"
        return PlutoPlotly.plot(scatter(x=[0.0], y=[0.0], text=["No valid high-correlation picks"]))
    end

    y_min = 0.9 * minimum(all_vels)
    y_max = 1.1 * maximum(all_vels)

    layout = Layout(
        title=attr(text="$(title) (threshold=$(correlation_threshold), high-corr samples=$(total_hi))", font=attr(size=font_size + 2, family=font_family)),
        xaxis=attr(title="Period (s)", type="linear", range=[minimum(visible_periods), maximum(visible_periods)], showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        yaxis=attr(title="Group Velocity (km/s)", range=[y_min, y_max], showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        width=width,
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=attr(l=80, r=120, t=80, b=80),
        showlegend=true,
        legend=attr(x=1.02, y=1.0, font=attr(size=font_size - 2), bgcolor="rgba(255,255,255,0.8)", borderwidth=1)
    )

    return _plotly_plot(traces, layout)
end

function _resolve_period_index(periods::AbstractVector{<:Real};
                               period::Union{Nothing,Real}=nothing,
                               period_index::Union{Nothing,Integer}=nothing)
    if !isnothing(period) && !isnothing(period_index)
        throw(ArgumentError("Provide either `period` or `period_index`, not both"))
    end

    if !isnothing(period_index)
        idx = Int(period_index)
        (idx < 1 || idx > length(periods)) && throw(BoundsError(periods, idx))
        return idx, Float64(periods[idx])
    end

    if isnothing(period)
        idx = max(1, Int(cld(length(periods), 2)))
        return idx, Float64(periods[idx])
    end

    idx = argmin(abs.(Float64.(periods) .- Float64(period)))
    return idx, Float64(periods[idx])
end

# ╔═╡ 88284a1d-00ff-4f69-ac6e-a39dbec70863
"""
    plot_filtered_traces_by_period(result::BranchBatchAnalysisResult;
                                   period=nothing,
                                   period_index=nothing,
                                   correlation_threshold=nothing,
                                   normalize_each=true,
                                   scale=0.7,
                                   spacing=2.2,
                                   colorscale="Viridis",
                                   width=1000,
                                   height=900,
                                   font_family="Arial, sans-serif",
                                   font_size=12,
                                   title="Filtered Traces Across Source States")

Plot filtered causal+acausal traces for all source states at one selected period.

For each source state, the displayed trace is:
- negative lag side: reversed acausal filtered trace
- positive lag side: causal filtered trace

Selection can be provided as a period value (`period`) or index (`period_index`).
If `correlation_threshold` is provided, only source states with correlation
at the selected period above threshold are shown.
"""
function plot_filtered_traces_by_period(result::BranchBatchAnalysisResult;
                                        period::Union{Nothing,Real}=nothing,
                                        period_index::Union{Nothing,Integer}=nothing,
                                        correlation_threshold::Union{Nothing,Float64}=nothing,
                                        normalize_each::Bool=true,
                                        scale::Float64=0.7,
                                        spacing::Float64=2.2,
                                        colorscale::String="Viridis",
                                        width::Int=1000,
                                        height::Int=900,
                                        font_family::String="Arial, sans-serif",
                                        font_size::Int=12,
                                        title::String="Filtered Traces Across Source States")
    nstates = length(result.state_results)
    nstates == 0 && return PlutoPlotly.plot(scatter(x=[0.0], y=[0.0], text=["No source states"]))

    idx, period_used = _resolve_period_index(result.periods; period=period, period_index=period_index)
    state_colors = _sample_scheme_colors(colorscale, max(nstates, 2))

    traces = [scatter()]
    plotted = 0
    max_abs_y = 0.0

    for i in 1:nstates
        st = result.state_results[i]
        label = result.state_labels[i]
        corr = st.branch_correlation[idx]

        if !isnothing(correlation_threshold)
            if !(isfinite(corr) && corr >= correlation_threshold)
                continue
            end
        end

        if !mft_has_full_storage(st.acausal_result) || !mft_has_full_storage(st.causal_result)
            continue
        end
        ac = st.acausal_result.filtered_traces[:, idx]
        ca = st.causal_result.filtered_traces[:, idx]
        if isempty(ac) || isempty(ca)
            continue
        end

        t_neg = -reverse(st.acausal_result.time)
        t_pos = st.causal_result.time
        t_full = vcat(t_neg, t_pos)
        x_full = vcat(reverse(ac), ca)

        if normalize_each
            amp = maximum(abs, x_full)
            x_full = amp > 0 ? x_full ./ amp : x_full
        end

        plotted += 1
        offset = (plotted - 1) * spacing
        y = scale .* x_full .+ offset
        max_abs_y = max(max_abs_y, maximum(abs, y))

        push!(traces, scatter(
            x=t_full,
            y=y,
            mode="lines",
            line=attr(color=state_colors[i], width=1.7),
            name=label,
            hovertemplate="State: $(label)<br>Period: $(round(period_used; digits=3)) s<br>Corr: $(isfinite(corr) ? round(corr; digits=3) : NaN)<br>Lag: %{x:.2f} s<extra></extra>"
        ))
    end

    if plotted == 0
        msg = isnothing(correlation_threshold) ? "No valid traces at selected period" : "No states pass threshold at selected period"
        @warn msg
        return PlutoPlotly.plot(scatter(x=[0.0], y=[0.0], text=[msg]))
    end

    subtitle = "period=$(round(period_used; digits=3)) s (index=$(idx))"
    if !isnothing(correlation_threshold)
        subtitle *= ", threshold=$(correlation_threshold)"
    end

    layout = Layout(
        title=attr(text="$(title) [$(subtitle)]", font=attr(size=font_size + 2, family=font_family)),
        xaxis=attr(title="Lag Time (s)", showgrid=true, gridcolor="rgba(128,128,128,0.2)", zeroline=true, zerolinecolor="rgba(100,100,100,0.45)"),
        yaxis=attr(title="Source State (offset filtered traces)", showgrid=false, zeroline=false, range=[-scale * 1.3, max_abs_y + scale]),
        width=width,
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=attr(l=85, r=120, t=85, b=75),
        showlegend=true,
        legend=attr(x=1.02, y=1.0, font=attr(size=font_size - 2), bgcolor="rgba(255,255,255,0.85)", borderwidth=1)
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
    overlay_title = isempty(overlay_labels) ? "" :
        " | overlay: $(join(first(overlay_labels, min(2, length(overlay_labels))), " + "))$(length(overlay_labels) > 2 ? " + ..." : "")"
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
        hover = ["$(label)$(isempty(branch) ? "" : "<br>$(branch)")<br>T=$(round(r.period; digits=3)) s<br>U=$(round(r.group_velocity; digits=3)) km/s<br>c=$(round(r.phase_velocity; digits=3)) km/s<br>U_pred=$(round(r.predicted_group_velocity; digits=3)) km/s<br>rel err=$(round(r.relative_error; digits=3))<br>quality=$(round(r.quality; digits=2))" for r in rs]
        push!(traces, scatter(
            x=[r.period for r in rs], y=[r.group_velocity for r in rs],
            yaxis="y", mode="markers+lines", name="overlay$(branch_suffix) U",
            legendgroup="overlay $(igroup)",
            marker=attr(size=8.5, color=marker_color, symbol="star",
                line=attr(color=line_color, width=0.8)),
            line=attr(color=line_color, width=2.0), text=hover, hoverinfo="text"))
        push!(traces, scatter(
            x=[r.period for r in rs], y=[r.predicted_group_velocity for r in rs],
            yaxis="y", mode="lines+markers", name="overlay$(branch_suffix) U from c(T)",
            legendgroup="overlay $(igroup)",
            marker=attr(size=5, color=line_color, symbol="star-open"),
            line=attr(color=line_color, width=1.8, dash="dash"),
            hovertemplate="overlay$(branch_suffix) U_pred<br>T=%{x:.3f} s<br>U_pred=%{y:.3f} km/s<extra></extra>"))
        push!(traces, scatter(
            x=[r.period for r in rs], y=[r.phase_velocity for r in rs],
            yaxis="y", mode="lines+markers", name="overlay$(branch_suffix) c",
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

    all_vels = Float64[]
    for r in base_rows
        append!(all_vels, [Float64(r.group_velocity), Float64(r.phase_velocity), Float64(r.predicted_group_velocity)])
    end
    for r in overlay_all
        append!(all_vels, [Float64(r.group_velocity), Float64(r.phase_velocity), Float64(r.predicted_group_velocity)])
    end
    filter!(v -> isfinite(v) && v > 0.0, all_vels)
    y_range = isnothing(velocity_range) ?
        (isempty(all_vels) ? nothing : [0.9 * minimum(all_vels), 1.1 * maximum(all_vels)]) :
        [Float64(velocity_range[1]), Float64(velocity_range[2])]
    title_pair = isnothing(pair_label) ? "" : ": $(pair_label)"
    _plotly_plot(traces, Layout(
        title="$(title)$(title_pair)$(overlay_title)",
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
        margin=attr(l=78, r=145, t=78, b=72),
        legend=attr(orientation="h", x=0.0, y=1.10,
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor="rgba(0,0,0,0.15)", borderwidth=1)))
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
    if isnothing(velocity_range)
        vel_min = 0.9 * minimum(vals)
        vel_max = 1.1 * maximum(vals)
    else
        vel_min, vel_max = velocity_range
    end

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
        title=attr(text=title, font=attr(size=font_size + 2, family=font_family)),
        xaxis=attr(title="Period (s)", type="linear", showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        yaxis=attr(title="Group Velocity (km/s)", range=[vel_min, vel_max], showgrid=true, gridcolor="rgba(128,128,128,0.2)"),
        width=width, height=height,
        plot_bgcolor="white", paper_bgcolor="white",
        margin=attr(l=80, r=100, t=80, b=80),
        showlegend=true
    )
    return _plotly_plot(traces, layout)
end
