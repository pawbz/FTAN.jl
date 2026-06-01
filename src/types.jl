begin
	# Data structures for MFT analysis
	"""
	    MFTResult
	
	Structure to hold Multiple Filter Technique analysis results.
	"""
	struct MFTResult
	    periods::Vector{Float64}                    # Analysis periods [s]
	    frequencies::Vector{Float64}                # Analysis frequencies [Hz]
	    group_velocities::Vector{Float64}           # Group velocities [km/s]
	    phase_velocities::Vector{Float64}           # Phase velocities [km/s]
	    phase_branch_numbers::Vector{Int}           # Candidate integer cycle counts N
	    phase_velocity_branches::Matrix{Float64}    # Phase velocities [period, branch]
	    measured_phases::Vector{Float64}            # Wrapped analytic phase at group arrival [rad]
	    selected_phase_branches::Vector{Int}        # Selected integer cycle count N per period
	    phase_suspect::Vector{Bool}                 # Smoothness flags for manual review
	    u_predicted_from_phase::Vector{Float64}     # Group velocity implied by c(T)
	    amplitudes::Vector{Float64}                 # Spectral amplitudes
	    filtered_traces::AbstractMatrix{<:Real}     # Filtered time series [time, frequency]; empty for picks-only runs
	    envelopes::AbstractMatrix{<:Real}           # Amplitude envelopes [time, frequency]; empty for picks-only runs
	    arrival_times::Vector{Float64}              # Group arrival times [s]
	    distance::Float64                           # Source-receiver distance [km]
	    quality_factors::Vector{Float64}  # Quality assessment for each measurement
		all_peaks::Vector{Vector{Tuple{Float64,Float64}}}  # Up to 4 strongest envelope peaks per period: (time, amplitude)
		time::Vector{Float64 }
	    storage_mode::Symbol
	end
	
	mft_has_full_storage(res::MFTResult) =
	    res.storage_mode == :full && !isempty(res.filtered_traces) && !isempty(res.envelopes)
end

struct MultimodalDispersion
    periods::Vector{Float64}
    arrival_times::Vector{Float64}
    group_velocities::Vector{Float64}
    peak_amplitudes::Vector{Float64}
    mode_index::Int

    function MultimodalDispersion(periods, arrival_times, group_velocities, peak_amplitudes, mode_index)
        n = length(periods)
        @assert length(arrival_times) == n "arrival_times length must match periods"
        @assert length(group_velocities) == n "group_velocities length must match periods"
        @assert length(peak_amplitudes) == n "peak_amplitudes length must match periods"
        new(periods, arrival_times, group_velocities, peak_amplitudes, mode_index)
    end
end

struct BranchAnalysisResult
    causal_result::MFTResult
    acausal_result::MFTResult
    causal_modes::Vector{MultimodalDispersion}
    acausal_modes::Vector{MultimodalDispersion}
    branch_correlation::Vector{Float64}
    periods::Vector{Float64}
    distance::Float64
end

struct BranchBatchAnalysisResult
    state_results::Vector{BranchAnalysisResult}
    branch_correlation::Matrix{Float64}
    periods::Vector{Float64}
    state_labels::Vector{String}
end

struct PhaseRegressionObservation
    pair_label::String
    period::Float64
    frequency::Float64
    omega::Float64
    distance::Float64
    t_peak::Float64
    phi_meas::Float64
    base_phase::Float64
    sigma_phi::Float64
    quality::Float64
    azimuth_deg::Float64
end

# ╔═╡ a2000002-0000-0000-0000-000000000001
"""
    PhaseVelocityPeriodFit

Robust fit result for one period. `unwrapped_phase[i]` equals
`base_phase + 2pi * cycle_number[i]` for the corresponding observation.
"""
struct PhaseVelocityPeriodFit
    period::Float64
    frequency::Float64
    omega::Float64
    phase_velocity::Float64
    sigma_c::Float64
    slope::Float64
    sigma_slope::Float64
    observations::Vector{PhaseRegressionObservation}
    cycle_numbers::Vector{Int}
    unwrapped_phase::Vector{Float64}
    residuals::Vector{Float64}
    inlier::Vector{Bool}
    azimuth_outlier::Vector{Bool}
    weights::Vector{Float64}
    prior_c::Float64
    ransac_score::Float64
end

# ╔═╡ a2000003-0000-0000-0000-000000000001
struct PhaseVelocityCurveFit
    periods::Vector{Float64}
    phase_velocities::Vector{Float64}
    sigma_c::Vector{Float64}
    fits::Vector{PhaseVelocityPeriodFit}
end
