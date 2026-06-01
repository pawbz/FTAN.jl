function narrow_band_filter(data::AbstractArray, dt::Float64, center_freq::Float64, 
                           bandwidth_factor::Float64=0.1)
    npts = size(data, 1)
    
    # FFT parameters
    fft_data = fft(data, 1)
    freqs = cat(collect(fftfreq(npts, 1.0/dt)), dims=ndims(data))
    
    # Gaussian filter parameters
    sigma_freq = center_freq * bandwidth_factor / 2.0  # Standard deviation
    
    # Create Gaussian filter in frequency domain
    filter_response = exp.(-0.5 * ((freqs .- center_freq) ./ sigma_freq).^2)
    filter_response += exp.(-0.5 * ((freqs .+ center_freq) ./ sigma_freq).^2)  # Negative frequencies
    
    # Apply filter
    filtered_fft = fft_data .* filter_response
    
    # Transform back to time domain
    filtered_data = real.(ifft(filtered_fft, 1))
    
    return filtered_data
end

function narrow_band_filter_analytic(data::AbstractVector{<:Real}, dt::Float64,
                                     center_freq::Float64,
                                     bandwidth_factor::Float64=sqrt(2.0 / 25.0))
    npts  = length(data)
    spec  = fft(data)
    freqs = fftfreq(npts, 1.0 / dt)
    out   = zeros(ComplexF64, npts)
    sigma_freq = center_freq * bandwidth_factor / 2.0

    for k in eachindex(freqs)
        f = freqs[k]
        if f > 0.0
            H = exp(-0.5 * ((f - center_freq) / sigma_freq)^2)
            out[k] = 2.0 * H * spec[k]   # ×2: one-sided → analytic convention
        elseif f == 0.0
            # DC: near-zero contribution at any usable f₀
            H = exp(-0.5 * (center_freq / sigma_freq)^2)
            out[k] = H * spec[k]
        end
        # f < 0: out[k] stays zero (one-sided spectrum)
    end
    return ifft(out)
end

# ╔═╡ e9bfb21b-f29d-410f-827f-e66d1117020f
"""
    compute_phase_velocity(analytic_signal, time, t_group, freq, distance;
                           phvel_source_phase=0.0, branch=0) -> Float64

Estimate phase velocity from the instantaneous phase of the analytic signal at
the group arrival time for one specified integer cycle branch.

> **Note:** `perform_mft_analysis` now evaluates explicit integer `N`
> candidates and resolves the branch ambiguity with a smoothness pass. This
> utility is for single-point checks when `branch` is already known.

The total phase accumulated travelling distance Δ at phase velocity c is:

    φ_total = 2π f₀ Δ / c

The instantaneous phase of the analytic signal at the group arrival:

    φ(tᵍ) = angle( z(tᵍ) ) = 2π f₀ tᵍ − φ_total − φ_source + 2πN

Rearranging (n = 0 branch):

    c = 2π f₀ Δ / ( 2π f₀ tᵍ − φ(tᵍ) − φ_source + 2πN )

# Arguments
- `analytic_signal`  : complex analytic signal from `narrow_band_filter_analytic`
- `time`             : time vector [s] (same grid)
- `t_group`          : group arrival time [s]
- `freq`             : centre frequency f₀ [Hz]
- `distance`         : source–receiver distance [km]
- `phvel_source_phase` : source/CCF phase shift [rad], subtracted from the denominator
- `branch`             : integer cycle count N

# Returns
- Phase velocity [km/s], or `NaN` if undefined or implausible
"""
function compute_phase_velocity(analytic_signal::AbstractVector{<:Complex},
                                time::Vector{Float64},
                                t_group::Float64,
                                freq::Float64,
                                distance::Float64;
                                phvel_correction::Float64=0.0,
                                phvel_source_phase::Float64=phvel_correction,
                                branch::Int=0)
    (isnan(t_group) || t_group <= 0.0) && return NaN
    dt_local = time[2] - time[1]
    i_g = clamp(round(Int, (t_group - time[1]) / dt_local) + 1, 1, length(time))

    # Instantaneous phase at group arrival (wrapped to (−π, π])
    phi_measured = angle(analytic_signal[i_g])

    omega = 2π * freq
    denom = omega * t_group - phi_measured - phvel_source_phase + 2π * branch
    abs(denom) < 1e-6 && return NaN

    c = omega * distance / denom
    # Sanity check: plausible seismic phase velocity (0.5 – 20 km/s)
    (c < 0.5 || c > 20.0) && return NaN
    return c
end

# ╔═╡ f71ee0c9-1e34-49c2-9385-4589c561f4e4
"""
    multiple_narrow_band_filters(data, dt, periods::Vector{Float64}, bandwidth_percent::Float64)

Apply multiple narrowband filters at different central periods and sum them.
Useful for enhancing specific frequency bands before CSS.
"""
function multiple_narrow_band_filters(traces, dt, periods::Vector{Float64}, bandwidth_percent)
   
    filtered_components = []
    
    # Filter at each central period and accumulate
     filtered_components = map(periods) do T
        center_freq = 1.0 / T
        filtered_data = narrow_band_filter(traces, dt, center_freq, bandwidth_percent / 100.0)
    end
    
    return filtered_components
end

function compute_envelope(x)
    return abs.(hilbert(x));
end

begin
	begin
		"""
		    MFTFilterBank
		
		Pre-computed Gaussian filter bank for efficient repeated MFT analysis.
		Stores FFTW plans, workspace buffers, and filter coefficients so that
		`perform_mft_analysis_batch!` requires only one batch FFT + nfreq IFFTs
		for N waveforms, with no heap allocation after the first call.
		
		Construct once per notebook session (parameters fixed by dt, npts, periods).
		Set `N_initial` to the maximum waveform-column count the bang calls will use.
		"""
		mutable struct MFTFilterBank{T<:AbstractFloat}
		    # geometry
		    periods::Vector{Float64}
		    frequencies::Vector{Float64}
		    dt::Float64                 # upsampled = dt_original / upsample_factor
		    dt_original::Float64
		    npts_original::Int          # length after upsampling (before pad)
		    npts_padded::Int            # npts_original * zero_pad_factor
		    nfreq::Int
		    velocity_range::Tuple{Float64,Float64}
		    bandwidth_factor::Float64
		    zero_pad_factor::Int
		    upsample_factor::Float64
		    time::Vector{Float64}       # (npts_original,) time axis starting at dt
		
		    # pre-computed filter weights: one-sided Gaussians, ×2 on positive freqs, neg freqs = 0
		    H_full::Matrix{Complex{T}}  # (npts_padded × nfreq)
		
		    # workspace preallocated for N_buf waveform columns
		    N_buf::Int
		    W_pad_buf::Matrix{Complex{T}}    # (npts_padded × N_buf) input after upsample+pad
		    SPEC_buf::Matrix{Complex{T}}     # (npts_padded × N_buf) FFT of W_pad_buf
		    Z_buf::Matrix{Complex{T}}        # (npts_padded × N_buf) per-freq filtered spectrum
		    Z_time_buf::Matrix{Complex{T}}   # (npts_padded × N_buf) IFFT result (analytic signal)
		    envelopes_buf::Union{Nothing,Array{T,3}}  # Full-storage Hilbert envelopes
		    filtered_buf::Union{Nothing,Array{T,3}}   # Full-storage real filtered traces
		    arrivals_buf::Matrix{Float64}    # (nfreq × N_buf) group arrival times
		    phases_buf::Matrix{Float64}      # (nfreq × N_buf) wrapped phases at group arrival
	        storage_mode::Symbol
		
		    # FFTW out-of-place plans: mul!(dst, plan, src) writes into dst without allocating
		    plan_fwd::Any   # plan_fft(W_pad_buf, 1)  — mul!(SPEC_buf, plan_fwd, W_pad_buf)
		    plan_inv::Any   # plan_ifft(Z_buf, 1)     — mul!(Z_time_buf, plan_inv, Z_buf)
		end
		
		function MFTFilterBank(dt_original::Real, npts_raw::Int,
		                       periods::Vector{Float64};
		                       bandwidth_factor::Float64=sqrt(2.0/25.0),
		                       zero_pad_factor::Int=4,
		                       upsample_factor::Real=2.0,
		                       velocity_range::Tuple{Float64,Float64}=(2.0, 6.0),
	                           precision::Type{T}=Float32,
	                           storage_mode::Symbol=:picks_only,
		                       N_initial::Int=1) where {T<:AbstractFloat}
	        storage_mode in (:picks_only, :full) ||
	            throw(ArgumentError("storage_mode must be :picks_only or :full"))
			
		    upsample = Float64(upsample_factor)
		    upsample > 0.0 || throw(ArgumentError("upsample_factor must be positive"))
	
		    dt_original = Float64(dt_original)
		    dt = dt_original / upsample
		    # Derive exact upsampled length using the same matrix code path as perform_mft_analysis_batch!
		    npts_original = upsample == 1.0 ? npts_raw :
		        size(DSP.resample(zeros(T, npts_raw, 1), upsample; dims=1), 1)
		    npts_padded   = npts_original * zero_pad_factor
		    nfreq         = length(periods)
		    frequencies   = 1.0 ./ periods
		
		    # Build H_full: one-sided complex Gaussian filter matrix
		    freqs_full = fftfreq(npts_padded, 1.0 / dt)
		    H_full = zeros(Complex{T}, npts_padded, nfreq)
		    for ifreq in 1:nfreq
		        f0    = frequencies[ifreq]
		        sigma = f0 * bandwidth_factor / 2.0
		        for k in eachindex(freqs_full)
		            f = freqs_full[k]
		            if f > 0.0
		                H_full[k, ifreq] = T(2.0 * exp(-0.5 * ((f - f0) / sigma)^2))
		            elseif f == 0.0
		                H_full[k, ifreq] = T(exp(-0.5 * (f0 / sigma)^2))
		            end
		            # f < 0: stays zero — one-sided analytic spectrum
		        end
		    end
		
		    time = collect(range(dt, step=dt, length=npts_original))
		
		    W_pad_buf  = zeros(Complex{T}, npts_padded, N_initial)
		    SPEC_buf   = zeros(Complex{T}, npts_padded, N_initial)
		    Z_buf      = zeros(Complex{T}, npts_padded, N_initial)
		    Z_time_buf = zeros(Complex{T}, npts_padded, N_initial)
		    envelopes_buf = storage_mode == :full ? zeros(T, npts_original, nfreq, N_initial) : nothing
		    filtered_buf  = storage_mode == :full ? zeros(T, npts_original, nfreq, N_initial) : nothing
		    arrivals_buf  = fill(NaN, nfreq, N_initial)
		    phases_buf    = fill(NaN, nfreq, N_initial)
		
		    # Out-of-place plans measured against the initial buffer shapes
		    plan_fwd = FFTW.plan_fft(W_pad_buf, 1)
		    plan_inv = FFTW.plan_ifft(Z_buf, 1)
		
			    MFTFilterBank{T}(
			        periods, frequencies, dt, dt_original, npts_original, npts_padded, nfreq,
			        velocity_range, bandwidth_factor, zero_pad_factor, upsample, time, H_full,
		        N_initial, W_pad_buf, SPEC_buf, Z_buf, Z_time_buf,
		        envelopes_buf, filtered_buf, arrivals_buf, phases_buf, storage_mode,
		        plan_fwd, plan_inv,
		    )
		end
	end
	
	_bank_float_type(::MFTFilterBank{T}) where {T} = T
end
