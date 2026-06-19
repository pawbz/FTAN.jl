module FTAN

using ColorSchemes
using Colors
using DSP
using FFTW
using LinearAlgebra
using Peaks
using PlutoPlotly
using Printf
using Statistics

export MFTResult, MFTFilterBank, MultimodalDispersion
export DEFAULT_VELOCITY_RANGE
export mft_has_full_storage
export narrow_band_filter, narrow_band_filter_analytic, multiple_narrow_band_filters
export find_group_arrival, find_group_arrivals, compute_envelope
export compute_phase_velocity, compute_group_velocity_from_phase
export perform_mft_analysis, perform_mft_analysis!, perform_mft_analysis_batch, perform_mft_analysis_batch!
export wavelength_valid_period, haversine_distance_km, linear_array_order
export extract_all_peaks_matrix, match_peaks_continuous, extract_all_modes
export pdsurftomo_dispersion_rows, write_pdsurftomo_dispersion
export plot_envelopes, plot_dispersion_curve, plot_crosscorr_bundle, plot_phase_velocities
export plot_multipeak_dispersion
export plot_uc_consistency_comparison, plot_group_phase_consistency
export plot_source_state_waveforms, plot_cluster_histogram, plot_state_ncc_heatmap

include("types.jl")
include("phase_velocity.jl")
include("filtering.jl")
include("picking.jl")
include("analysis.jl")
include("utils.jl")
include("plotting.jl")
include("source_state_plots.jl")

end
