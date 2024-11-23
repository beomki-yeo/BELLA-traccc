/** BELLA experiment track reconstruction framework
 *
 * (c) 2024 Lawrence Berkeley National Laboratory
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/fitting/fitting_algorithm.hpp"
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"
#include "traccc/io/read_geometry.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/performance.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/track_propagation.hpp"
#include "traccc/utils/seed_generator.hpp"
#include "traccc/utils/event_data.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"
#include "detray/core/detector_metadata.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"

// Local include(s).
#include "src/field_options.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace traccc;
namespace po = boost::program_options;

// The main routine
//
int main(int argc, char *argv[])
{
    // Program options.
    traccc::opts::detector detector_opts;
    traccc::opts::input_data input_opts;
    traccc::opts::track_propagation propagation_opts;
    traccc::opts::field_options field_opts;
    traccc::opts::program_options program_opts{
        "Truth Track Fitting on the Host",
        {detector_opts, input_opts, propagation_opts, field_opts},
        argc,
        argv};

    /// Type declarations
    using host_detector_type = detray::detector<detray::default_metadata,
                                                detray::host_container_types>;

    using b_field_t = covfie::field<detray::bfield::inhom_bknd_t>;
    using rk_stepper_type =
        detray::rk_stepper<b_field_t::view_t, traccc::default_algebra,
                           detray::constrained_step<>>;

    using host_navigator_type = detray::navigator<const host_detector_type>;
    using host_fitter_type =
        traccc::kalman_fitter<rk_stepper_type, host_navigator_type>;

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;

    /*****************************
     * Build a geometry
     *****************************/

    // B field value and its type
    b_field_t field = detray::io::read_bfield<b_field_t>(field_opts.bfield_file);

    // Read the detector
    detray::io::detector_reader_config reader_cfg{};
    reader_cfg.add_file(detector_opts.detector_file);
    if (!detector_opts.material_file.empty())
    {
        reader_cfg.add_file(detector_opts.material_file);
    }
    if (!detector_opts.grid_file.empty())
    {
        reader_cfg.add_file(detector_opts.grid_file);
    }
    const auto [host_det, names] =
        detray::io::read_detector<host_detector_type>(host_mr, reader_cfg);

    /*****************************
     * Do the reconstruction
     *****************************/

    // Fitting algorithm object
    typename traccc::fitting_algorithm<host_fitter_type>::config_type fit_cfg;
    fit_cfg.propagation = propagation_opts;

    traccc::fitting_algorithm<host_fitter_type> host_fitting(fit_cfg);

    // Residual file
    std::ofstream residual_file;
    residual_file.open("residual.csv");
    residual_file << "fit_qop, fit_qopT, fit_qopz, ";
    residual_file << "truth_qop, truth_qopT, truth_qopz,";
    residual_file << "qop_residual, qopT_residual, qopz_residual";
    residual_file << std::endl;

    // Track state file
    std::ofstream state_file;
    state_file.open("state.csv");
    state_file << "event_id, fit_track_id, x, y, z";
    state_file << std::endl;

    // Iterate over events
    for (auto event = input_opts.skip;
         event < input_opts.events + input_opts.skip; ++event)
    {

        // Truth Track Candidates
        traccc::event_data evt_data(input_opts.directory, event, host_mr,
                                    input_opts.use_acts_geom_source, &host_det,
                                    input_opts.format, false);

        assert(evt_data.m_particle_map.size() > 0u);
        /// Assume that all particle momentum are the same
        const auto truth_mom = evt_data.m_particle_map.begin()->second.momentum;
        const auto charge = evt_data.m_particle_map.begin()->second.charge;

        const scalar qop_stddev = 0.05f * math::abs(charge) / getter::norm(truth_mom);

        /// Standard deviations for seed track parameters
        const std::array<scalar, e_bound_size> stddevs = {
            0.02f * detray::unit<scalar>::mm,
            0.02f * detray::unit<scalar>::mm,
            0.0085f,
            0.0085f,
            qop_stddev,
            1.f * detray::unit<scalar>::ns};

        // Seed generator
        traccc::seed_generator<host_detector_type> sg(host_det, stddevs);

        traccc::track_candidate_container_types::host truth_track_candidates =
            evt_data.generate_truth_candidates(sg, host_mr);

        // Run fitting
        auto track_states =
            host_fitting(host_det, field, truth_track_candidates);

        std::cout << "Number of fitted tracks: " << track_states.size()
                  << std::endl;

        const decltype(track_states)::size_type n_fitted_tracks =
            track_states.size();

        for (unsigned int i = 0; i < n_fitted_tracks; i++)
        {
            const auto &trk_states_per_track = track_states.at(i).items;

            if (trk_states_per_track.size() == 0u)
            {
                throw std::runtime_error("track states is empty");
            }
            /*
            if (trk_states_per_track.size() < 6u)
            {
                throw std::runtime_error(
                    "The number of track states per track (" +
                    std::to_string(trk_states_per_track.size()) +
                    ") is less than 6");
            }
            */

            const auto &fit_res = track_states[i].header;

            /************************************
             *  Write Residuals of qop
             * **********************************/

            // Fit qop
            const auto &fit_par = trk_states_per_track.at(0).smoothed();
            const scalar fit_qop = fit_par.qop();
            const scalar fit_qopT = fit_par.qopT();
            // @NOTE: qopz is a signed value
            const scalar fit_qopz = fit_par.qopz();

            // Truth qop
            const measurement meas = trk_states_per_track.at(0).get_measurement();
            const auto global_mom = evt_data.m_meas_to_param_map.at(meas).second;
            const auto p = getter::norm(global_mom);
            const auto pT = getter::perp(global_mom);
            // @NOTE: pz is a signed value
            const auto pz = global_mom[2];

            std::map<particle, uint64_t> contributing_particles =
                evt_data.m_meas_to_ptc_map.at(meas);
            const particle ptc = contributing_particles.begin()->first;
            const auto q = ptc.charge;

            const scalar truth_qop = q / p;
            const scalar truth_qopT = q / pT;
            const scalar truth_qopz = q / pz;

            residual_file << fit_qop << "," << fit_qopT << "," << fit_qopz << ",";
            residual_file << truth_qop << "," << truth_qopT << "," << truth_qopz << ",";
            residual_file << fit_qop - truth_qop << ",";
            residual_file << fit_qopT - truth_qopT << ",";
            residual_file << fit_qopz - truth_qopz << " \n";

            for (const auto &st : trk_states_per_track)
            {
                const detray::tracking_surface sf{host_det, st.surface_link()};
                const auto xyz = sf.bound_to_global({}, st.smoothed().bound_local(), st.smoothed().dir());
                state_file << event << "," << i << "," 
                           << xyz[0] << "," << xyz[1] << "," << xyz[2] << "\n";
            }
        }
    }

    residual_file.close();
    state_file.close();

    return EXIT_SUCCESS;
}
