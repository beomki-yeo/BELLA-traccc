/** BELLA experiment track reconstruction framework
 *
 * (c) 2024 Lawrence Berkeley National Laboratory
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/options/generation.hpp"
#include "traccc/options/output_data.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/telescope_detector.hpp"
#include "traccc/options/track_propagation.hpp"
#include "traccc/simulation/measurement_smearer.hpp"
#include "traccc/simulation/simulator.hpp"
#include "traccc/simulation/smearing_writer.hpp"

// detray include(s).
#include "detray/detectors/bfield.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/io/frontend/detector_writer.hpp"
#include "detray/materials/material.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/test/utils/detectors/build_telescope_detector.hpp"
#include "detray/test/utils/simulation/event_generator/track_generators.hpp"

// Local include(s).
#include "src/field_options.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// Boost include(s).
#include <boost/filesystem.hpp>

using namespace traccc;

// The main routine
//
int main(int argc, char *argv[])
{
    // Program options.
    traccc::opts::generation generation_opts;
    traccc::opts::output_data output_opts;
    traccc::opts::track_propagation propagation_opts;
    traccc::opts::field_options field_opts;
    traccc::opts::program_options program_opts{
        "Telescope-Detector Simulation",
        {generation_opts, output_opts, propagation_opts, field_opts},
        argc,
        argv};

    // Use deterministic random number generator for testing
    using uniform_gen_t =
        detray::detail::random_numbers<scalar,
                                       std::uniform_real_distribution<scalar>>;

    // Memory resource
    vecmem::host_memory_resource host_mr;

    /*****************************
     * Build the Bella Detector
     *****************************/

    // Plane alignment direction (Planes are aligned along the x-axis)
    const vector3 align_axis{1.f, 0.f, 0.f};
    detray::detail::ray<traccc::default_algebra> pilot_track{
        {0, 0, 0}, 0, align_axis, -1};

    // Positions of 12 sensitive planes (in mm unit)
    std::vector<scalar> sensitive_positions{
        10.f, 20.f, 30.f,    // 3 SCC
        60.f, 70.f, 80.f,    // 3 SCC
        180.f, 190.f, 200.f, // 3 SCC
        230.f, 240.f, 250.f  // 3 SCC
    };

    // Set sensitive planes material, thickness and its size
    detray::material<scalar> sensitive_mat = detray::silicon<scalar>();
    const scalar sensitive_thickness = 5.f * traccc::unit<scalar>::mm;
    detray::mask<detray::rectangle2D> sensitive_rect{0u,
                                                     100.f * traccc::unit<scalar>::mm,
                                                     100.f * traccc::unit<scalar>::mm};

    // Create the telescope geometry with 12 sensitive planes
    detray::tel_det_config<detray::rectangle2D,
                           detray::detail::ray<traccc::default_algebra>>
        tel_cfg{sensitive_rect, pilot_track};
    tel_cfg.positions(sensitive_positions);
    tel_cfg.module_material(sensitive_mat);
    tel_cfg.mat_thickness(sensitive_thickness);
    tel_cfg.envelope(100.f * traccc::unit<scalar>::mm);

    const auto [det, name_map] = build_telescope_detector(host_mr, tel_cfg);

    // (WIP) Add the Magnet

    // (WIP) Add any structure needed

    // (WIP) Add the attenuator to the telescope geometry

    /*
    // Set attenuator material and thickness
    detray::material<scalar> attenuator_mat = detray::iron<scalar>();
    const scalar attenuator_thickness = 30.f * traccc::unit<scalar>::mm;
    detray::mask<detray::rectangle2D> attenuator_rect{0u,
                                                      100.f * traccc::unit<scalar>::mm,
                                                      100.f * traccc::unit<scalar>::mm};
    auto &surfaces = det.surfaces();
    auto &masks = det.mask_store();
    auto &materials = det.material_store();
    auto &transforms = det.transform_store();
    surfaces.push_back({trf_index, mask_link, material_link, volume_idx,
                        surface_id::e_sensitive},
                       invalid_src_link);
    masks.template emplace_back<mask_id>(empty_context{}, m_boundaries,
                                         mask_volume_link);
    materials.template emplace_back<material_id>(empty_context{}, m_boundaries,
                                                 mask_volume_link);
    transforms.emplace_back(ctx, mod_placement.pos, m_local_z,
                            m_local_x);
    */

    /***************************
     * Run the muon simulation
     ***************************/

    // Origin of particles
    using generator_type =
        detray::random_track_generator<traccc::free_track_parameters,
                                       uniform_gen_t>;
    generator_type::configuration gen_cfg{};
    gen_cfg.n_tracks(generation_opts.gen_nparticles);
    gen_cfg.origin(traccc::point3{generation_opts.vertex[0],
                                  generation_opts.vertex[1],
                                  generation_opts.vertex[2]});
    gen_cfg.origin_stddev(traccc::point3{generation_opts.vertex_stddev[0],
                                         generation_opts.vertex_stddev[1],
                                         generation_opts.vertex_stddev[2]});
    gen_cfg.phi_range(generation_opts.phi_range);
    gen_cfg.theta_range(generation_opts.theta_range);
    gen_cfg.mom_range(generation_opts.mom_range);
    gen_cfg.charge(generation_opts.ptc_type.charge());
    generator_type generator(gen_cfg);

    // Smearing value for measurements (measurement noise)
    traccc::measurement_smearer<traccc::default_algebra> meas_smearer(
        50.f * traccc::unit<scalar>::mm, 50.f * traccc::unit<scalar>::mm);

    // Type declarations
    using detector_type = decltype(det);
    using writer_type = traccc::smearing_writer<
        traccc::measurement_smearer<traccc::default_algebra>>;

    // Writer config
    typename writer_type::config smearer_writer_cfg{meas_smearer};

    // Run simulator
    const std::string full_path = output_opts.directory;
    boost::filesystem::create_directories(full_path);

    // B field value and its type
    using b_field_t = covfie::field<detray::bfield::inhom_bknd_t>;
    const auto field = detray::io::read_bfield<b_field_t>(field_opts.bfield_file);

    auto sim = traccc::simulator<detector_type, b_field_t, generator_type,
                                 writer_type>(
        generation_opts.ptc_type, generation_opts.events, det, field,
        std::move(generator), std::move(smearer_writer_cfg), full_path);
    sim.get_config().propagation = propagation_opts;

    sim.run();

    // Create detector file
    auto writer_cfg = detray::io::detector_writer_config{}
                          .format(detray::io::format::json)
                          .replace_files(true);
    detray::io::write_detector(det, name_map, writer_cfg);

    return 1;
}