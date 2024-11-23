/** BELLA experiment track reconstruction framework
 *
 * (c) 2024 Lawrence Berkeley National Laboratory
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/options/details/interface.hpp"

namespace traccc::opts
{

    /// Convenience namespace shorthand
    namespace po = boost::program_options;

    /// Command line options used in the propagation tests
    class field_options : public interface
    {

    public:
        /// Constructor
        field_options() : interface("B Field Options")
        {

            m_desc.add_options()("bfield-file",
                                 po::value(&(bfield_file))
                                     ->default_value(""),
                                 "B field file name");
        }

        std::string bfield_file;

    }; // class track_propagation

} // namespace traccc::opts
