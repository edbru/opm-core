#ifndef ECLIPSEREADER_HPP
#define ECLIPSEREADER_HPP

#include <iostream>
#include <opm/core/simulator/WellState.hpp>

namespace Opm
{
    void restore_kw( const std::string restart_filename, const char * kw, int report_step, WellState& wellState );
}

#endif // ECLIPSEREADER_HPP
