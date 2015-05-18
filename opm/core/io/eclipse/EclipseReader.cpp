/*
  Copyright 2015 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM. If not, see <http://www.gnu.org/licenses/>.
*/

#include "EclipseReader.hpp"
#include <iostream>
#include <opm/core/simulator/WellState.hpp>

#include <ert/ecl/ecl_file.h>

namespace Opm
{
    void restore_kw( const std::string restart_filename, const char * kw, int reportStep, WellState& wellState ) {
        const char * filename = restart_filename.c_str();

        ecl_file_type * file_type = ecl_file_open(filename, 0);
        bool block_selected = ecl_file_select_rstblock_report_step(file_type , reportStep);

        if (block_selected){
            ecl_kw_type * xwel = ecl_file_iget_named_kw(file_type , kw, 0);
            const double* xwel_data = ecl_kw_get_double_ptr(xwel);
            size_t offset = 0;

            for(size_t i = 0; i < wellState.bhp().size(); i++){
                wellState.bhp()[i] = xwel_data[offset];
                std::cout << "WellState.bhp : " << wellState.bhp()[i] << std::endl;
                offset++;
            }

            for(size_t i = 0; i < wellState.perfPress().size(); i++){
                wellState.perfPress()[i] = xwel_data[offset];
                std::cout << "WellState.perfPress : " << wellState.perfPress()[i] << std::endl;
                offset++;
            }

            for(size_t i = 0; i < wellState.perfRates().size(); i++){
                wellState.perfRates()[i] = xwel_data[offset];
                std::cout << "WellState.perfRates : " << wellState.perfRates()[i] << std::endl;
                offset++;
            }

            for(size_t i = 0; i < wellState.temperature().size(); i++){
                wellState.temperature()[i] = xwel_data[offset];
                std::cout << "WellState.temperature : " << wellState.temperature()[i] << std::endl;
                offset++;
            }

            for(size_t i = 0; i < wellState.wellRates().size(); i++){
                wellState.wellRates()[i] = xwel_data[offset];
                std::cout << "WellState.wellRates : " << wellState.wellRates()[i] << std::endl;
                offset++;
            }
        }
    }
} // namespace Opm
