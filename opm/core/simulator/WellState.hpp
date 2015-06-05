/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_WELLSTATE_HEADER_INCLUDED
#define OPM_WELLSTATE_HEADER_INCLUDED

#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <vector>
#include <cassert>
#include <utility>
#include <iostream>
#include <string>

namespace Opm
{
    struct WellIndex{
        size_t wellNumber;
        size_t wellCompletionOffset;
        std::map<std::size_t, size_t> completionMap;
    };

    /// The state of a set of wells.
    class WellState
    {
    public:
        /// Allocate and initialize if wells is non-null.
        /// Also tries to give useful initial values to the bhp() and
        /// wellRates() fields, depending on controls.  The
        /// perfRates() field is filled with zero, and perfPress()
        /// with -1e100.
        template <class State>
        void init(const Wells* wells, const State& state)
        {
            if (wells) {
                const int nw = wells->number_of_wells;
                const int np = wells->number_of_phases;
                bhp_.resize(nw);
                temperature_.resize(nw, 273.15 + 20); // standard temperature for now
                wellrates_.resize(nw * np, 0.0);
                for (int w = 0; w < nw; ++w) {
                    assert((wells->type[w] == INJECTOR) || (wells->type[w] == PRODUCER));
                    const WellControls* ctrl = wells->ctrls[w];
                    if (well_controls_well_is_stopped(ctrl)) {
                        // Stopped well:
                        // 1. Assign zero well rates.
                        for (int p = 0; p < np; ++p) {
                            wellrates_[np*w + p] = 0.0;
                        }
                        // 2. Assign bhp equal to bhp control, if
                        //    applicable, otherwise assign equal to
                        //    first perforation cell pressure.
                        if (well_controls_get_current_type(ctrl) == BHP) {
                            bhp_[w] = well_controls_get_current_target( ctrl );
                        } else {
                            const int first_cell = wells->well_cells[wells->well_connpos[w]];
                            bhp_[w] = state.pressure()[first_cell];
                        }
                    } else {
                        // Open well:
                        // 1. Initialize well rates to match controls
                        //    if type is SURFACE_RATE.  Otherwise, we
                        //    cannot set the correct value here, so we
                        //    assign a small rate with the correct
                        //    sign so that any logic depending on that
                        //    sign will work as expected.
                        if (well_controls_get_current_type(ctrl) == SURFACE_RATE) {
                            const double rate_target = well_controls_get_current_target(ctrl);
                            const double * distr = well_controls_get_current_distr( ctrl );
                            for (int p = 0; p < np; ++p) {
                                wellrates_[np*w + p] = rate_target * distr[p];
                            }
                        } else {
                            const double small_rate = 1e-14;
                            const double sign = (wells->type[w] == INJECTOR) ? 1.0 : -1.0;
                            for (int p = 0; p < np; ++p) {
                                wellrates_[np*w + p] = small_rate * sign;
                            }
                        }
                        // 2. Initialize bhp to be target pressure if
                        //    bhp-controlled well, otherwise set to a
                        //    little above or below (depending on if
                        //    the well is an injector or producer)
                        //    pressure in first perforation cell.
                        if (well_controls_get_current_type(ctrl) == BHP) {
                            bhp_[w] = well_controls_get_current_target( ctrl );
                        } else {
                            const int first_cell = wells->well_cells[wells->well_connpos[w]];
                            const double safety_factor = (wells->type[w] == INJECTOR) ? 1.01 : 0.99;
                            bhp_[w] = safety_factor*state.pressure()[first_cell];
                        }
                    }
                }
                // The perforation rates and perforation pressures are
                // not expected to be consistent with bhp_ and wellrates_
                // after init().
                perfrates_.resize(wells->well_connpos[nw], 0.0);
                perfpress_.resize(wells->well_connpos[nw], -1e100);
            }
        }

        void resize(const EclipseState& state, size_t reportStep)
        {
            std::vector<WellPtr> wells = state.getSchedule()->getOpenWells(reportStep);
            int numPhases = state.getNumPhases();
            int numCompl = getNumCompletions(wells, reportStep);
            EclipseGridConstPtr grid = state.getEclipseGrid();

            size_t numWells = wells.size();
            std::vector<double> bhpTemp(bhp_.begin(), bhp_.end());
            std::vector<double> temperatureTemp(temperature_.begin(), temperature_.end());
            std::vector<double> wellRatesTemp(wellrates_.begin(), wellrates_.end());
            std::vector<double> perfratesTemp(perfrates_.begin(), perfrates_.end());
            std::vector<double> perfpressTemp(perfpress_.begin(), perfpress_.end());

            std::map<std::string, WellIndex> newMap;

            size_t totalWellCompletionOffset = 0;
            for (size_t i = 0; i < numWells; i++){
                WellIndex wi;
                wi.wellNumber = i;
                wi.wellCompletionOffset = totalWellCompletionOffset;

                totalWellCompletionOffset += getNumCompletions(wells[i], reportStep);
                fillCompletionsMap(grid, wells[i], wi, reportStep);
                std::pair <std::string, WellIndex> p(wells[i]->name(), wi);
                newMap.insert(p);
            }

            bhp_.resize(numWells);
            temperature_.resize(numWells, 273.15 + 20);
            wellrates_.resize(numWells * numPhases, 0.0);
            perfrates_.resize(numCompl, 0.0);
            perfpress_.resize(numCompl, -1e100);

            for (auto oldIter = wellIndexMap_.begin(); oldIter != wellIndexMap_.end(); ++oldIter){
                const std::string& wellName = (*oldIter).first;
                auto newIter = newMap.find(wellName);

                if(newIter != newMap.end()){
                    const WellIndex& oldWi = (*oldIter).second;
                    size_t oldIndex = oldWi.wellNumber;

                    const WellIndex& newWi = (*newIter).second;
                    size_t newIndex = newWi.wellNumber;

                    bhp_[newIndex] = bhpTemp[oldIndex];
                    temperature_[newIndex] = temperatureTemp[oldIndex];

                    for (int i = 0; i < numPhases; i++){
                        wellrates_[newIndex * numPhases + i] = wellRatesTemp[oldIndex * numPhases + i];
                    }

                    for (auto oldComplIter = oldWi.completionMap.begin(); oldComplIter != oldWi.completionMap.end(); ++oldComplIter){
                        size_t completionId = (*oldComplIter).first;
                        auto newComplIter = newWi.completionMap.find(completionId);

                        if(newComplIter != newWi.completionMap.end()){
                            size_t oldComplIndex = (*oldComplIter).second;
                            size_t newComplIndex = (*newComplIter).second;

                            size_t oldWellCompletionOffset = oldWi.wellCompletionOffset;
                            size_t newWellCompletionOffset = newWi.wellCompletionOffset;

                            perfrates_[newWellCompletionOffset + newComplIndex] = perfratesTemp[oldWellCompletionOffset + oldComplIndex];
                            perfpress_[newWellCompletionOffset + newComplIndex] = perfpressTemp[oldWellCompletionOffset + oldComplIndex];
                        }
                    }
                }
            }

            wellIndexMap_ = newMap;
        }

        void setBhpValue(size_t index, double value){
            if (index < bhp_.size()){
                bhp_[index] = value;
            }
        }

        void setBhpValue(std::string wellName, double value){
            int index = getWellMapIndex(wellName);

            if(index > -1){
                bhp_[index] = value;
            }
        }

        double getBhpValue(size_t index) const{
            double retVal;

            if (index < bhp_.size()){
                retVal = bhp_[index];
            }

            return retVal;
        }

        double getBhpValue(std::string wellName) const{
            double retVal;
            int index = getWellMapIndex(wellName);

            if(index > -1){
                retVal = bhp_[index];
            }

            return retVal;
        }

        void setTemperatureValue(size_t index, double value){
            if (index < temperature_.size()){
                temperature_[index] = value;
            }
        }

        void setTemperatureValue(std::string wellName, double value){
            int index = getWellMapIndex(wellName);

            if(index > -1){
                temperature_[index] = value;
            }
        }

        double getTemperatureValue(size_t index) const{
            double retVal;

            if (index < temperature_.size()){
                retVal = temperature_[index];
            }

            return retVal;
        }

        double getTemperatureValue(std::string wellName) const{
            double retVal;
            int index = getWellMapIndex(wellName);

            if(index > -1){
                retVal = temperature_[index];
            }

            return retVal;
        }

        void setWellRatesValue(size_t index, double value){
            if (index < wellrates_.size()){
                wellrates_[index] = value;
            }
        }

        void setWellRatesValue(std::string wellName, size_t phase, size_t numPhases, double value){
            int index = getWellMapIndex(wellName);

            if((index > -1) && (phase <= numPhases) && (phase > 0)){
                wellrates_[(index * numPhases) + (phase - 1)] = value;
            }
        }

        double getWellRatesValue(size_t index) const{
            double retVal;

            if (index < wellrates_.size()){
                retVal = wellrates_[index];
            }

            return retVal;
        }

        double getWellRatesValue(std::string wellName, size_t phase, size_t numPhases) const{
            double retVal;
            int index = getWellMapIndex(wellName);

            if((index > -1) && (phase <= numPhases) && (phase > 0)){
                retVal = wellrates_[(index * numPhases) + (phase - 1)];
            }

            return retVal;
        }

        void setPerfRatesValue(size_t index, double value){
            if (index < perfrates_.size()){
                perfrates_[index] = value;
            }
        }

        void setPerfRatesValue(std::string wellName, double value){
            int index = getWellMapIndex(wellName);

            if(index > -1){
                perfrates_[index] = value;
            }
        }

        double getPerfRatesValue(size_t index) const{
            double retVal;

            if (index < perfrates_.size()){
                retVal = perfrates_[index];
            }

            return retVal;
        }

        double getPerfRatesValue(std::string wellName) const{
            double retVal;
            int index = getWellMapIndex(wellName);

            if(index > -1){
                retVal = perfrates_[index];
            }

            return retVal;
        }

        void setPerfPressValue(size_t index, double value){
            if (index < perfpress_.size()){
                perfpress_[index] = value;
            }
        }

        void setPerfPressValue(std::string wellName, double value){
            int index = getWellMapIndex(wellName);

            if(index > -1){
                perfpress_[index] = value;
            }
        }

        double getPerfPressValue(size_t index) const{
            double retVal;

            if (index < perfpress_.size()){
                retVal = perfpress_[index];
            }

            return retVal;
        }

        double getPerfPressValue(std::string wellName) const{
            double retVal;
            int index = getWellMapIndex(wellName);

            if(index > -1){
                retVal = perfpress_[index];
            }

            return retVal;
        }

        /// One bhp pressure per well.
        std::vector<double>& bhp() { return bhp_; }
        const std::vector<double>& bhp() const { return bhp_; }

        /// One temperature per well.
        std::vector<double>& temperature() { return temperature_; }
        const std::vector<double>& temperature() const { return temperature_; }

        /// One rate per well and phase.
        std::vector<double>& wellRates() { return wellrates_; }
        const std::vector<double>& wellRates() const { return wellrates_; }

        /// One rate per well connection.
        std::vector<double>& perfRates() { return perfrates_; }
        const std::vector<double>& perfRates() const { return perfrates_; }

        /// One pressure per well connection.
        std::vector<double>& perfPress() { return perfpress_; }
        const std::vector<double>& perfPress() const { return perfpress_; }

    private:
        std::vector<double> bhp_;
        std::vector<double> temperature_;
        std::vector<double> wellrates_;
        std::vector<double> perfrates_;
        std::vector<double> perfpress_;

        std::map<std::string, WellIndex> wellIndexMap_;

        int getWellMapIndex(std::string wellName) const{
            int retVal = -1;
            auto iter = wellIndexMap_.find(wellName);

            if(iter != wellIndexMap_.end()){
                const WellIndex& wi = (*iter).second;
                retVal = wi.wellNumber;
            }

            return retVal;
        }

        int getNumCompletions(std::vector<WellPtr>& wells, int reportStep) const{
            int numCompl = 0;

            for(size_t i = 0; i < wells.size(); i++){
                int numComplForWell = getNumCompletions(wells[i], reportStep);
                numCompl += numComplForWell;
            }

            return numCompl;
        }

        int getNumCompletions(WellPtr& well, int reportStep) const{
            int numComplForWell = 0;
            CompletionSetConstPtr complSet = well->getCompletions(reportStep);

            for (size_t i = 0; i < complSet->size(); i++){
                WellCompletion::StateEnum complState = complSet->get(i)->getState();

                if(complState == WellCompletion::StateEnum::OPEN){
                    numComplForWell++;
                }
            }

            return numComplForWell;
        }   

        void fillCompletionsMap(EclipseGridConstPtr& grid, WellPtr& well, WellIndex& wi, int reportStep) const{
            CompletionSetConstPtr complSet = well->getCompletions(reportStep);
            size_t open_completions_count = 0;

            for (size_t idx = 0; idx < complSet->size(); idx++){
                CompletionConstPtr completion = complSet->get(idx);

                if(completion->getState() == WellCompletion::StateEnum::OPEN){
                    int i = completion->getI();
                    int j = completion->getJ();
                    int k = completion->getK();
                    size_t globalIndex = grid->getGlobalIndex(i, j, k);
                    std::pair <size_t, size_t> p(globalIndex, open_completions_count);
                    wi.completionMap.insert(p);
                    ++open_completions_count;
                }
            }
        }
    };

} // namespace Opm

#endif // OPM_WELLSTATE_HEADER_INCLUDED
