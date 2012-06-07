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

#ifndef OPM_INCOMPTPFA_HEADER_INCLUDED
#define OPM_INCOMPTPFA_HEADER_INCLUDED


#include <vector>

struct UnstructuredGrid;
struct ifs_tpfa_data;
struct Wells;
struct FlowBoundaryConditions;

namespace Opm
{

    class LinearSolverInterface;

    /// Encapsulating a tpfa pressure solver for the incompressible-fluid case.
    /// Supports gravity, wells controlled by bhp or reservoir rates,
    /// boundary conditions and simple sources as driving forces.
    /// Rock compressibility can be included, but any nonlinear iterations
    /// are not handled in this class.
    /// Below we use the shortcuts D for the number of dimensions, N
    /// for the number of cells and F for the number of faces.
    class IncompTpfa
    {
    public:
	/// Construct solver.
	/// \param[in] g             A 2d or 3d grid.
	/// \param[in] permeability  Array of permeability tensors, the array
	///                          should have size N*D^2, if D == g.dimensions
	///                          and N == g.number_of_cells.
	/// \param[in] gravity       Gravity vector. If nonzero, the array should
	///                          have D elements.
	/// \param[in] linsolver     A linear solver.
        /// \param[in] wells         The wells argument. Will be used in solution, 
        ///                          is ignored if NULL
	IncompTpfa(const UnstructuredGrid& g,
		   const double* permeability,
		   const double* gravity,
                   LinearSolverInterface& linsolver,
                   const Wells* wells);

	/// Destructor.
	~IncompTpfa();

	/// Assemble and solve incompressible pressure system.
	/// \param[in]  totmob     Must contain N total mobility values (one per cell).
	///                        \f$totmob = \sum_{p} kr_p/mu_p\f$.
	/// \param[in]  omega      Must be empty if constructor gravity argument was null.
	///                        Otherwise must contain N mobility-weighted density values (one per cell).
	///                        \f$omega = \frac{\sum_{p} mob_p rho_p}{\sum_p rho_p}\f$.
	/// \param[in]  src        Must contain N source rates (one per cell).
	///                        Positive values represent total inflow rates,
	///                        negative values represent total outflow rates.
	/// \param[in]  wdp        Should contain the differences between
        ///                        well BHP and perforation pressures.
	///                        May be empty if there are no wells.
	/// \param[in]  bcs        If non-null, specifies boundary conditions.
	///                        If null, noflow conditions are assumed.
	/// \param[out] pressure   Will contain N cell-pressure values.
	/// \param[out] faceflux   Will contain F signed face flux values.
        /// \param[out] well_bhp   Will contain bhp values for each well passed
        ///                        in the constructor
        /// \param[out] well_rate  Will contain rate values for each well passed
        ///                        in the constructor
	void solve(const std::vector<double>& totmob,
		   const std::vector<double>& omega,
		   const std::vector<double>& src,
                   const std::vector<double>& wdp,
		   const FlowBoundaryConditions* bcs,
		   std::vector<double>& pressure,
		   std::vector<double>& faceflux, 
                   std::vector<double>& well_bhp,
                   std::vector<double>& well_rate);

        /// Assemble and solve pressure system with rock compressibility (assumed constant per cell).
        /// \param[in]  totmob     Must contain N total mobility values (one per cell).
        ///                        totmob = \sum_{p} kr_p/mu_p.
        /// \param[in]  omega      Must be empty if constructor gravity argument was null.
        ///                        Otherwise must contain N fractional-flow-weighted density
        ///                        values (one per cell).
        ///                        omega = \frac{\sum_{p} mob_p rho_p}{\sum_p mob_p}.
        /// \param[in]  src        Must contain N source rates (one per cell).
        ///                        Positive values represent total inflow rates,
        ///                        negative values represent total outflow rates.
	/// \param[in]  wdp        Should contain the differences between
        ///                        well BHP and perforation pressures.
	///                        May be empty if there are no wells.
        /// \param[in]  bcs        If non-null, specifies boundary conditions.
        ///                        If null, noflow conditions are assumed.
        /// \param[in]  porevol    Must contain N pore volumes.
        /// \param[in]  rock_comp  Must contain N rock compressibilities.
        ///                        rock_comp = (d poro / d p)*(1/poro).
        /// \param[in]  dt         Timestep.
        /// \param[out] pressure   Will contain N cell-pressure values.
        /// \param[out] faceflux   Will contain F signed face flux values.
        /// \param[out] well_bhp   Will contain bhp values for each well passed
        ///                        in the constructor
        /// \param[out] well_rate  Will contain rate values for each well passed
        ///                        in the constructor
        void solve(const std::vector<double>& totmob,
                   const std::vector<double>& omega,
                   const std::vector<double>& src,
                   const std::vector<double>& wdp,
                   const FlowBoundaryConditions* bcs,
                   const std::vector<double>& porevol,
                   const std::vector<double>& rock_comp,
                   const double dt,
                   std::vector<double>& pressure,
                   std::vector<double>& faceflux,
                   std::vector<double>& well_bhp,
                   std::vector<double>& well_rate);

        void solveIncrement(const std::vector<double>& totmob,
			    const std::vector<double>& omega,
			    const std::vector<double>& src,
			    const std::vector<double>& wdp,
			    const FlowBoundaryConditions* bcs,
			    const std::vector<double>& porevol,
			    const std::vector<double>& rock_comp,
			    const std::vector<double>& prev_pressure,
			    const std::vector<double>& initial_porevol,
			    const double dt,
			    std::vector<double>& pressure_increment);

	void computeFaceFlux(const std::vector<double>& totmob,
					 const std::vector<double>& omega,
					 const std::vector<double>& src,
					 const std::vector<double>& wdp,
					 const FlowBoundaryConditions* bcs,
					 std::vector<double>& pressure,
					 std::vector<double>& faceflux,
					 std::vector<double>& well_bhp,
					 std::vector<double>& well_rate);

        /// Expose read-only reference to internal half-transmissibility.
        const ::std::vector<double>& getHalfTrans() const { return htrans_; }

        /// Set tolerance for the linear solver.
        /// \param[in] tol         tolerance value
        void setTolerance(const double tol);

        /// Get tolerance of the linear solver.
        /// \param[out] tolerance value
        double getTolerance() const;


    private:
	const UnstructuredGrid& grid_;
        LinearSolverInterface& linsolver_;
	::std::vector<double> htrans_;
	::std::vector<double> trans_ ;
	::std::vector<double> gpress_;
	::std::vector<double> gpress_omegaweighted_;
        
        const struct Wells* wells_;
	struct ifs_tpfa_data* h_;
    };

} // namespace Opm

#endif // OPM_INCOMPTPFA_HEADER_INCLUDED
