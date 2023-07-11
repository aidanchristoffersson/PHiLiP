#include <cmath>
#include <vector>

#include "ADTypes.hpp"

#include "euler.h"

#include "euler_source_term.h"

namespace PHiLiP {
namespace Physics {

template <int dim, int nstate, typename real>
EulerPotentialFlow<dim, nstate, real>::EulerPotentialFlow( 
    const Parameters::AllParameters *const                    parameters_input,
    const double                                              ref_length,
    const double                                              gamma_gas,
    const double                                              mach_inf,
    const double                                              angle_of_attack,
    const double                                              side_slip_angle,
    const two_point_num_flux_enum                             two_point_num_flux_type)
    : Euler<dim,nstate,real>(parameters_input,
                             ref_length, 
                             gamma_gas, 
                             mach_inf, 
                             angle_of_attack, 
                             side_slip_angle, 
                             manufactured_solution_function,
                             two_point_num_flux_type,
                             false,  // has_nonzero_diffusion = false
                             true) 	// has_nonzero_physical_source = true
{
    static_assert(nstate==dim+2, "Physics::EulerPotentialFlow() should be created with nstate=dim+2");
    // Nothing to do here so far
}

bool EulerPotentialFlow<dim,nstate,real>
::within_physical_object (
		const dealii::Point<dim,real> &pos) const
{
	// extracting coordinate
	real x = pos[0];
	if constexpr(dim>1) real y = pos[1];
	if constexpr(dim>2) real z = pos[2];

	// radius
	real r = sqrt(x*x + y*y);

	return (r < 5);
};


// double EulerPotentialFlow<dim,nstate,real>
// ::lift_calc () {};

// double EulerPotentialFlow<dim,nstate,real>
// ::drag_calc () {};


// returns the body force at a specific quadrature node and position
template <int dim, int nstate, typename real>
std::array<real,nstate> EulerPotentialFlow<dim,nstate,real>
::physical_source_term (
        const dealii::Point<dim,real> &pos,
        const std::array<real,nstate> &/*conservative_soln*/,
        const std::array<dealii::Tensor<1,dim,real>,nstate> &/*solution_gradient*/,
        const dealii::types::global_dof_index /*cell_index*/) const
{
	
    std::array<real,nstate> physical_source;
   	// for density [0]
	// for momentums do [i+1] loop d to 0,dim
	// energy [nstate-1]
    
    for (int i=0; i<nstate; i++) {
        physical_source[i] = 0;
    }

    return physical_source;
}

// Instantiate explicitly
template class EulerPotentialFlow < PHILIP_DIM, PHILIP_DIM+2, double     >;
template class EulerPotentialFlow < PHILIP_DIM, PHILIP_DIM+2, FadType    >;
template class EulerPotentialFlow < PHILIP_DIM, PHILIP_DIM+2, RadType    >;
template class EulerPotentialFlow < PHILIP_DIM, PHILIP_DIM+2, FadFadType >;
template class EulerPotentialFlow < PHILIP_DIM, PHILIP_DIM+2, RadFadType >;


} // Physics namespace
} // PHiLiP namespace

