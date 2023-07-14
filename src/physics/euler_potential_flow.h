#ifndef __EULER_POTENTIAL_FLOW__
#define __EULER_POTENTIAL_FLOW__

#include "euler.h"
#include "parameters/parameters_navier_stokes.h"

namespace PHiLiP {
namespace Physics {

/// EulerPotentialFlow equations. Derived from Euler, which is derived from PhysicsBase. 
template <int dim, int nstate, typename real>
class EulerPotentialFlow : public Euler <dim, nstate, real>
{
public:
    using two_point_num_flux_enum = Parameters::AllParameters::TwoPointNumericalFlux;
    
    /// Constructor
    EulerPotentialFlow ( 
        const Parameters::AllParameters *const                    parameters_input,
        const double                                              ref_length,
        const double                                              gamma_gas,
        const double                                              mach_inf,
        const double                                              angle_of_attack,
        const double                                              side_slip_angle,
        std::shared_ptr< ManufacturedSolutionFunction<dim,real> > manufactured_solution_function = nullptr,
        const two_point_num_flux_enum                             two_point_num_flux_type = two_point_num_flux_enum::KG,
        const bool                                                has_nonzero_diffusion = false,
        const bool                                                has_nonzero_physical_source = true);

    /// Physical source term
    std::array<real,nstate> physical_source_term (
        const dealii::Point<dim,real> &pos,
        const std::array<real,nstate> &/*conservative_solution*/,
        const std::array<dealii::Tensor<1,dim,real>,nstate> &/*solution_gradient*/,
        const dealii::types::global_dof_index /*cell_index*/) const override;

private:
    // checking if node is within object 
    bool within_physical_object (const dealii::Point<dim,real> &pos) const;

    // returns boundaries of airfoil
    std::array<real,(2 * dim)> NACA_series (const dealii::Point<dim,real> &pos) const;

    // returns source term for drag (x)
    template<typename real2> 
    real2 drag_calc () const;

    // returns source term for lift (y)
    template<typename real2> 
    real2 lift_calc () const;

};

} // Physics namespace
} // PHiLiP namespace

#endif
