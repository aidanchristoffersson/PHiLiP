#include <cmath>
#include <vector>

#include "ADTypes.hpp"

#include "euler.h"

#include "euler_potential_flow.h"

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
    std::shared_ptr< ManufacturedSolutionFunction<dim,real> > manufactured_solution_function,
    const two_point_num_flux_enum                             two_point_num_flux_type,
    const bool                                                has_nonzero_diffusion,
    const bool                                                has_nonzero_physical_source)
    : Euler<dim,nstate,real>(parameters_input,
                             ref_length, 
                             gamma_gas, 
                             mach_inf, 
                             angle_of_attack, 
                             side_slip_angle, 
                             manufactured_solution_function,
                             two_point_num_flux_type,
                             has_nonzero_diffusion,  // has_nonzero_diffusion = false
                             has_nonzero_physical_source), 	// has_nonzero_physical_source = true
    NACA_code(parameters_input->potential_source_param.NACA_code)
{
    static_assert(nstate==dim+2, "Physics::EulerPotentialFlow() should be created with nstate=dim+2");
    // no tasks to do here yet
}

// template <int dim, int nstate, typename real>
// std::array<real, 4> EulerPotentialFlow<dim, nstate, real>
// ::NACA_parameters() const
// {
//     const std::array<real, 4> NACA_code { {4, 4, 12, 1.0} };     // NACA {a, b, cd, chord}

//     return NACA_code;
// }


template <int dim, int nstate, typename real>
real EulerPotentialFlow<dim, nstate, real>
::NACA_volume () const
{
    // // extract NACA code values with: NACA {a, b, cd, chord}
    // const std::array<real, 4> NACA_code = NACA_parameters();
    const real chord {this->NACA_code[3]};

    const real max_thickness {this->NACA_code[2] / 100};    // max thickness as a fraction of chord (cd = 100t)

    real airfoil_volume {chord};    // initialize volume of airfoil

    if constexpr(dim > 1)    // 2D or greater -> initialize area
    {
        // 2D - volume (unit depth)
        airfoil_volume = 0.685083333 * chord * chord * max_thickness;
    }
    if constexpr(dim > 2)    // 3D -> extrude to form volume
    {
        real zMin {-10};
        real zMax {10};

        airfoil_volume = airfoil_volume * (zMax - zMin);
    }

    return airfoil_volume;
}

template <int dim, int nstate, typename real>
std::array<real,(2 * dim)> EulerPotentialFlow<dim, nstate, real>
::NACA_series (const dealii::Point<dim,real> &pos) const // airfoil situated from x_offset -> x_offset + c (default 0 -> 1)
{
    // NACA parameters defined here (or taken in from parameters file): NACA {a, b, cd, chord}
    // const std::array<real, 4> NACA_code = NACA_parameters();
    const real chord {this->NACA_code[3]};

    const real x_offset {0};

    std::array<real,(2 * dim)> airfoil_boundaries {};   // [0]: xMin, [1]: xMax, [2]: yMin, [3]: yMax, [4]: zMin, [5]: zMax, [6]: volume

    // 1D - maximum and minimum x within chord
    airfoil_boundaries[0] = x_offset;
    airfoil_boundaries[1] = x_offset + chord;


    if constexpr(dim>1) // 2D 
    {
        // change of coordinates
        real s = (pos[0] - x_offset); // change of coordinates

        // define shape variables
        const real max_camber_val {this->NACA_code[0] / 100};   // value of maximum camber as fraction of cord (a = 100m)
        const real p_max_camber {this->NACA_code[1] / 10};      // location of maximum camber as fraction of cord (b = 10p)
        const real max_thickness {this->NACA_code[2] / 100};    // max thickness as a fraction of chord (cd = 100t)
        real c_r {s / chord};   // chord_ratio
        real y_c {};    // camber line
        
        // -----------------------------------------------------------
        if (s < 0 || s > chord) // if outside domain (I functionally end up checking this twice)
        {
            return airfoil_boundaries;
        }
        else if (s <= p_max_camber)
        {
            y_c = max_camber_val * chord * ((2 * p_max_camber * c_r) - (c_r * c_r)) / (p_max_camber * p_max_camber);
        }
        else if (s > p_max_camber)
        {
            y_c = max_camber_val * chord * ((1 - (2 * p_max_camber)) + (2 * p_max_camber * c_r) - (c_r * c_r)) / ((1 - p_max_camber) * (1 - p_max_camber));
        }

        real y_t = 5 * chord * max_thickness * ((0.2969 * sqrt(c_r)) - (0.1260 * c_r) - (0.3516 * c_r * c_r) + (0.2843 * pow(c_r, 3)) - (0.1015 * pow(c_r, 4)));

        // -----------------------------------------------------------

        // 2D - maximum and minimum y at specific chord location
        airfoil_boundaries[2] = y_c - y_t;
        airfoil_boundaries[3] = y_c + y_t;

        // check if y within bounds, if not, do not check 3rd dimension
        real y = pos[1];
        if (y < airfoil_boundaries[2] || y > airfoil_boundaries[3])
        {
            return airfoil_boundaries;
        }

    }

    // 3D - check z within extrusion
    if constexpr(dim>2) 
    {
        real zMin {-10};
        real zMax {10};

        airfoil_boundaries[4] = zMin;
        airfoil_boundaries[5] = zMax;
    }

    return airfoil_boundaries;
}


// evaluates if node within potential body
template <int dim, int nstate, typename real>
bool EulerPotentialFlow<dim,nstate,real>
::within_physical_object (
		const dealii::Point<dim,real> &pos) const
{
    bool locationStatus {false};

	// extracting boundary to airfoil at those coordinates
    std::array<real,2 * dim> airfoil_boundaries = NACA_series(pos);

	real x = pos[0];
    locationStatus = (x >= airfoil_boundaries[0]) && (x <= airfoil_boundaries[1]);    // x in [xMin, xMax]

	if constexpr(dim>1) // use 3 different NACA functions, depending on number of arguments. x -> within chord? x, y -> within area x,y,z -> within volume
    {
        if (locationStatus) // I don't like this but I also can't put in constexpr, is there a better way?
        {
            real y = pos[1];
            locationStatus = (y >= airfoil_boundaries[2]) && (y <= airfoil_boundaries[3]);  // y in [yMin, yMax]
        }
    }

	if constexpr(dim>2) 
    {
        if (locationStatus)
        {
            real z = pos[2]; 
            locationStatus = (z >= airfoil_boundaries[4]) && (z <= airfoil_boundaries[5]);  // z in [zMin, zMax]
        }
    }

	return locationStatus;
}

template <int dim, int nstate, typename real>
template <typename real2>
real2 EulerPotentialFlow<dim,nstate,real>
::lift_calc (
        const dealii::Point<dim,real> &/*pos*/,
        const std::array<real,nstate> &/*conservative_soln*/) const     // conservative_soln: [0]: density, [1]: momentum x, [2]: momentum y, [nstate - 1]: energy
{
    real2 coeff_lift {};
    real2 lift_force;
    // const std::array<real2, 4> NACA_code = NACA_parameters();
    const real2 chord {this->NACA_code[3]};

    const double pi = 4 * atan(1);
    const double angle_of_attack = this->angle_of_attack;

    // freestream properties
    const double mach_inf = this->mach_inf;
    const double density_inf = this->density_inf;
    const double pressure_inf = this->pressure_inf;
    const double gamma_gas = this->gam;
    const double freestream_velocity {mach_inf * sqrt(gamma_gas * pressure_inf / density_inf)};

    

    #if 0   ///// 2D flat plate
    coeff_lift = 2 * pi * angle_of_attack;
    #endif

    #if 1   ///// 2D airfoil

    // define shape variables
    const real2 max_camber {this->NACA_code[0] / 100};   // value of maximum camber as fraction of cord (a = 100m)
    const real2 p_max_camber {this->NACA_code[1] / 10};      // location of maximum camber as fraction of cord (b = 10p)
    const real2 theta_max_p { acos(1 - 2 * p_max_camber) };  // location of maximum camber in radians

    if (max_camber < 0.01) { // case of 0 camber is treated as a flat plate
        coeff_lift = 2 * pi * angle_of_attack;
    }
    else {
        // geometric fourier series representation (first two terms are satisfactory for overall force calculations)
        const real2 A0 = ((max_camber / (pi * p_max_camber * p_max_camber)) * ((2 * p_max_camber - 1) * theta_max_p + sin(theta_max_p))
                + ((max_camber / (pi * (1 - p_max_camber) * (1 - p_max_camber) )) * ((2 * p_max_camber - 1) * (pi - theta_max_p) - sin(theta_max_p))) );

        const real2 A1 = ( ((2 * max_camber / (pi * p_max_camber * p_max_camber)) * ((2 * p_max_camber - 1) * sin(theta_max_p) + 
            (sin(2 * theta_max_p) / 4) + (theta_max_p / 2))) - ((2 * max_camber / (pi * (1 - p_max_camber) * (1 - p_max_camber) )) * 
            ((2 * p_max_camber - 1) * sin(theta_max_p) + (sin(theta_max_p) / 4) - ((pi - theta_max_p) / 2))) );

        coeff_lift = pi * (A1 - 2 * A0) + (2 * pi * angle_of_attack);
    }
    #endif

    // lift force = cL * rho * U^2 * area -> unit depth
    lift_force = coeff_lift * density_inf * freestream_velocity * freestream_velocity * chord;

    return lift_force;
}

template <int dim, int nstate, typename real>
template <typename real2>
real2 EulerPotentialFlow<dim,nstate,real>
::drag_calc (
        const dealii::Point<dim,real> &/*pos*/,
        const std::array<real,nstate> &/*conservative_soln*/) const
{
    real2 drag {0};

    return drag;
}



// returns the body force at a specific quadrature node and position
template <int dim, int nstate, typename real>
std::array<real,nstate> EulerPotentialFlow<dim,nstate,real>
::physical_source_term (
        const dealii::Point<dim,real> &pos,
        const std::array<real,nstate> &conservative_soln,
        const std::array<dealii::Tensor<1,dim,real>,nstate> &/*solution_gradient*/,
        const dealii::types::global_dof_index /*cell_index*/) const
{
	
    std::array<real,nstate> physical_source;
    real airfoil_volume = NACA_volume();

	if (within_physical_object(pos)) 
    {
        // density
        physical_source[0] = 0;    

        // momentum
        physical_source[1] = - drag_calc<real>(pos, conservative_soln) / airfoil_volume;   // x

        if constexpr(dim>1)
	    {
	        physical_source[2] = - lift_calc<real>(pos, conservative_soln) / airfoil_volume;   // y
        }

        if constexpr(dim>2)
        {
            physical_source[3] = 0;     // z
        }

        // energy
        physical_source[nstate - 1] = 0;
        
    }

    else 
    {
        for (int i=0; i<nstate; i++) 
        {
            physical_source[i] = 0;
        }
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

