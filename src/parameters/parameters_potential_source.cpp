#include "parameters/parameters_potential_source.h"

namespace PHiLiP {
namespace Parameters {

// Potential source geometry inputs
PotentialSourceParam::PotentialSourceParam () {}

void PotentialSourceParam::declare_parameters (dealii::ParameterHandler &prm)
{
    prm.enter_subsection("potential_source");
    {
        prm.declare_entry("trailing_edge_serration_frequency", "0.015",
                          dealii::Patterns::Double(1e-15, 10000000),
                          "Distance between consecutive trailing edge serrations.");
        prm.declare_entry("half_trailing_edge_serration_length", "0.015",
                          dealii::Patterns::Double(1e-15, 10000000),
                          "Half the length of trailing edge serrations.");
        prm.declare_entry("trailing_edge_serration_thickness", "0.001",
                  dealii::Patterns::Double(1e-15, 1000),
                  "Thickness of trailing edge serrations.");
        prm.declare_entry("trailing_edge_serration_effective_length_factor", "0.5",
                          dealii::Patterns::Double(1e-15, 1),
                          "Effective length factor for trailing edge serrations.");
        prm.declare_entry("angle_of_serration", "0.0",
                          dealii::Patterns::Double(-180, 180),
                          "Angle of trailing edge flap in degrees (gamma).");
        prm.declare_entry("viscous_drag","true",
                              dealii::Patterns::Bool(),
                              "Set as true by default (i.e. apply viscous drag term). " 
                              "If true, includes the contribution of the viscous drag.");

        prm.declare_entry("NACA", "[0, 0, 12, 1.0]",
                              dealii::Patterns::Anything(),
                              "Stores a vector containing the NACA code. ",
                              "Formatted in square brackets and seperated by commas, eg. \"[a, b, cd, chord]\"");
}
    prm.leave_subsection();
}

void PotentialSourceParam ::parse_parameters (dealii::ParameterHandler &prm)
{
    prm.enter_subsection("potential_source");
    {
        TES_frequency   = prm.get_double("trailing_edge_serration_frequency");
        TES_h           = prm.get_double("half_trailing_edge_serration_length");
        TES_thickness   = prm.get_double("trailing_edge_serration_thickness");
        TES_effective_length_factor = prm.get_double("trailing_edge_serration_effective_length_factor");
        const double pi = atan(1.0) * 4.0;
        TES_gamma       = prm.get_double("angle_of_serration") * pi/180.0;
        use_viscous_drag = prm.get_bool("viscous_drag");

        std::string NACA_string = prm.get("NACA");

        // removing formatting
        std::string removeChars = "[]";
        for(unsigned int i = 0; i < removeChars.length(); ++i)
            NACA_string.erase(std::remove(NACA_string.begin(), NACA_string.end(), removeChars.at(i)), NACA_string.end());
        std::vector<std::string> NACA_string_vector = dealii::Utilities::split_string_list(NACA_string);

        // pushing into a vector
        for(unsigned int i = 0; i < 4; i++)
            NACA_code[i] = dealii::Utilities::string_to_double(NACA_string_vector[i]);
    }
    prm.leave_subsection();
}

} // Parameters namespace
} // PHiLiP namespace
