/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "ScalarField.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates a constant scalar field given params for initial
//! matter config
class InitialScalarData
{
  public:
    struct params_t
    {
        double eta;
        double lambda;
        double R0;
        std::array<double, CH_SPACEDIM> center;
    };

    //! The constructor for the class
    InitialScalarData(const params_t a_params, const double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

    Coordinates<data_t> coords(current_cell, m_dx, m_params.center);


    current_cell.store_vars(compute_phi(coords),c_phi);
    current_cell.store_vars(0.0 ,c_Pi);
    }

  

  protected:
    const double m_dx;
    const params_t m_params;


  // Compute the value of phi at the current point
template <class data_t>
data_t compute_phi(Coordinates<data_t> coords) const
{
    data_t x = coords.x;
    data_t y = coords.y;
    data_t rho = sqrt(x*x+y*y);
    data_t out_phi = m_params.eta* tanh(sqrt(2.0*m_params.lambda)*m_params.eta*(rho-m_params.R0)/2.0);

    return out_phi;
}
};



#endif /* INITIALSCALARDATA_HPP_ */
