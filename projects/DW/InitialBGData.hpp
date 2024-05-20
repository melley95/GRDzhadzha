/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALBGDATA_HPP_
#define INITIALBGDATA_HPP_

#include "ADMFixedBGVars.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "Cell.hpp"
#include "TensorAlgebra.hpp"

//! Class which creates an initial flat FLRW background

class InitialBGData
{
  public:
    struct params_t
    {
        double K0;
        std::array<double, CH_SPACEDIM> center;
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    //! The constructor for the class
    InitialBGData(const params_t a_params, const double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        Vars<data_t> vars;
    //    VarsTools::assign(vars, 0.);

  

        FOR(i) vars.gamma[i][i] = 1.;
        
  

    
        vars.lapse = 1.;

        
        data_t chi = TensorAlgebra::compute_determinant_sym(vars.gamma);
        chi = pow(chi, -1.0 / 2.0);
        current_cell.store_vars(chi, c_chi);
    
   
     
        vars.K = m_params.K0;
        FOR(i) vars.K_tensor[i][i] = (1.0/2.0)*vars.gamma[i][i]*m_params.K0; 
        
        
        data_t K = vars.K;
        current_cell.store_vars(K, c_K);

    //    current_cell.store_vars(vars);
    }

  protected:
    const params_t m_params;
    const double m_dx;
};

#endif /* INITIALBGDATA_HPP_ */
