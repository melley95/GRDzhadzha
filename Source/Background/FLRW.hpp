/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FLRW_HPP_
#define FLRW_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the initial conditions for a Kerr Schild BH
//! https://arxiv.org/pdf/gr-qc/9805023.pdf
//! https://arxiv.org/pdf/2011.07870.pdf

class FLRW
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double K0; //t0 = 3/(2K0)
        std::array<double, CH_SPACEDIM> center; 
                           

    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    const params_t m_params;

    const double m_time;
    const double m_dx;
 

   FLRW(params_t a_params, double a_dx, double a_time) : m_params(a_params), m_dx(a_dx), m_time(a_time){}

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

       

        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        Vars<data_t> metric_vars;
        compute_metric_background(metric_vars, coords);

        
        data_t chi = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        chi = pow(chi, -1.0 / 3.0);
        current_cell.store_vars(chi, c_chi);
        

        data_t K = metric_vars.K;
        current_cell.store_vars(K, c_K);
    }

     template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Coordinates<data_t> &coords) const
    {

        vars.lapse = 1.0;
        vars.K = 0.0;
        FOR(i) vars.shift[i] = 0.0;
        FOR(i,j) vars.gamma[i][j] = 0.0;
        FOR(i,j) vars.K_tensor[i][j] = 0.0;
        
        FOR(i,j, k) vars.d1_gamma[i][j][k] = 0.0;
        FOR(i) vars.d1_lapse[i] = 0.0;
        FOR(i,j) vars.d1_shift[i][j] = 0.0;

        

        if (m_params.K0 == 0.0){

        FOR(i) vars.gamma[i][i] = 1.0; //assuming a0 = 1

       
        }
        else{
            
        double t0 = -3.0/(2.0*m_params.K0);
        vars.K = -1.5*pow((m_time+t0), -1);
        FOR(i) vars.gamma[i][i] = pow(((m_time+t0)/t0), 1.0); //assuming a0 = 1

        FOR(i) vars.K_tensor[i][i] = (1.0/3.0)*vars.K*vars.gamma[i][i]; 
        
        }
       

    }

};

#endif /* FLRW_HPP_ */
