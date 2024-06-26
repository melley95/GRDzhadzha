/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCALARPOTENTIAL_HPP_
#define SCALARPOTENTIAL_HPP_

#include "simd.hpp"

class ScalarPotential
{
  protected:
    //    const double m_mu;
    const InitialScalarData::params_t m_initial_params;

  public:
    //! The constructor
    ScalarPotential(const InitialScalarData::params_t a_initial_params)
        : m_initial_params(a_initial_params)
    {
    }

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        
        V_of_phi = m_initial_params.lambda*0.25*pow(pow(vars.phi, 2.0)-pow(m_initial_params.eta, 2.0), 2.0);

      
        dVdphi = m_initial_params.lambda*vars.phi*(pow(vars.phi, 2.0)-pow(m_initial_params.eta, 2.0));
    }
};

#endif /* SCALARPOTENTIAL_HPP_ */
