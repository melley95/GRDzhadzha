/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "FixedBGSimulationParametersBase.hpp"
#include "GRParmParse.hpp"

// Problem specific includes:
#include "InitialScalarData.hpp"
#include "InitialBGData.hpp"
#include "FLRW.hpp"


class SimulationParameters : public FixedBGSimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : FixedBGSimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {

        // Initial SF data
        pp.load("eta", initial_params.eta, 1.0);
        pp.load("lambda", initial_params.lambda, 1.0);
        pp.load("R0", initial_params.R0, 10.0);
        pp.load("scalar_center", initial_params.center, center);

        // Initial BG data
        pp.load("K0", initial_bg_params.K0, 0.0);
        pp.load("center", initial_bg_params.center, center);
        pp.load("center", flrw_params.center, center);
        flrw_params.K0 = initial_bg_params.K0;

 
    }

    void check_params()
    {
        
        
    }

    // Problem specific parameters
 
    // Collection of parameters necessary for the initial conditions
    InitialScalarData::params_t initial_params;
    InitialBGData::params_t initial_bg_params;
    FLRW::params_t flrw_params;
    // Collection of parameters necessary for the background metric
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
