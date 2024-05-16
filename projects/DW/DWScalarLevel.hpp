/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DWSCALARLEVEL_HPP_
#define DWSCALARLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "ScalarField.hpp"
#include "ScalarPotential.hpp"



class DWScalarLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<DWScalarLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // Typedef for scalar field
    typedef ScalarField<ScalarPotential> ScalarFieldWithPotential;

    //! Initialize data for the field and metric variables
    virtual void initialData();

    //! RHS routines used at each RK4 step
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time);

    //! To do after each timestep
    virtual void specificPostTimeStep();

    //! Tell Chombo how to tag cells for regridding
    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state);
};

#endif /* DWSCALARLEVEL_HPP_ */
