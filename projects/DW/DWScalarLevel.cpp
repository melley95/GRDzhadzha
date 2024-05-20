/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "DWScalarLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// For tag cells
#include "ChiTaggingCriterion.hpp"

// Problem specific includes
#include "InitialScalarData.hpp"
#include "InitialBGData.hpp"
#include "MatterEvolution.hpp"
#include "ScalarField.hpp"
#include "ScalarPotential.hpp"
#include "FLRW.hpp"

// Initial data for field and metric variables
void DWScalarLevel::initialData()
{
    CH_TIME("DWScalarLevel::initialData");
    if (m_verbosity)
        pout() << "DWScalarLevel::initialData " << m_level << endl;

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then initial conditions for fields
    SetValue set_zero(0.0);
    InitialBGData initial_bg(m_p.initial_bg_params, m_dx);
    auto compute_pack = make_compute_pack(set_zero, initial_bg);
    BoxLoops::loop(compute_pack, m_state_diagnostics, m_state_diagnostics,
                   SKIP_GHOST_CELLS);

    InitialScalarData initial_sf(m_p.initial_params, m_dx);
    BoxLoops::loop(initial_sf, m_state_new, m_state_new, FILL_GHOST_CELLS);   
}
// Things to do in RHS update, at each RK4 step
void DWScalarLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                        const double a_time)
{
    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    // We don't want undefined values floating around in the constraints so
    // zero these
    ScalarPotential potential(m_p.initial_params);
    ScalarFieldWithPotential scalar_field(potential);
    FLRW flrw(m_p.flrw_params, m_dx, m_time);
 //   BoxLoops::loop(flrw, m_state_new, m_state_new, INCLUDE_GHOST_CELLS); //Do I need this?
    MatterEvolution<ScalarFieldWithPotential, FLRW> my_matter(
        scalar_field, flrw, m_p.sigma, m_dx, m_p.center);
    BoxLoops::loop(my_matter, a_soln, a_rhs, SKIP_GHOST_CELLS);

  
}

void DWScalarLevel::specificPostTimeStep()
{
    // Check for nans on every level
    if (m_p.nan_check){
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());
    }
    FLRW flrw(m_p.flrw_params, m_dx, m_time);
    BoxLoops::loop(flrw, m_state_diagnostics, m_state_diagnostics,
                   SKIP_GHOST_CELLS);          //m_state_new?       


}

void DWScalarLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                                const FArrayBox &current_state)
{
    BoxLoops::loop(ChiTaggingCriterion(m_dx), current_state, tagging_criterion);
}
