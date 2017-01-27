/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MContraction/WeakHamiltonianNonEye.cc

Copyright (C) 2017

Author: Andrew Lawson    <andrew.lawson1991@gmail.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#include <Grid/Hadrons/Modules/MContraction/WeakHamiltonianNonEye.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MContraction;

/*
 * Weak Hamiltonian current-current contractions, Non-Eye-type.
 * 
 * These contractions are generated by the Q1 and Q2 operators in the physical
 * basis (see e.g. Fig 3 of arXiv:1507.03094).
 * 
 * Schematic:     
 *            q2             q3          |           q2              q3
 *          /--<--¬       /--<--¬        |        /--<--¬         /--<--¬       
 *         /       \     /       \       |       /       \       /       \      
 *        /         \   /         \      |      /         \     /         \     
 *       /           \ /           \     |     /           \   /           \    
 *    i *             * H_W         *  f |  i *             * * H_W         * f 
 *      \             *             |    |     \           /   \           /
 *       \           / \           /     |      \         /     \         /    
 *        \         /   \         /      |       \       /       \       /  
 *         \       /     \       /       |        \-->--/         \-->--/      
 *          \-->--/       \-->--/        |          q1               q4 
 *            q1             q4          |
 *                Connected (C)          |                 Wing (W)
 *
 * C: trace(q1*adj(q2)*g5*gL[mu]*q3*adj(q4)*g5*gL[mu])
 * W: trace(q1*adj(q2)*g5*gL[mu])*trace(q3*adj(q4)*g5*gL[mu])
 * 
 */

/******************************************************************************
 *                  TWeakHamiltonianNonEye implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TWeakHamiltonianNonEye::TWeakHamiltonianNonEye(const std::string name)
: Module<WeakHamiltonianPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TWeakHamiltonianNonEye::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

std::vector<std::string> TWeakHamiltonianNonEye::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TWeakHamiltonianNonEye::setup(void)
{

}

// execution ///////////////////////////////////////////////////////////////////
void TWeakHamiltonianNonEye::execute(void)
{
    XmlWriter             writer(par().output);
    PropagatorField &q1 = *env().template getObject<PropagatorField>(par().q1);
    PropagatorField &q2 = *env().template getObject<PropagatorField>(par().q2);
    PropagatorField &q3 = *env().template getObject<PropagatorField>(par().q3);
    PropagatorField &q4 = *env().template getObject<PropagatorField>(par().q4);
    SpinMatrix g5       = makeGammaProd(Ns*Ns - 1); // Soon to be deprecated.
    LatticeComplex        expbuf(env().getGrid());
    std::vector<TComplex> corrbuf;
    std::vector<Result>   result;

    PropagatorField              tmp1(env().getGrid());
    LatticeComplex               tmp2(env().getGrid());
    std::vector<PropagatorField> C_i_side_loop(Nd, tmp1);
    std::vector<PropagatorField> C_f_side_loop(Nd, tmp1);
    std::vector<LatticeComplex>  W_i_side_loop(Nd, tmp2);
    std::vector<LatticeComplex>  W_f_side_loop(Nd, tmp2);

    // Soon to be deprecated. Keeping V and A parts distinct to take advantage
    // of zero-flop gamma products, when implemented.
    std::vector<std::vector<SpinMatrix>> gL;
    gL.resize(n_i);
    gL[i_V].resize(Nd);
    gL[i_A].resize(Nd);
    for (unsigned int mu = 0; mu < Nd; ++mu)
    {
        gL[i_V][mu] = makeGammaProd(mu);
        gL[i_A][mu] = g5*makeGammaProd(mu);
    }

    // Setup for C-type contractions.
    for (int mu = 0; mu < Nd; ++mu)
    {
        C_i_side_loop[mu] = MAKE_CW_SUBDIAG(q1, q2, gL[i_V][mu]) -
                            MAKE_CW_SUBDIAG(q1, q2, gL[i_A][mu]);
        C_f_side_loop[mu] = MAKE_CW_SUBDIAG(q3, q4, gL[i_V][mu]) -
                            MAKE_CW_SUBDIAG(q3, q4, gL[i_A][mu]);
    }

    // Perform C-type contractions.    
    SUM_MU(expbuf, trace(C_i_side_loop[mu]*C_f_side_loop[mu]))
    MAKE_DIAG(expbuf, corrbuf, result[C_diag], "3pt_HW_C")

    // Recycle sub-expressions for W-type contractions.
    for (unsigned int mu = 0; mu < Nd; ++mu)
    {
        W_i_side_loop[mu] = trace(C_i_side_loop[mu]);
        W_f_side_loop[mu] = trace(C_f_side_loop[mu]);
    }

    // Perform W-type contractions.
    SUM_MU(expbuf, W_i_side_loop[mu]*W_f_side_loop[mu])
    MAKE_DIAG(expbuf, corrbuf, result[W_diag], "3pt_HW_W")

    write(writer, "HW_3pt_NonEye", result[C_diag]);
    write(writer, "HW_3pt_NonEye", result[W_diag]);
}
