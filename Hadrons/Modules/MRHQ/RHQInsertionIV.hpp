/*
 * RHQInsertionIV.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
 * Author: Alessandro Barone <barone1618@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */
 
/*  END LEGAL */

#ifndef Hadrons_MRHQ_RHQInsertionIV_hpp_
#define Hadrons_MRHQ_RHQInsertionIV_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                            RHQInsertionIV                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MRHQ)

GRID_SERIALIZABLE_ENUM(OpFlag, undef, Chroma, 0, LeftRight, 1);

class RHQInsertionIVPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RHQInsertionIVPar,
                                    std::string,    q,
                                    unsigned int,   index,
                                    Gamma::Algebra, gamma5,
                                    std::string,    gauge,
                                    OpFlag,         flag);
};

template <typename FImpl, typename GImpl>
class TRHQInsertionIV: public Module<RHQInsertionIVPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,);
public:
    // constructor
    TRHQInsertionIV(const std::string name);
    // destructor
    virtual ~TRHQInsertionIV(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(RHQInsertionIV, ARG(TRHQInsertionIV<FIMPL, GIMPL>), MRHQ);

/******************************************************************************
 *                       TRHQInsertionIV implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
TRHQInsertionIV<FImpl, GImpl>::TRHQInsertionIV(const std::string name)
: Module<RHQInsertionIVPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQInsertionIV<FImpl, GImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q, par().gauge};
    
    return in;
}

template <typename FImpl, typename GImpl>
std::vector<std::string> TRHQInsertionIV<FImpl, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionIV<FImpl, GImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename GImpl>
void TRHQInsertionIV<FImpl, GImpl>::execute(void)
{
    // add flag parameter
    LOG(Message) << "Applying Improvement term IV with index " << par().index
                 << " and gamma5=" << par().gamma5 
                 << " to '" << par().q 
                 << "' with flag '" << par().flag << "'"
                 << std::endl;
    
    if (par().flag == OpFlag::Chroma)
    {
        Gamma g5(par().gamma5); 
        Gamma::Algebra gi; // gamma index
        switch(par().index){
            case 0:
                gi = Gamma::Algebra::GammaX;
                break;
            case 1:
                gi = Gamma::Algebra::GammaY;
                break;
            case 2:
                gi = Gamma::Algebra::GammaZ;
                break;
            case 3:
                gi = Gamma::Algebra::GammaT;
                break;
            // not sure if I should put a default as an error, in case the index is not valid;
            // in other modules there is no such check (ex in index for RHQ action)
            default:
                HADRONS_ERROR(Argument, "Index must be in {'0', '1', '2', '3'}."); 
        }

        auto &field = envGet(PropagatorField, par().q);
        const auto &gaugefield  = envGet(GaugeField, par().gauge);
        const auto gauge_x = peekLorentz(gaugefield, 0);
        const auto gauge_y = peekLorentz(gaugefield, 1);
        const auto gauge_z = peekLorentz(gaugefield, 2);

        Gamma gx(Gamma::Algebra::GammaX);
        Gamma gy(Gamma::Algebra::GammaY);
        Gamma gz(Gamma::Algebra::GammaZ);
        PropagatorField insertion = 
            gx*g5*gi * (GImpl::CovShiftForward(gauge_x,0,field) - GImpl::CovShiftBackward(gauge_x,0,field))
          + gy*g5*gi * (GImpl::CovShiftForward(gauge_y,1,field) - GImpl::CovShiftBackward(gauge_y,1,field))
          + gz*g5*gi * (GImpl::CovShiftForward(gauge_z,2,field) - GImpl::CovShiftBackward(gauge_z,2,field));
        
        auto &out = envGet(PropagatorField, getName());
        out = insertion;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MRHQ_RHQInsertionIV_hpp_
