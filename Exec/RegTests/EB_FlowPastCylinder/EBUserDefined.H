#ifndef _EBUSERDEFINED_H_
#define _EBUSERDEFINED_H_

using namespace amrex;

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
void
EBUserDefined(const Geometry& geom,
              const int required_coarsening_level,
              const int max_coarsening_level)
{
    Array<Real, AMREX_SPACEDIM> CylCenter1;
    Array<Real, AMREX_SPACEDIM> CylCenter2;
    Real CylRadius1 = 0.0;
    Real CylRadius2 = 0.0;
    // ParmParse your geometry parameters
    ParmParse pp("TwoCyl");
    pp.query("Cyl1_center", CylCenter1); 
    pp.query("Cyl2_center", CylCenter2); 
    pp.query("CylRadius1", CylRadius1); 
    pp.query("CylRadius2", CylRadius2); 

    // Build geometry pieces using EB2::* methods
    EB2::SphereIF Sph1(CylRadius1,CylCenter1,false);
    EB2::SphereIF Sph2(CylRadius2,CylCenter2,false);
    auto twoSpheres = EB2::makeUnion(Sph1, Sph2);

    // Build your geometry shop using EB2::makeShop
    auto gshop = EB2::makeShop(twoSpheres);

    // Build geom using EB2::Build
    EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
}
#endif
#endif
