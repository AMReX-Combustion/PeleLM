#ifndef _Efield_Constant_H_
#define _Efield_Constant_H_
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
namespace EFConst {
   extern AMREX_GPU_DEVICE_MANAGED amrex::Real eps0;          //Free space permittivity (C/(V.m))
   extern AMREX_GPU_DEVICE_MANAGED amrex::Real epsr;
   extern AMREX_GPU_DEVICE_MANAGED amrex::Real elemCharge;    //Coulomb per charge
   extern AMREX_GPU_DEVICE_MANAGED amrex::Real Na;            //Avogadro's number
}
static constexpr amrex::Real small_vel = 1.e-10;
#endif
