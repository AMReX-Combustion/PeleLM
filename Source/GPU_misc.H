#ifndef GPU_MISC_H_
#define GPU_MISC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE
inline
void
gpu_flatten(int icell, int i, int j, int k,  
		amrex::Array4<amrex::Real> const& rhoY,
		amrex::Array4<amrex::Real> const& frcExt,
                amrex::Real * tmp_vect, amrex::Real * tmp_src_vect,
                amrex::Real * tmp_vect_energy, amrex::Real * tmp_src_vect_energy) noexcept
{
        amrex::Real * TV_cell  = tmp_vect + icell * (NUM_SPECIES+1);
        amrex::Real * TVS_cell = tmp_src_vect + icell * (NUM_SPECIES);
	for (int n = 0; n < NUM_SPECIES; n++) {
		TV_cell[n]   = rhoY(i,j,k,n) * 1.e-3;
		TVS_cell[n]  = frcExt(i,j,k,n) * 1.e-3;
	}
        TV_cell[NUM_SPECIES]       = rhoY(i,j,k,NUM_SPECIES+1);  
        tmp_vect_energy[icell]     = rhoY(i,j,k,NUM_SPECIES) * 10.0;
        tmp_src_vect_energy[icell] = frcExt(i,j,k,NUM_SPECIES) * 10.0; 
}


AMREX_GPU_HOST_DEVICE
inline
void
gpu_unflatten(int icell, int i, int j, int k,  
		amrex::Array4<amrex::Real> const& rhoY,
                amrex::Real * tmp_vect,
                amrex::Real * tmp_vect_energy) noexcept
{
        amrex::Real * TV_cell  = tmp_vect + icell * (NUM_SPECIES+1);
	for (int n = 0; n < NUM_SPECIES; n++) {
		rhoY(i,j,k,n) = TV_cell[n] * 1.e+3;
	}
        rhoY(i,j,k,NUM_SPECIES+1) = TV_cell[NUM_SPECIES];  
        rhoY(i,j,k,NUM_SPECIES) = tmp_vect_energy[icell] * 1.e-01;
}




#endif
