#ifndef PELELM_DERIVE_H_
#define PELELM_DERIVE_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

void pelelm_dermprho (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                      const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                      amrex::Real time, const int* bcrec, int level);

void pelelm_derdvrho (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                      const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                      amrex::Real time, const int* bcrec, int level);

void pelelm_derRhoY (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                     const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                     amrex::Real time, const int* bcrec, int level);

void pelelm_dermassfrac (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                         amrex::Real time, const int* bcrec, int level);

void pelelm_dermolefrac (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                         amrex::Real time, const int* bcrec, int level);

void pelelm_dermolweight (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                         amrex::Real time, const int* bcrec, int level);

void pelelm_dercpmix (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                         amrex::Real time, const int* bcrec, int level);

void pelelm_drhomry (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                     const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                     amrex::Real time, const int* bcrec, int level);

void pelelm_dsrhoydot (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                       const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                       amrex::Real time, const int* bcrec, int level);

void pelelm_deravgpres (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                        amrex::Real time, const int* bcrec, int level);

void pelelm_dergrdpx (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                        amrex::Real time, const int* bcrec, int level);

void pelelm_dergrdpy (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                        amrex::Real time, const int* bcrec, int level);

#if (AMREX_SPACEDIM == 3)
void pelelm_dergrdpz (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                        amrex::Real time, const int* bcrec, int level);
#endif

void pelelm_dertransportcoeff (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                        amrex::Real time, const int* bcrec, int level);

void pelelm_dermixanddiss (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                        amrex::Real time, const int* bcrec, int level);

void pelelm_dermixfrac (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                        amrex::Real time, const int* bcrec, int level);

void pelelm_dhrr (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                        amrex::Real time, const int* bcrec, int level);

void pelelm_derconcentration (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                        amrex::Real time, const int* bcrec, int level);

void pelelm_dcma (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& datafab, const amrex::Geometry& geomdata,
                        amrex::Real time, const int* bcrec, int level);

#endif
