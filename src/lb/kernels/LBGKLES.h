
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_LBGKLES_H
#define HEMELB_LB_KERNELS_LBGKLES_H

#include <cstdlib>
#include "lb/HFunction.h"
#include "util/utilityFunctions.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
	namespace lb
	{
		namespace kernels
		{
			/**
			 * LBGKLES: This class implements the LBGK single-relaxation time with LES turbulent model kernel.
			 */
			template<class LatticeType>
				class LBGKLES : public BaseKernel<LBGKLES<LatticeType>, LatticeType>
			{
				public:
					LBGKLES(InitParams& initParams)
					{
					}

					inline void DoCalculateDensityMomentumFeq(HydroVars<LBGKLES<LatticeType> >& hydroVars, site_t index)
					{
						LatticeType::CalculateDensityMomentumFEq(hydroVars.f,
								hydroVars.density,
								hydroVars.momentum.x,
								hydroVars.momentum.y,
								hydroVars.momentum.z,
								hydroVars.velocity.x,
								hydroVars.velocity.y,
								hydroVars.velocity.z,
								hydroVars.f_eq.f);

						for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
						{
							hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
						}

						// Compute the local relaxation time that will be used in the next time step
						hydroVars.tau = compute_tau_smagorinsky(hydroVars);
					}

					inline void DoCalculateFeq(HydroVars<LBGKLES>& hydroVars, site_t index)
					{
						LatticeType::CalculateFeq(hydroVars.density,
								hydroVars.momentum.x,
								hydroVars.momentum.y,
								hydroVars.momentum.z,
								hydroVars.f_eq.f);

						for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
						{
							hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
						}
					}

					inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<LBGKLES>& hydroVars)
					{
						double omega = -1.0 / hydroVars.tau;
						for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
						{
							hydroVars.SetFPostCollision(direction, hydroVars.f[direction] + hydroVars.f_neq.f[direction] * omega);
						}
					}
				private:
					/** 
						* compute_tau_smagorinsky update the relaxation time tau based on the LES Smagorinsky model
						* @param  {HydroVars<LBGKLES>} hydro variables from LBM simulation : 
						*/
						double compute_tau_smagorinsky(HydroVars<LBGKLES>& hydroVars) const
						{
							// calculate tau using smagorinsky local correction
							double dx = 1.0;
							double dt = 1.0;
							double C = dx / dt;
							double rho1 = 1.0;
							double localTau;
							double C_smag = 0.1;
							// Compute non-equilibrium values
							for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
							{
								hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
							}
							double Q_12 = 0.0;
							// Calculate diagonal and upper diagonal of the non equilibrium stress tensor
							for (int i = 0; i < 3; ++i) {
								for (int j = 0; j < 3; ++j) {
									double qij = 0.0;
									for (int v = 0; v < LatticeType::NUMVECTORS; ++v) {
										qij += LatticeType::discreteVelocityVectors[i][v] * LatticeType::discreteVelocityVectors[j][v] * hydroVars.f_neq.f[v];
									}
									Q_12 += qij * qij;
								}
							}
							Q_12 = sqrt(Q_12);
							// eq 36 Koda 2015, csmag is smagorinsky constant here is c_smag^2 in the paper is c_smag
							localTau = 1. / 2. *
										(hydroVars.tau + sqrt((hydroVars.tau * rho1 * C)*(hydroVars.tau * rho1 * C) + 
										18.0 * 1.4142135623730950488016887242097 * rho1 * C_smag * C_smag * Q_12) / (rho1 * C));

							return localTau;
						}
						
				
			};

		}
	}
}

#endif /* HEMELB_LB_KERNELS_LBGKLES_H */
