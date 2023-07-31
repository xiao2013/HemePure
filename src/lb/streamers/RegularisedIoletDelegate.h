
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_REGULARISEDDELEGATE_H
#define HEMELB_LB_STREAMERS_REGULARISEDDELEGATE_H

#include "util/utilityFunctions.h"
#include "lb/streamers/BaseStreamerDelegate.h"

namespace hemelb
{
	namespace lb
	{
		namespace streamers
		{
			template<typename CollisionImpl>
				class RegularisedIoletDelegate : public BaseStreamerDelegate<CollisionImpl>
			{
				public:
					typedef CollisionImpl CollisionType;
					typedef typename CollisionType::CKernel::LatticeType LatticeType;

					RegularisedIoletDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
						collider(delegatorCollider), iolet(*initParams.boundaryObject)
				{
					// Copied from YangPressureDelegate
					for (int i = 0; i < iolet.GetLocalIoletCount(); ++i)
					{
						SortDirectionsCloseToIOletNormal(iolet.GetLocalIolet(i));
					}
				}

					inline void StreamLink(const LbmParameters* lbmParams,
							geometry::LatticeData* const latticeData,
							const geometry::Site<geometry::LatticeData>& site,
							kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
							const Direction& direction)
					{
						int boundaryId = site.GetIoletId();
						iolets::InOutLet* localIOlet = iolet.GetIolets()[boundaryId];
						Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];

						// Couple with an external system if there is
						if (unstreamed == localIOlet->GetDirectionCloseToNormal(0))
						{
							localIOlet->DoPreStreamCoupling(site.GetIndex(), iolet.GetTimeStep(), site.GetGlobalSiteCoords(),
															hydroVars.density, hydroVars.velocity);
						}
						// hydroVars.density -> current cell
						// how to get current cell/current location
						// Set the density at the "ghost" site to be the density of the iolet.
						distribn_t ghostDensity = hydroVars.density;

						// Calculate the velocity at the ghost site, as the component normal to the iolet.
						util::Vector3D<float> ioletNormal = iolet.GetIolets()[boundaryId]->GetNormal();

						// Note that the division by density compensates for the fact that v_x etc have momentum
						// not velocity.
						distribn_t component = (hydroVars.momentum / hydroVars.density).Dot(ioletNormal);

						// TODO it's ugly that we have to do this.
						// TODO having to give 0 as an argument is also ugly.
						// TODO it's ugly that we have to give hydroVars a nonsense distribution vector
						// that doesn't get used.
						// this is the cell next to the fluid neighbour?
						kernels::HydroVars<typename CollisionType::CKernel> ghostHydrovars(site);

						ghostHydrovars.density = ghostDensity;
						ghostHydrovars.momentum = ioletNormal * component * ghostDensity;
						// hydrovars.GetFE
						collider.kernel.CalculateFeq(ghostHydrovars, 0);
						// Q: find neighbour density
						// Q: find neighbour velocity
						// Q: store last timestep variable
						// Q: reverse direction?
						// Q:
            // double Sxx = 0.0, Syy = 0.0, Szz = 0.0;
            // double Sxy = 0.0, Sxz = 0.0, Syz = 0.0; 
            // TODO: add for loop for local Sxx, Syy, Sxy, Sxz, Syz part
            // for (i = 0; i < Q; i++) {
            //   // fneq may be obtained from fluid inside
            //   fneq = collideField[currentIndex(xpos, y, z) * Q + i] - feq[i];
            //   // bug! missing c_s * c_s in Sxx Syy and Szz component
            //   Sxx += (LatticeType::CX[i] * LatticeType::CX[i] - C_S * C_S) * fneq;
            //   Syy += (LatticeType::CY[i] * LatticeType::CY[i] - C_S * C_S) * fneq;
            //   Szz += (LatticeType::CZ[i] * LatticeType::CZ[i] - C_S * C_S) * fneq;
            //   Sxy += LatticeType::CX[i] * LatticeType::CY[i] * fneq;
            //   Sxz += LatticeType::CX[i] * LatticeType::CZ[i] * fneq;
            //   Syz += LatticeType::CY[i] * LatticeType::CZ[i] * fneq;
            // }
            // for (i = 0; i < Q; i++) {
            //   // computeQPopTensor(Q_tensor, i);
            //   // // double Q_dd_P = doubleDotProduct(Q_tensor, p_hat_bc_1);
            //   double Q_dd_P = ((LatticeType::CX[i] * LatticeType::CX[i] - C_S * C_S) * Sxx) +
            //                   ((LatticeType::CY[i] * LatticeType::CY[i] - C_S * C_S) * Syy) +
            //                   ((LatticeType::CZ[i] * LatticeType::CZ[i] - C_S * C_S) * Szz) +
            //                   (2.0 * LatticeType::CX[i] * LatticeType::CY[i] * Sxy) +
            //                   (2.0 * LatticeType::CX[i] * LatticeType::CZ[i] * Sxz) +
            //                   (2.0 * LatticeType::CY[i] * LatticeType::CZ[i] * Syz);

            //   // collideField[currentIndex(xpos, y, z) * OPTLBM_Q + i] =
            //   //     feq[i] + LatticeType::EQMWEIGHTS[i] / (2. * C_S * C_S * C_S * C_S) * Q_dd_P;
            // }
						// get non-equilibrium part of the population eq(11)
						// std::cout << "hydroVars.GetFEq()[unstreamed]: " << hydroVars.GetFEq()[unstreamed] << std::endl;
						// if(std::abs(hydroVars.velocity[2]) > 0.01){
						// 	std::cout << "hydroVars.GetFNeq()[unstreamed]: " << hydroVars.GetFNeq()[unstreamed] << std::endl;
						// }
						*latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS + unstreamed)
							= ghostHydrovars.GetFEq()[unstreamed] + (1.0 - 1.0 / hydroVars.tau) * hydroVars.GetFNeq()[unstreamed];
					}
					// doesn't do anything for now
					inline void PostStepLink(geometry::LatticeData* const latticeData,
							const geometry::Site<geometry::LatticeData>& site,
							const Direction& direction)
					{
						int boundaryId = site.GetIoletId();
						iolets::InOutLet* localIOlet = iolet.GetIolets()[boundaryId];
						Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];

						// Finalise the coupling with the external system
						if (unstreamed == localIOlet->GetDirectionCloseToNormal(0))
						{
							localIOlet->DoPostStreamCoupling(site.GetIndex(), iolet.GetTimeStep(), site.GetGlobalSiteCoords());
						}
					}

				protected:
					// Copied from YangPressureDelegate
					void SortDirectionsCloseToIOletNormal(iolets::InOutLet* localIOlet)
					{
						const LatticePosition& ioletNormal = localIOlet->GetNormal();
						std::array<Direction, LatticeType::NUMVECTORS> dirs;
        				std::array<LatticeDistance, LatticeType::NUMVECTORS> dist;

	        			for (Direction k = 0; k < LatticeType::NUMVECTORS; ++k)
						{
							const LatticePosition ck = LatticePosition(LatticeType::CXD[k],
												                       LatticeType::CYD[k],
															           LatticeType::CZD[k]);
							const LatticePosition unitCk = ck.GetNormalised();
							dist[k] = LatticePosition(unitCk - ioletNormal).GetMagnitudeSquared();
        	  				dirs[k] = k;
        				}

						// Sort dirs by comparing any two elements of dist.
       	 				std::sort(dirs.begin(), dirs.end(), [&dist](Direction i, Direction j) {return dist[i] < dist[j];});

						// Store the results in the iolet object.
						localIOlet->SetDirectionsCloseToNormal(dirs.begin(), dirs.end());
      				}

					CollisionType& collider;
					iolets::BoundaryValues& iolet;
			};
		}
	}
}

#endif // HEMELB_LB_STREAMERS_REGULARISEDDELEGATE_H
