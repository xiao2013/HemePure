
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_GUOZHENGSHIIOLETDELEGATE_H
#define HEMELB_LB_STREAMERS_GUOZHENGSHIIOLETDELEGATE_H

#include "util/utilityFunctions.h"
#include "lb/streamers/BaseStreamerDelegate.h"

namespace hemelb
{
	namespace lb
	{
		namespace streamers
		{
			template<typename CollisionImpl>
				class GuoZhengShiIoletDelegate : public BaseStreamerDelegate<CollisionImpl>
			{
				public:
					typedef CollisionImpl CollisionType;
					typedef typename CollisionType::CKernel::LatticeType LatticeType;

					GuoZhengShiIoletDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
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
            
						kernels::HydroVars<typename CollisionType::CKernel> ghostHydrovars(site);

						ghostHydrovars.density = ghostDensity;
						ghostHydrovars.momentum = ioletNormal * component * ghostDensity;
						// hydrovars.GetFE
						collider.kernel.CalculateFeq(ghostHydrovars, 0);

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

#endif // HEMELB_LB_STREAMERS_GUOZHENGSHIIOLETDELEGATE_H
