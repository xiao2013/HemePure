
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_REGULARISEDIOLET_H
#define HEMELB_LB_STREAMERS_REGULARISEDIOLET_H

#include "lb/streamers/StreamerTypeFactory.h"
#include "lb/streamers/RegularisedIoletDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"
#include "lb/streamers/GuoZhengShiDelegate.h"
#include "lb/streamers/JunkYangFactory.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<class CollisionType>
      struct RegularisedIolet
      {
          typedef IoletStreamerTypeFactory<CollisionType, RegularisedIoletDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct RegularisedIoletSBB
      {
          typedef WallIoletStreamerTypeFactory<CollisionType,
              SimpleBounceBackDelegate<CollisionType>, RegularisedIoletDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct RegularisedIoletBFL
      {
          typedef WallIoletStreamerTypeFactory<CollisionType,
              BouzidiFirdaousLallemandDelegate<CollisionType>, RegularisedIoletDelegate<CollisionType> > Type;
      };

      // template<class CollisionType>
      // struct RegularisedIoletGZS
      // {
      //     typedef WallIoletStreamerTypeFactory<CollisionType, GuoZhengShiDelegate<CollisionType>,
      //         RegularisedIoletDelegate<CollisionType> > Type;
      // };
      
      // template<class CollisionType>
      // struct RegularisedIoletGZSE
      // {
      //     typedef WallIoletStreamerTypeFactory<CollisionType, GuoZhengShiElasticWallDelegate<CollisionType>,
      //         RegularisedIoletDelegate<CollisionType> > Type;
      // };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_REGULARISEDIOLET_H */
