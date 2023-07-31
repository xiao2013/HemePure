
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_GUOZHENGSHIIOLET_H
#define HEMELB_LB_STREAMERS_GUOZHENGSHIIOLET_H

#include "lb/streamers/StreamerTypeFactory.h"
#include "lb/streamers/GuoZhengShiDelegate.h"
#include "lb/streamers/SimpleBounceBackDelegate.h"
#include "lb/streamers/BouzidiFirdaousLallemandDelegate.h"
#include "lb/streamers/GuoZhengShiIoletDelegate.h"
#include "lb/streamers/JunkYangFactory.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<class CollisionType>
      struct GuoZhengShiIolet
      {
          typedef IoletStreamerTypeFactory<CollisionType, GuoZhengShiIoletDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct GuoZhengShiIoletSBB
      {
          typedef WallIoletStreamerTypeFactory<CollisionType,
              SimpleBounceBackDelegate<CollisionType>, GuoZhengShiIoletDelegate<CollisionType> > Type;
      };

      template<class CollisionType>
      struct GuoZhengShiIoletBFL
      {
          typedef WallIoletStreamerTypeFactory<CollisionType,
              BouzidiFirdaousLallemandDelegate<CollisionType>, GuoZhengShiIoletDelegate<CollisionType> > Type;
      };

      // template<class CollisionType>
      // struct GuoZhengShiIoletGZS
      // {
      //     typedef WallIoletStreamerTypeFactory<CollisionType, GuoZhengShiDelegate<CollisionType>,
      //         GuoZhengShiIoletDelegate<CollisionType> > Type;
      // };
      
      // template<class CollisionType>
      // struct GuoZhengShiIoletGZSE
      // {
      //     typedef WallIoletStreamerTypeFactory<CollisionType, GuoZhengShiElasticWallDelegate<CollisionType>,
      //         GuoZhengShiIoletDelegate<CollisionType> > Type;
      // };
    }
  }
}

#endif /* HEMELB_LB_STREAMERS_GUOZHENGSHIIOLET_H */
