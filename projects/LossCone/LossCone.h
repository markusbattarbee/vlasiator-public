/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef LOSSCONE_H
#define LOSSCONE_H

#include <stdlib.h>

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {

   struct LossConeSpeciesParameters {
      Real DENSITY;
      Real TEMPERATUREX;
      Real TEMPERATUREY;
      Real TEMPERATUREZ;
      Real densityPertRelAmp;
      Real velocityPertAbsAmp;
      Real maxwCutoff;
      uint nSpaceSamples;
      uint nVelocitySamples;

      // Parameters for instead implementing monotonic population
      Real monotonic_maxv;
      Real monotonic_amp;
      Real monotonic_base;
      bool monotonic_subtract;
   };

   class LossCone: public TriAxisSearch {
   public:
      LossCone();
      virtual ~LossCone();
      
      virtual bool initialize(void);
      static void addParameters(void);
      virtual void getParameters(void);
      virtual void setProjectBField(
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, 2>& perBGrid,
         FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, 2>& BgBGrid,
         FsGrid< fsgrids::technical, 2>& technicalGrid
      );
      virtual std::vector<std::array<Real, 3> > getV0(
         creal x,
         creal y,
         creal z,
         const uint popID
      ) const;
      virtual void hook(
         cuint& stage,
         const dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
         FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
         FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid
      ) const;
   protected:
      Real getDistribValue(creal& vx, creal& vy, creal& vz, const uint popID) const;
      virtual void calcCellParameters(spatial_cell::SpatialCell* cell,creal& t);
      virtual Real calcPhaseSpaceDensity(
         creal& x, creal& y, creal& z,
         creal& dx, creal& dy, creal& dz,
         creal& vx, creal& vy, creal& vz,
         creal& dvx, creal& dvy, creal& dvz,const uint popID
      ) const;

      bool SaveDispersionData;
      Real BX0;
      Real BY0;
      Real BZ0;
      Real magXPertAbsAmp;
      Real magYPertAbsAmp;
      Real magZPertAbsAmp;
      uint seed;
      std::vector<LossConeSpeciesParameters> speciesParams;

      static Real rndRho, rndVel[3];
      #pragma omp threadprivate(rndRho,rndVel)
   } ; // class LossCone
} // namespace projects
#endif
