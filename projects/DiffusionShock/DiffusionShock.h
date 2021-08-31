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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

//#pragma once

#ifndef DIFFUSIONSHOCK_H
#define DIFFUSIONSHOCK_H

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {

   struct DiffusionShockSpeciesParameters {
      Real V0u[3];
      Real DENSITYu;
      Real TEMPERATUREu;

      Real V0d[3];
      Real DENSITYd;
      Real TEMPERATUREd;

      Real mushape_minv,mushape_maxv;
      Real mushape_A,mushape_B;

      Real maxwCutoff;
      uint nSpaceSamples;
      uint nVelocitySamples;
   };

   class DiffusionShock: public TriAxisSearch {
      public:
         DiffusionShock();
         virtual ~DiffusionShock();

         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);

         virtual void setProjectBField(
            FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
            FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
            FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
         );
         virtual Real calcPhaseSpaceDensity(
               creal& x, creal& y, creal& z,
               creal& dx, creal& dy, creal& dz,
               creal& vx, creal& vy, creal& vz,
               creal& dvx, creal& dvy, creal& dvz,
               const uint popID
               ) const;



      protected:
         Real getDistribValue(
               creal& x,creal& y, creal& z,
               creal& vx, creal& vy, creal& vz,
               creal& dvx, creal& dvy, creal& dvz,
               const uint popID
               ) const;
         virtual std::vector<std::array<Real, 3>> getV0(creal x, creal y, creal z, const uint popID) const;
         //virtual void calcCellParameters(Real* cellParams,creal& t);
         virtual void calcCellParameters(spatial_cell::SpatialCell* cell, creal& t);
	 bool refineSpatialCells( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const;
         // Interpolate between up- and downstream quantities
         // based on position
         Real interpolate(Real u, Real d, Real x) const;

         // Upstream bulk values
         Real B0u[3];
         // Downstream bulk values
         Real B0d[3];

         // Flow direction definitions
         Real Bucosphi;
         Real Bdcosphi;
         int Byusign;
         int Bydsign;
         int Bzusign;
         int Bzdsign;

         Real Shockwidth;
	 Real AMR_L1width;
	 Real AMR_L2width;
	 Real AMR_L3width;
	 Real AMR_L4width;

         std::vector<DiffusionShockSpeciesParameters> speciesParams;

   } ; //class DiffusionShock
} // namespace projects

#endif
