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

%
Interplanetary shock project by Markus Battarbee (markus.battarbee@gmail.com)
Based on SilvaShock project by Urs Ganse
Previous development version name was UtuShock
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <vector>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"
#include "../../object_wrapper.h"

#include "Diffusionshock.h"

using namespace std;
using namespace spatial_cell;

namespace projects {
  Diffusionshock::Diffusionshock(): TriAxisSearch() { }
  Diffusionshock::~Diffusionshock() { }
  
  bool Diffusionshock::initialize() {
    return Project::initialize();
  }
  
  void Diffusionshock::addParameters() {
    typedef Readparameters RP;
    // Common (field / etc.) parameters
    RP::add("Diffusionshock.BX0u", "Upstream mag. field value (T)", 1.0e-9);
    RP::add("Diffusionshock.BY0u", "Upstream mag. field value (T)", 2.0e-9);
    RP::add("Diffusionshock.BZ0u", "Upstream mag. field value (T)", 3.0e-9);
    RP::add("Diffusionshock.BX0d", "Downstream mag. field value (T)", 1.0e-9);
    RP::add("Diffusionshock.BY0d", "Downstream mag. field value (T)", 2.0e-9);
    RP::add("Diffusionshock.BZ0d", "Downstream mag. field value (T)", 3.0e-9);
    RP::add("Diffusionshock.Width", "Shock Width (m)", 50000);

    RP::add("Diffusionshock.AMR_L1width", "L1 AMR region width (m)", 0);
    RP::add("Diffusionshock.AMR_L2width", "L2 AMR region width (m)", 0);
    RP::add("Diffusionshock.AMR_L3width", "L3 AMR region width (m)", 0);
    RP::add("Diffusionshock.AMR_L4width", "L4 AMR region width (m)", 0);

    // Per-population parameters
    for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
       const std::string& pop = getObjectWrapper().particleSpecies[i].name;
       RP::add(pop + "_Diffusionshock.VX0u", "Upstream Bulk velocity in x", 0.0);
       RP::add(pop + "_Diffusionshock.VY0u", "Upstream Bulk velocity in y", 0.0);
       RP::add(pop + "_Diffusionshock.VZ0u", "Upstream Bulk velocuty in z", 0.0);
       RP::add(pop + "_Diffusionshock.rhou", "Upstream Number density (m^-3)", 1.0e7);
       RP::add(pop + "_Diffusionshock.Temperatureu", "Upstream Temperature (K)", 2.0e6);

       RP::add(pop + "_Diffusionshock.VX0d", "Downstream Bulk velocity in x", 0.0);
       RP::add(pop + "_Diffusionshock.VY0d", "Downstream Bulk velocity in y", 0.0);
       RP::add(pop + "_Diffusionshock.VZ0d", "Downstream Bulk velocuty in z", 0.0);
       RP::add(pop + "_Diffusionshock.rhod", "Downstream Number density (m^-3)", 1.0e7);
       RP::add(pop + "_Diffusionshock.Temperatured", "Downstream Temperature (K)", 2.0e6);

       RP::add(pop + "_Diffusionshock.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
       RP::add(pop + "_Diffusionshock.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
       RP::add(pop + "_Diffusionshock.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);

       RP::add(pop + "_Diffusionshock.mushape_minv", "Minimum initialised velocity for parabolic mu-dist (m^-6 s^3)", 0.0e6);
       RP::add(pop + "_Diffusionshock.mushape_maxv", "Maximum initialised velocity for parabolic mu-dist (m^-6 s^3)", 0.0e6);
       RP::add(pop + "_Diffusionshock.mushape_A", "Shape coefficient A for parabolic mu-dist (m^-6 s^3)", 2.0e-12);
       RP::add(pop + "_Diffusionshock.mushape_B", "Base phase-space density B for parabolic mu-dist (m^-6 s^3)", 1.0e-12);
    }

  }

  void Diffusionshock::getParameters() {
    Project::getParameters();

    typedef Readparameters RP;
    RP::get("Diffusionshock.BX0u", this->B0u[0]);
    RP::get("Diffusionshock.BY0u", this->B0u[1]);
    RP::get("Diffusionshock.BZ0u", this->B0u[2]);
    RP::get("Diffusionshock.BX0d", this->B0d[0]);
    RP::get("Diffusionshock.BY0d", this->B0d[1]);
    RP::get("Diffusionshock.BZ0d", this->B0d[2]);
    RP::get("Diffusionshock.Width", this->Shockwidth);

    if (abs(this->B0u[0]-this->B0d[0]) > std::numeric_limits<double>::min()) {
       std::cerr<<"Warning! Bx component changes between upstream and downstream. Divergence of B!"<<std::endl;
    }

    RP::get("Diffusionshock.AMR_L1width", this->AMR_L1width);
    RP::get("Diffusionshock.AMR_L2width", this->AMR_L2width);
    RP::get("Diffusionshock.AMR_L3width", this->AMR_L3width);
    RP::get("Diffusionshock.AMR_L4width", this->AMR_L4width);

    // Per-population parameters
    for(uint i=0; i< getObjectWrapper().particleSpecies.size(); i++) {
       const std::string& pop = getObjectWrapper().particleSpecies[i].name;
       DiffusionshockSpeciesParameters sP;

       RP::get(pop + "_Diffusionshock.VX0u", sP.V0u[0]);
       RP::get(pop + "_Diffusionshock.VY0u", sP.V0u[1]);
       RP::get(pop + "_Diffusionshock.VZ0u", sP.V0u[2]);
       RP::get(pop + "_Diffusionshock.rhou", sP.DENSITYu);
       RP::get(pop + "_Diffusionshock.Temperatureu", sP.TEMPERATUREu);

       RP::get(pop + "_Diffusionshock.VX0d", sP.V0d[0]);
       RP::get(pop + "_Diffusionshock.VY0d", sP.V0d[1]);
       RP::get(pop + "_Diffusionshock.VZ0d", sP.V0d[2]);
       RP::get(pop + "_Diffusionshock.rhod", sP.DENSITYd);
       RP::get(pop + "_Diffusionshock.Temperatured", sP.TEMPERATUREd);

       RP::get(pop + "_Diffusionshock.nSpaceSamples", sP.nSpaceSamples);
       RP::get(pop + "_Diffusionshock.nVelocitySamples", sP.nVelocitySamples);
       RP::get(pop + "_Diffusionshock.maxwCutoff", sP.maxwCutoff);

       RP::get(pop + "_Diffusionshock.mushape_minv", sP.mushape_minv);
       RP::get(pop + "_Diffusionshock.mushape_maxv", sP.mushape_maxv);
       RP::get(pop + "_Diffusionshock.mushape_A", sP.mushape_A);
       RP::get(pop + "_Diffusionshock.mushape_B", sP.mushape_B);
       speciesParams.push_back(sP);
    }
  }


  std::vector<std::array<Real, 3>> Diffusionshock::getV0(creal x, creal y, creal z, const uint popID) const {

    Real mass = getObjectWrapper().particleSpecies[popID].mass;
    Real mu0 = physicalconstants::MU_0;
    const DiffusionshockSpeciesParameters& sP = this->speciesParams[popID];

    // Interpolate all values between upstream and downstream
    Real DENSITY = interpolate(sP.DENSITYu,sP.DENSITYd, x);
    if (DENSITY < 1e-20) {
      std::cout<<"density too low! "<<DENSITY<<" x "<<x<<" y "<<y<<" z "<<z<<std::endl;
    }
    
    // Solve tangential components for B and V
    Real VX = interpolate(sP.V0u[0],sP.V0d[0], x);
    Real VY = interpolate(sP.V0u[1],sP.V0d[1], x);
    Real VZ = interpolate(sP.V0u[2],sP.V0d[2], x);

    // Disable compiler warnings: (unused variables but the function is inherited)
    (void)y;
    (void)z;
    
    std::array<Real, 3> V0 {{VX, VY, VZ}};
    std::vector<std::array<Real, 3>> retval;
    retval.push_back(V0);

    return retval;
  }

  Real Diffusionshock::getDistribValue(creal& x, creal& y, creal& z, 
                                       creal& vx, creal& vy, creal& vz, 
                                       creal& dvx, creal& dvy, creal& dvz,
                                       const uint popID) const {

    Real mass = getObjectWrapper().particleSpecies[popID].mass;
    Real KB = physicalconstants::K_B;
    Real mu0 = physicalconstants::MU_0;
    Real adiab = 5./3.;
    const DiffusionshockSpeciesParameters& sP = this->speciesParams[popID];

    // Special case: distribution parabolic in mu
    if (sP.mushape_maxv>std::numeric_limits<double>::min()) {
       Real v = sqrt(vx*vx+vy*vy+vz*vz);
       if (v<sP.mushape_minv) return 0.0;
       if (v>sP.mushape_maxv) return 0.0;
       // find mu
       // For now, assume B=Bx
       Real mu = vx/sqrt(vy*vy+vz*vz);
       return sP.mushape_A*(mu+1.)*(mu+1.) + sP.mushape_B;
    }

    // Interpolate between upstream and downstream
    Real DENSITY = interpolate(sP.DENSITYu,sP.DENSITYd, x);
    if (DENSITY < 1e-20) {
      std::cout<<"density too low! "<<DENSITY<<" x "<<x<<" y "<<y<<" z "<<z<<std::endl;
    }
    Real TEMPERATURE = interpolate(sP.TEMPERATUREu,sP.TEMPERATUREd, x);
    Real hereVX = interpolate(sP.V0u[0],sP.V0d[0], x);
    Real hereVY = interpolate(sP.V0u[1],sP.V0d[1], x);
    Real hereVZ = interpolate(sP.V0u[2],sP.V0d[2], x);
    std::array<Real, 3> pertV0 {{hereVX, hereVY, hereVZ}};

    Real result = 0.0;

    result = DENSITY * std::pow(mass / (2.0 * M_PI * KB * TEMPERATURE), 1.5) *
      exp(- mass * ((vx-pertV0[0])*(vx-pertV0[0]) + (vy-pertV0[1])*(vy-pertV0[1]) + (vz-pertV0[2])*(vz-pertV0[2])) / (2.0 * KB * TEMPERATURE));

    return result;
  }

  Real Diffusionshock::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz, const uint popID) const {

    const DiffusionshockSpeciesParameters& sP = this->speciesParams[popID];
    Real result = 0.0;
    if((sP.nSpaceSamples > 1) && (sP.nVelocitySamples > 1)) {
      creal d_x = dx / (sP.nSpaceSamples-1);
      creal d_y = dy / (sP.nSpaceSamples-1);
      creal d_z = dz / (sP.nSpaceSamples-1);
      creal d_vx = dvx / (sP.nVelocitySamples-1);
      creal d_vy = dvy / (sP.nVelocitySamples-1);
      creal d_vz = dvz / (sP.nVelocitySamples-1);

      Real avg = 0.0;
      
      for (uint i=0; i<sP.nSpaceSamples; ++i)
         for (uint j=0; j<sP.nSpaceSamples; ++j)
            for (uint k=0; k<sP.nSpaceSamples; ++k)      
               for (uint vi=0; vi<sP.nVelocitySamples; ++vi)
                  for (uint vj=0; vj<sP.nVelocitySamples; ++vj)
                     for (uint vk=0; vk<sP.nVelocitySamples; ++vk)
                     {
                        avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz, popID);
                     }
      
      result = avg /
	(sP.nSpaceSamples*sP.nSpaceSamples*sP.nSpaceSamples) / 
	(sP.nVelocitySamples*sP.nVelocitySamples*sP.nVelocitySamples);
    } else {
      result = getDistribValue(x+0.5*dx, y+0.5*dy, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, dvx, dvy, dvz, popID);
    }               

    if(result < sP.maxwCutoff) {
      return 0.0;
    } else {
      return result;
    }
  }
  
  void Diffusionshock::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) { }

  Real Diffusionshock::interpolate(Real upstream, Real downstream, Real x) const {
    Real coord = 0.5 + x/this->Shockwidth; //Now shock will be from 0 to 1
    //x /= 0.5 * this->Shockwidth;
    Real a = 0.0;
    if (coord <= 0.0) a = downstream;
    if (coord >= 1.0) a = upstream;
    if ((coord > 0.0) && (coord < 1.0)) {
      // Ken Perlin Smootherstep
      Real interpolation = ( 6.0 * coord * coord - 15.0 * coord +10. ) * coord * coord * coord;
      a = upstream * interpolation + downstream * (1. - interpolation);
    }
    return a;
  }

  void Diffusionshock::setProjectBField(
     FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
     FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
     FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid
  ) {
      setBackgroundFieldToZero(BgBGrid);
      
      if(!P::isRestart) {
         auto localSize = perBGrid.getLocalSize().data();
      
#pragma omp parallel for collapse(3)
         for (int x = 0; x < localSize[0]; ++x) {
            for (int y = 0; y < localSize[1]; ++y) {
               for (int z = 0; z < localSize[2]; ++z) {
                  const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(x, y, z);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);
                  
                  Real BX = interpolate(this->B0u[0],this->B0d[0], xyz[0]);
                  Real BY = interpolate(this->B0u[1],this->B0d[1], xyz[0]);
                  Real BZ = interpolate(this->B0u[2],this->B0d[2], xyz[0]);
                  cell->at(fsgrids::bfield::PERBX) = BX;
                  cell->at(fsgrids::bfield::PERBY) = BY;
                  cell->at(fsgrids::bfield::PERBZ) = BZ;
               }
            }
         }
      }
   }


   bool Diffusionshock::refineSpatialCells( dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) const {
 
     int myRank;       
     MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

     std::vector<CellID> refinedCells;

     if(myRank == MASTER_RANK) std::cout << "Maximum refinement level is " << mpiGrid.mapping.get_maximum_refinement_level() << std::endl;
      
     // Leave boundary cells and a bit of safety margin
//      const int bw = 2* VLASOV_STENCIL_WIDTH;
//      const int bw2 = 2*(bw + VLASOV_STENCIL_WIDTH);
//      const int bw3 = 2*(bw2 + VLASOV_STENCIL_WIDTH);

     // Calculate regions for refinement
     if (P::amrMaxSpatialRefLevel > 0) {
	// L1 refinement.
	for (uint i = 0; i < P::xcells_ini; ++i) {
	   for (uint j = 0; j < P::ycells_ini; ++j) {
	      for (uint k = 0; k < P::zcells_ini; ++k) {

		 std::array<double,3> xyz;
		 xyz[0] = P::xmin + (i+0.5)*P::dx_ini;
		 xyz[1] = P::ymin + (j+0.5)*P::dy_ini;
		 xyz[2] = P::zmin + (k+0.5)*P::dz_ini;

		 if (abs(xyz[0]) < AMR_L1width)
		    {
		       CellID myCell = mpiGrid.get_existing_cell(xyz);
		       mpiGrid.refine_completely(myCell);
		    }
	      }
	   }
	}
	refinedCells = mpiGrid.stop_refining(true);      
	if(myRank == MASTER_RANK) std::cout << "Finished first level of refinement" << endl;
	mpiGrid.balance_load();
     }

     if (P::amrMaxSpatialRefLevel > 1) {
	// L2 refinement.
	for (uint i = 0; i < 2*P::xcells_ini; ++i) {
	   for (uint j = 0; j < 2*P::ycells_ini; ++j) {
	      for (uint k = 0; k < 2*P::zcells_ini; ++k) {

		 std::array<double,3> xyz;
		 xyz[0] = P::xmin + (i+0.5)*0.5*P::dx_ini;
		 xyz[1] = P::ymin + (j+0.5)*0.5*P::dy_ini;
		 xyz[2] = P::zmin + (k+0.5)*0.5*P::dz_ini;

		 if (abs(xyz[0]) < AMR_L2width)
		    {
		       CellID myCell = mpiGrid.get_existing_cell(xyz);
		       mpiGrid.refine_completely(myCell);
		    }
	      }
	   }
	}
	refinedCells = mpiGrid.stop_refining(true);      
	if(myRank == MASTER_RANK) std::cout << "Finished second level of refinement" << endl;
	mpiGrid.balance_load();
     }

     if (P::amrMaxSpatialRefLevel > 2) {
	// L3 refinement.
	for (uint i = 0; i < 4*P::xcells_ini; ++i) {
	   for (uint j = 0; j < 4*P::ycells_ini; ++j) {
	      for (uint k = 0; k < 4*P::zcells_ini; ++k) {

		 std::array<double,3> xyz;
		 xyz[0] = P::xmin + (i+0.5)*0.25*P::dx_ini;
		 xyz[1] = P::ymin + (j+0.5)*0.25*P::dy_ini;
		 xyz[2] = P::zmin + (k+0.5)*0.25*P::dz_ini;

		 if (abs(xyz[0]) < AMR_L3width)
		    {
		       CellID myCell = mpiGrid.get_existing_cell(xyz);
		       mpiGrid.refine_completely(myCell);
		    }
	      }
	   }
	}
	refinedCells = mpiGrid.stop_refining(true);      
	if(myRank == MASTER_RANK) std::cout << "Finished third level of refinement" << endl;
	mpiGrid.balance_load();
     }

     if (P::amrMaxSpatialRefLevel > 3) {
	// L4 refinement.
	for (uint i = 0; i < 8*P::xcells_ini; ++i) {
	   for (uint j = 0; j < 8*P::ycells_ini; ++j) {
	      for (uint k = 0; k < 8*P::zcells_ini; ++k) {

		 std::array<double,3> xyz;
		 xyz[0] = P::xmin + (i+0.5)*0.125*P::dx_ini;
		 xyz[1] = P::ymin + (j+0.5)*0.125*P::dy_ini;
		 xyz[2] = P::zmin + (k+0.5)*0.125*P::dz_ini;

		 if (abs(xyz[0]) < AMR_L4width)
		    {
		       CellID myCell = mpiGrid.get_existing_cell(xyz);
		       mpiGrid.refine_completely(myCell);
		    }
	      }
	   }
	}
	refinedCells = mpiGrid.stop_refining(true);      
	if(myRank == MASTER_RANK) std::cout << "Finished fourth level of refinement" << endl;
	mpiGrid.balance_load();
     }

     
     return true;
   }

}//namespace projects
