/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 * Copyright 2017-2019 University of Helsinki
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
/*
Background magnetic field class of Vlasiator.
*/

#ifndef VECTORDIPOLE_HPP
#define VECTORDIPOLE_HPP
#include "fieldfunction.hpp"



class VectorDipole {
private:
   bool initialized = false;
   double q[3];      // Dipole moment; set to (0,0,moment) for z-aligned
   double center[3]; // Coordinates where the dipole sits; set to (0,0,0)
   double xlimit[2]; // X-coordinate extents of full and zero dipole
   double rlimit[2]; // Radial extents for scaling in IMF field
   double IMF[3];    // IMF value to scale to
   /**
      IMF field is added either based on X-coordinate:
        only in upstream
        scales in as dipole is scaled out, with zero IMF at X<xlimit[0] and full IMF at X>xlimit[1]
      or based on radial coordinate:
        everywhere outside inner magnetosphere
        scales in with zero IMT at R<rlimit[0] and full IMF at R>rlimit[1]

        Note: dipole itself always only follows xlimits, not rlimits
   **/
   
public:
   
   VectorDipole(){};
   void initialize(const double moment,const double center_x, const double center_y, const double center_z, const double tilt_angle_phi, const double tilt_angle_theta, const double xlimit_f, const double xlimit_z, const double IMF_rz, const double IMF_rf, const double IMF_Bx, const double IMF_By, const double IMF_Bz);
   double operator()(double x, double y, double z, coordinate component, unsigned int derivative=0, coordinate dcomponent=X) const;
};

#endif

