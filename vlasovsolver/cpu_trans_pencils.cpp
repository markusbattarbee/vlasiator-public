// use DCCRG version Nov 8th 2018 01482cfba8
#include "../grid.h"

using namespace std;
using namespace spatial_cell;







/* Get the one-dimensional neighborhood index for a given direction and neighborhood size.
 * 
 * @param dimension spatial dimension of neighborhood
 * @param stencil neighborhood size in cells
 * @return neighborhood index that can be passed to DCCRG functions
 */
int getNeighborhood(const uint dimension, const uint stencil) {

   int neighborhood = 0;

   if (stencil == 1) {
      switch (dimension) {
      case 0:
         neighborhood = VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID;
         break;
      case 1:
         neighborhood = VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID;
         break;
      case 2:
         neighborhood = VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID;
         break;
      }
   }

   if (stencil > 1) {
      switch (dimension) {
      case 0:
         neighborhood = VLASOV_SOLVER_X_NEIGHBORHOOD_ID;
         break;
      case 1:
         neighborhood = VLASOV_SOLVER_Y_NEIGHBORHOOD_ID;
         break;
      case 2:
         neighborhood = VLASOV_SOLVER_Z_NEIGHBORHOOD_ID;
         break;
      }
   }
   
   return neighborhood;
   
}

void flagSpatialCellsForAmrCommunication(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                         const vector<CellID>& localPropagatedCells) {

   // Only flag/unflag cells if AMR is active
   if (mpiGrid.get_maximum_refinement_level()==0) return;

   // return if there's no cells to flag
   if(localPropagatedCells.size() == 0) {
      std::cerr<<"No cells!"<<std::endl;
      return;
   }

   for (int dimension=0; dimension<3; dimension++) {
      // These neighborhoods now include the AMR addition beyond the regular vlasov stencil
      int neighborhood = getNeighborhood(dimension,VLASOV_STENCIL_WIDTH);

      // Set flags: loop over local cells
#pragma omp parallel for
      for (uint i=0; i<localPropagatedCells.size(); i++) {
         CellID c = localPropagatedCells[i];
         SpatialCell *ccell = mpiGrid[c];
         if (!ccell) continue;

         // Translated cells also need to be included in order to communicate boundary cell VDFs.
         // Attempting to leave these out for the x or y dimensions also resulted in diffs.
         // if (!do_translate_cell(ccell)) continue;

         // Start with false
         ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_X+dimension] = false;

         // In dimension, check iteratively if any neighbors up to VLASOV_STENCIL_WIDTH distance away are on a different process
         const auto* NbrPairs = mpiGrid.get_neighbors_of(c, neighborhood);

         // Create list of unique distances
         std::set< int > distancesplus;
         std::set< int > distancesminus;
         std::set<CellID> foundNeighborsP;
         std::set<CellID> foundNeighborsM;
         /** Using sets of cells as well, we should only get one distance per
             (potentially less refined) cell. This should result in safe behaviour
             as long as the neighborhood of a cell does not contain cells with a
             refinement level more than 1 level apart from the cell itself.
         */
         for (const auto& nbrPair : *NbrPairs) {
            if(nbrPair.second[dimension] > 0) {
               if (foundNeighborsP.find(nbrPair.first) == foundNeighborsP.end()) {
                  distancesplus.insert(nbrPair.second[dimension]);
                  foundNeighborsP.insert(nbrPair.first);
               }
            }
            if(nbrPair.second[dimension] < 0) {
               if (foundNeighborsM.find(nbrPair.first) == foundNeighborsM.end()) {
                  distancesminus.insert(-nbrPair.second[dimension]);
                  foundNeighborsM.insert(nbrPair.first);
               }
            }
         }

         foundNeighborsP.clear();
         foundNeighborsM.clear();

         int iSrc = VLASOV_STENCIL_WIDTH - 1;
         // Iterate through positive distances for VLASOV_STENCIL_WIDTH elements starting from the smallest distance.
         for (auto it = distancesplus.begin(); it != distancesplus.end(); ++it) {
            if (ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_X+dimension] == true) iSrc = -1;
            if (iSrc < 0) break; // found enough elements
            // Check all neighbors at distance *it
            for (const auto& nbrPair : *NbrPairs) {
               SpatialCell *ncell = mpiGrid[nbrPair.first];
               if (!ncell) continue;
               int distanceInRefinedCells = nbrPair.second[dimension];
               if (distanceInRefinedCells == *it) {
                  if (foundNeighborsP.find(nbrPair.first) != foundNeighborsP.end()) continue;
                  foundNeighborsP.insert(nbrPair.first);
                  if (!mpiGrid.is_local(nbrPair.first)) {
                     ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_X+dimension] = true;
                     break;
                  }
               }
            } // end loop over neighbors
            iSrc--;
         } // end loop over positive distances

         iSrc = VLASOV_STENCIL_WIDTH - 1;
         // Iterate through negtive distances for VLASOV_STENCIL_WIDTH elements starting from the smallest distance.
         for (auto it = distancesminus.begin(); it != distancesminus.end(); ++it) {
            if (ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_X+dimension] == true) iSrc = -1;
            if (iSrc < 0) break; // found enough elements
            // Check all neighbors at distance *it
            for (const auto& nbrPair : *NbrPairs) {
               SpatialCell *ncell = mpiGrid[nbrPair.first];
               if (!ncell) continue;
               int distanceInRefinedCells = -nbrPair.second[dimension];
               if (distanceInRefinedCells == *it) {
                  if (foundNeighborsM.find(nbrPair.first) != foundNeighborsM.end()) continue;
                  foundNeighborsM.insert(nbrPair.first);
                  if (!mpiGrid.is_local(nbrPair.first)) {
                     ccell->SpatialCell::parameters[CellParams::AMR_TRANSLATE_COMM_X+dimension] = true;
                     break;
                  }
               }
            } // end loop over neighbors
            iSrc--;
         } // end loop over negative distances
      } // end loop over local propagated cells
   } // end loop over dimensions
   return;
}

/* Get pointers to spatial cells that are considered source cells for a pencil.
 * Source cells are cells that the pencil reads data from to compute polynomial
 * fits that are used for propagation in the vlasov solver. All cells included
 * in the pencil + VLASOV_STENCIL_WIDTH cells on both ends are source cells.
 * Invalid cells are replaced by closest good cells.
 * Boundary cells are included.
 *
 * @param [in] mpiGrid DCCRG grid object
 * @param [in] pencils pencil data struct
 * @param [in] ipencil index of a pencil in the pencils data struct
 * @param [in] dimension spatial dimension
 * @param [out] sourceCells pointer to an array of pointers to SpatialCell objects for the source cells
 */
void computeSpatialSourceCellsForPencil(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                                        setOfPencils& pencils,
                                        const uint iPencil,
                                        const uint dimension,
                                        SpatialCell **sourceCells){

   // L = length of the pencil iPencil
   int L = pencils.lengthOfPencils[iPencil];
   vector<CellID> ids = pencils.getIds(iPencil);

   // These neighborhoods now include the AMR addition beyond the regular vlasov stencil
   int neighborhood = getNeighborhood(dimension,VLASOV_STENCIL_WIDTH);

   // Get pointers for each cell id of the pencil
   for (int i = 0; i < L; ++i) {
      sourceCells[i + VLASOV_STENCIL_WIDTH] = mpiGrid[ids[i]];
   }

   // Insert pointers for neighbors of ids.front() and ids.back()
   const auto* frontNbrPairs = mpiGrid.get_neighbors_of(ids.front(), neighborhood);
   const auto* backNbrPairs  = mpiGrid.get_neighbors_of(ids.back(),  neighborhood);

   // Create list of unique distances in the negative direction from the first cell in pencil
   std::set< int > distances;
   for (const auto& nbrPair : *frontNbrPairs) {
      if(nbrPair.second[dimension] < 0) {
         // gather positive distance values
         distances.insert(-nbrPair.second[dimension]);
      }
   }

   int iSrc = VLASOV_STENCIL_WIDTH - 1;

   // Iterate through distances for VLASOV_STENCIL_WIDTH elements starting from the smallest distance.
   for (auto it = distances.begin(); it != distances.end(); ++it) {
      if (iSrc < 0) break; // found enough elements

      // Collect all neighbors at distance *it to a vector
      std::vector< CellID > neighbors;
      for (const auto& nbrPair : *frontNbrPairs) {
         int distanceInRefinedCells = -nbrPair.second[dimension];
         if(distanceInRefinedCells == *it) neighbors.push_back(nbrPair.first);
      }
      // Get rid of duplicate neighbor cells at single distance
      neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());

      int refLvl = mpiGrid.get_refinement_level(ids.front());
      if (neighbors.size() == 1) {
         if (sourceCells[iSrc+1] == mpiGrid[neighbors.at(0)]) continue; // already found this cell for different distance         
         sourceCells[iSrc--] = mpiGrid[neighbors.at(0)];
      } else if ( pencils.path[iPencil][refLvl] < neighbors.size() ) {
         if (sourceCells[iSrc+1] == mpiGrid[neighbors.at(pencils.path[iPencil][refLvl])]) continue; // already found this cell for different distance (should not happen)
         sourceCells[iSrc--] = mpiGrid[neighbors.at(pencils.path[iPencil][refLvl])];
         // Code for alternate approach to verify that multiple neighbors are in correct ordering (z-y-x)
         // int ix=0,iy=0;
         // switch(dimension) {
         //    case 0:
         //       ix = 1;
         //       iy = 2;
         //       break;
         //    case 1:
         //       ix = 0;
         //       iy = 2;
         //       break;
         //    case 2:
         //       ix = 0;
         //       iy = 1;
         //       break;
         // }
         // bool accept = false;
         // std::array<double, 3> parentCoords = mpiGrid.get_center(ids.front());
         // for (CellID n : neighbors) {
         //    std::array<double, 3> myCoords = mpiGrid.get_center(n);
         //    switch (pencils.path[iPencil][refLvl]) {
         //       case 0:
         //          if (myCoords[ix] < parentCoords[ix] && myCoords[iy] < parentCoords[iy]) accept=true;
         //          break;
         //       case 1:
         //          if (myCoords[ix] > parentCoords[ix] && myCoords[iy] < parentCoords[iy]) accept=true;
         //             break;
         //       case 2:
         //          if (myCoords[ix] < parentCoords[ix] && myCoords[iy] > parentCoords[iy]) accept=true;
         //          break;
         //       case 3:
         //          if (myCoords[ix] > parentCoords[ix] && myCoords[iy] > parentCoords[iy]) accept=true;
         //          break;
         //    }
         //    if (accept) {
         //       sourceCells[iSrc--] = mpiGrid[n];
         //       break;
         //    }
         // }
      } else {
         std::cerr<<"error too few neighbors for path! "<<std::endl; 
      }
   }

   iSrc = L + VLASOV_STENCIL_WIDTH;
   distances.clear();
   // Create list of unique distances in the positive direction from the last cell in pencil
   for (const auto& nbrPair : *backNbrPairs) {
      if(nbrPair.second[dimension] > 0) {
         distances.insert(nbrPair.second[dimension]);
      }
   }

   // Iterate through distances for VLASOV_STENCIL_WIDTH elements starting from the smallest distance.
   // Distances are positive here so smallest distance has smallest value.
   for (auto it = distances.begin(); it != distances.end(); ++it) {
      if (iSrc >= L+2*VLASOV_STENCIL_WIDTH) break; // Found enough cells

      // Collect all neighbors at distance *it to a vector
      std::vector< CellID > neighbors;
      for (const auto& nbrPair : *backNbrPairs) {
         int distanceInRefinedCells = nbrPair.second[dimension];
         if(distanceInRefinedCells == *it) neighbors.push_back(nbrPair.first);
      }
      // Get rid of duplicate neighbor cells at single distance
      neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());

      int refLvl = mpiGrid.get_refinement_level(ids.back());
      if (neighbors.size() == 1) {
         if (sourceCells[iSrc-1] == mpiGrid[neighbors.at(0)]) continue; // already found this cell for different distance
         sourceCells[iSrc++] = mpiGrid[neighbors.at(0)];
      } else if ( pencils.path[iPencil][refLvl] < neighbors.size() ) {
         if (sourceCells[iSrc-1] == mpiGrid[neighbors.at(pencils.path[iPencil][refLvl])]) continue; // already found this cell for different distance (should not happen)
         sourceCells[iSrc++] = mpiGrid[neighbors.at(pencils.path[iPencil][refLvl])];
      } else {
         std::cerr<<"error too few neighbors for path!"<<std::endl;
      }
   }

   /*loop to negative side and replace all invalid cells with the closest good cell*/
   SpatialCell* lastGoodCell = mpiGrid[ids.front()];
   for(int i = VLASOV_STENCIL_WIDTH - 1; i >= 0 ;--i){
      if(sourceCells[i] == NULL || sourceCells[i]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE)
         sourceCells[i] = lastGoodCell;
      else
         lastGoodCell = sourceCells[i];
   }
   /*loop to positive side and replace all invalid cells with the closest good cell*/
   lastGoodCell = mpiGrid[ids.back()];
   for(int i = VLASOV_STENCIL_WIDTH + L; i < (L + 2*VLASOV_STENCIL_WIDTH); ++i){
      if(sourceCells[i] == NULL || sourceCells[i]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE)
         sourceCells[i] = lastGoodCell;
      else
         lastGoodCell = sourceCells[i];
   }
}

/* Select one nearest neighbor of a cell on the + side in a given dimension. If the neighbor
 * has a higher level of refinement, a path variable is needed to make the selection.
 * Returns INVALID_CELLID if the nearest neighbor is not local to this process.
 * 
 * @param grid DCCRG grid object
 * @param id DCCRG cell id
 * @param dimension spatial dimension
 * @param path index of the desired face neighbor
 * @return neighbor DCCRG cell id of the neighbor
 */
CellID selectNeighbor(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> &grid,
                      CellID id, int dimension = 0, uint path = 0) {

   //int neighborhood = getNeighborhood(dimension,1);
   //const auto* nbrPairs = grid.get_neighbors_of(id, neighborhood);
   
   vector < CellID > myNeighbors;
   CellID neighbor = INVALID_CELLID;
   
   // Iterate through neighbor ids in the positive direction of the chosen dimension,
   // select the neighbor indicated by path, if it is local to this process.
   for (const auto& [neighbor, dir] : grid.get_face_neighbors_of(id)) {
      if (dir == ((int)dimension + 1)) {
         myNeighbors.push_back(neighbor);
      }
   }
   
   if( myNeighbors.size() == 0 ) {
      return neighbor;
   }
   
   int neighborIndex = 0;
   if (myNeighbors.size() > 1) {
      neighborIndex = path;
   }
   
   if (grid.is_local(myNeighbors[neighborIndex])) {
      neighbor = myNeighbors[neighborIndex];
   }
   
   return neighbor;
}


/* Recursive function for building one-dimensional pencils to cover local DCCRG cells.
 * Starts from a given seedID and proceeds finding the nearest neighbor in the given dimension
 * and adding it to the pencil until no neighbors are found or an endId is met. When a higher
 * refinement level (ie. multiple nearest neighbors) is met, the pencil splits into four
 * copies to remain at a width of 1 cell. This is done by the function calling itself recursively
 * and passing as inputs the cells added so far. The cell selected by each copy of the function
 * at a split is stored in the path variable, the same path has to be followed if a refinement
 * level is encoutered multiple times.
 *
 * @param [in] grid DCCRG grid object
 * @param [out] pencils Pencil data struct
 * @param [in] seedId DCCRG cell id where we start building the pencil. 
 *             The pencil will continue in the + direction in the given dimension until an end condition is met
 * @param [in] dimension Spatial dimension
 * @param [in] path Integer value that determines which neighbor is added to the pencil when a higher refinement level is met
 * @param [in] endIds Prescribed end conditions for the pencil. If any of these cell ids is about to be added to the pencil,
 *             the builder terminates.
 */
setOfPencils buildPencilsWithNeighbors( const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry> &grid, 
					setOfPencils &pencils, const CellID seedId,
					vector<CellID> ids, const uint dimension, 
					vector<uint> path, const vector<CellID> &endIds) {

   const bool debug = false;
   CellID nextNeighbor;
   CellID id = seedId;
   int startingRefLvl = grid.get_refinement_level(id);
   bool periodic = false;
   // If this is a new pencil (instead of being a result of a pencil being split
   if( ids.size() == 0 )
      ids.push_back(seedId);

   // If the cell where we start is refined, we need to figure out which path
   // to follow in future refined cells. This is a bit hacky but we have to
   // use the order or the children of the parent cell to figure out which
   // corner we are in.

   std::array<double, 3> coordinates = grid.get_center(seedId);
   int startingPathSize = path.size();
   auto it = path.end();
   if( startingRefLvl > startingPathSize ) {

      CellID myId = seedId;
      
      for ( int i = path.size(); i < startingRefLvl; ++i) {

         //CellID parentId = grid.mapping.get_parent(myId);
         CellID parentId = grid.get_parent(myId);
         
         auto myCoords = grid.get_center(myId);
         auto parentCoords = grid.get_center(parentId);

         int ix=0,iy=0;

         switch(dimension) {
         case 0:
            ix = 1;
            iy = 2;
            break;
         case 1:
            ix = 0;
            iy = 2;
            break;
         case 2:
            ix = 0;
            iy = 1;
            break;
         }
         //int ix = (dimension + 1) % 3; // incorrect for DCCRG
         //int iy = (dimension + 2) % 3;

         int step = -1;
         
         if        (myCoords[ix] < parentCoords[ix] && myCoords[iy] < parentCoords[iy]) {
            step = 0;
         } else if (myCoords[ix] > parentCoords[ix] && myCoords[iy] < parentCoords[iy]) {
            step = 1;
         } else if (myCoords[ix] < parentCoords[ix] && myCoords[iy] > parentCoords[iy]) {
            step = 2;
         } else if (myCoords[ix] > parentCoords[ix] && myCoords[iy] > parentCoords[iy]) {
            step = 3;
         }

         it = path.insert(it, step);

         myId = parentId;
      }
   }
   
   while (id != INVALID_CELLID) {

      periodic = false;
      bool neighborExists = false;
      int refLvl = 0;
      
      // Find the refinement level in the neighboring cell. Check all possible neighbors
      // in case some of them are remote.
      for (int tmpPath = 0; tmpPath < 4; ++tmpPath) {
         nextNeighbor = selectNeighbor(grid,id,dimension,tmpPath);
         if(nextNeighbor != INVALID_CELLID) {
            refLvl = max(refLvl,grid.get_refinement_level(nextNeighbor));
            neighborExists = true;
         }
      }
         
      // If there are no neighbors, we can stop.
      if (!neighborExists)
         break;   

      if (refLvl > 0) {
    
         // If we have encountered this refinement level before and stored
         // the path this builder follows, we will just take the same path
         // again.
         if ( static_cast<int>(path.size()) >= refLvl ) {
      
            if(debug) {
               std::cout << "I am cell " << id << ". ";
               std::cout << "I have seen refinement level " << refLvl << " before. Path is ";
               for (auto k = path.begin(); k != path.end(); ++k)
                  std::cout << *k << " ";
               std::cout << std::endl;
            }
	
            nextNeighbor = selectNeighbor(grid,id,dimension,path[refLvl - 1]);      
	    coordinates = grid.get_center(nextNeighbor);

         } else {
	
            if(debug) {
               std::cout << "I am cell " << id << ". ";
               std::cout << "I have NOT seen refinement level " << refLvl << " before. Path is ";
               for (auto k = path.begin(); k != path.end(); ++k)
                  std::cout << *k << ' ';
               std::cout << std::endl;
            }
	
            // New refinement level, create a path through each neighbor cell
            for ( uint i : {0,1,2,3} ) {
	  
               vector < uint > myPath = path;
               myPath.push_back(i);
	  
               nextNeighbor = selectNeighbor(grid,id,dimension,myPath.back());
	  
               if ( i == 3 ) {
	    
                  // This builder continues with neighbor 3
                  path = myPath;
		  coordinates = grid.get_center(nextNeighbor);

               } else {
	    
                  // Spawn new builders for neighbors 0,1,2
                  buildPencilsWithNeighbors(grid,pencils,id,ids,dimension,myPath,endIds);
	    
               }
	  
            }
	
         }

      } else {
         if(debug) {
            std::cout << "I am cell " << id << ". ";
            std::cout << " This pencil has reached refinement level 0." << std::endl;
         }
      }// Closes if (refLvl == 0)

      // If we found a neighbor, add it to the list of ids for this pencil.
      if(nextNeighbor != INVALID_CELLID) {
         if (debug) {
            std::cout << " Next neighbor is " << nextNeighbor << "." << std::endl;
         }

         if ( std::any_of(endIds.begin(), endIds.end(), [nextNeighbor](uint i){return i == nextNeighbor;}) ||
              !do_translate_cell(grid[nextNeighbor])) {
            
            nextNeighbor = INVALID_CELLID;
         } else {
            ids.push_back(nextNeighbor);
         }
      }
      
      id = nextNeighbor;
   } // Closes while loop

   // Get the x,y - coordinates of the pencil (in the direction perpendicular to the pencil)
   double x,y;
   int ix=0,iy=0;

   switch(dimension) {
   case 0:
      ix = 1;
      iy = 2;
      break;
   case 1:
      ix = 0;
      iy = 2;
      break;
   case 2:
      ix = 0;
      iy = 1;
      break;
   }
   //ix = (dimension + 1) % 3; // incorrect for DCCRG
   //iy = (dimension + 2) % 3;
      
   x = coordinates[ix];
   y = coordinates[iy];

   pencils.addPencil(ids,x,y,periodic,path);
   
   // TODO why do we have both return value and the argument modified in place? Could be made consistent.
   return pencils;
  
}

/* Determine which cells in the local DCCRG mesh should be starting points for pencils.
 * If a neighbor cell is non-local, across a periodic boundary, or in non-periodic boundary layer 1
 * then we use this cell as a seed for pencils
 *
 * @param [in] mpiGrid DCCRG grid object
 * @param [in] localPropagatedCells List of local cells that get propagated
 * ie. not boundary or DO_NOT_COMPUTE
 * @param [in] dimension Spatial dimension
 * @param [out] seedIds list of cell ids that will be starting points for pencils
 */
void getSeedIds(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                const vector<CellID> &localPropagatedCells,
                const uint dimension,
                vector<CellID> &seedIds) {

   const bool debug = false;
   int myRank;
   if (debug) MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   
   // These neighborhoods now include the AMR addition beyond the regular vlasov stencil
   int neighborhood = getNeighborhood(dimension,VLASOV_STENCIL_WIDTH);

#pragma omp parallel for
   for (uint i=0; i<localPropagatedCells.size(); i++) {
      CellID celli = localPropagatedCells[i];

      bool addToSeedIds = P::amrTransShortPencils;
      if (addToSeedIds) {
#pragma omp critical
         seedIds.push_back(celli);
         continue;
      }
      auto myIndices = mpiGrid.mapping.get_indices(celli);
      int myRefLevel;

      /* -----------------------------------------
         | A |   | B |   |_|_|_|_|   |   | C |   |
         |   |   |   |   | | | | |   |   |   |   |
         -----------------------------------------
         For optimal pencil generation, we need seedids at A, B, and C.
         A Is triggered in the first if-clause. Pencils starting from B
         will be split (but won't cause A to split), and pencils from
         C will be able to remain again un-split. These checks need to be done
         only if we aren't already at the maximum refinement level.
      */

      // First check negative face neighbors (A)
      // Returns all neighbors as (id, direction-dimension) pair pointers.
      for (const auto& [neighbor, dir] : mpiGrid.get_face_neighbors_of(celli) ) {
         if ( dir == -((int)dimension + 1) ) {
            // Check that the neighbor is not across a periodic boundary by calculating
            // the distance in indices between this cell and its neighbor.
            auto nbrIndices = mpiGrid.mapping.get_indices(neighbor);

            // If a neighbor is non-local, across a periodic boundary, or in non-periodic boundary layer 1
            // then we use this cell as a seed for pencils
            if (abs ( (int64_t)(myIndices[dimension] - nbrIndices[dimension]) ) > pow(2,mpiGrid.get_maximum_refinement_level()) ||
               !mpiGrid.is_local(neighbor) ||
               !do_translate_cell(mpiGrid[neighbor]) ) 
            {
               addToSeedIds = true;
               break;
            }
               }
      } // finish check A
      if ( addToSeedIds ) {
#pragma omp critical
         seedIds.push_back(celli);
         continue;
      }
      myRefLevel = mpiGrid.get_refinement_level(celli);
      if (mpiGrid.get_maximum_refinement_level() == myRefLevel) continue;

      // Gather neighbours in neighbourhood stencil
      const auto* nbrPairs  = mpiGrid.get_neighbors_of(celli, neighborhood);
      // Create list of unique neighbour distances in both directions
      std::set< int > distancesplus;
      std::set< int > distancesminus;
      for (const auto& nbrPair : *nbrPairs) {
         if(nbrPair.second[dimension] > 0) {
            distancesplus.insert(nbrPair.second[dimension]);
         }
         if(nbrPair.second[dimension] < 0) {
            // gather positive distance values
            distancesminus.insert(-nbrPair.second[dimension]);
         }
      }
      /* Proceed with B, checking if the next positive neighbour has the same refinement level as ccell, but the
         second neighbour a higher one. Iterate through positive distances for VLASOV_STENCIL_WIDTH elements
         starting from the smallest distance. */
      int iSrc = VLASOV_STENCIL_WIDTH-1;
      for (auto it = distancesplus.begin(); it != distancesplus.end(); ++it) {
         if (iSrc < 0) break; // found enough elements
         for (const auto& nbrPair : *nbrPairs) {
            int distanceInRefinedCells = nbrPair.second[dimension];
            if(distanceInRefinedCells == *it) {
               // Break search if we are not at the final entry, and have different refinement level
               if (iSrc!=0 && mpiGrid.get_refinement_level(nbrPair.first)!=myRefLevel) {
                  iSrc = -1;
                  break;
               }
               // Flag as seed id if VLASOV_STENCIL_WIDTH positive neighbour is at higher refinement level
               if (iSrc==0 && mpiGrid.get_refinement_level(nbrPair.first)>myRefLevel) {
                  addToSeedIds = true;
                  break;
               }
            }
         }
         iSrc--;
      } // Finish B check

      if ( addToSeedIds ) {
         #pragma omp critical
         seedIds.push_back(celli);
         continue;
      }
      /* Proceed with C, checking if the next two negative neighbours have the same refinement level as ccell, but the
         third neighbour a higher one. Iterate through negative distances for VLASOV_STENCIL_WIDTH+1 elements
         starting from the smallest distance. */
      iSrc = VLASOV_STENCIL_WIDTH;
      for (auto it = distancesminus.begin(); it != distancesminus.end(); ++it) {
         if (iSrc < 0) break; // found enough elements
         for (const auto& nbrPair : *nbrPairs) {
            int distanceInRefinedCells = -nbrPair.second[dimension];
            if(distanceInRefinedCells == *it) {
               // Break search if we are not at the final entry, and have different refinement level
               if (iSrc!=0 && mpiGrid.get_refinement_level(nbrPair.first)!=myRefLevel) {
                  iSrc = -1;
                  break;
               }
               // Flag as seed id if VLASOV_STENCIL_WIDTH+1 positive neighbour is at higher refinement level
               if (iSrc==0 && mpiGrid.get_refinement_level(nbrPair.first)>myRefLevel) {
                  addToSeedIds = true;
                  break;
               }
            }
         }
         iSrc--;
      } // Finish C check

      if ( addToSeedIds ) {
#pragma omp critical
         seedIds.push_back(celli);
      }
   }

   if(debug) {
      cout << "Rank " << myRank << ", Seed ids are: ";
      for (const auto seedId : seedIds) {
         cout << seedId << " ";
      }
      cout << endl;
   }   
}

/* Check whether the ghost cells around the pencil contain higher refinement than the pencil does.
 * If they do, the pencil must be split to match the finest refined ghost cell.
 *
 * @param mpiGrid DCCRG grid object
 * @param pencils Pencil data struct
 * @param dimension Spatial dimension
 */
void check_ghost_cells(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                       setOfPencils& pencils,
                       uint dimension) {

   const bool debug = false;
   int neighborhoodId = getNeighborhood(dimension,VLASOV_STENCIL_WIDTH);

   int myRank;
   if(debug) {
      MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   }
   
   std::vector<CellID> idsToSplit;

// Thread this loop here
#pragma omp parallel for   
   for (uint pencili = 0; pencili < pencils.N; ++pencili) {

      if(pencils.periodic[pencili]) continue;
         
      auto ids = pencils.getIds(pencili);

      // It is possible that the pencil has already been refined by the pencil building algorithm
      // and is on a higher refinement level than the refinement level of any of the cells it contains
      // due to e.g. process boundaries.
      int maxPencilRefLvl = pencils.path[pencili].size();
      int maxNbrRefLvl = 0;

      const auto* frontNeighbors = mpiGrid.get_neighbors_of(ids.front(),neighborhoodId);
      const auto* backNeighbors  = mpiGrid.get_neighbors_of(ids.back() ,neighborhoodId);


      // Create list of unique distances in the negative direction from the first cell in pencil
      std::set< int > distances;
      for (const auto& nbrPair : *frontNeighbors) {
         if(nbrPair.second[dimension] < 0) {
            // gather positive values
            distances.insert(-nbrPair.second[dimension]);
         }
      }
      int foundcells = 0;
      CellID lastcell = INVALID_CELLID;
      // Iterate through distances for VLASOV_STENCIL_WIDTH elements starting from the smallest distance.
      for (auto it = distances.begin(); it != distances.end(); ++it) {
         for (const auto& nbrPair : *frontNeighbors) {
            if (nbrPair.first==lastcell) continue;
            int distanceInRefinedCells = -nbrPair.second[dimension];
            if(distanceInRefinedCells == *it) {
               maxNbrRefLvl = max(maxNbrRefLvl,mpiGrid.get_refinement_level(nbrPair.first));
               lastcell = nbrPair.first;
               foundcells++;
               continue;
            }
         }
         if (foundcells >= VLASOV_STENCIL_WIDTH) break; // checked enough distances
      }

      // Create list of unique distances in the positive direction from the last cell in pencil
      distances.clear();
      for (const auto& nbrPair : *backNeighbors) {
         if(nbrPair.second[dimension] > 0) {
            distances.insert(nbrPair.second[dimension]);
         }
      }
      foundcells = 0;
      lastcell = INVALID_CELLID;
      for (auto it = distances.begin(); it != distances.end(); ++it) {
         for (const auto& nbrPair : *backNeighbors) {
            if (nbrPair.first==lastcell) continue;
            int distanceInRefinedCells = nbrPair.second[dimension];
            if(distanceInRefinedCells == *it) {
               maxNbrRefLvl = max(maxNbrRefLvl,mpiGrid.get_refinement_level(nbrPair.first));
               lastcell = nbrPair.first;
               foundcells++;
               continue;
            }
         }
         if (foundcells >= VLASOV_STENCIL_WIDTH) break; // checked enough distances
      }

      // Old version which can check needlessly far
      // for (const auto nbrPair: *frontNeighbors) {
      //    maxNbrRefLvl = max(maxNbrRefLvl,mpiGrid.get_refinement_level(nbrPair.first));
      // }
      // for (const auto nbrPair: *backNeighbors) {
      //    maxNbrRefLvl = max(maxNbrRefLvl,mpiGrid.get_refinement_level(nbrPair.first));
      // }


      if (maxNbrRefLvl > maxPencilRefLvl) {
         if(debug) {
            std::cout << "I am rank " << myRank << ". ";
            std::cout << "Found refinement level " << maxNbrRefLvl << " in one of the ghost cells of pencil " << pencili << ". ";
            std::cout << "Highest refinement level in this pencil is " << maxPencilRefLvl;
            std::cout << ". Splitting pencil " << pencili << endl;
         }
         // Let's avoid modifying pencils while we are looping over it. Write down the indices of pencils
         // that need to be split and split them later.
#pragma omp critical
         {
            idsToSplit.push_back(pencili);
         }
      }
   }

// No threading here, probably more efficient to thread inside the splitting
   for (auto pencili: idsToSplit) {

      Realv dx = 0.0;
      Realv dy = 0.0;
      // TODO: Double-check that this gives you the right dimensions!
      auto ids = pencils.getIds(pencili);
      switch(dimension) {
      case 0:
         dx = mpiGrid[ids[0]]->SpatialCell::parameters[CellParams::DY];
         dy = mpiGrid[ids[0]]->SpatialCell::parameters[CellParams::DZ];
         break;
      case 1:
         dx = mpiGrid[ids[0]]->SpatialCell::parameters[CellParams::DX];
         dy = mpiGrid[ids[0]]->SpatialCell::parameters[CellParams::DZ];
         break;
      case 2:
         dx = mpiGrid[ids[0]]->SpatialCell::parameters[CellParams::DX];
         dy = mpiGrid[ids[0]]->SpatialCell::parameters[CellParams::DY];
         break;
      }

// WARNING threading inside this function
      pencils.split(pencili,dx,dy);
         
   }
}

/* Checks that each local spatial cell appears in pencils at least 1 time.
 *
 * @param mpiGrid DCCRG grid object
 * @param cells Local spatial cells
 * @param pencils Pencil data struct
 */
bool checkPencils(
   const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const std::vector<CellID> &cells,
   const setOfPencils& pencils
) {

   bool correct = true;

   for (auto id : cells) {

      if (mpiGrid[id]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY )  {
      
         int myCount = std::count(pencils.ids.begin(), pencils.ids.end(), id);
         
         if( myCount == 0) {
            
            std::cerr << "ERROR: Cell ID " << id << " Appears in pencils " << myCount << " times!"<< std::endl;            
            correct = false;
         }

      }
      
   }

   for (uint ipencil = 0; ipencil < pencils.N; ++ipencil) {
      cint nPencilsThroughThisCell = pow(pow(2,pencils.path[ipencil].size()),2);
      auto ids = pencils.getIds(ipencil);
      
      for (auto id : ids) {

         cint myCount = std::count(pencils.ids.begin(), pencils.ids.end(), id);

         if (myCount > nPencilsThroughThisCell) {

            std::cerr << "ERROR: Cell ID " << id << " Appears in pencils " << myCount << " times!"<< std::endl;
            std::cerr << "       It should not appear more than " << nPencilsThroughThisCell << " times." << std::endl;
            correct = false;

         }

      }

   }

   return correct;
   
}

/* Debugging function, prints the list of cells in each pencil
 *
 * @param pencils Pencil data struct
 * @param dimension Spatial dimension
 * @param myRank MPI rank
 */
void printPencilsFunc(const setOfPencils& pencils, const uint dimension, const int myRank) {
   
   // Print out ids of pencils (if needed for debugging)
   uint ibeg = 0;
   uint iend = 0;
   std::cout << "I am rank " << myRank << ", I have " << pencils.N << " pencils along dimension " << dimension << ":\n";
   MPI_Barrier(MPI_COMM_WORLD);
   if(myRank == MASTER_RANK) {
      std::cout << "t, N, mpirank, dimension (x, y): indices {path} " << std::endl;
      std::cout << "-----------------------------------------------------------------" << std::endl;
   }
   MPI_Barrier(MPI_COMM_WORLD);
   for (uint i = 0; i < pencils.N; i++) {
      iend += pencils.lengthOfPencils[i];
      std::cout << P::t << ", ";
      std::cout << i << ", ";
      std::cout << myRank << ", ";
      std::cout << dimension << ", ";
      std::cout << "(" << pencils.x[i] << ", " << pencils.y[i] << "): ";
      for (auto j = pencils.ids.begin() + ibeg; j != pencils.ids.begin() + iend; ++j) {
         std::cout << *j << " ";
      }
      ibeg  = iend;
      
      std::cout << "{";         
      for (auto step : pencils.path[i]) {
         std::cout << step << ", ";
      }
      std::cout << "}";
      
      std::cout << std::endl;
   }

   MPI_Barrier(MPI_COMM_WORLD);
}

/* Wrapper function for calling seed ID selection and pencil generation, per dimension.
 * Includes threading and gathering of pencils into thread-containers.
 *
 * @param [in] mpiGrid DCCRG grid object
 * @param [in] dimension Spatial dimension
 */
void prepareSeedIdsAndPencils(const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                              const uint dimension) {

   const bool printPencils = false;
   int myRank;
   if(printPencils) MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

   const vector<CellID>& localCells = getLocalCells();
   vector<CellID> localPropagatedCells;
   // Figure out which spatial cells are translated,
   // result independent of particle species.
   for (size_t c=0; c<localCells.size(); ++c) {
      if (do_translate_cell(mpiGrid[localCells[c]])) {
         localPropagatedCells.push_back(localCells[c]);
      }
   }

   phiprof::Timer getSeedsTimer {"getSeedIds"};
   vector<CellID> seedIds;
   getSeedIds(mpiGrid, localPropagatedCells, dimension, seedIds);
   getSeedsTimer.stop();

   phiprof::Timer buildPencilsTimer {"buildPencils"};
   // Output vectors for ready pencils
   //setOfPencils pencils;

   // Clear previous set
   DimensionPencils[dimension].removeAllPencils();
   
#pragma omp parallel
   {
      // Empty vectors for internal use of buildPencilsWithNeighbors. Could be default values but
      // default vectors are complicated. Should overload buildPencilsWithNeighbors like suggested here
      // https://stackoverflow.com/questions/3147274/c-default-argument-for-vectorint
      vector<CellID> ids;
      vector<uint> path;
      // thread-internal pencil set to be accumulated at the end
      setOfPencils thread_pencils;
      // iterators used in the accumulation
      std::vector<CellID>::iterator ibeg, iend;

#pragma omp for schedule(guided)
      for (uint i=0; i<seedIds.size(); i++) {
         cuint seedId = seedIds[i];
         // Construct pencils from the seedIds into a set of pencils.
         thread_pencils = buildPencilsWithNeighbors(mpiGrid, thread_pencils, seedId, ids, dimension, path, seedIds);
      }

      // accumulate thread results in global set of pencils
#pragma omp critical
      {
         for (uint i=0; i<thread_pencils.N; i++) {
            // Use vector range constructor
            ibeg = thread_pencils.ids.begin() + thread_pencils.idsStart[i];
            iend = ibeg + thread_pencils.lengthOfPencils[i];
            std::vector<CellID> pencilIds(ibeg, iend);
            DimensionPencils[dimension].addPencil(pencilIds,thread_pencils.x[i],thread_pencils.y[i],thread_pencils.periodic[i],thread_pencils.path[i]);
         }
      }
   }

   phiprof::Timer checkGhostsTimer {"check_ghost_cells"};
   // Check refinement of two ghost cells on each end of each pencil
   check_ghost_cells(mpiGrid,DimensionPencils[dimension],dimension);
   checkGhostsTimer.stop();

   // ****************************************************************************

   if(printPencils) {
      printPencilsFunc(DimensionPencils[dimension],dimension,myRank);
   }
   buildPencilsTimer.stop();
}

