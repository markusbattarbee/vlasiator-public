/*
 * This file is part of Vlasiator.
 * Copyright 2010-2024 Finnish Meteorological Institute and University of Helsinki
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

#include "spatial_batch_gpu.hpp"
#include "arch/gpu_base.hpp"
#include "object_wrapper.h"
#include "velocity_mesh_parameters.h"

using namespace std;

   // void adjust_velocity_blocks(
   //    dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   //    const std::vector<SpatialCell*>& spatial_neighbors,
   //    const uint popID,
   //    bool doDeleteEmptyBlocks=true);
   // vmesh::LocalID adjust_velocity_blocks_caller(
   //    dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   //    const uint popID);

   // void update_velocity_block_content_lists(
   //    dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   //    const uint popID=0);

/** GPU kernel for identifying which blocks have relevant content */
__global__ void batch_update_velocity_block_content_lists_kernel (
   vmesh::VelocityMesh **vmeshes,
   vmesh::VelocityBlockContainer **blockContainers,
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** allMaps,
   Real* velocity_block_min_values
   ) {
   // launch griddim3 grid(nCells,launchBlocks,1);
   const int cellIndex = blockIdx.x;
   const int blocki = blockIdx.y;
   const uint ti = threadIdx.x;

   vmesh::VelocityMesh* vmesh = vmeshes[cellIndex];
   vmesh::VelocityBlockContainer* blockContainer = blockContainers[cellIndex];
   Real velocity_block_min_value = velocity_block_min_values[cellIndex];
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwcl_map = allMaps[2*cellIndex];
   Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>* vbwncl_map = allMaps[2*cellIndex +1];

   // Each GPU block / workunit can theoretically manage several Vlasiator velocity blocks at once.
   const uint vlasiBlocksPerWorkUnit = 1;
   const uint workUnitIndex = 0; // [0,vlasiBlocksPerWorkUnit)
   // const uint vlasiBlocksPerWorkUnit = WARPSPERBLOCK * GPUTHREADS / WID3;
   // const uint workUnitIndex = ti / WID3; // [0,vlasiBlocksPerWorkUnit)
   const uint b_tid = ti % WID3; // [0,WID3)
   const uint blockLID = blocki * vlasiBlocksPerWorkUnit + workUnitIndex; // [0,nBlocksToChange)

   __shared__ int has_content[WARPSPERBLOCK * GPUTHREADS];
   const uint nBlocks = vmesh->size();
   //for (uint blockLID=blocki; blockLID<nBlocks; blockLID += gpuBlocks) {
   if (blockLID < nBlocks) {
      const vmesh::GlobalID blockGID = vmesh->getGlobalID(blockLID);
      #ifdef DEBUG_SPATIAL_CELL
      if (blockGID == vmesh->invalidGlobalID()) {
         if (b_tid==0) printf("Invalid GID encountered in update_velocity_block_content_lists_kernel!\n");
         return;
      }
      if (blockLID == vmesh->invalidLocalID()) {
         if (b_tid==0) printf("Invalid LID encountered in update_velocity_block_content_lists_kernel!\n");
         return;
      }
      #endif
      // Check each velocity cell if it is above the threshold
      const Realf* avgs = blockContainer->getData(blockLID);
      has_content[ti] = avgs[b_tid] >= velocity_block_min_value ? 1 : 0;
      __syncthreads(); // THIS SYNC IS CRUCIAL!
      // Implemented just a simple non-optimized thread OR
      // GPUTODO reductions via warp voting
      for (unsigned int s=WID3/2; s>0; s>>=1) {
         if (b_tid < s) {
            has_content[ti] = has_content[ti] || has_content[ti + s];
         }
         __syncthreads();
      }
      // Insert into map only from threads 0...WARPSIZE
      if (b_tid < GPUTHREADS) {
         if (has_content[0]) {
            vbwcl_map->warpInsert(blockGID,blockLID,b_tid);
         } else {
            vbwncl_map->warpInsert(blockGID,blockLID,b_tid);
         }
      }
      __syncthreads();
   }
}

/** Bulk call over listed cells of spatial grid
    Prepares the content / no-content velocity block lists
    for all requested cells, for the requested popID
**/
   void update_velocity_block_content_lists(
      dccrg::Dccrg<spatial_cell::SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
      const vector<CellID>& cells,
      const uint popID) {

      const uint nCells = cells.size();
#ifdef DEBUG_SPATIAL_CELL
      if (popID >= populations.size()) {
         std::cerr << "ERROR, popID " << popID << " exceeds populations.size() " << populations.size() << " in ";
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         exit(1);
      }
#endif

      // Allocate buffers for GPU operations
      Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>** host_allMaps, **dev_allMaps;
      //size_t* host_sizePowers, *dev_sizePowers;
      Real* host_minValues, *dev_minValues;
      vmesh::VelocityMesh** host_vmeshes, **dev_vmeshes;
      vmesh::VelocityBlockContainer** host_VBCs, **dev_VBCs;
      //vmesh::LocalID* host_contentSizes, *dev_contentSizes;
      CHK_ERR( gpuMalloc((void**)&dev_allMaps, 2*nCells*sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*)) );
      CHK_ERR( gpuMallocHost((void**)&host_allMaps, 2*nCells*sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*)) );
      // CHK_ERR( gpuMalloc((void**)&dev_sizePowers, 2*nCells*sizeof(size_t)) );
      // CHK_ERR( gpuMallocHost((void**)&host_sizePowers, 2*nCells*sizeof(size_t)) );
      CHK_ERR( gpuMalloc((void**)&dev_minValues,nCells*sizeof(Real)) );
      CHK_ERR( gpuMallocHost((void**)&host_minValues, nCells*sizeof(Real)) );
      CHK_ERR( gpuMalloc((void**)&dev_vmeshes,nCells*sizeof(vmesh::VelocityMesh*)) );
      CHK_ERR( gpuMallocHost((void**)&host_vmeshes,nCells*sizeof(vmesh::VelocityMesh*)) );
      CHK_ERR( gpuMalloc((void**)&dev_VBCs,nCells*sizeof(vmesh::VelocityBlockContainer*)) );
      CHK_ERR( gpuMallocHost((void**)&host_VBCs,nCells*sizeof(vmesh::VelocityBlockContainer*)) );
      // CHK_ERR( gpuMalloc((void**)&dev_contentSizes,nCells*sizeof(vmesh::LocalID)) );
      // CHK_ERR( gpuMallocHost((void**)&host_contentSizes,nCells*sizeof(vmesh::LocalID)) );

      phiprof::Timer sparsityTimer {"update Sparsity values and clear maps"};
      size_t largestSizePower = 0;
      size_t largestVelMesh = 0;
#pragma omp parallel
      {
         const gpuStream_t stream = gpu_getStream();
         size_t threadLargestVelMesh = 0;
         size_t threadLargestSizePower = 0;
         SpatialCell *SC;
#pragma omp for
         for (uint i=0; i<nCells; ++i) {
            SC = mpiGrid[cells[i]];

            // Clear hashmaps here, as batch operation isn't yet operational
            SC->velocity_block_with_content_map->clear<false>(Hashinator::targets::device,stream,std::pow(2,SC->vbwcl_sizePower));
            SC->velocity_block_with_no_content_map->clear<false>(Hashinator::targets::device,stream,std::pow(2,SC->vbwncl_sizePower));
            SC->velocity_block_with_content_list_size = 0;
            // Note: clear before applying reservation

            SC->updateSparseMinValue(popID);
            SC->applyReservation(popID);

            // Store values and pointers
            host_vmeshes[i] = SC->dev_get_velocity_mesh(popID);
            host_VBCs[i] = SC->dev_get_velocity_blocks(popID);
            host_minValues[i] = SC->getVelocityBlockMinValue(popID);
            host_allMaps[2*i  ] = SC->dev_velocity_block_with_content_map;
            host_allMaps[2*i+1] = SC->dev_velocity_block_with_no_content_map;
            // host_sizePowers[2*i  ] = SC->vbwcl_sizePower;
            // host_sizePowers[2*i+1] = SC->vbwncl_sizePower;

            // Gather largest values
            vmesh::VelocityMesh* vmesh = SC->get_velocity_mesh(popID);
            threadLargestVelMesh = threadLargestVelMesh > vmesh->size() ? threadLargestVelMesh : vmesh->size();
            threadLargestSizePower = threadLargestSizePower > SC->vbwcl_sizePower ? threadLargestSizePower : SC->vbwcl_sizePower;
            threadLargestSizePower = threadLargestSizePower > SC->vbwncl_sizePower ? threadLargestSizePower : SC->vbwncl_sizePower;
         }
         #pragma omp critical
         {
            largestVelMesh = threadLargestVelMesh > largestVelMesh ? threadLargestVelMesh : largestVelMesh;
            largestSizePower = threadLargestSizePower > largestSizePower ? threadLargestSizePower : largestSizePower;
         }
      }
      sparsityTimer.stop();

      // Copy pointers and counters over to device
      CHK_ERR( gpuMemcpy(dev_allMaps, host_allMaps, 2*nCells*sizeof(Hashinator::Hashmap<vmesh::GlobalID,vmesh::LocalID>*), gpuMemcpyHostToDevice) );
      //CHK_ERR( gpuMemcpy(dev_sizePowers, host_sizePowers, 2*nCells*sizeof(size_t), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(dev_minValues, host_minValues, nCells*sizeof(Real), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(dev_vmeshes, host_vmeshes, nCells*sizeof(vmesh::VelocityMesh*), gpuMemcpyHostToDevice) );
      CHK_ERR( gpuMemcpy(dev_VBCs, host_VBCs, nCells*sizeof(vmesh::VelocityBlockContainer*), gpuMemcpyHostToDevice) );

      // First task: bulk clear of hash maps
      // velocity_block_with_content_map->clear<false>(Hashinator::targets::device,stream,std::pow(2,vbwcl_sizePower));
      // velocity_block_with_no_content_map->clear<false>(Hashinator::targets::device,stream,std::pow(2,vbwncl_sizePower));
      //const gpuStream_t stream = gpu_getStream();
      // Note: the clear call requires a DeviceHasher instance so is called as a method.
      //mpiGrid[cells[0]]->velocity_block_with_content_map->batch_clear((size_t)(2*nCells), dev_allMaps, dev_sizePowers, largestSizePower, stream);

      // Batch gather GID-LID-pairs into two maps
      phiprof::Timer blockKernelTimer {"update content lists kernel 1"};
      const uint vlasiBlocksPerWorkUnit = 1;
      // ceil int division
      const uint launchBlocks = 1 + ((largestVelMesh - 1) / vlasiBlocksPerWorkUnit);
      dim3 grid(nCells,launchBlocks,1);

      batch_update_velocity_block_content_lists_kernel<<<grid, (vlasiBlocksPerWorkUnit * WID3), 0, 0>>> (
         dev_vmeshes,
         dev_VBCs,
         dev_allMaps,
         dev_minValues
         );
      CHK_ERR( gpuPeekAtLastError() );
      CHK_ERR( gpuDeviceSynchronize() );
      blockKernelTimer.stop();

      phiprof::Timer blocklistTimer {"update content lists extract"};
#pragma omp parallel
      {
         const gpuStream_t stream = gpu_getStream();
         split::SplitInfo info;
         SpatialCell *SC;
#pragma omp for
         for (uint i=0; i<nCells; ++i) {
            SC = mpiGrid[cells[i]];
            //SC->update_velocity_block_content_lists(popID);

            // Now extract values from the map
            SC->velocity_block_with_content_map->extractAllKeysLoop(*(SC->dev_velocity_block_with_content_list),stream);
            // split::SplitInfo info;
            // Hashinator::MapInfo info_m;
            // TODO: in batch extraction loop, store end size in batch array
            SC->velocity_block_with_content_list->copyMetadata(&info, stream);
            CHK_ERR( gpuStreamSynchronize(stream) );
            SC->velocity_block_with_content_list_size = info.size;
            //velocity_block_with_content_list_size = velocity_block_with_content_list->size();
         }
      }
      blocklistTimer.stop();
   }

   //CHK_ERR( gpuMemcpy(host_contentSizes, dev_contentSizes, nCells*sizeof(vmesh::LocalID), gpuMemcpyDeviceToHost) );

   //    const gpuStream_t stream = gpu_getStream();
   //    // phiprof::Timer updateListsTimer {"GPU update spatial cell block lists"};

   //    phiprof::Timer reservationTimer {"GPU apply reservation"};
   //    applyReservation(popID);
   //    reservationTimer.stop();

   //    phiprof::Timer clearTimer {"GPU clear maps"};
   //    velocity_block_with_content_list_size = 0;
   //    velocity_block_with_content_map->clear<false>(Hashinator::targets::device,stream,std::pow(2,vbwcl_sizePower));
   //    velocity_block_with_no_content_map->clear<false>(Hashinator::targets::device,stream,std::pow(2,vbwncl_sizePower));
   //    CHK_ERR( gpuStreamSynchronize(stream) );
   //    clearTimer.stop();
   //    phiprof::Timer sizeTimer {"GPU size"};
   //    const uint nBlocks = populations[popID].vmesh->size();
   //    if (nBlocks==0) {
   //       return;
   //    }
   //    CHK_ERR( gpuStreamSynchronize(stream) );
   //    sizeTimer.stop();
   //    phiprof::Timer updateListsTimer {"GPU update spatial cell block lists"};
   //    const Real velocity_block_min_value = getVelocityBlockMinValue(popID);
   //    // Each GPU block / workunit can manage several Vlasiator velocity blocks at once. (TODO FIX)
   //    //const uint vlasiBlocksPerWorkUnit = WARPSPERBLOCK * GPUTHREADS / WID3;
   //    const uint vlasiBlocksPerWorkUnit = 1;
   //    // ceil int division
   //    const uint launchBlocks = 1 + ((nBlocks - 1) / vlasiBlocksPerWorkUnit);
   //    //std::cerr<<launchBlocks<<" "<<vlasiBlocksPerWorkUnit<<" vmesh "<<populations[popID].dev_vmesh<<" vbc "<<populations[popID].dev_blockContainer<<" wcm "<<dev_velocity_block_with_content_map<<" wncm "<<dev_velocity_block_with_no_content_map<<" minval "<<velocity_block_min_value<<std::endl;

   //    // Third argument specifies the number of bytes in *shared memory* that is
   //    // dynamically allocated per block for this call in addition to the statically allocated memory.
   //    //update_velocity_block_content_lists_kernel<<<launchBlocks, WID3, WID3*sizeof(int), stream>>> (
   //    update_velocity_block_content_lists_kernel<<<launchBlocks, (vlasiBlocksPerWorkUnit * WID3), 0, stream>>> (
   //       populations[popID].dev_vmesh,
   //       populations[popID].dev_blockContainer,
   //       dev_velocity_block_with_content_map,
   //       dev_velocity_block_with_no_content_map,
   //       velocity_block_min_value
   //       );
   //    CHK_ERR( gpuPeekAtLastError() );

   //    // Now extract values from the map
   //    velocity_block_with_content_map->extractAllKeysLoop(*dev_velocity_block_with_content_list,stream);
   //    // split::SplitInfo info;
   //    // Hashinator::MapInfo info_m;
   //    split::SplitInfo info;
   //    velocity_block_with_content_list->copyMetadata(&info, stream);
   //    CHK_ERR( gpuStreamSynchronize(stream) );
   //    velocity_block_with_content_list_size = info.size;
   //    //velocity_block_with_content_list_size = velocity_block_with_content_list->size();
   // }
