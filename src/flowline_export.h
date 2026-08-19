/*
    File: tmp/tools/instant_meshes_fork/src/flowline_export.h
    Purpose: Export Instant Meshes orientation flow lines for batch / Forge View.
*/

#pragma once

#include "common.h"
#include "hierarchy.h"
#include "bvh.h"
#include "meshstats.h"

void export_orientation_flowlines_obj(
    const MultiResolutionHierarchy &mRes,
    BVH *bvh,
    const MeshStats &stats,
    int rosy,
    Float density_scale,
    const std::string &output_obj);
