/*
    batch.h -- command line interface to Instant Meshes

    This file is part of the implementation of

        Instant Field-Aligned Meshes
        Wenzel Jakob, Daniele Panozzo, Marco Tarini, and Olga Sorkine-Hornung
        In ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2015)

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE.txt file.
*/

#pragma once

#include "common.h"

extern void batch_process(const std::string &input, const std::string &output,
                          int rosy, int posy, Float scale, int face_count,
                          int vertex_count, Float creaseAngle, bool extrinsic,
                          bool align_to_boundaries, int smooth_iter,
                          int knn_points, bool dominant, bool deterministic,
                          const std::string &field_export = std::string(),
                          bool field_only = false,
                          const std::string &flowline_export = std::string(),
                          Float flowline_density = 1.0f,
                          const std::string &orientation_constraints = std::string(),
                          const std::string &position_constraints = std::string());
