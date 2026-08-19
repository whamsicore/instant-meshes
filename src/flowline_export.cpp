/*
    File: tmp/tools/instant_meshes_fork/src/flowline_export.cpp
    Purpose: Batch export of Instant Meshes orientation flow lines (same tracer as GUI).
             BLDR: used for Forge View without Blender inventing glyphs.
*/

#include "flowline_export.h"
#include "field.h"
#include "meshio.h"
#include "pcg32.h"
#include <fstream>
#include <limits>
#include <set>
#include <vector>

void export_orientation_flowlines_obj(
    const MultiResolutionHierarchy &mRes,
    BVH *bvh,
    const MeshStats &stats,
    int rosy,
    Float density_scale,
    const std::string &output_obj) {

    if (!bvh)
        throw std::runtime_error("export_orientation_flowlines_obj: BVH required");
    if (mRes.F().size() == 0)
        throw std::runtime_error("export_orientation_flowlines_obj: point clouds unsupported");

    auto rotate = rosy == 2 ? rotate180 : (rosy == 4 ? rotate90 : rotate60);

    const MatrixXf &V = mRes.V();
    const MatrixXf &N = mRes.N();
    const MatrixXf &Q = mRes.Q();
    const MatrixXu &F = mRes.F();

    Float edgeLength = stats.mAverageEdgeLength * 2;
    const Float targetLength = edgeLength * 20;
    const Float stepSize = edgeLength / 2.0f;
    const Float thickness = edgeLength / 3.0f;
    const Float eps = edgeLength / 10;
    const uint32_t nSteps = (uint32_t) (targetLength / stepSize);
    uint32_t nLines = (uint32_t) (
        stats.mSurfaceArea / (edgeLength * thickness) * 0.02f * density_scale);
    if (nLines < 64)
        nLines = 64;
    if (nLines > 8000)
        nLines = 8000;

    const size_t nTriangles = (size_t) nLines * nSteps * 2;
    const size_t nVertices = (size_t) nLines * (nSteps + 1) * 2;
    MatrixXu indices(3, (uint32_t) nTriangles);
    MatrixXf position(3, (uint32_t) nVertices);
    indices.setZero();
    position.setZero();

    std::vector<tbb::spin_mutex> locks((size_t) F.cols() * (size_t) (rosy / 2));
    cout << "Tracing " << nLines << " orientation flow lines for export ..";
    cout.flush();
    Timer<> timer;

    tbb::parallel_for(
        tbb::blocked_range<uint32_t>(0u, nLines, GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
            std::set<uint32_t> locked;
            pcg32 rng;
            for (uint32_t k = range.begin(); k != range.end(); ++k) {
                rng.seed(1, k);
                int nTries = 0;
                Float height = rng.nextFloat() * eps + eps;
                do {
                    size_t vertexIdx = (size_t) k * 2 * (nSteps + 1);
                    size_t triangleIdx = (size_t) k * nSteps * 2;
                    uint32_t startIdx = rng.nextUInt((uint32_t) V.cols());

                    Vector3f p = V.col(startIdx), n = N.col(startIdx), q = Q.col(startIdx);
                    for (int i = 0, nrot = (int) rng.nextUInt((uint32_t) rosy); i < nrot; ++i)
                        q = rotate(q, n);

                    Vector3f t = n.cross(q);
                    position.col((uint32_t) vertexIdx++) = p + n * height;
                    position.col((uint32_t) vertexIdx++) = p + n * height;

                    bool fail = false;
                    locked.clear();
                    for (uint32_t step = 0; step < nSteps; ++step) {
                        p += q * stepSize;
                        Vector3f o = p;
                        Ray ray1(o, n, 0, edgeLength * 2);
                        Ray ray2(o, -n, 0, edgeLength * 2);
                        uint32_t idx1 = 0, idx2 = 0;
                        Float t1 = 0, t2 = 0;
                        Vector2f uv1, uv2;
                        bool hit1 = bvh->rayIntersect(ray1, idx1, t1, &uv1);
                        bool hit2 = bvh->rayIntersect(ray2, idx2, t2, &uv2);
                        if (!hit1 && !hit2) {
                            fail = true;
                            break;
                        }
                        if (t2 < t1) {
                            p = ray2(t2);
                            uv1 = uv2;
                            t1 = t2;
                            idx1 = idx2;
                        } else {
                            p = ray1(t1);
                        }

                        n = ((1 - uv1.sum()) * N.col(F(0, idx1)) + uv1.x() * N.col(F(1, idx1)) +
                             uv1.y() * N.col(F(2, idx1)))
                                .normalized();
                        Vector3f qp = Q.col(F(0, idx1));
                        q = (q - n.dot(q) * n).normalized();
                        qp = (qp - n.dot(qp) * n).normalized();

                        Vector3f best = Vector3f::Zero();
                        Float best_dp = -std::numeric_limits<Float>::infinity();
                        int best_index = -1;
                        for (int j = 0; j < rosy; ++j) {
                            Float dp = qp.dot(q);
                            if (dp > best_dp) {
                                best_dp = dp;
                                best = qp;
                                best_index = j;
                            }
                            qp = rotate(qp, n);
                        }

                        uint32_t lock_idx = idx1 * (uint32_t) (rosy / 2) + (uint32_t) (best_index % (rosy / 2));
                        if (locked.find(lock_idx) == locked.end()) {
                            if (!locks[lock_idx].try_lock()) {
                                fail = true;
                                break;
                            }
                            locked.insert(lock_idx);
                        }

                        q = (best - n.dot(best) * n).normalized();
                        t = n.cross(q);
                        indices.col((uint32_t) triangleIdx++) << (uint32_t) vertexIdx,
                            (uint32_t) vertexIdx - 2, (uint32_t) vertexIdx - 1;
                        indices.col((uint32_t) triangleIdx++) << (uint32_t) vertexIdx,
                            (uint32_t) vertexIdx - 1, (uint32_t) vertexIdx + 1;

                        Float thicknessP =
                            (1 - std::pow((2 * ((step + 1) / (Float) nSteps) - 1), 2.0f)) * thickness;
                        position.col((uint32_t) vertexIdx++) = p + t * thicknessP + n * height;
                        position.col((uint32_t) vertexIdx++) = p - t * thicknessP + n * height;
                    }
                    if (fail) {
                        for (auto l : locked)
                            locks[l].unlock();
                        if (++nTries > 10) {
                            indices.block(0, (int) (k * nSteps * 2), 3, (int) (nSteps * 2)).setConstant(0);
                            break;
                        }
                        continue;
                    } else {
                        break;
                    }
                } while (true);
            }
        });

    cout << "done. (took " << timeString(timer.value()) << ")" << endl;

    // Drop failed (zero) triangles.
    std::vector<uint32_t> keep;
    keep.reserve(nTriangles);
    for (uint32_t i = 0; i < (uint32_t) nTriangles; ++i) {
        if (indices(0, i) != 0 || indices(1, i) != 0 || indices(2, i) != 0)
            keep.push_back(i);
    }
    MatrixXu F_out(3, (uint32_t) keep.size());
    for (uint32_t i = 0; i < (uint32_t) keep.size(); ++i)
        F_out.col(i) = indices.col(keep[i]);

    write_obj(output_obj, F_out, position);
    cout << "Wrote flow lines OBJ faces=" << keep.size() << " → " << output_obj << endl;
}
