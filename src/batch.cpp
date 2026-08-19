/*
    batch.cpp -- command line interface to Instant Meshes

    This file is part of the implementation of

        Instant Field-Aligned Meshes
        Wenzel Jakob, Daniele Panozzo, Marco Tarini, and Olga Sorkine-Hornung
        In ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2015)

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE.txt file.
*/

#include "batch.h"
#include "meshio.h"
#include "dedge.h"
#include "subdivide.h"
#include "meshstats.h"
#include "hierarchy.h"
#include "field.h"
#include "normal.h"
#include "extract.h"
#include "bvh.h"
#include "flowline_export.h"
#include <fstream>
#include <map>
#include <sstream>
#include <tuple>
#include <limits>

static void apply_orientation_constraints_csv(MultiResolutionHierarchy &mRes,
                                              const std::string &path,
                                              int rosy, int posy) {
    if (path.empty())
        return;
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Failed to open orientation constraints: " + path);

    mRes.clearConstraints();
    const MatrixXf &V = mRes.V();
    const MatrixXf &N = mRes.N();
    MatrixXf &CQ = mRes.CQ();
    VectorXf &CQw = mRes.CQw();

    /* Build a coarse hash grid for nearest vertex lookup. */
    Vector3f vmin = V.rowwise().minCoeff();
    Vector3f vmax = V.rowwise().maxCoeff();
    Float extent = (vmax - vmin).norm();
    Float cell = std::max(extent / 64.0f, (Float) 1e-4);
    std::map<std::tuple<int,int,int>, std::vector<uint32_t>> buckets;
    for (uint32_t i = 0; i < (uint32_t) V.cols(); ++i) {
        Vector3f p = V.col(i);
        int ix = (int) std::floor((p.x() - vmin.x()) / cell);
        int iy = (int) std::floor((p.y() - vmin.y()) / cell);
        int iz = (int) std::floor((p.z() - vmin.z()) / cell);
        buckets[std::make_tuple(ix, iy, iz)].push_back(i);
    }

    auto nearest = [&](const Vector3f &p) -> uint32_t {
        int ix = (int) std::floor((p.x() - vmin.x()) / cell);
        int iy = (int) std::floor((p.y() - vmin.y()) / cell);
        int iz = (int) std::floor((p.z() - vmin.z()) / cell);
        Float best_d = std::numeric_limits<Float>::infinity();
        uint32_t best_i = 0;
        bool found = false;
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    auto it = buckets.find(std::make_tuple(ix+dx, iy+dy, iz+dz));
                    if (it == buckets.end())
                        continue;
                    for (uint32_t j : it->second) {
                        Float d = (V.col(j) - p).squaredNorm();
                        if (d < best_d) {
                            best_d = d;
                            best_i = j;
                            found = true;
                        }
                    }
                }
        if (!found) {
            for (uint32_t j = 0; j < (uint32_t) V.cols(); ++j) {
                Float d = (V.col(j) - p).squaredNorm();
                if (d < best_d) {
                    best_d = d;
                    best_i = j;
                }
            }
        }
        return best_i;
    };

    uint32_t applied = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        for (char &c : line)
            if (c == ',')
                c = ' ';
        std::istringstream iss(line);
        Float px, py, pz, qx, qy, qz, w;
        if (!(iss >> px >> py >> pz >> qx >> qy >> qz >> w))
            continue;
        if (w <= 0)
            continue;
        Vector3f p(px, py, pz);
        Vector3f q(qx, qy, qz);
        uint32_t v = nearest(p);
        Vector3f n = N.col(v);
        q -= n * n.dot(q);
        Float qn = q.norm();
        if (qn < RCPOVERFLOW)
            continue;
        q /= qn;
        CQ.col(v) = q;
        CQw[v] = std::min((Float) 1.0, std::max(CQw[v], w));
        ++applied;
    }
    mRes.propagateConstraints(rosy, posy);
    cout << "Applied orientation constraints from " << path
         << " (" << applied << " samples mapped to vertices)." << endl;
}

static void apply_position_constraints_csv(MultiResolutionHierarchy &mRes,
                                           const std::string &path,
                                           int rosy, int posy) {
    /* Edge-brush: same CSV as CQ, also writes CO/COw. Does not clearConstraints
       so TokenRig orientation samples survive. */
    if (path.empty())
        return;
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Failed to open position constraints: " + path);

    const MatrixXf &V = mRes.V();
    const MatrixXf &N = mRes.N();
    MatrixXf &CQ = mRes.CQ();
    VectorXf &CQw = mRes.CQw();
    MatrixXf &CO = mRes.CO();
    VectorXf &COw = mRes.COw();

    Vector3f vmin = V.rowwise().minCoeff();
    Vector3f vmax = V.rowwise().maxCoeff();
    Float extent = (vmax - vmin).norm();
    Float cell = std::max(extent / 64.0f, (Float) 1e-4);
    std::map<std::tuple<int,int,int>, std::vector<uint32_t>> buckets;
    for (uint32_t i = 0; i < (uint32_t) V.cols(); ++i) {
        Vector3f p = V.col(i);
        int ix = (int) std::floor((p.x() - vmin.x()) / cell);
        int iy = (int) std::floor((p.y() - vmin.y()) / cell);
        int iz = (int) std::floor((p.z() - vmin.z()) / cell);
        buckets[std::make_tuple(ix, iy, iz)].push_back(i);
    }

    auto nearest = [&](const Vector3f &p) -> uint32_t {
        int ix = (int) std::floor((p.x() - vmin.x()) / cell);
        int iy = (int) std::floor((p.y() - vmin.y()) / cell);
        int iz = (int) std::floor((p.z() - vmin.z()) / cell);
        Float best_d = std::numeric_limits<Float>::infinity();
        uint32_t best_i = 0;
        bool found = false;
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    auto it = buckets.find(std::make_tuple(ix+dx, iy+dy, iz+dz));
                    if (it == buckets.end())
                        continue;
                    for (uint32_t j : it->second) {
                        Float d = (V.col(j) - p).squaredNorm();
                        if (d < best_d) {
                            best_d = d;
                            best_i = j;
                            found = true;
                        }
                    }
                }
        if (!found) {
            for (uint32_t j = 0; j < (uint32_t) V.cols(); ++j) {
                Float d = (V.col(j) - p).squaredNorm();
                if (d < best_d) {
                    best_d = d;
                    best_i = j;
                }
            }
        }
        return best_i;
    };

    uint32_t applied = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        for (char &c : line)
            if (c == ',')
                c = ' ';
        std::istringstream iss(line);
        Float px, py, pz, qx, qy, qz, w;
        if (!(iss >> px >> py >> pz >> qx >> qy >> qz >> w))
            continue;
        if (w <= 0)
            continue;
        Vector3f p(px, py, pz);
        Vector3f q(qx, qy, qz);
        uint32_t v = nearest(p);
        Vector3f n = N.col(v);
        q -= n * n.dot(q);
        Float qn = q.norm();
        if (qn < RCPOVERFLOW)
            continue;
        q /= qn;
        CQ.col(v) = q;
        CQw[v] = std::min((Float) 1.0, std::max(CQw[v], w));
        CO.col(v) = p;
        COw[v] = std::min((Float) 1.0, std::max(COw[v], w));
        ++applied;
    }
    mRes.propagateConstraints(rosy, posy);
    cout << "Applied position (edge-brush) constraints from " << path
         << " (" << applied << " samples mapped to vertices)." << endl;
}

static void write_orientation_field_csv(const std::string &path,
                                        const MultiResolutionHierarchy &mRes,
                                        const std::map<uint32_t, uint32_t> &sing,
                                        int rosy) {
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Failed to open orientation field export: " + path);
    const MatrixXf &V = mRes.V();
    const MatrixXf &N = mRes.N();
    const MatrixXf &Q = mRes.Q();
    out << "# instant_meshes_fork_orientation_field\n";
    out << "# rosy=" << rosy << "\n";
    out << "# vertex_count=" << V.cols() << "\n";
    out << "# singularity_count=" << sing.size() << "\n";
    out << "# columns: i,vx,vy,vz,nx,ny,nz,qx,qy,qz,is_singularity\n";
    for (uint32_t i = 0; i < (uint32_t) V.cols(); ++i) {
        Vector3f v = V.col(i);
        Vector3f n = N.col(i);
        Vector3f q = Q.col(i);
        int is_sing = sing.count(i) ? 1 : 0;
        out << i << ','
            << v.x() << ',' << v.y() << ',' << v.z() << ','
            << n.x() << ',' << n.y() << ',' << n.z() << ','
            << q.x() << ',' << q.y() << ',' << q.z() << ','
            << is_sing << '\n';
    }
    cout << "Wrote orientation field CSV (" << V.cols() << " verts, "
         << sing.size() << " singularities): " << path << endl;
}

void batch_process(const std::string &input, const std::string &output,
                   int rosy, int posy, Float scale, int face_count,
                   int vertex_count, Float creaseAngle, bool extrinsic,
                   bool align_to_boundaries, int smooth_iter, int knn_points,
                   bool pure_quad, bool deterministic,
                   const std::string &field_export, bool field_only,
                   const std::string &flowline_export, Float flowline_density,
                   const std::string &orientation_constraints,
                   const std::string &position_constraints) {
    cout << endl;
    cout << "Running in batch mode:" << endl;
    cout << "   Input file             = " << input << endl;
    cout << "   Output file            = " << (output.empty() ? "(none)" : output) << endl;
    cout << "   Field export CSV       = " << (field_export.empty() ? "(none)" : field_export) << endl;
    cout << "   Flowline export OBJ    = " << (flowline_export.empty() ? "(none)" : flowline_export) << endl;
    cout << "   Flowline density scale = " << flowline_density << endl;
    cout << "   Orientation CQ CSV     = " << (orientation_constraints.empty() ? "(none)" : orientation_constraints) << endl;
    cout << "   Position CO CSV        = " << (position_constraints.empty() ? "(none)" : position_constraints) << endl;
    cout << "   Field-only (no extrude)= " << (field_only ? "yes" : "no") << endl;
    cout << "   Rotation symmetry type = " << rosy << endl;
    cout << "   Position symmetry type = " << (posy==3?6:posy) << endl;
    cout << "   Crease angle threshold = ";
    if (creaseAngle > 0)
        cout << creaseAngle << endl;
    else
        cout << "disabled" << endl;
    cout << "   Extrinsic mode         = " << (extrinsic ? "enabled" : "disabled") << endl;
    cout << "   Align to boundaries    = " << (align_to_boundaries ? "yes" : "no") << endl;
    cout << "   kNN points             = " << knn_points << " (only applies to point clouds)"<< endl;
    cout << "   Fully deterministic    = " << (deterministic ? "yes" : "no") << endl;
    if (posy == 4)
        cout << "   Output mode            = " << (pure_quad ? "pure quad mesh" : "quad-dominant mesh") << endl;
    cout << endl;

    MatrixXu F;
    MatrixXf V, N;
    VectorXf A;
    std::set<uint32_t> crease_in, crease_out;
    BVH *bvh = nullptr;
    AdjacencyMatrix adj = nullptr;

    /* Load the input mesh */
    load_mesh_or_pointcloud(input, F, V, N);

    bool pointcloud = F.size() == 0;

    Timer<> timer;
    MeshStats stats = compute_mesh_stats(F, V, deterministic);

    if (pointcloud) {
        bvh = new BVH(&F, &V, &N, stats.mAABB);
        bvh->build();
        adj = generate_adjacency_matrix_pointcloud(V, N, bvh, stats, knn_points, deterministic);
        A.resize(V.cols());
        A.setConstant(1.0f);
    }

    if (scale < 0 && vertex_count < 0 && face_count < 0) {
        if (field_only) {
            /* Keep input average edge length so field preview does not mega-subdivide
               (old face_count=input faces → tiny scale → millions of verts). */
            scale = stats.mAverageEdgeLength;
            cout << "Field-only: defaulting scale to average input edge length "
                 << scale << endl;
        } else {
            cout << "No target vertex count/face count/scale argument provided. "
                    "Setting to the default of 1/16 * input vertex count." << endl;
            vertex_count = V.cols() / 16;
        }
    }

    if (scale > 0 && face_count < 0 && vertex_count < 0) {
        Float face_area = posy == 4 ? (scale*scale) : (std::sqrt(3.f)/4.f*scale*scale);
        face_count = stats.mSurfaceArea / face_area;
        vertex_count = posy == 4 ? face_count : (face_count / 2);
    } else if (face_count > 0) {
        Float face_area = stats.mSurfaceArea / face_count;
        vertex_count = posy == 4 ? face_count : (face_count / 2);
        scale = posy == 4 ? std::sqrt(face_area) : (2*std::sqrt(face_area * std::sqrt(1.f/3.f)));
    } else if (vertex_count > 0) {
        face_count = posy == 4 ? vertex_count : (vertex_count * 2);
        Float face_area = stats.mSurfaceArea / face_count;
        scale = posy == 4 ? std::sqrt(face_area) : (2*std::sqrt(face_area * std::sqrt(1.f/3.f)));
    }

    cout << "Output mesh goals (approximate)" << endl;
    cout << "   Vertex count           = " << vertex_count << endl;
    cout << "   Face count             = " << face_count << endl;
    cout << "   Edge length            = " << scale << endl;

    MultiResolutionHierarchy mRes;

    if (!pointcloud) {
        /* Subdivide the mesh if necessary — skip for field-only previews to avoid
           exploding Trellis characters to millions of verts. */
        VectorXu V2E, E2E;
        VectorXb boundary, nonManifold;
        bool need_subdiv =
            !field_only &&
            (stats.mMaximumEdgeLength*2 > scale ||
             stats.mMaximumEdgeLength > stats.mAverageEdgeLength * 2);
        if (need_subdiv) {
            cout << "Input mesh is too coarse for the desired output edge length "
                    "(max input mesh edge length=" << stats.mMaximumEdgeLength
                 << "), subdividing .." << endl;
            build_dedge(F, V, V2E, E2E, boundary, nonManifold);
            subdivide(F, V, V2E, E2E, boundary, nonManifold, std::min(scale/2, (Float) stats.mAverageEdgeLength*2), deterministic);
        }

        /* Compute a directed edge data structure */
        build_dedge(F, V, V2E, E2E, boundary, nonManifold);

        /* Compute adjacency matrix */
        adj = generate_adjacency_matrix_uniform(F, V2E, E2E, nonManifold);

        /* Compute vertex/crease normals */
        if (creaseAngle >= 0)
            generate_crease_normals(F, V, V2E, E2E, boundary, nonManifold, creaseAngle, N, crease_in);
        else
            generate_smooth_normals(F, V, V2E, E2E, nonManifold, N);

        /* Compute dual vertex areas */
        compute_dual_vertex_areas(F, V, V2E, E2E, nonManifold, A);

        mRes.setE2E(std::move(E2E));
    }

    /* Build multi-resolution hierarrchy */
    mRes.setAdj(std::move(adj));
    mRes.setF(std::move(F));
    mRes.setV(std::move(V));
    mRes.setA(std::move(A));
    mRes.setN(std::move(N));
    mRes.setScale(scale);
    mRes.build(deterministic);
    mRes.resetSolution();

    if (align_to_boundaries && !pointcloud) {
        mRes.clearConstraints();
        for (uint32_t i=0; i<3*mRes.F().cols(); ++i) {
            if (mRes.E2E()[i] == INVALID) {
                uint32_t i0 = mRes.F()(i%3, i/3);
                uint32_t i1 = mRes.F()((i+1)%3, i/3);
                Vector3f p0 = mRes.V().col(i0), p1 = mRes.V().col(i1);
                Vector3f edge = p1-p0;
                if (edge.squaredNorm() > 0) {
                    edge.normalize();
                    mRes.CO().col(i0) = p0;
                    mRes.CO().col(i1) = p1;
                    mRes.CQ().col(i0) = mRes.CQ().col(i1) = edge;
                    mRes.CQw()[i0] = mRes.CQw()[i1] = mRes.COw()[i0] =
                        mRes.COw()[i1] = 1.0f;
                }
            }
        }
        mRes.propagateConstraints(rosy, posy);
    }

    if (!orientation_constraints.empty()) {
        if (align_to_boundaries)
            cout << "Note: orientation-constraints CSV applied after boundary constraints." << endl;
        apply_orientation_constraints_csv(mRes, orientation_constraints, rosy, posy);
    }

    if (!position_constraints.empty()) {
        if (orientation_constraints.empty() && !align_to_boundaries)
            mRes.clearConstraints();
        else if (align_to_boundaries && orientation_constraints.empty())
            cout << "Note: position-constraints CSV applied after boundary constraints." << endl;
        apply_position_constraints_csv(mRes, position_constraints, rosy, posy);
    }

    if (bvh) {
        bvh->setData(&mRes.F(), &mRes.V(), &mRes.N());
    } else if (smooth_iter > 0 || !flowline_export.empty()) {
        bvh = new BVH(&mRes.F(), &mRes.V(), &mRes.N(), stats.mAABB);
        bvh->build();
    }

    cout << "Preprocessing is done. (total time excluding file I/O: "
         << timeString(timer.reset()) << ")" << endl;

    Optimizer optimizer(mRes, false);
    optimizer.setRoSy(rosy);
    optimizer.setPoSy(posy);
    optimizer.setExtrinsic(extrinsic);

    cout << "Optimizing orientation field .. ";
    cout.flush();
    optimizer.optimizeOrientations(-1);
    optimizer.notify();
    optimizer.wait();
    cout << "done. (took " << timeString(timer.reset()) << ")" << endl;

    std::map<uint32_t, uint32_t> sing;
    compute_orientation_singularities(mRes, sing, extrinsic, rosy);
    cout << "Orientation field has " << sing.size() << " singularities." << endl;
    timer.reset();

    if (!field_export.empty())
        write_orientation_field_csv(field_export, mRes, sing, rosy);

    if (!flowline_export.empty()) {
        if (!bvh) {
            bvh = new BVH(&mRes.F(), &mRes.V(), &mRes.N(), stats.mAABB);
            bvh->build();
        }
        export_orientation_flowlines_obj(mRes, bvh, stats, rosy, flowline_density, flowline_export);
    }

    if (field_only) {
        cout << "Field-only mode: skipping position field + mesh extraction." << endl;
        optimizer.shutdown();
        if (bvh)
            delete bvh;
        return;
    }

    if (output.empty())
        throw std::runtime_error("Batch remesh requires -o/--output (or use --field-only)");

    cout << "Optimizing position field .. ";
    cout.flush();
    optimizer.optimizePositions(-1);
    optimizer.notify();
    optimizer.wait();
    cout << "done. (took " << timeString(timer.reset()) << ")" << endl;

    optimizer.shutdown();

    MatrixXf O_extr, N_extr, Nf_extr;
    std::vector<std::vector<TaggedLink>> adj_extr;
    extract_graph(mRes, extrinsic, rosy, posy, adj_extr, O_extr, N_extr,
                  crease_in, crease_out, deterministic);

    MatrixXu F_extr;
    extract_faces(adj_extr, O_extr, N_extr, Nf_extr, F_extr, posy,
            mRes.scale(), crease_out, true, pure_quad, bvh, smooth_iter);
    cout << "Extraction is done. (total time: " << timeString(timer.reset()) << ")" << endl;

    write_mesh(output, F_extr, O_extr, MatrixXf(), Nf_extr);
    if (bvh)
        delete bvh;
}
