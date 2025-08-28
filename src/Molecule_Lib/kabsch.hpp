#ifndef MOLECULE_LIB_KABSCH_H_
#define MOLECULE_LIB_KABSCH_H_

#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/src/Geometry/Transform.h"

namespace kabsch {

// Kabsch based superimposition. The original implementation was with u3b.f
// That is highly efficient, no temporary arrays, but since it is Fortran
// I want an alternative if needed.

// Input matrices must be as below for molecules `m1` and `m2`.
//  Eigen::MatrixXd from(3, matoms), to(3, matoms);
//  Eigen::MatrixXd to(matoms, 3);
//  for (int i = 0; i < matoms; ++i) {
//    const Atom& a = m1[i];
//    from(0, i) = a.x();
//    from(1, i) = a.y();
//    from(2, i) = a.z();

//  same thing for the `to` vector.
// 
//  Eigen::Affine3d tx;
//  kabsch::kabsch(from, to, tx);
// 
//  for (int i = 0; i < matoms; i++) {
//    const Atom& a2 = m2[i];
//    Eigen::Vector3d f(a2.x(), a2.y(), a2.z());
//    auto rotated = tx * f;
//    m2.setxyz(i, rotated(0), rotated(1), rotated(2));

//    If computing an RMS.
//    const Atom& a1 = m1[i];
//    Eigen::Vector3d t(a1.x(), a1.y(), a1.z());
//    double err = (t - rotated).norm();
//    residual_rms += (err * err);


bool kabsch(const std::vector<Eigen::Vector3d>& _tgt, const std::vector<Eigen::Vector3d>& _src, Eigen::Affine3d& _tx, double& _residual_rms);
bool kabsch(const std::vector<Eigen::Vector3d>& _tgt, const std::vector<Eigen::Vector3d>& _src, Eigen::Affine3d& _tx);
bool kabsch(Eigen::MatrixXd& tgt, Eigen::MatrixXd& src, Eigen::Affine3d& tx);

}  // namespace kabsch

#endif // MOLECULE_LIB_KABSCH_H_
