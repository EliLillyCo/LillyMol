// Copied from https://gist.github.com/sjmduncan/2a11d818ca019b379f2b8d8e3646413c

#include <iostream>
#include <vector>

#include "kabsch.hpp"

namespace kabsch {

using std::cerr;
using std::vector;
using Eigen::Affine3d;
using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::Vector3d;

/*
 * Estimate _tx such that:  _tgt = _tx * _src
 * Requires at least four pairs of points.
 * See:
 *   https://en.wikipedia.org/wiki/Kabsch_algorithm
 *   https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp
 *   Closed-form solution of absolute orientation using unit quaternions - Horn 1987
 */
bool kabsch(const vector<Vector3d>& _src, const vector<Vector3d>& _tgt, Affine3d& _tx, double& _residual_rms) {
  if (_tgt.size() != _src.size()) {
    cerr << "Kabsch: _tgt.size() != _src.size() " << _tgt.size() << " and " << _src.size() <<
            " check input data\n";
    return false;
  }
  if (_tgt.size() < 4) {
    cerr << "Kabsch: need at least four point pairs (more is better). Only got " << _tgt.size() << " points\n";
    return false;
  }

  const size_t n = _tgt.size();
  MatrixXd src(3, n), tgt(3, n);
  for (size_t i = 0; i < n; i++) {
    tgt.col(i) = _tgt[i];
    src.col(i) = _src[i];
  }

  Vector3d centroid_src(0, 0, 0), centroid_tgt(0, 0, 0);
  for (size_t c = 0; c < n; c++) {
    centroid_src += src.col(c);
    centroid_tgt += tgt.col(c);
  }
  centroid_src /= n;
  centroid_tgt /= n;
  for (size_t c = 0; c < n; c++) {
    src.col(c) -= centroid_src;
    tgt.col(c) -= centroid_tgt;
  }

  // cerr << "centroid_src\n" << centroid_src << "\ncentroid_tgt\n" << centroid_tgt << '\n';

  double scale = 1.0;
  double mean_dist_src = 0, mean_dist_tgt = 0;
  {
    for (size_t c = 0; c < n; c++) {
      mean_dist_src += src.col(c).norm();
      mean_dist_tgt += tgt.col(c).norm();
    }
    scale = mean_dist_tgt / mean_dist_src;
    centroid_tgt /= scale;
    cerr << "Kabsch: scale=" << scale << '\n';
  }

  /*
   * Kabsch estimates the rotation matrix which minimizes the RMSD.
   * The error is a function of both the angle and the length of the
   * vector - that is the distance from the centroid.
   *
   * This selective normalization avoids biasing the result towards 
   * points which are further than the average distance from the centroid.
   */
  for (size_t c = 0; c < n; c++) {
    if(src.col(c).norm() > mean_dist_src)
      src.col(c).normalize();
    else
      src.col(c) /= mean_dist_src;

    if(tgt.col(c).norm() > mean_dist_tgt)
      tgt.col(c).normalize();
    else
      tgt.col(c) /= mean_dist_tgt;
  }

  Matrix3d rotation = Matrix3d::Identity();
  {
    MatrixXd cov = src * tgt.transpose();
    Eigen::JacobiSVD<MatrixXd> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (d > 0) d = 1.0;
    else d = -1.0;
    Matrix3d I = Matrix3d::Identity();
    I(2, 2) = d;
    rotation = svd.matrixV() * I * svd.matrixU().transpose();
  }

  _tx.setIdentity();
  _tx.linear() = scale * rotation;
  _tx.translation() = scale * (centroid_tgt - rotation * centroid_src);

  _residual_rms = 0;
  for (size_t c = 0; c < n; c++) {
    auto err = (_tgt[c] - _tx * _src[c]).norm();
    _residual_rms += err * err;
//    cerr << " row " << c << " err " << err << '\n';
  }
  _residual_rms /= n;
  _residual_rms = sqrt(_residual_rms);

  cerr << "Kabsch: RMSD=" << _residual_rms << '\n';

  return true;
}

bool kabsch(const vector<Vector3d>& _src, const vector<Vector3d>& _tgt, Affine3d& _tx) {
  if (_tgt.size() != _src.size()) {
    cerr << "Kabsch: _tgt.size() != _src.size() " << _tgt.size() << " and " << _src.size() <<
            " check input data\n";
    return false;
  }
  if (_tgt.size() < 4) {
    cerr << "Kabsch: need at least four point pairs (more is better). Only got " << _tgt.size() << " points\n";
    return false;
  }

  const size_t n = _tgt.size();
  MatrixXd src(3, n), tgt(3, n);
  for (size_t i = 0; i < n; i++) {
    tgt.col(i) = _tgt[i];
    src.col(i) = _src[i];
  }

  Vector3d centroid_src(0, 0, 0), centroid_tgt(0, 0, 0);
  for (size_t c = 0; c < n; c++) {
    centroid_src += src.col(c);
    centroid_tgt += tgt.col(c);
  }
  centroid_src /= n;
  centroid_tgt /= n;
  for (size_t c = 0; c < n; c++) {
    src.col(c) -= centroid_src;
    tgt.col(c) -= centroid_tgt;
  }

  double scale = 1.0;
  double mean_dist_src = 0, mean_dist_tgt = 0;
  {
    for (size_t c = 0; c < n; c++) {
      mean_dist_src += src.col(c).norm();
      mean_dist_tgt += tgt.col(c).norm();
    }
    scale = mean_dist_tgt / mean_dist_src;
    centroid_tgt /= scale;
    // cerr << "Kabsch: scale=" << scale << '\n';
  }

  /*
   * Kabsch estimates the rotation matrix which minimizes the RMSD.
   * The error is a function of both the angle and the length of the
   * vector - that is the distance from the centroid.
   *
   * This selective normalization avoids biasing the result towards 
   * points which are further than the average distance from the centroid.
   */
  for (size_t c = 0; c < n; c++) {
    if(src.col(c).norm() > mean_dist_src)
      src.col(c).normalize();
    else
      src.col(c) /= mean_dist_src;

    if(tgt.col(c).norm() > mean_dist_tgt)
      tgt.col(c).normalize();
    else
      tgt.col(c) /= mean_dist_tgt;
  }

  Matrix3d rotation = Matrix3d::Identity();
  {
    MatrixXd cov = src * tgt.transpose();
    Eigen::JacobiSVD<MatrixXd> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (d > 0) d = 1.0;
    else d = -1.0;
    Matrix3d I = Matrix3d::Identity();
    I(2, 2) = d;
    rotation = svd.matrixV() * I * svd.matrixU().transpose();
  }

  _tx.setIdentity();
  _tx.linear() = scale * rotation;
  _tx.translation() = scale * (centroid_tgt - rotation * centroid_src);

  return true;
}

bool kabsch(MatrixXd& src, MatrixXd& tgt, Affine3d& _tx) {
  if (tgt.size() != src.size()) {
    cerr << "Kabsch: tgt.size() != src.size() " << tgt.size() << " and " << src.size() <<
            " check input data\n";
    return false;
  }
  if (tgt.size() < 4) {
    cerr << "Kabsch: need at least four point pairs (more is better). Only got " << tgt.size() << " points\n";
    return false;
  }

  const size_t n = tgt.cols();

  Vector3d centroid_src(0, 0, 0), centroid_tgt(0, 0, 0);
  for (size_t c = 0; c < n; c++) {
    centroid_src += src.col(c);
    centroid_tgt += tgt.col(c);
  }
  centroid_src /= n;
  centroid_tgt /= n;
  for (size_t c = 0; c < n; c++) {
    src.col(c) -= centroid_src;
    tgt.col(c) -= centroid_tgt;
  }

  double scale = 1.0;
  double mean_dist_src = 0, mean_dist_tgt = 0;
  {
    for (size_t c = 0; c < n; c++) {
      mean_dist_src += src.col(c).norm();
      mean_dist_tgt += tgt.col(c).norm();
    }
    scale = mean_dist_tgt / mean_dist_src;
    centroid_tgt /= scale;
    // cerr << "Kabsch: scale=" << scale << '\n';
  }

  /*
   * Kabsch estimates the rotation matrix which minimizes the RMSD.
   * The error is a function of both the angle and the length of the
   * vector - that is the distance from the centroid.
   *
   * This selective normalization avoids biasing the result towards 
   * points which are further than the average distance from the centroid.
   */
  for (size_t c = 0; c < n; c++) {
    if(src.col(c).norm() > mean_dist_src)
      src.col(c).normalize();
    else
      src.col(c) /= mean_dist_src;

    if(tgt.col(c).norm() > mean_dist_tgt)
      tgt.col(c).normalize();
    else
      tgt.col(c) /= mean_dist_tgt;
  }

  Matrix3d rotation = Matrix3d::Identity();
  {
    MatrixXd cov = src * tgt.transpose();
    Eigen::JacobiSVD<MatrixXd> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (d > 0) d = 1.0;
    else d = -1.0;
    Matrix3d I = Matrix3d::Identity();
    I(2, 2) = d;
    rotation = svd.matrixV() * I * svd.matrixU().transpose();
  }

  _tx.setIdentity();
  _tx.linear() = scale * rotation;
  _tx.translation() = scale * (centroid_tgt - rotation * centroid_src);

  return true;
}

}  // namespace kabsch
