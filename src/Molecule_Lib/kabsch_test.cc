// Tests for kabsch

#include <random>
#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "kabsch.hpp"
#include "molecule.h"

namespace {

using Eigen::Vector3d;
using Eigen::MatrixXd;

using std::cerr;

double
AverageInterAtomicDistance(const Molecule& m) {
  double tot = 0.0;
  int n = 0;
  const int matoms = m.natoms();
  for (int i = 0; i < matoms; ++i) {
    for (int j = i + 1; j < matoms; ++j) {
      distance_t d = m.distance_between_atoms(i, j);
      tot += d;
      n += 1;
    }
  }

  return tot / n;
}

TEST(TestSameAtoms, NoChange) {
  Molecule m1;
  ASSERT_TRUE(m1.build_from_smiles("C(C)(C)C"));
  m1.setxyz(0, -1.0, 0.0, 0.0);
  m1.setxyz(1,  0.0, 0.0, 0.0);
  m1.setxyz(2,  0.0, 1.0, 0.0);
  m1.setxyz(3,  0.0, 0.0, 1.0);

  Molecule m2(m1);

  const int matoms = m1.natoms();

  std::vector<Vector3d> from(matoms), to(matoms);
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m1[i];
    from[i][0] = a.x();
    from[i][1] = a.y();
    from[i][2] = a.z();

    to[i] = from[i];
  }

  Eigen::Affine3d tx;
  double rms;
  ASSERT_TRUE(kabsch::kabsch(from, to, tx, rms)) << "kabsch failed";
  EXPECT_NEAR(rms, 0.0, 1.0e-07);

#ifdef CHECK_TRANSFORM
  double residual_rms = 0.0;
  for (int i = 0; i < matoms; i++) {
    const Atom& a2 = m2[i];
    Eigen::Vector3d f(a2.x(), a2.y(), a2.z());
    auto rotated = tx * f;
    m2.setxyz(i, rotated(0), rotated(1), rotated(2));

    const Atom& a1 = m1[i];
    Eigen::Vector3d t(a1.x(), a1.y(), a1.z());
    double err = (t - tx * f).norm();
    residual_rms += (err * err);
  }
  residual_rms /= matoms;
  residual_rms = sqrt(residual_rms);

  EXPECT_NEAR(residual_rms, 0.0, 1.0e-12);
#endif
}

TEST(TestSameAtoms, Rot90) {
  Molecule m1;
  ASSERT_TRUE(m1.build_from_smiles("C(N)(O)C"));
  m1.setxyz(0, -1.0, 0.0, 0.0);
  m1.setxyz(1,  0.0, 0.0, 0.0);
  m1.setxyz(2,  0.5, 1.0, -0.5);
  m1.setxyz(3,  0.5, 1.0,  0.5);

  Molecule m2(m1);
  m2.setxyz(2, 0.5, -1.0,  0.5);
  m2.setxyz(3, 0.5, -1.0, -0.5);

  ASSERT_TRUE(m1.discern_chirality_from_3d_structure());
  ASSERT_TRUE(m2.discern_chirality_from_3d_structure());
  ASSERT_EQ(m1.unique_smiles(), m2.unique_smiles());

  EXPECT_NEAR(AverageInterAtomicDistance(m1), AverageInterAtomicDistance(m2), 1.0e-05);

#ifdef DEBUG_KABSCH_TEST
  m1.debug_print(std::cerr);
  m2.debug_print(std::cerr);
#endif

  const int matoms = m1.natoms();

  std::vector<Vector3d> from(matoms), to(matoms);
  for (int i = 0; i < matoms; ++i) {
    const Atom& a1 = m1[i];
    to[i](0) = a1.x();
    to[i](1) = a1.y();
    to[i](2) = a1.z();

    const Atom& a2 = m2[i];
    from[i](0) = a2.x();
    from[i](1) = a2.y();
    from[i](2) = a2.z();
  }
#ifdef DEBUG_KABSCH_TEST
  std::cerr << "FROM\n" << from << '\n';
  std::cerr << "TO  \n" << to << '\n';
#endif

  Eigen::Affine3d tx;
  double rms;
  ASSERT_TRUE(kabsch::kabsch(from, to, tx, rms)) << "kabsch failed";

#ifdef DEBUG_KABSCH_TEST
  std::cerr << "translation\n" << tx.translation() << '\n';
  std::cerr << "linear\n" << tx.linear() << '\n';
#endif
  EXPECT_NEAR(rms, 0.0, 1.0e-8);

  double residual_rms = 0.0;
  for (int i = 0; i < matoms; i++) {
    const Atom& a1 = m1[i];
    const Atom& a2 = m2[i];
    Eigen::Vector3d f(a2.x(), a2.y(), a2.z());
    Eigen::Vector3d t(a1.x(), a1.y(), a1.z());
#ifdef DEBUG_KABSCH_TEST
    std::cerr << "f\n" << f << '\n';
    std::cerr << "t\n" << f << '\n';
    std::cerr << "diff\n" << (t - tx * f) << '\n';
#endif
    auto rotated = tx * f;
    m2.setxyz(i, rotated(0), rotated(1), rotated(2));

    double err = (t - rotated).norm();
    residual_rms += (err * err);
  }
  residual_rms /= matoms;
  residual_rms = sqrt(residual_rms);

  EXPECT_NEAR(residual_rms, 0.0, 1.0e-12);
}

TEST(TestSameAtoms, NoChangeV2) {
  Molecule m1;
  ASSERT_TRUE(m1.build_from_smiles("C(C)(C)C"));
  m1.setxyz(0, -1.0, 0.0, 0.0);
  m1.setxyz(1,  0.0, 0.0, 0.0);
  m1.setxyz(2,  0.0, 1.0, 0.0);
  m1.setxyz(3,  0.0, 0.0, 1.0);

  Molecule m2(m1);

  const int matoms = m1.natoms();

  MatrixXd from(3, matoms), to(3, matoms);
  for (int i = 0; i < matoms; ++i) {
    const Atom& a = m1[i];
    from(0, i) = a.x();
    from(1, i) = a.y();
    from(2, i) = a.z();

    to(0, i) = a.x();
    to(1, i) = a.y();
    to(2, i) = a.z();
  }

  Eigen::Affine3d tx;
  ASSERT_TRUE(kabsch::kabsch(from, to, tx)) << "kabsch failed";
  return;

  double residual_rms = 0.0;
  for (int i = 0; i < matoms; i++) {
    const Atom& a2 = m2[i];

    Eigen::Vector3d f(a2.x(), a2.y(), a2.z());
    auto rotated = tx * f;
    m2.setxyz(i, rotated(0), rotated(1), rotated(2));

    const Atom& a1 = m1[i];
    Eigen::Vector3d t(a1.x(), a1.y(), a1.z());
    double err = (t - tx * f).norm();
    residual_rms += (err * err);
  }
  residual_rms /= matoms;
  residual_rms = sqrt(residual_rms);

  EXPECT_NEAR(residual_rms, 0.0, 1.0e-12);
}


TEST(TestKabsch, TestManyRotations) {
  // Acknowledge that random numbers in tests are not a great idea.
  using std::cerr;

  Molecule m1;
  ASSERT_TRUE(m1.build_from_smiles("C{{-0.73,4.777,0.282}}1(=N{{-2.064,4.703,0.191}}C{{-2.521,3.46,0.387}}(=C{{-1.502,2.535,0.641}}(S{{0.064,3.279,0.624}}1)C{{-1.656,1.155,0.89}}1=C{{-0.533,0.341,1.131}}C{{-0.744,-1.024,1.372}}=N{{-2,-1.529,1.371}}C{{-3.072,-0.739,1.139}}(=N{{-2.887,0.58,0.903}}1)N{{-4.304,-1.255,1.142}})C{{-3.959,3.173,0.327}})C{{-0.003,6.032,0.107}} 2-Anilino-4-_thiazol-5-yl_pyrimidine"));

  const int matoms = m1.natoms();

  const double initial_ave_dist = AverageInterAtomicDistance(m1);

  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<float> translate(-10.0, 10.0);
  std::uniform_real_distribution<float> u01(0.0, 1.0);
  std::uniform_real_distribution<float> uhalfpi(0.0f, M_PI * 0.50);

  int constexpr kNexpts = 10;
  for (int i = 0; i < kNexpts; ++i) {
    Molecule m2(m1);

    float cx = translate(rng);
    float cy = translate(rng);
    float cz = translate(rng);
    Coordinates rotation_axis;
    rotation_axis.setxyz(cx, cy, cz);
    rotation_axis.normalise();
  
    const angle_t theta = uhalfpi(rng);
  
    m2.rotate_atoms(rotation_axis, theta);

    float dx = translate(rng);
    float dy = translate(rng);
    float dz = translate(rng);
    m2.translate_atoms(dx, dy, dz);
    EXPECT_NEAR(AverageInterAtomicDistance(m2), initial_ave_dist, 1.0e-05);

    std::vector<Vector3d> to(matoms);
    for (int j = 0; j < matoms; ++j) {
      const Atom& a1 = m1[j];
      to[j](0) = a1.x();
      to[j](1) = a1.y();
      to[j](2) = a1.z();
    }

    std::vector<Vector3d> from(matoms);
    for (int j = 0; j < matoms; ++j) {
      const Atom& a2 = m2[j];
      from[j](0) = a2.x();
      from[j](1) = a2.y();
      from[j](2) = a2.z();
    }

    Eigen::Affine3d tx;
    ASSERT_TRUE(kabsch::kabsch(from, to, tx)) << "kabsch failed " << i;
#ifdef DEBUG_KABSCH_TEST
    std::cerr << "tx.transate\n" << tx.translation() << '\n';
    std::cerr << "tx.linear\n" << tx.linear() << '\n';
#endif

#define CHECK_TRANSFORM
#ifdef CHECK_TRANSFORM
    double residual_rms = 0.0;
    for (int j = 0; j < matoms; j++) {
      const Atom& a2 = m2[j];
      Eigen::Vector3d f(a2.x(), a2.y(), a2.z());
      auto rotated = tx * f;
      m2.setxyz(j, rotated(0), rotated(1), rotated(2));

      const Atom& a1 = m1[j];
      Eigen::Vector3d t(a1.x(), a1.y(), a1.z());
      double err = (to[j] - tx * from[j]).norm();
      residual_rms += (err * err);
    }
    residual_rms /= matoms;
    residual_rms = sqrt(residual_rms);
#ifdef DEBUG_KABSCH_TEST
    std::cerr << "Translated back\n";
    m2.debug_print(std::cerr);
    std::cerr << "RMS " << residual_rms << '\n';
#endif

    EXPECT_NEAR(residual_rms, 0.0, 1.0e-4);
#endif
  }
}

TEST(TestKabsch, TestManyRotationsV2) {
  // Acknowledge that random numbers in tests are not a great idea.
  using std::cerr;

  Molecule m1;
  ASSERT_TRUE(m1.build_from_smiles("C{{-0.73,4.777,0.282}}1(=N{{-2.064,4.703,0.191}}C{{-2.521,3.46,0.387}}(=C{{-1.502,2.535,0.641}}(S{{0.064,3.279,0.624}}1)C{{-1.656,1.155,0.89}}1=C{{-0.533,0.341,1.131}}C{{-0.744,-1.024,1.372}}=N{{-2,-1.529,1.371}}C{{-3.072,-0.739,1.139}}(=N{{-2.887,0.58,0.903}}1)N{{-4.304,-1.255,1.142}})C{{-3.959,3.173,0.327}})C{{-0.003,6.032,0.107}} 2-Anilino-4-_thiazol-5-yl_pyrimidine"));

  const int matoms = m1.natoms();

  const double initial_ave_dist = AverageInterAtomicDistance(m1);

  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<float> translate(-10.0, 10.0);
  std::uniform_real_distribution<float> u01(0.0, 1.0);
  std::uniform_real_distribution<float> uhalfpi(0.0f, M_PI * 0.50);

  int constexpr kNexpts = 10;
  for (int i = 0; i < kNexpts; ++i) {
    Molecule m2(m1);

    float cx = translate(rng);
    float cy = translate(rng);
    float cz = translate(rng);
    Coordinates rotation_axis;
    rotation_axis.setxyz(cx, cy, cz);
    rotation_axis.normalise();
  
    const angle_t theta = uhalfpi(rng);
  
    m2.rotate_atoms(rotation_axis, theta);

    float dx = translate(rng);
    float dy = translate(rng);
    float dz = translate(rng);
    m2.translate_atoms(dx, dy, dz);
    EXPECT_NEAR(AverageInterAtomicDistance(m2), initial_ave_dist, 1.0e-05);

    Eigen::MatrixXd from(3, matoms), to(3, matoms);
    for (int j = 0; j < matoms; ++j) {
      const Atom& a1 = m1[j];
      to(0, j) = a1.x();
      to(1, j) = a1.y();
      to(2, j) = a1.z();
    }

    for (int j = 0; j < matoms; ++j) {
      const Atom& a2 = m2[j];
      from(0, j) = a2.x();
      from(1, j) = a2.y();
      from(2, j) = a2.z();
    }

    Eigen::Affine3d tx;
    ASSERT_TRUE(kabsch::kabsch(from, to, tx)) << "kabsch failed " << i;
#ifdef DEBUG_KABSCH_TEST
    std::cerr << "tx.transate\n" << tx.translation() << '\n';
    std::cerr << "tx.linear\n" << tx.linear() << '\n';
#endif

#define CHECK_TRANSFORM
#ifdef CHECK_TRANSFORM
    double residual_rms = 0.0;
    for (int j = 0; j < matoms; j++) {
      const Atom& a2 = m2[j];
      Eigen::Vector3d f(a2.x(), a2.y(), a2.z());
      auto rotated = tx * f;
      m2.setxyz(j, rotated(0), rotated(1), rotated(2));

      const Atom& a1 = m1[j];
      Eigen::Vector3d t(a1.x(), a1.y(), a1.z());
      double err = (t - tx * f).norm();
      residual_rms += (err * err);
    }
    residual_rms /= matoms;
    residual_rms = sqrt(residual_rms);
#ifdef DEBUG_KABSCH_TEST
    std::cerr << "Translated back\n";
    m2.debug_print(std::cerr);
    std::cerr << "RMS " << residual_rms << '\n';
#endif

    EXPECT_NEAR(residual_rms, 0.0, 1.0e-4);
#endif
  }
}

}  // namespace
