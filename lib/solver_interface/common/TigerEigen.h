#ifndef TIGEREIGEN_H
#define TIGEREIGEN_H
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"
typedef Eigen::VectorXd VX;
typedef Eigen::RowVectorXd rVX;
typedef Eigen::MatrixXd MX;
typedef Eigen::VectorXi VXi;
typedef Eigen::Map<VX> MapV;
typedef Eigen::Map<MX> MapM;
typedef Eigen::Map<VXi> MapVi;
typedef Eigen::Map<const MX> cMapM;
typedef Eigen::Map<const VX> cMapV;
typedef const Eigen::Ref<const VX> cRefV;
typedef const Eigen::Ref<const MX> cRefM;
typedef const Eigen::Ref<const VXi> cRefVi;
typedef Eigen::Ref<VX> RefV;
typedef Eigen::Ref<MX> RefM;
typedef Eigen::Ref<VXi> RefVi;
//For row matrices
typedef Eigen::Matrix<double, -1, -1, Eigen::RowMajor> rMX;
typedef Eigen::Map<rVX> rMapV;
typedef Eigen::Map<rMX> rMapM;
typedef const Eigen::Ref<const rVX> crRefV;
typedef const Eigen::Ref<const rMX> crRefM;
typedef Eigen::Ref<rVX> rRefV;
typedef Eigen::Ref<rMX> rRefM;
//For sparse matrices
typedef Eigen::SparseMatrix<double> SpMX; // declares a column-major sparse matrix type of double
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> rSpMX; // declares a row-major sparse matrix type of double
typedef Eigen::Ref<SpMX> RefSpMX;
#endif
