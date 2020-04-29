/*
 * A python wrapper for sqopt solver
 */
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <string>
#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"

#include "sqoptEZ.h"


typedef Eigen::VectorXd VX;
typedef Eigen::VectorXi VXi;
typedef Eigen::Ref<VX> RefVX;
typedef Eigen::Ref<const VX> cRefVX;
typedef Eigen::Ref<const VXi> cRefVXi;


int PRINT_LEVEL = 0;


class PySqopt : public SqoptEZ {
public:
    PySqopt(int n, int m) : SqoptEZ(n, m) {}

    void set_H(cRefVX val, cRefVXi row, cRefVXi col) {
        setMatrixCSC('H', val.data(), row.data(), col.data());
    }

    void set_q_zero() {
        setLinearCostToZero();
    }

    void set_q(cRefVX q){
        setVector('c', q.data());
    }

    void set_A(cRefVX val, cRefVXi row, cRefVXi col) {
        setMatrixCSC('A', val.data(), row.data(), col.data());
    }

    void set_bound(cRefVX xlb, cRefVX xub, cRefVX Alb, cRefVX Aub) {
        VX lb(n + m - 1), ub(n + m - 1);
        lb << xlb, Alb;
        ub << xub, Aub;
        if (PRINT_LEVEL > 0)
            std::cout << "lb" << lb << " ub " << ub << std::endl;
        setVector('l', lb.data());
        setVector('u', ub.data());
    }

    int get_info() {
        return inform;
    }

    double get_obj() {
        return getObjective();
    }

    VX get_solution() {
        VX sol(n);
        getPrimalSolution(sol.data());
        return sol;
    }

    VX get_lambda() {
        VX lmd(m + n - 1);
        if (PRINT_LEVEL > 0)
            std::cout << lmd << std::endl;
        getDualSolution(lmd.data());
        if (PRINT_LEVEL > 0)
            std::cout << lmd << std::endl;
        return lmd;
    }

    VXi get_hs() {
        VXi hs(m + n - 1);
        getHsSolution(hs.data());
        return hs;
    }
};


namespace py = pybind11;
PYBIND11_MODULE(libsqopt, m) {
    m.doc() = R"pbdoc(
        A Python interface for SQOPT
    )pbdoc";

    m.def("set_print_level", [](int level) {PRINT_LEVEL = level;});

    py::class_<PySqopt>(m, "Sqopt", R"pbdoc(
        A class that provides a pipeline for QP problem solving.

        Args:
            n (int): dimension of variables to be optimized
            m (int): dimension of linear constraints
    )pbdoc")
    .def(py::init<int, int>())
    .def("set_H", &PySqopt::set_H, R"pbdoc(
            Set H by a CSC sparse matrix. It has to be full, not triangular half.

            Args:
                val (ndarray): the value part for CSC matrix
                row (ndarray): the row part for CSC matrix
                col (ndarray): the column part for CSC matrix
        )pbdoc")
    .def("set_q_zero", &PySqopt::set_q_zero, R"pbdoc(
            Set linear cost part to be zero.
        )pbdoc")
    .def("set_q", &PySqopt::set_q, R"pbdoc(
            Set linear cost part by a vector.

            Args:
                q (ndarray): the q vector
        )pbdoc")
    .def("set_A", &PySqopt::set_A, R"pbdoc(
            Set A by a CSC sparse matrix.

            Args:
                val (ndarray): the value part for CSC matrix
                row (ndarray): the row part for CSC matrix
                col (ndarray): the column part for CSC matrix
        )pbdoc")
    .def("set_bound", &PySqopt::set_bound, R"pbdoc(
            Set bounds on x itself and linear constraints.

            Args:
                xlb (ndarray): the lower bound of x
                xub (ndarray): the upper bound of x
                Alb (ndarray): the lower bound of linear constraints
                Aub (ndarray): the upper bound of linear constraints
        )pbdoc")
    .def("solve", &PySqopt::solve, R"pbdoc(
            Solve the problem using Cold or Warm start specified by input type.

            Args:
                type (int): 0 for cold start, 1 for basic, and 2 for warm start.
            Returns:
                type (int): the inform of calling solver. 1 means okay.
        )pbdoc")
    .def("get_info", &PySqopt::get_info, R"pbdoc(
            Get solution status of sqopt.

            Returns:
                int: the inform member.
        )pbdoc")
    .def("get_obj", &PySqopt::get_obj, R"pbdoc(
            Get the cost function.

            Returns:
                float: the cost function.
        )pbdoc")
    .def("get_solution", &PySqopt::get_solution, R"pbdoc(
            Get the primal solution.

            Returns:
                ndarray: the solution.
        )pbdoc")
    .def("get_lambda", &PySqopt::get_lambda, R"pbdoc(
            Get the dual lambda.

            Returns:
                ndarray: the Lagrangian multipliers.
        )pbdoc")
    .def_readonly("n", &PySqopt::n)
    .def_readonly("m", &PySqopt::m)
    ;
}
