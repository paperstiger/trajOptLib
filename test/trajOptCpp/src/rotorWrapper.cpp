/*
 * solverwrapper.cpp
 * Copyright (C) 2017 Gao <gao.tang@duke.edu>
 *
 * Distributed under terms of the  license.
 */

/* Write a python wrapper for the problem */

#include "pybind11/eigen.h"
#include "pybind11/functional.h"
#include "TigerTools/TigerEigen.h"


class Rotor{
private:
    int dimx = 12, dimu = 4;
public:
    double m = 0.5, g = 9.81, kF = 1., kM = 0.0245, L = 0.175;
    VX inertia;
    double cg0[204] = {0.0};
    double *In;
    Rotor(){ 
        inertia.resize(3); 
        inertia << 0.0023, 0.0023, 0.004; 
        In = inertia.data();};
    void dyn(const double t, cRefV x, cRefV u, RefV f, RefM df){
        //First we call dronedyn to flush cg0
        dronedyn(t, x.data(), u.data());
        df.setZero();
        MapV dx(cg0, dimx);
        f = dx;
        MapM J(cg0 + dimx, dimx, dimx + dimu);
        df.middleCols(1, dimx + dimu) = J;
    }
    void dronedyn(const double t, const double *x, const double *u){
        //extract variables
        double phi = x[3], theta = x[4], psi = x[5], xd = x[6], yd = x[7], zd = x[8], p = x[9], q = x[10], r = x[11];
        double t1 = cos(theta);
        double t2 = sin(theta);
        double t3 = p * t1 + r * t2;
        double t4 = sin(phi);
        double t5 = cos(phi);
        double t6 = 0.1e1 / t5;
        double t7 = t1 * r;
        double t8 = t2 * p;
        double t9 = t8 - t7;
        double t10 = t6 * t9;
        double t11 = cos(psi);
        double t12 = sin(psi);
        double t13 = t1 * t12;
        double t14 = t11 * t2;
        double t15 = (u[0] + u[1] + u[2] + u[3]) * kF;
        t11 = t11 * t1;
        t12 = t12 * t2;
        double t16 = -t11 * t4 + t12;
        double t17 = 0.1e1 / m;
        t5 = t17 * t5;
        double t18 = t5 * t1;
        double t19 = -In[1] + In[2];
        double t20 = q * t19;
        double t21 = In[0] - In[2];
        double t22 = p * t21;
        double t23 = L * kF * (u[0] - u[2]) + t22 * r;
        double t24 = In[0] - In[1];
        double t25 = p * t24;
        double t26 = (u[0] - u[1] + u[2] - u[3]) * kM + t25 * q;
        double t27 = pow(t6, 0.2e1);
        double t28 = t27 * pow(t4, 0.2e1) + 0.1e1;
        t7 = -t7 * t28 + t8 * t28;
        t8 = t6 * t3;
        t28 = 0.1e1 / In[1];
        double t29 = 0.1e1 / In[0];
        double t30 = 0.1e1 / In[2];
        double t31 = t18 * kF;
        double t32 = t17 * (t13 * t4 + t14);
        double t33 = t32 * t15;
        t32 = t32 * kF;
        double t34 = t17 * t16;
        double t35 = t34 * kF;
        double t36 = t28 * L * kF;
        double t37 = t29 * L * kF;
        double t38 = t30 * kM;
        double t39 = t1 * t6;
        t6 = t2 * t6;
        cg0[0] = xd;
        cg0[1] = yd;
        cg0[2] = zd;
        cg0[3] = t3;
        cg0[4] = t10 * t4 + q;
        cg0[5] = -t10;
        cg0[6] = t33;
        cg0[7] = t34 * t15;
        cg0[8] = t18 * t15 - g;
        cg0[9] = -t29 * (-L * kF * (u[1] - u[3]) + t20 * r);
        cg0[10] = -t28 * t23;
        cg0[11] = t30 * t26;
        cg0[52] = t7;
        cg0[53] = -t27 * t4 * t9;
        cg0[54] = t5 * t13 * t15;
        cg0[55] = -t5 * t11 * t15;
        cg0[56] = -t17 * t4 * t1 * t15;
        cg0[63] = -t9;
        cg0[64] = t8 * t4;
        cg0[65] = -t8;
        cg0[66] = t17 * (-t12 * t4 + t11) * t15;
        cg0[67] = t17 * (t14 * t4 + t13) * t15;
        cg0[68] = -t5 * t2 * t15;
        cg0[78] = -t17 * t16 * t15;
        cg0[79] = t33;
        cg0[84] = 1;
        cg0[97] = 1;
        cg0[110] = 1;
        cg0[123] = t1;
        cg0[124] = t6 * t4;
        cg0[125] = -t6;
        cg0[130] = -t28 * r * t21;
        cg0[131] = t30 * q * t24;
        cg0[136] = 1;
        cg0[141] = -t29 * r * t19;
        cg0[143] = t30 * t25;
        cg0[147] = t2;
        cg0[148] = -t39 * t4;
        cg0[149] = t39;
        cg0[153] = -t29 * t20;
        cg0[154] = -t28 * t22;
        cg0[162] = t32;
        cg0[163] = t35;
        cg0[164] = t31;
        cg0[166] = -t36;
        cg0[167] = t38;
        cg0[174] = t32;
        cg0[175] = t35;
        cg0[176] = t31;
        cg0[177] = t37;
        cg0[179] = -t38;
        cg0[186] = t32;
        cg0[187] = t35;
        cg0[188] = t31;
        cg0[190] = t36;
        cg0[191] = t38;
        cg0[198] = t32;
        cg0[199] = t35;
        cg0[200] = t31;
        cg0[201] = -t37;
        cg0[203] = -t38;
    }
};


namespace py = pybind11;
PYBIND11_MODULE(libRotor, m){
    py::class_<Rotor>(m, "Rotor")
        .def(py::init<>())
        .def_readwrite("m", &Rotor::m)
        .def_readwrite("g", &Rotor::g)
        .def_readwrite("kF", &Rotor::kF)
        .def_readwrite("kM", &Rotor::kM)
        .def_readwrite("L", &Rotor::L)
        .def_readwrite("inertia", &Rotor::inertia)
        .def("dyn", &Rotor::dyn);
}
