#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
testManifoldFirstOrder.py

Implement mposa's paper with first order thing.
I temporarily give up analytic gradient now and only subclass probFun for quick implementation


It turns out that Posa's approach does not work well for such a simple problem. I have to think more.
Is this problem simply too simple for this approach to handle? Or the guess is simply too bad.
I can test it by providing a guess from collocation approach without special care on vel/acce.
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
sys.path.append('../')
from trajOptLib.io import getOnOffArgs
from trajOptLib import daeSystem, trajOptCollocProblem
from trajOptLib import lqrObj, nonLinearObj
from trajOptLib import system
from trajOptLib import manifoldConstr, nonLinearPointConstr
from trajOptLib import snoptConfig, solver, probFun
from trajOptLib import trajOptCollocProblem
from trajOptLib.utility import showSol

from carCommon import OmniCar, FirstOrderOmniCar, CircleConstr, CollocCircleConstr


class penObj(nonLinearObj):
    """Penalty function."""
    def __init__(self, prob):
        self.prob = prob
        nonLinearObj.__init__(self, prob.nx, nG=0)
        self.weight_u = 1
        self.weight_lmd = 1

    def __callg__(self, x, y, ):
        psx = self.prob.parseX(x)
        obj = self.weight_u * np.sum(psx['U']) + self.weight_lmd * np.sum(psx['P'])
        y[0] = obj


class CarProb(probFun):
    """This class is an implementation of mPosa's approach on a simple car problem.
    It calculates no gradient and does not inherent from those complex classes.
    This is only used for prototype."""
    def __init__(self, sys, con, x0, xf, N, h):
        """We assume a fixed time stuff, a first order system.
        We use simplified version, we only optimize state (no derivative) at knots, impose dyn constr at collocation.
        We only impose manifold constraint on knots, we correct dyn constr at collocation points

        Parameters
        ----------
        sys : the system instance, it gives information on a few dimensions
        con : the manifold constraint
        x0 : initial state
        xf : final state
        N : discretization size
        h : grid size

        """
        self.N = N
        self.h = h
        self.x0 = x0
        self.xf = xf
        self.sys = sys
        self.dimx = sys.nx
        self.dimu = sys.nu
        self.dimp = sys.np

        self.con = con
        self.nc = con.nc  # it equals the correction gamma
        self.nc_man = con.nf
        # construct problem
        numSol = N * (sys.nx + sys.nu + sys.np) + (N - 1) * (self.dimp + self.nc)
        numF = 1 + (N - 1) * sys.nx + (N - 2) * self.nc_man + 2 * self.nc  # 2*dimx means we release initial and final one

        # possible update, make xc be on manifold, and let xm (from integration) be closest to the manifold
        # numF += (N - 1) * self.nc_man

        probFun.__init__(self, numSol, numF)

        # set bounds for them
        self.setBounds()

    def __callf__(self, x, y):
        psx = self.parseX(x)
        psf = self.parseF(y)
        X = psx['X']
        U = psx['U']
        P = psx['P']
        Lmd = psx['Lmd']
        Gamma = psx['Gamma']
        obj = psf['obj']
        dyn = psf['dyn']
        man_mid = psf['man_mid']
        man_acce_0 = psf['man_acce_0']
        man_acce_f = psf['man_acce_f']
        # calculate obj, it is lqr cost
        obj[0] = np.sum(U**2) + (np.sum(P**2))
        # impose dyn constr
        for i in range(self.N - 1):
            if i == 0:
                dx0 = self.sys.dyn(0, X[i], U[i], P[i])
            else:
                dx0 = dx1
            dx1 = self.sys.dyn(0, X[i+1], U[i+1], P[i+1])
            if i < self.N - 2:  # impose manifold constraints for knot points
                catX = np.concatenate((X[i + 1], dx1))
                self.con.__callf__(catX, man_mid[i])
            xmid = 0.5*(X[i]+X[i+1])+self.h/8*(dx0 - dx1)
            xmiddot = 1.5/self.h*(X[i+1] - X[i]) - 0.25*(dx0 + dx1)
            umid = (U[i] + U[i+1])/2
            dxm = self.sys.dyn(0, xmid, umid, Lmd[i])
            corr = self.con.__return_correction__(xmid, Gamma[i])
            dxm[:self.dimx/2] += corr
            dyn[i] = xmiddot - dxm
        # impose acce constr on
        dx0 = self.sys.dyn(0, X[0], U[0], P[0])
        dxf = self.sys.dyn(0, X[-1], U[-1], P[-1])
        self.con.__callf__(np.concatenate((X[0], dx0)), man_acce_0, acce=True)
        self.con.__callf__(np.concatenate((X[-1], dxf)), man_acce_f, acce=True)

    def setBounds(self):
        xlb = -1e20*np.ones(self.nx)
        xub = -xlb
        lb = np.zeros(self.nf)
        ub = np.zeros(self.nf)
        psxlb = self.parseX(xlb)
        psxub = self.parseX(xub)
        psxlb['X'][0, :self.dimx] = self.x0
        psxub['X'][0, :self.dimx] = self.x0
        psxlb['X'][-1, :self.dimx] = self.xf
        psxub['X'][-1, :self.dimx] = self.xf
        # psxub['Lmd'][:] = 0  # lmd should be negative
        # psxub['P'][:] = 0  # lmd should be negative
        self.lb = lb
        self.ub = ub
        self.xlb = xlb
        self.xub = xub

    def parseX(self, x):
        """Parse a long vector x into parts"""
        n0, n1 = 0, self.N * (self.dimx + self.dimu + self.dimp)
        XUP = np.reshape(x[:n1], (self.N, self.dimx+self.dimu+self.dimp))  # state part
        X = XUP[:, :self.dimx]
        U = XUP[:, self.dimx:self.dimx+self.dimu]
        P = XUP[:, self.dimx+self.dimu:self.dimx+self.dimu+self.dimp]
        n0 = n1
        n1 = n0 + (self.N - 1) * self.dimp  # support force
        Lmd = np.reshape(x[n0:n1], (self.N - 1, self.dimp))
        n0 = n1
        # Gamma term
        n1 = n0 + (self.N - 1) * self.nc
        Gamma = np.reshape(x[n0:n1], (self.N - 1, self.nc))
        assert n1 == self.nx
        return {'X': X, 'U': U, 'P': P, 'Lmd': Lmd, 'Gamma': Gamma}

    def parseF(self, f):
        """Parse f"""
        obj = f[0:1]
        n0 = 1
        n1 = n0 + self.dimx * (self.N - 1)
        dyn = np.reshape(f[n0:n1], (self.N - 1, self.dimx))  # dynamics at mid
        n0 = n1
        n1 = n0 + (self.N - 2) * self.nc_man
        man_mid = np.reshape(f[n0:n1], (self.N - 2, self.nc_man))
        n0 = n1
        n1 = n0 + self.nc  # the acce constr at t0
        man_acce_0 = f[n0:n1]
        n0 = n1
        n1 = n0 + self.nc
        man_acce_f = f[n0:n1]
        assert n1 == self.nf
        return {'obj': obj, 'dyn': dyn, 'man_mid': man_mid, 'man_acce_0': man_acce_0, 'man_acce_f': man_acce_f}

    def genTraj(self, x):
        """Generate the trajectory in the dense format, same with redundant one"""
        psx = self.parseX(x)
        X = psx['X']
        U = psx['U']
        P = psx['P']
        Lmd = psx['Lmd']
        Gamma = psx['Gamma']
        # create space for xmid
        nqdqddq = 3*self.dimx // 2
        nq = self.dimx // 2
        Xout = np.zeros((2*self.N - 1, nqdqddq + self.dimu + self.dimp))
        xvout = Xout[:, :self.dimx]
        aout = Xout[:, self.dimx: nqdqddq]
        uout = Xout[:, nqdqddq: nqdqddq + self.dimu]
        pout = Xout[:, nqdqddq + self.dimu: nqdqddq + self.dimu + self.dimp]
        # use cubic to calc Xmid and ddq
        for i in range(self.N - 1):
            if i == 0:
                dx0 = self.sys.dyn(0, X[i], U[i], P[i])
                xvout[0, :] = X[i]
                aout[0, :] = dx0[nq:]
                uout[0, :] = U[i]
                pout[0, :] = P[i]
            else:
                dx0 = dx1
            dx1 = self.sys.dyn(0, X[i+1], U[i+1], P[i+1])
            xvout[2*(i+1), :] = X[i+1]
            aout[2*(i+1), :] = dx1[nq:]
            uout[2*(i+1), :] = U[i+1]
            pout[2*(i+1), :] = P[i+1]
            xmid = 0.5*(X[i] + X[i+1]) + self.h / 8 * (dx0 - dx1)
            xmiddot = 1.5/self.h*(X[i+1] - X[i]) - 0.25*(dx0 + dx1)
            xvout[2*i+1] = xmid
            aout[2*i + 1] = xmiddot[nq:]
            uout[2*i + 1] = (U[i] + U[i + 1]) / 2
            pout[2*i + 1] = Lmd[i]
        return {'X': Xout[:, :nqdqddq], 'U': uout, 'P': pout, 'Gamma': Gamma}


def main():
    args = getOnOffArgs('test', 'debug', 'drake', 'theone', 'colloc', 'initcolloc')
    if args.debug:
        sys.path.append('/home/motion/GAO/ThirdPartyLibs/pycharm-debug.egg')
        import pydevd
        pydevd.settrace('10.194.148.47', port=10000, stdoutToServer=True, stderrToServer=True)
    if args.test:
        runTest(init=args.initcolloc)
    if args.drake:
        checkDrake()
    if args.colloc:
        runColloc()


def runInitColloc():
    """Solve colloc to get initial guess for ultimate approach."""
    pass


def runColloc(show=True, skipfile=True):
    """The idea is, use direct collocation approach to solve position-only problem, use it to initialize."""
    fnm = 'run_colloc_solution.npz'
    if not skipfile:
        if os.path.exists(fnm):
            tmp = np.load(fnm)
            return tmp['flag'], tmp['sol'].item()
    x0 = np.array([1, 0, 0, 0])
    xf = np.array([0, 1, 0, 0])
    sys = FirstOrderOmniCar()
    N = 10
    t0 = 0
    tf = 2
    prob = trajOptCollocProblem(sys, N, t0, tf)
    setCollocBound(prob)

    # lqr = lqrObj(R=np.ones(2), P=np.ones(1))
    lqr = lqrObj(R=np.ones(2), P=np.ones(1))
    prob.addLQRObj(lqr)
    # obj = penObj(prob)
    # prob.addNonLinearObj(obj)

    constr = CollocCircleConstr(with_vel=True, with_acce=False, order=1)
    prob.addConstr(constr, path=True)

    # we have a few choices here:
    # 1. defect both p and u, not impose acce
    # 2. not defect p or u, impose vel and acce. It turns out to be saw-type solution

    prob.preProcess(colloc_constr_is_on=False, defect_p=True, defect_u=True)
    # construct solver
    cfg = snoptConfig()
    cfg.printLevel = 1
    cfg.printFile = 'test.out'
    cfg.verifyLevel = 3
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    sol = prob.parseSol(rst.sol)
    if rst.flag == 1 and show:
        con = np.array([constr.fast_call(x) for x in sol['x']])
        fig, ax = plt.subplots()
        ax.plot(con[:, 0], label='pos')
        ax.plot(con[:, 1], label='vel')
        ax.plot(con[:, 2], label='acce')
        ax.legend()
        ax.set_title('constraint vio')
        showSol(sol)
    np.savez(fnm, flag=rst.flag, sol=sol)
    del slv
    del prob
    return rst.flag, sol


def checkDrake():
    """Build the same problem, but read the solution output from drake, we use it to check"""
    x0 = np.array([1, 0, 0, 0])
    xf = np.array([0, 1, 0, 0])
    N = 10
    h = 0.3
    prob = CarProb(OmniCar(), CircleConstr(), x0, xf, N, h)
    sol_dict = dict()
    with open('drake_sol.txt', 'r') as f:
        for line in f:
            # we maintain a few lists of names
            ind = line.index('(')
            nm = line[:ind]
            # get the number
            ind2 = line.index('=')
            number = float(line[ind2+1:])
            if nm not in sol_dict:
                sol_dict[nm] = [number]
            else:
                sol_dict[nm].append(number)
    sol_dict = {key: np.array(sol_dict[key]) for key in sol_dict.keys()}
    # construct that guy
    x0 = np.zeros(prob.nx)
    y = np.zeros(prob.nf)
    psx = prob.parseX(x0)
    psx['X'][:].flat = sol_dict['x']
    psx['U'][:].flat = sol_dict['u']
    psx['P'][:].flat = sol_dict['lambda']
    psx['Lmd'][:].flat = sol_dict['lambda_c']
    psx['Gamma'][:].flat = sol_dict['v_c']
    prob.__callf__(x0, y)
    rst = prob.parseF(y)
    fig, ax = plt.subplots()
    X = psx['X']
    ax.plot(X[:, 0], X[:, 1])
    plt.show()
    pass


def runTest(init=False, show=True):
    """Test the simple car problem with finite difference approximation."""
    x0 = np.array([1, 0, 0, 0])
    xf = np.array([0, 1, 0, 0])
    N = 10
    h = 2.0 / (N - 1)  # to make time meet
    prob = CarProb(OmniCar(), CircleConstr(), x0, xf, N, h)
    cfg = snoptConfig()
    cfg.printLevel = 1
    cfg.optTol = 1e-4
    cfg.feaTol = 1e-5
    cfg.printFile = 'first_order.out'
    slver = solver(prob, cfg)
    if init:
        x0 = np.zeros(prob.nx)
        psx = prob.parseX(x0)
        # solve for a solution
        flag, sol = runColloc(show=False, skipfile=False)
        psx['X'][:] = sol['x'][::2, :4]
        psx['U'][:] = sol['u'][::2]
        psx['P'][:] = sol['p'][::2]
        rst = slver.solveGuess(x0)
    else:
        rst = slver.solveRand()
    print(rst.flag)
    rstF = np.zeros(prob.nf)
    prob.__callf__(rst.sol, rstF)
    psx = prob.parseX(rst.sol)
    traj = prob.genTraj(rst.sol)
    if show:
        fig, ax = plt.subplots()
        for i in range(traj['X'].shape[1]):
            ax.plot(traj['X'][:, i], label='x%d' % i)
        ax.legend()
        fig, ax = plt.subplots()
        for i in range(traj['U'].shape[1]):
            ax.plot(traj['U'][:, i], label='u%d' % i)
        ax.legend()
        fig, ax = plt.subplots()
        for i in range(traj['P'].shape[1]):
            ax.plot(traj['P'][:, i], label='p%d' % i)
        ax.legend()
        fig, ax = plt.subplots()
        for i in range(traj['Gamma'].shape[1]):
            ax.plot(traj['Gamma'][:, i], label='gamma%d' % i)
        ax.legend()
        plt.show()
    np.savez('first_order_traj.npz', **traj)


def setCollocBound(prob):
    x0lb = np.array([1, 0, 0, 0, -1e20, -1e20, -1e20, -1e20])
    x0ub = np.array([1, 0, 0, 0, 1e20, 1e20, 1e20, 1e20])
    prob.x0bd = [x0lb, x0ub]
    xflb = np.array([0, 1, 0, 0, -1e20, -1e20, -1e20, -1e20])
    xfub = np.array([0, 1, 0, 0, 1e20, 1e20, 1e20, 1e20])
    prob.xfbd = [xflb, xfub]
    prob.xbd = [-1e20*np.ones(8), 1e20*np.ones(8)]
    prob.ubd = [-1e20*np.ones(2), 1e20*np.ones(2)]
    if prob.numP > 0:
        prob.pbd = [-1e20*np.ones(1), 1e20*np.ones(1)]


if __name__ == '__main__':
    main()
