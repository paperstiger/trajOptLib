# import from other directories
import pyoptsolver
from pyoptsolver import OptProblem as probFun, OptResult as result
from pyoptsolver import OptProblem as SnoptProblem
if pyoptsolver.__with_snopt__:
    from pyoptsolver import SnoptConfig, SnoptSolver
    from pyoptsolver import SnoptConfig as snoptConfig, SnoptSolver as solver
else:
    from pyoptsolver import IpoptSolver as SnoptSolver
    from pyoptsolver import IpoptConfig as IpoptSolver
if pyoptsolver.__with_ipopt__:
    from pyoptsolver import IpoptSolver as ipSolver, IpoptConfig as ipOption
    from pyoptsolver import IpoptSolver as IpSolver, IpoptConfig as IpOption
else:
    from pyoptsolver import SnoptSolver as ipSolver, SnoptConfig as ipOption
    from pyoptsolver import SnoptSolver as IpSolver, SnoptConfig as IpOption

from .trajOptProblem import trajOptProblem
from .trajOptBase import system, daeSystem, baseFun, addX
from .trajOptBase import lqrObj, linearPointObj, nonLinearObj, nonLinearPointObj, linearObj
from .trajOptBase import nonLinearPointConstr, nonLinearConstr
from .trajOptBase import linearPointConstr, linearConstr
# the ipopt solver
# the collocation version
from .trajOptCollocationProblem import trajOptCollocProblem
# the multi-phase version
from .trajOptMultiPhaseCollocationProblem import NonLinearConnectConstr
from .trajOptMultiPhaseCollocationProblem import LinearConnectConstr
from .trajOptMultiPhaseCollocationProblem import TrajOptMultiPhaseCollocProblem
from .trajOptManifoldCollocationProblem import manifoldConstr, trajOptManifoldCollocProblem

from .classBuilder import systemWrapper, daeSystemWrapper, nonLinearPointConstrWrapper, blockIndex
from .systemBuilder import GeometricRobot, GeometricRobotSystem

from .oopInterface import AbstractSolver

# for upper case alias
from .trajOptProblem import trajOptProblem as TrajOptProblem
from .trajOptBase import system as System, daeSystem as DaeSystem, baseFun as BaseFun, addX as AddX
from .trajOptBase import lqrObj as LQRObj, linearPointObj as LinearPointObj, nonLinearObj as NonLinearObj, nonLinearPointObj as NonLinearPointObj, linearObj as LinearObj
from .trajOptBase import nonLinearPointConstr as NonLinearPointConstr, nonLinearConstr as NonLinearConstr
from .trajOptBase import linearPointConstr as LinearPointConstr, linearConstr as LinearConstr
from .utility import parseX, showSol, getInf, InfBuilder
# import from other directories
from .snoptWrapper import directSolve, inDirectSolve, gradSolve, inGradSolve, spGradSolve, inSpGradSolve
# the collocation version
from .trajOptCollocationProblem import trajOptCollocProblem as TrajOptCollocProblem
# the multi-phase version
from .trajOptManifoldCollocationProblem import manifoldConstr as ManifoldConstr, trajOptManifoldCollocProblem as TrajOptManifoldCollocProblem

from .utility import OneBuilder, ZeroBuilder
inf_ = InfBuilder()
one_ = OneBuilder()
zero_ = ZeroBuilder()
