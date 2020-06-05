# import from other directories
__all__ = ['System', 'DaeSystem', 'BaseFun', 'AddX',
           'LqrObj', 'LinearPointObj', 'NonLinearObj', 'NonLinearPointObj', 'LinearObj',
           'NonLinearPointConstr', 'NonLinearConstr',
           'LinearPointConstr', 'LinearConstr',
           'TrajOptProblem', 'TrajOptCollocProblem',
           'NonLinearConnectConstr', 'LinearConnectConstr',
           'TrajOptMultiPhaseCollocProblem',
           'ManifoldConstr', 'TrajOptManifoldCollocProblem',
           'SystemWrapper', 'DaeSystemWrapper', 'NonLinearPointConstrWrapper', 'block_index',
           'GeometricRobot', 'GeometricRobotSystem',
           'parse_X', 'show_sol', 'get_inf',
           'OptSolver', 'OptConfig',
           '__version__']
__version__ = '1.0.0'

from pyoptsolver import OptProblem, OptSolver, OptConfig
# import basic elements
from .trajOptBase import System, DaeSystem, BaseFun, AddX
from .trajOptBase import LqrObj, LinearPointObj, NonLinearObj, NonLinearPointObj, LinearObj
from .trajOptBase import NonLinearPointConstr, NonLinearConstr
from .trajOptBase import LinearPointConstr, LinearConstr
# import basic solvers
from .trajOptProblem import TrajOptProblem
from .trajOptDisProblem import TrajOptDisProblem
from .trajOptCollocationProblem import TrajOptCollocProblem
# the multi-phase version
from .trajOptMultiPhaseCollocationProblem import NonLinearConnectConstr, LinearConnectConstr
from .trajOptMultiPhaseCollocationProblem import TrajOptMultiPhaseCollocProblem
from .trajOptManifoldCollocationProblem import ManifoldConstr, TrajOptManifoldCollocProblem

# those wrappers, might be useful but who knows
from .classBuilder import SystemWrapper, DaeSystemWrapper, NonLinearPointConstrWrapper, blockIndex as block_index
from .systemBuilder import GeometricRobot, GeometricRobotSystem

from . import io
from . import plot

# utilities functions
from .utility import parse_X, show_sol, get_inf, InfBuilder, OneBuilder, ZeroBuilder
inf_ = InfBuilder()
one_ = OneBuilder()
zero_ = ZeroBuilder()
# legacy code
