# store the legacy code
# for backward compacity, just import trajoptlib.legacy as trajOptLib
__all__ = ['probFun', 'solver', 'snoptConfig', 'ipSolver', 'ipOption',
           'trajOptProblem', 'system', 'daeSystem', 'baseFun', 'addX',
           'lqrObj', 'linearPointObj', 'nonLinearObj', 'nonLinearPointObj',
           'linearObj', 'nonLinearPointConstr', 'nonLinearConstr', 'linearPointConstr',
           'linearConstr',
           'trajOptProblem', 'trajOptCollocProblem',
           'NonLinearConnectConstr', 'LinearConnectConstr',
           'TrajOptMultiPhaseCollocProblem',
           'manifoldConstr', 'trajOptManifoldCollocProblem',
           'systemWrapper', 'daeSystemWrapper', 'nonLinearPointConstrWrapper',
           'blockIndex', 'GeometricRobot', 'GeometricRobotSystem',
           'parseX', 'showSol', 'getInf',
           'inf_', 'one_', 'zero_', '__version__']
import pyoptsolver
from pyoptsolver import OptProblem as probFun
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
from . import (System as system,
               DaeSystem as daeSystem,
               BaseFun as baseFun,
               AddX as addX)

from . import (LqrObj as lqrObj,
               LinearPointObj as linearPointObj,
               NonLinearObj as nonLinearObj,
               NonLinearPointObj as nonLinearPointObj,
               LinearObj as linearObj)
from . import (NonLinearPointConstr as nonLinearPointConstr,
               NonLinearConstr as nonLinearConstr,
               LinearPointConstr as linearPointConstr,
               LinearConstr as linearConstr)
from . import TrajOptProblem as trajOptProblem
from . import TrajOptCollocProblem as trajOptCollocProblem
from . import NonLinearConnectConstr, LinearConnectConstr
from . import TrajOptMultiPhaseCollocProblem
from . import (ManifoldConstr as manifoldConstr,
               TrajOptManifoldCollocProblem as trajOptManifoldCollocProblem)
from .classBuilder import (SystemWrapper as systemWrapper,
                           DaeSystemWrapper as daeSystemWrapper,
                           NonLinearPointConstrWrapper as nonLinearPointConstrWrapper,
                           blockIndex)
from .systemBuilder import GeometricRobot, GeometricRobotSystem
from .utility import parse_X as parseX, show_sol as showSol, get_inf as getInf
from . import inf_, one_, zero_
from . import __version__
