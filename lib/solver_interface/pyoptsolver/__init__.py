__all__ = ['IpoptConfig', 'IpoptSolver', 'SnoptConfig', 'SnoptSolver', 'OptProblem', 'OptResult',
           'OptSolver', 'OptConfig']
from .pyoptsolvercpp import __with_snopt__, __with_ipopt__, __version__
if __with_snopt__:
    from .pyoptsolvercpp import IpoptConfig, IpoptSolver
else:
    print("Cannot import IpoptWrapper")
if __with_snopt__:
    from .pyoptsolvercpp import SnoptConfig, SnoptSolver
else:
    print("Cannot import SnoptWrapper")
from .pyoptsolvercpp import OptProblem, OptResult
from .pyoptsolver import OptSolver, OptConfig
