"""
Solve one-dimensional bounce problem with shooting method
==========================================================

export LD_PRELOAD=/usr/lib/libginac.so.2.1.0:/usr/lib/x86_64-linux-gnu/libboost_log.so.1.61.0:/usr/lib/x86_64-linux-gnu/libnlopt.so.0.8.1:/usr/lib/x86_64-linux-gnu/libalglib.so.3.8.0

"""

import ctypes
import os


PATH = os.path.dirname(os.path.abspath(__file__))
LIB_NAME = '{}/../lib/libbubbler.so'.format(PATH)

lib = ctypes.CDLL(LIB_NAME)
c_action = lib.action
c_action.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.Structure]
c_action.restype = ctypes.c_double


class Settings(ctypes.Structure):
    """
    Struct for Shooting::default.
    """
    _fields_ = [('shoot_bisect_bits', ctypes.c_int),
                ('action_arrived_rel', ctypes.c_double),
                ('shoot_ode_abs', ctypes.c_double),
                ('shoot_ode_rel', ctypes.c_double),
                ('action_ode_abs', ctypes.c_double),
                ('action_ode_rel', ctypes.c_double),
                ('drho_frac', ctypes.c_double),
                ('evolve_change_rel', ctypes.c_double),
                ('evolve_rho_min', ctypes.c_double),
                ('bisect_lambda_max', ctypes.c_double),
                ('iter_max', ctypes.c_int),
                ('periods_max', ctypes.c_double),
                ('f_y_max', ctypes.c_double),
                ('f_y_min', ctypes.c_double)]

# Make struct of default arguments

default = Settings()
default.shoot_bisect_bits = 5
default.action_arrived_rel = 1.e-3
default.shoot_ode_abs = 1.e-4
default.shoot_ode_rel = 1.e-4
default.action_ode_abs = 1.e-6
default.action_ode_rel = 1.e-6
default.drho_frac = 1.e-3
default.evolve_change_rel = 1.e-2
default.evolve_rho_min = 1.e-3
default.bisect_lambda_max = 5
default.iter_max = 100000
default.periods_max = 1.e2
default.f_y_max = 1.e6
default.f_y_min = 1.e-3


def make_struct(settings):
  """
  Populate struct for Shooting::default from a dictionary.
  """
  struct = default
  
  for k, v in settings.iteritems():
      struct.__setattr__(k, v)

  return struct

def action(E, alpha, dim=3, **kwargs):
    """
    Wrapper for C function that calculates action.
    """
    return c_action(E, alpha, dim, make_struct(kwargs))
    

if __name__ == "__main__":
    
    E = 1.
    alpha = 0.65
    print "alpha = {}".format(alpha)
    print "E = {}".format(E)
    print "action = {}".format(action(E, alpha))
