# Shared parameter loader; channel drives are configured in channels.py.

import os
import sys

import mcell as m

import shared


MODEL_PATH = os.path.dirname(os.path.abspath(__file__))


# Apply command-line/environment overrides before loading the BNGL parameters.
for parameter_name, value in shared.parameter_overrides.items():
    setattr(sys.modules[__name__], parameter_name, value)


def not_defined(parameter_name):
    return parameter_name not in globals()


bngl_params = m.bngl_utils.load_bngl_parameters(
    os.path.join(MODEL_PATH, 'model.bngl'), shared.parameter_overrides
)
for parameter_name, value in bngl_params.items():
    globals()[parameter_name] = value

# The BNGL loader also returns override values as floating-point numbers.
# Restore the explicitly typed runtime values (notably integer SEED and
# ITERATIONS) supplied by customization.py.
for parameter_name, value in shared.parameter_overrides.items():
    globals()[parameter_name] = value


if not_defined('ITERATIONS'):
    ITERATIONS = 10000000  # 1 s at TIME_STEP = 1e-7 s

if not_defined('TIME_STEP'):
    TIME_STEP = 1e-7

if not_defined('DUMP'):
    DUMP = False

if not_defined('EXPORT_DATA_MODEL'):
    EXPORT_DATA_MODEL = True

if not_defined('SEED'):
    SEED = 1
