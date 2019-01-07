from __future__ import absolute_import
import importlib
import os
import pkgutil

# Import all submodules of this module
modules = vars()
for _, module, _ in pkgutil.iter_modules([os.path.dirname(__file__)]):
    modules[module] = importlib.import_module('.' + module, __name__)
    del module

# Clean up
del importlib, modules, os, pkgutil
