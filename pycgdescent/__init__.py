__copyright__ = "Copyright (C) 2020 Alexandru Fikl"

__license__ = """
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

try:
    # python >=3.8 only
    from importlib import metadata
except ImportError:
    import importlib_metadata as metadata

__version__ = metadata.version("pycgdescent")

from _cg_descent import (
        cg_stats,
        cg_parameter,
        cg_default,
        cg_descent)


def cg_stats_repr(self):
    return (f"cg_stats<f={self.f:.16e}, gnorm={self.gnorm:.16e}, "
            f"iter={self.iter}, IterSub={self.IterSub}, NumSum={self.NumSub}, "
            f"nfunc={self.nfunc}, ngrad={self.ngrad}>")


cg_stats.__repr__ = cg_stats_repr
