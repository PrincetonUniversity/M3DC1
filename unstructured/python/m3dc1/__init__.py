# -*- coding: utf-8 -*-
"""
Initialization file for the M3D-C1 package

Author: Andreas Kleiner, Chris Smiet, Ralf Mackenbach
Date created: Oct 1 2019
Date updated: Jul 28 2020
"""
import fpy

#Stable modules
from m3dc1.eval_field           import eval_field
from m3dc1.plot_diagnostics     import plot_diagnostics
from m3dc1.plot_field           import plot_field
from m3dc1.plot_field_basic     import plot_field_basic
from m3dc1.plot_field_mesh      import plot_field_mesh
from m3dc1.plot_line            import plot_line
from m3dc1.plot_mesh            import plot_mesh
from m3dc1.plot_shape           import plot_shape
from m3dc1.plot_time_trace      import plot_time_trace
from m3dc1.time_trace_fast      import get_timetrace
from m3dc1.time_trace_fast      import plot_time_trace_fast
from m3dc1.time_trace_fast      import avg_time_trace
from m3dc1.time_trace_fast      import growth_rate
from m3dc1.time_trace_fast      import scan_n
from m3dc1.time_trace_fast      import plot_gamma_n
from m3dc1.time_trace_fast      import write_gamma_n
from m3dc1.time_trace_fast      import eval_growth_n
from m3dc1.time_trace_fast      import compare_gamma_n
from m3dc1.time_trace_fast      import omegastari

from m3dc1.unit_conv            import unit_conv
from m3dc1.compensate_renorm    import compensate_renorm
from m3dc1.plot_signal          import plot_signal
from m3dc1.flux_coordinates     import flux_coordinates
from m3dc1.flux_average         import flux_average
from m3dc1.plot_flux_average    import plot_flux_average
from m3dc1.eigenfunction        import eigenfunction
from m3dc1.eigenfunction        import mode_type
from m3dc1.extend_profile       import extend_profile
from m3dc1.mesh_size            import mesh_size
from m3dc1                      import read_h5

from m3dc1.extract_profiles     import extract_profiles
from m3dc1.extract_profiles     import convert_p

from m3dc1.gfile                import read_gfile
from m3dc1.gfile                import plot_gfile

#Modules in development
#from m3dc1.plot_vector_field    import plot_vector_field
