from distutils.core import setup, Extension, os

root = os.environ['M3DC1_ROOT']
arch = os.environ['M3DC1_ARCH']
srcdir= os.environ['SRCDIR']
copt = os.environ['CFLAGS']

fio_module = Extension('fio_py',
                       sources = [srcdir+'/python_interface.cpp'],
                       include_dirs = [root+'/m3dc1_lib'],
                       libraries = ['fusionio', 'm3dc1'],
                       library_dirs = [root+'/m3dc1_lib/_'+arch,
                                       root+'/fusion_io/_'+arch],
                       extra_linker_args = [copt])

setup (name = 'fio_py',
       author = 'N.M. Ferraro', 
       version = '1.0',
       description = 'Fusion IO Package',
       ext_modules = [fio_module])
