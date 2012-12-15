import fio_py

filename = "/u/ferraro/Share/data/145117/15_iterdb_ExB/C1.h5"

isrc = fio_py.open_source(fio_py.FIO_M3DC1_SOURCE,filename)

fio_py.get_options(isrc)
imag = fio_py.get_field(isrc, fio_py.FIO_VECTOR_POTENTIAL)
ipres = fio_py.get_field(isrc, fio_py.FIO_TOTAL_PRESSURE)
ivel = fio_py.get_field(isrc, fio_py.FIO_FLUID_VELOCITY)

fio_py.set_int_option(fio_py.FIO_SPECIES, fio_py.FIO_ELECTRON)
idens = fio_py.get_field(isrc, fio_py.FIO_DENSITY)

fio_py.set_int_option(fio_py.FIO_SPECIES, fio_py.FIO_MAIN_ION)
ipi = fio_py.get_field(isrc, fio_py.FIO_PRESSURE)

list = fio_py.get_available_fields(isrc)

print 'Available fields:'
for field in list:
    print ' ', fio_py.get_field_name(field)

x = (1.6, 0., 0.)
(br, bphi, bz) = fio_py.eval_vector_field(imag, x)
p = fio_py.eval_scalar_field(ipres, x)
pi = fio_py.eval_scalar_field(ipi, x)
ne = fio_py.eval_scalar_field(idens, x)
(vr, vphi, vz) = fio_py.eval_vector_field(ivel, x)

print 'At x = ', x
print ' magnetic field = ', (br, bphi, bz), ' T'
print ' fluid velocity = ', (vr, vphi, vz), ' m/s'
print ' toroidal rotation = ', vphi/x[0], ' rad / s' 
print ' pressure = ', p, ' Pa'
print ' ion pressure = ', pi, ' Pa'
print ' density = ', ne, ' '

fio_py.close_field(imag)
fio_py.close_field(ipres)
fio_py.close_source(isrc)


