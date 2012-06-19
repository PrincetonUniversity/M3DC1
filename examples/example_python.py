import fio_py

filename = "/Users/ferraro/data/DIII-D/126006/mesh21a_kap6_amu6_n=1/C1.h5"

isrc = fio_py.open_source(fio_py.FIO_M3DC1_SOURCE,filename)

fio_py.get_options(isrc)
imag = fio_py.get_field(isrc, fio_py.FIO_MAGNETIC_FIELD)
ipres = fio_py.get_field(isrc, fio_py.FIO_TOTAL_PRESSURE)

x = (1.6, 0., 0.)
(br, bphi, bz) = fio_py.eval_vector_field(imag, x)
p = fio_py.eval_scalar_field(ipres, x)

print 'At x = ', x
print ' magnetic field = ', (br, bphi, bz), 'T'
print ' pressure = ', p, ' Pa'

fio_py.close_field(imag)
fio_py.close_field(ipres)
fio_py.close_source(isrc)


