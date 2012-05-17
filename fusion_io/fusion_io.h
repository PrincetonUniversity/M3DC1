#ifndef FUSION_IO_H
#define FUSION_IO_H


#define FIO_GEQDSK_SOURCE    1
#define FIO_M3DC1_SOURCE     2


#define FIO_ALPHA     1
#define FIO_DENSITY   2
#define FIO_PRESSURE  3


#define FIO_MAGNETIC_FIELD  101
#define FIO_VELOCITY        102
#define FIO_ELECTRIC_FIELD  103


#define FIO_SUCCESS         0
#define FIO_UNSUPPORTED     10001
#define FIO_OUT_OF_BOUNDS   10002
#define FIO_FILE_ERROR      10003


#define FIO_INT_OPT_START    0
#define FIO_TIMESLICE        1
#define FIO_PERTURBED_ONLY   2
#define FIO_INT_OPT_END      3
#define FIO_IS_INT_OPT(x)    ((x)>FIO_INT_OPT_START && (x)<FIO_INT_OPT_END)


#define FIO_REAL_OPT_START   100
#define FIO_LINEAR_SCALE     101
#define FIO_REAL_OPT_END     102
#define FIO_IS_REAL_OPT(x)   ((x)>FIO_REAL_OPT_START && (x)<FIO_REAL_OPT_END)


#define FIO_STR_OPT_START    10000
#define FIO_STR_OPT_END      10000
#define FIO_IS_STR_OPT(x)    ((x)>FIO_STR_OPT_START && (x)<FIO_STR_OPT_END)


#define FIO_INTEGER          1
#define FIO_DOUBLE           2
#define FIO_STRING           3

#endif
