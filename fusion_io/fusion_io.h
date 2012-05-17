#ifndef FUSION_IO_H
#define FUSION_IO_H

// Scalar fields
#define FIO_ALPHA     1
#define FIO_DENSITY   2
#define FIO_PRESSURE  3

// Vector fields
#define FIO_MAGNETIC_FIELD  101
#define FIO_VELOCITY        102
#define FIO_ELECTRIC_FIELD  103

// Error codes
#define FIO_SUCCESS         0
#define FIO_UNSUPPORTED     10001
#define FIO_OUT_OF_BOUNDS   10002
#define FIO_FILE_ERROR      10003

#endif
