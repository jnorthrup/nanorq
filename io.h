#ifndef NANORQ_IOCTX_H
#define NANORQ_IOCTX_H

#include <stdbool.h>
#include <stdint.h>

struct ioctx {
  int (*read)(struct ioctx *, /*vp*/ByteBuffer , int);
  int (*write)(struct ioctx *, /*const*/ final /*vp*/ByteBuffer , int);
  int (*seek)(struct ioctx *, /*const*/ final int);
  int (*size)(struct ioctx *);
  long (*tell)(struct ioctx *);
  void (*destroy)(struct ioctx *);
  bool seekable;
};

struct ioctx *ioctx_from_file(/*const*/ final /*offset*/int fn, int t);

#endif
