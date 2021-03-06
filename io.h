#pragma once

struct ioctx {
  size_t (*read)(struct ioctx *, void *, int);
  size_t (*write)(struct ioctx *, const void *, int);
  int (*seek)(struct ioctx *, const int);
  size_t (*size)(struct ioctx *);

    void (*destroy)(struct ioctx *);
};

struct ioctx *ioctx_from_file(const char *fn, int t);

