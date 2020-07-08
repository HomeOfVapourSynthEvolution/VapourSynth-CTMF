Description
================

Constant-time median filtering.

Ported from https://nomis80.org/ctmf.html


Usage
=====

    ctmf.CTMF(clip clip[, int radius=2, int memsize=1048576, int[] planes=[0, 1, 2], int opt=0])

* clip: Clip to process. Any planar format with integer sample type of 8, 10, 12, 14, and 16 bit depth is supported. Float sample type of 32 bit depth is only supported for `radius=2`.

* radius: Median filter radius. The kernel will be a 2\*radius+1 by 2\*radius+1 square. Must be between 1 and 127 (inclusive).

* memsize: Maximum amount of memory to use, in bytes. Set this to the size of the L2 or L3 cache, then vary it slightly and measure the processing time to find the optimal value. For example, a 512 KB L2 cache would have memsize=512*1024 initially.

* planes: Sets which planes will be processed. Any unprocessed planes will be simply copied.

* opt: Sets which cpu optimizations to use.
  * 0 = auto detect
  * 1 = use c
  * 2 = use sse2
  * 3 = use avx2
  * 4 = use avx512


Compilation
===========

```
meson build
ninja -C build
ninja -C build install
```
