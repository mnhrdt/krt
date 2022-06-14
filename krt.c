#include <stdlib.h> // malloc, free

static float getpixel(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j*w+i];
}

static float *build_kernel_from_string(char *s, int *w, int *h)
{
	(void)s;
	*w = 7;
	*h = 7;
	int wh = *w * *h;
	float *k = malloc(wh * sizeof*k);
	for (int i = 0; i < wh; i++)
		k[i] = 1.0/wh;
	return k;
}

void kernel_rank_transform_bruteforce(
		char *kernel_string,   // textual description of the kernel
		float *y,              // output image data
		float *x,              // input image data
		int w,                 // image width
		int h                  // image height
		)
{
	int W, H; // kernel width, height
	float *k = build_kernel_from_string(kernel_string, &W, &H);

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float xij = getpixel(x, w, h, i, j);
		double a = 0;
		for (int q = 0; q < H; q++)
		for (int p = 0; p < W; p++)
		{
			// TODO: verify signs, symmetry, center, etc
			float xpq = getpixel(x, w, h, i+p-W/2, j+q-H/2);
			float kpq = getpixel(k, W, H, p, q);
			if (xpq > xij)
				a += kpq;
		}
		y[j*w + i] = a;
	}

	free(k);
}

void kernel_rank_transform_bruteforce_split(
		char *kernel_string,   // textual description of the kernel
		float *y,              // output image data
		float *x,              // input image data
		int w,                 // image width
		int h,                 // image height
		int pd                 // pixel dimension
		)
{
	for (int i = 0; i < pd; i++)
		kernel_rank_transform_bruteforce(kernel_string,
				y + i*w*h,
				x + i*w*h,
				w, h);
}

#define KRT_MAIN

#ifdef KRT_MAIN
#include <stdio.h>
#include "iio.h"
int main(int c, char *v[])
{
	// process input arguments
	if (c != 2 && c != 3 && c != 4)
		return fprintf(stderr,
			"usage:\n\t%s KERNEL [in.png [out.png]]\n", *v);
			//          0 1       2       3
	char *kernel_string = v[1];
	char *filename_in   = c > 2 ? v[2] : "-";
	char *filename_out  = c > 3 ? v[3] : "-";

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);

	// allocate space for output image
	float *y = malloc(w * h * pd * sizeof*y);

	// run the algorithm
	kernel_rank_transform_bruteforce_split(kernel_string, y, x, w, h, pd);

	// write result and exit
	iio_write_image_float_split(filename_out, y, w, h, pd);
}
#endif//KRT_MAIN
