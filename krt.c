#include <stdlib.h> // malloc, free

static float pixel(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j*w+i];
}

static float *build_kernel_from_string(char *s, int *w, int *h)
{
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

	for (int i = 0; i < w*h; i++)
		y[i] = 0;

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int q = 0; q < H; q++)
	for (int p = 0; p < W; p++)
		if (pixel(x, w, h, i+p-W/2, j+q-H/2) > pixel(x, w, h, i, j))
			y[j*w+i] += pixel(k, W, H, p, q);

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
#include <stdio.h>  // puts, fprintf
#include <string.h> // strcmp
#include "iio.h"
char *version = "krt 1.0\n\nWritten by RGvG and EML";
char *help = ""
"Krt computes the kernel rank transform of an image\n"
"\n"
"The kernel rank transform is like the classical rank transform, but\n"
"using weighted neighborhoods, defined by a user-supplied kernel.\n"
"The output image has the same size as the input and its values are between\n"
"0 and 1.  If the pixels have several components, each component image is\n"
"treated independently.\n"
"\n"
"Usage: krt KERNEL in out\n"
"   or: krt KERNEL in > out\n"
"   or: cat in | krt KERNEL > out\n"
"\n"
"Kernels:\n"
" squareN     square of size (2N+1) x (2N+1)\n"
" diskR       discrete disk of radius R\n"
" gaussS      gaussian kernel of size S\n"
" file.npy    read the kernel weights from an image file\n"
"\n"
"Examples:\n"
" krt square7 i.png o.png    Classic rank transform of 49-pixel neighborhood\n"
" krt gauss2.7 i.png o.png   Gaussian-weighted kernel rank transform\n"
"\n"
"Report bugs to <grompone@gmail.com>"
;
int main(int c, char *v[])
{
	// process help options
	if (c == 2 && 0 == strcmp(v[1], "--help"   )) return 0*puts(help);
	if (c == 2 && 0 == strcmp(v[1], "--version")) return 0*puts(version);

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
