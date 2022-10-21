#include <math.h>   // fabs, ceil, exp
#include <stdio.h>  // sscanf, fprintf
#include <stdlib.h> // malloc, free

static float pixel(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[j*w+i];
}

static float pixel_e(float *x, int w, int h, int i, int j)
{
	if (i < 0 || j < 0 || i >= w || j >= h)
	{
		fprintf(stderr, "bad pixel [%d,%d](%d,%d)!\n", w, h, i, j);
		return 0;
	}
	return x[j*w+i];
}

static float *build_kernel_from_string(char *s, int *w, int *h)
{
	float p; // kernel parameter

	if (1 == sscanf(s, "gauss%g", &p))
	{
		int n = 2*ceil(2*fabs(p)) + 1; // TODO: add this as an option?
		fprintf(stderr, "gaussian of σ = %g (n = %d)\n", fabs(p), n);
		*w = *h = n;
		float *k = malloc(n * n * sizeof*k);

		// fill-in gaussian kernel
		for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			float x = i - n/2;
			float y = j - n/2;
			k[n*j + i] = exp(-(x*x + y*y)/(2*p*p));
		}

		// set central pixel to 0 (cf. article)
		k[n*(n/2) + n/2] = 0;

		// normalize so that the sum of k is 1
		float K = 0;
		for (int i = 0; i < n*n; i++)
			K += k[i];
		for (int i = 0; i < n*n; i++)
			k[i] /= K;
		return k;
	}

	if (1 == sscanf(s, "laplace%g", &p))
	{
		int n = 2*ceil(3*fabs(p)) + 1; // TODO: add this as an option?
		fprintf(stderr, "laplacian of σ = %g (n = %d)\n", fabs(p), n);
		*w = *h = n;
		float *k = malloc(n * n * sizeof*k);

		// fill-in laplacian kernel
		for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			k[n*j + i] = exp(-hypot(i-n/2, j-n/2)/p);

		// set central pixel to 0 (cf. article)
		k[n*(n/2) + n/2] = 0;

		// normalize so that the sum of k is 1
		float K = 0;
		for (int i = 0; i < n*n; i++)
			K += k[i];
		for (int i = 0; i < n*n; i++)
			k[i] /= K;
		return k;
	}

	if (1 == sscanf(s, "cauchy%g", &p))
	{
		int n = 2*ceil(16*sqrt(p)) + 1;
		fprintf(stderr, "cauchy of σ = %g (n = %d)\n", fabs(p), n);
		*w = *h = n;
		float *k = malloc(n * n * sizeof*k);

		// fill-in cauchy kernel
		for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			k[n*j + i] = p/(p + pow(hypot(i-n/2, j-n/2),2) );

		// set central pixel to 0 (cf. article)
		k[n*(n/2) + n/2] = 0;

		// normalize so that the sum of k is 1
		float K = 0;
		for (int i = 0; i < n*n; i++)
			K += k[i];
		for (int i = 0; i < n*n; i++)
			k[i] /= K;
		return k;
	}

	if (1 == sscanf(s, "landc%g", &p)) // cut land (parameter=discrete size)
	{
		int n = 2*ceil(p) + 1;
		fprintf(stderr, "land of p = %g (n = %d)\n", fabs(p), n);
		*w = *h = n;
		float *k = malloc(n * n * sizeof*k);

		// fill-in land kernel
		for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			k[n*j + i] = 1/hypot(i-n/2, j-n/2);

		// set central pixel to 0 (cf. article)
		k[n*(n/2) + n/2] = 0;

		// normalize so that the sum of k is 1
		float K = 0;
		for (int i = 0; i < n*n; i++)
			K += k[i];
		for (int i = 0; i < n*n; i++)
			k[i] /= K;
		return k;
	}

	if (1 == sscanf(s, "square%g", &p))
	{
		int n = p;
		fprintf(stderr, "square of side n = %d\n", n);
		if (n >= 1 && n < 1000)
		{
			*w = *h = n;
			float *k = malloc(n * n * sizeof*k);
			for (int i = 0; i < n*n; i++)
				k[i] = 1.0/(n*n - 1);
			return k;
		}
	}

	if (1 == sscanf(s, "disk%g", &p))
	{
		;
	}

	*w = 7;
	*h = 7;
	int wh = *w * *h;
	float *k = malloc(wh * sizeof*k);
	for (int i = 0; i < wh; i++)
		k[i] = 1.0/(wh - 1);
	return k;
}

static float discrete_heaviside(float x) // like Heaviside, but gives 0.5 at 0
{
	if (x < 0) return 0;
	if (x > 0) return 1;
	return 0.5;
}

void kernel_rank_transform_bruteforce(
		char *kernel_string,  // textual description of the kernel
		float *v,             // output image data
		float *u,             // input image data
		int w,                // image width
		int h                 // image height
		)
{
	int W, H; // kernel width, height
	float *k = build_kernel_from_string(kernel_string, &W, &H);

	for (int i = 0; i < w*h; i++)
		v[i] = 0;

	for (int j = 0; j < h; j++)   // image line
	for (int i = 0; i < w; i++)   // image column
	for (int q = 0; q < H; q++)   // kernel line
	for (int p = 0; p < W; p++)   // kernel column
	if (p!=W/2 || q!=H/2)
	{
		float ux = pixel(u, w, h, i, j);
		float uy = pixel(u, w, h, i + p - W/2, j + q - H/2);
		float kxy = pixel(k, W, H, p, q);
		v[j*w+i] += kxy * discrete_heaviside(ux - uy);
		// TODO: verify that this is correct for asymmetric kernels
	}

	free(k);
}

void kernel_rank_transform_bruteforce_split(
		char *kernel_string,  // textual description of the kernel
		float *y,             // output image data
		float *x,             // input image data
		int w,                // image width
		int h,                // image height
		int pd                // pixel dimension
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
" gaussS      gaussian kernel of sigma S\n"
" file.npy    read the kernel weights from an image file\n"
"\n"
"Examples:\n"
" krt square7 i.png o.tiff   Classic rank transform of 49-pixel neighborhood\n"
" krt gauss2.7 i.png o.tiff  Gaussian-weighted kernel rank transform\n"
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
