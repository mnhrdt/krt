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

// find the truncation size so that the cropped land kernel has σ at least s
static int land_truncation_for_sigma(float s)
{
	int n = 1;
	float ss;
	float m;
	do {
		n = n + 2;
		ss = 0;
		m = 0;
		for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		if (i-n/2 != 0 || j-n/2 != 0)
		{
			m += 1 / hypot(i - n/2, j - n/2);
			ss += (i - n/2)*(i - n/2) / hypot(i - n/2, j - n/2);
		}
	} while (sqrt(ss/m) < s);
	return n;
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
		//int n = 2*ceil(p) + 1;
		int n = land_truncation_for_sigma(p);
		fprintf(stderr, "land of p = %g (n = %d)\n", fabs(p), n);
		*w = *h = n;
		float *k = malloc(n * n * sizeof*k);

		// fill-in land kernel
		for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
			k[n*j + i] = 1/hypot(i-n/2, j-n/2); // inf at center

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

	if (1 == sscanf(s, "actualsquare%g", &p))
	{
		int n = p;
		fprintf(stderr, "actual square of side n = %d\n", n);
		if (n >= 1 && n < 1000)
		{
			*w = *h = n;
			float *k = malloc(n * n * sizeof*k);
			for (int i = 0; i < n*n; i++)
				k[i] = 1.0/(n*n - 1);
			return k;
		}
	}

	if (1 == sscanf(s, "square%g", &p))
	{
		int n = 1 + 2*round(sqrt(3)*p);
		fprintf(stderr, "square of σ = %g, side=%d\n",
				(n-1)/sqrt(12), n);
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

static float gap_heaviside_parameter;
static float gap_heaviside(float x)
{
	if (x < -gap_heaviside_parameter) return 0;
	if (x >  gap_heaviside_parameter) return 1;
	if (gap_heaviside_parameter <= 0.0) return 0.5;
	return 0.5 * (x + gap_heaviside_parameter) / gap_heaviside_parameter;
}

static float logistic_h_parameter;
static float logistic_h(float x)
{
	return 1.0 / (1.0 + exp(-x/logistic_h_parameter));
}

static float arctan_h_parameter;
static float arctan_h(float x)
{
	return 0.5 + atan(x/arctan_h_parameter) / M_PI;
}

static float erf_h_parameter;
static float erf_h(float x)
{
	/* error function erf approximation based on the formula given in
	   "A handy approximation for the error function and its inverse"
	   by Sergei Winitzki, February 6, 2008,
	*/
	float erf;
	float a = 8.0 / 3.0 / M_PI * (M_PI - 3.0) / (4.0 - M_PI);
	x /= erf_h_parameter;
	erf = sqrt( 1.0 - exp(-x*x * (4.0/M_PI + a*x*x) / (1.0 + a*x*x)) );
	if( x < 0.0 ) erf = -erf;
	return 0.5 + 0.5 * erf;
}

static float (*build_heaviside_from_string(char *s))(float)
{
	float p; // parameter

	if (1 == sscanf(s, "gap%g", &p))
	{
		fprintf(stderr, "Heaviside with gap of radius = %g\n", fabs(p));
		gap_heaviside_parameter = fabs(p);
		return gap_heaviside;
	}

	if (1 == sscanf(s, "logistic%g", &p))
	{
		fprintf(stderr, "logistic 'Heaviside' scaled by %g\n", fabs(p));
		logistic_h_parameter = fabs(p);
		return logistic_h;
	}

	if (1 == sscanf(s, "arctan%g", &p))
	{
		fprintf(stderr, "arctan 'Heaviside' scaled by %g\n", fabs(p));
		arctan_h_parameter = fabs(p);
		return arctan_h;
	}

	if (1 == sscanf(s, "erf%g", &p))
	{
		fprintf(stderr, "erf 'Heaviside' scaled by %g\n", fabs(p));
		erf_h_parameter = fabs(p);
		return erf_h;
	}

	fprintf(stderr, "discrete Heaviside function\n");
	return discrete_heaviside;
}

void kernel_rank_transform_bruteforce(
		char *kernel_string,    // textual description of the kernel
		char *heaviside_string, // textual description of Heaviside f.
		float *v,               // output image data
		float *u,               // input image data
		int w,                  // image width
		int h                   // image height
		)
{
	int W, H; // kernel width, height
	float *k = build_kernel_from_string(kernel_string, &W, &H);
	float (*HH)(float) = build_heaviside_from_string(heaviside_string);

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
		v[j*w+i] += kxy * HH(ux - uy);
		//v[j*w+i] += kxy * discrete_heaviside(ux - uy);
		// TODO: verify that this is correct for asymmetric kernels
	}

	free(k);
}

void kernel_rank_transform_bruteforce_split(
		char *kernel_string,    // textual description of the kernel
		char *heaviside_string, // textual description of Heaviside f.
		float *y,               // output image data
		float *x,               // input image data
		int w,                  // image width
		int h,                  // image height
		int pd                  // pixel dimension
		)
{
	for (int i = 0; i < pd; i++)
		kernel_rank_transform_bruteforce(kernel_string,
				heaviside_string,
				y + i*w*h,
				x + i*w*h,
				w, h);
}

#define KRT_MAIN

#ifdef KRT_MAIN
#include <stdio.h>  // puts, fprintf
#include <string.h> // strcmp
#include <stdlib.h> // setenv
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
"Usage: krt KERNEL HEAVISIDE in out\n"
"   or: krt KERNEL HEAVISIDE in > out\n"
"   or: cat in | krt KERNEL HEAVISIDE > out\n"
"\n"
"Kernels:\n"
" squareN     square of size (2N+1) x (2N+1)\n"
" diskR       discrete disk of radius R\n"
" gaussS      gaussian kernel of sigma S\n"
" file.npy    read the kernel weights from an image file\n"
"\n"
"Heaviside:\n"
" h           Discrete Heaviside function\n"
" gapR        Heaviside function with a gap of radius R\n"
" logisticK   logistic function scaled by a factor K\n"
" arctanK     arc tangent function scaled by a factor K\n"
" erfK        error function scaled by a factor K\n"
"\n"
"Examples:\n"
" krt square7 h i.png o.tif   Classic rank transform of 49-pixel neighborhood\n"
" krt gauss2.7 h i.png o.tif  Gaussian-weighted kernel rank transform\n"
"\n"
"Report bugs to <grompone@gmail.com>"
;

// @c pointer to original argc
// @v pointer to original argv
// @o option name (after hyphen)
// @d default value
static char *pick_option(int *c, char ***v, char *o, char *d)
{
	int argc = *c;
	char **argv = *v;
	int id = d ? 1 : 0;
	for (int i = 0; i < argc - id; i++)
		if (argv[i][0] == '-' && 0 == strcmp(argv[i]+1, o))
		{
			char *r = argv[i+id]+1-id;
			*c -= id+1;
			for (int j = i; j < argc - id; j++)
				(*v)[j] = (*v)[j+id+1];
			return r;
		}
	return d;
}

int main(int c, char *v[])
{
	// process help options
	if (c == 2 && 0 == strcmp(v[1], "--help"   )) return 0*puts(help);
	if (c == 2 && 0 == strcmp(v[1], "--version")) return 0*puts(version);

	// extract_named_parameters
	char *heaviside_string = pick_option(&c, &v, "h", "h");

	// process input arguments
	if (c != 2 && c != 3 && c != 4)
		return fprintf(stderr,
			"usage:\n\t%s KERNEL [-h HEAVISIDE] [in [out]]\n", *v);
			//          0 1                      2   3
	char *kernel_string = v[1];
	char *filename_in   = c > 2 ? v[2] : "-";
	char *filename_out  = c > 3 ? v[3] : "-";

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);

	// allocate space for output image
	float *y = malloc(w * h * pd * sizeof*y);

	// run the algorithm
	kernel_rank_transform_bruteforce_split(kernel_string,
			heaviside_string, y, x, w, h, pd);

	// write result and exit
	setenv("IIO_REM", version, 0);
	iio_write_image_float_split(filename_out, y, w, h, pd);
}
#endif//KRT_MAIN
