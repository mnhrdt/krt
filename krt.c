#include <assert.h>
#include <math.h>   // fabs, ceil, exp
#include <stdio.h>  // sscanf, fprintf
#include <stdlib.h> // malloc, free

// pixel management                                                         {{{1

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

// construciton of the kernels                                              {{{1

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

static float kernel_square(float x, float y, float p)
{
	return fmax(fabs(x), fabs(y)) < p ? 1 : 0;
}

static float kernel_disk(float x, float y, float p)
{
	return hypot(x, y) < p ? 1 : 0;
}

static float kernel_gauss(float x, float y, float p)
{
	return exp(-(x*x + y*y)/(2*p*p));
}

static float kernel_laplace(float x, float y, float p)
{
	return exp(-hypot(x, y)/p);
}

static float kernel_cauchy(float x, float y, float p)
{
	return p/(p + pow(hypot(x,y),2));
}

static float kernel_land(float x, float y, float p)
{
	return 1/hypot(x, y);
}

static float kernel_riesz(float x, float y, float p)
{
	return pow(hypot(x, y), -2-p);
}

typedef float (*kernel_t)(float,float,float);

static float *fill_kernel_into_square(kernel_t f, float p, int d)
{
}

// input s: string describing the kernel, like "disk3"
// output n: number of data points in the kernel
// output k: array of size 3*n with [x,y,k(x,y)] for each data point
static float *build_kernel_from_string(char *s, int *n)
{
	float p; // kernel parameter

	if (1 == sscanf(s, "actualsquare%g", &p))
	{
		int d = p;
		fprintf(stderr, "actual square of side d = %d\n", d);
		if (d >= 1)
		{
			*n = d*d;
			float *k = malloc(3 * d * d * sizeof*k);
			int l = 0;
			for (int i = 0; i < d; i++)
			for (int j = 0; j < d; j++)
			{
				k[3*l+0] = i - d/2;
				k[3*l+1] = j - d/2;
				k[3*l+2] = 1.0;
				l += 1;
			}
			assert(l == *n);
			return k;
		}
	}

	if (1 == sscanf(s, "square%g", &p))
	{
		int d = 1 + 2*round(sqrt(3)*p);
		fprintf(stderr, "square of σ = %g, side=%d\n",
				(d-1)/sqrt(12), d);
		if (d >= 1 && d < 1000)
		{
			*n = d*d;
			float *k = malloc(3 * d * d * sizeof*k);
			int l = 0;
			for (int i = 0; i < d; i++)
			for (int j = 0; j < d; j++)
			{
				k[3*l+0] = i - d/2;
				k[3*l+1] = j - d/2;
				k[3*l+2] = 1.0;
				l += 1;
			}
			assert(l == *n);
			return k;
		}
	}

	if (1 == sscanf(s, "disk%g", &p))
	{
		int d = 2*ceil(2*fabs(p)) + 1;
		fprintf(stderr, "disk of p = %g (d = %d)\n", p, d);
		float *k = malloc(3 * d * d * sizeof*k);

		// fill-in disk kernel
		int l = 0;
		for (int j = 0; j < d; j++)
		for (int i = 0; i < d; i++)
		{
			int x = i - d/2;
			int y = j - d/2;
			if (hypot(x, y) >= p) continue;
			k[3*l + 0] = x;
			k[3*l + 1] = y;
			k[3*l + 2] = 1;
			l += 1;
		}
		*n = l;
		return k;
	}

	if (1 == sscanf(s, "gauss%g", &p))
	{
		int d = 2*ceil(2*fabs(p)) + 1; // TODO: add this as an option?
		fprintf(stderr, "gaussian of σ = %g (d = %d)\n", p, d);
		float *k = malloc(3 * d * d * sizeof*k);

		// fill-in gauss kernel
		int l = 0;
		for (int j = 0; j < d; j++)
		for (int i = 0; i < d; i++)
		{
			int x = i - d/2;
			int y = j - d/2;
			// TODO add this test back
			//if (hypot(x, y) > d/2.0) continue;
			k[3*l + 0] = x;
			k[3*l + 1] = y;
			k[3*l + 2] = exp(-(x*x + y*y)/(2*p*p));
			l += 1;
		}
		*n = l;
		return k;
	}

	if (1 == sscanf(s, "laplace%g", &p))
	{
		int d = 2*ceil(3*fabs(p)) + 1;
		fprintf(stderr, "laplacian of σ = %g (d = %d)\n", p, d);
		float *k = malloc(3 * d * d * sizeof*k);

		// fill-in laplace kernel
		int l = 0;
		for (int j = 0; j < d; j++)
		for (int i = 0; i < d; i++)
		{
			int x = i - d/2;
			int y = j - d/2;
			//TODO add this test back
			//if (hypot(x, y) > d/2.0) continue;
			k[3*l + 0] = x;
			k[3*l + 1] = y;
			k[3*l + 2] = exp(-hypot(x, y)/p);
			l += 1;
		}
		*n = l;
		return k;
	}

	if (1 == sscanf(s, "cauchy%g", &p))
	{
		int d = 2*ceil(16*sqrt(p)) + 1;
		fprintf(stderr, "cauchy of σ = %g (d = %d)\n", p, d);
		float *k = malloc(3 * d * d * sizeof*k);

		// fill-in cauchy kernel
		int l = 0;
		for (int j = 0; j < d; j++)
		for (int i = 0; i < d; i++)
		{
			int x = i - d/2;
			int y = j - d/2;
			//TODO add this test back
			//if (hypot(x, y) > d/2.0) continue;
			k[3*l + 0] = x;
			k[3*l + 1] = y;
			k[3*l + 2] = p/(p + pow(hypot(x, y),2) );
			l += 1;
		}
		*n = l;
		return k;
	}

	if (1 == sscanf(s, "landc%g", &p)) // cut land (parameter=discrete size)
	{
		int d = land_truncation_for_sigma(p);
		fprintf(stderr, "land of σ = %g (d = %d)\n", p, d);
		float *k = malloc(3 * d * d * sizeof*k);

		// fill-in land kernel
		int l = 0;
		for (int j = 0; j < d; j++)
		for (int i = 0; i < d; i++)
		{
			int x = i - d/2;
			int y = j - d/2;
			// TODO
			//if (hypot(x, y) > d/2.0) continue;
			k[3*l + 0] = x;
			k[3*l + 1] = y;
			k[3*l + 2] = 1/hypot(x, y); // inf at center
			l += 1;
		}
		*n = l;
		return k;
	}

	exit(fprintf(stderr, "unrecognized kernel string \"%s\"\n", s));
}


// comparison functions                                                     {{{1

static float discrete_heaviside(float x) // like Heaviside, but gives 0.5 at 0
{
	if (x < 0) return 0;
	if (x > 0) return 1;
	return 0.5;
}

static float heaviside_H0(float x)
{
	if (x < 0) return 0;
	if (x > 0) return 1;
	return 0;
}

static float heaviside_H1(float x)
{
	if (x < 0) return 0;
	if (x > 0) return 1;
	return 1;
}

// NOTE: the "gap" parameters are normalized so that for
// parameter=1 the slope of the function equals 1

static float gap_heaviside_parameter;
static float gap_heaviside(float x)
{
	if (2*x < -gap_heaviside_parameter) return 0;
	if (2*x >  gap_heaviside_parameter) return 1;
	if (gap_heaviside_parameter <= 0.0) return 0.5;
	return 0.5 * (2*x + gap_heaviside_parameter) / gap_heaviside_parameter;
}

static float logistic_h_parameter;
static float logistic_h(float x)
{
	return 1.0 / (1.0 + exp(-4*x/logistic_h_parameter));
}

static float arctan_h_parameter;
static float arctan_h(float x)
{
	return 0.5 + atan(M_PI*x/arctan_h_parameter) / M_PI;
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
	x *= sqrt(M_PI) / erf_h_parameter;
	erf = sqrt( 1.0 - exp(-x*x * (4.0/M_PI + a*x*x) / (1.0 + a*x*x)) );
	if( x < 0.0 ) erf = -erf;
	return 0.5 + 0.5 * erf;
}

static float (*get_heaviside_from_string(char *s))(float)
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

	if (1 == sscanf(s, "H%g", &p))
	{
		fprintf(stderr, "Heaviside H%g\n", p);
		if (p == 1) return heaviside_H1;
		else return heaviside_H0;
	}

	fprintf(stderr, "discrete Heaviside function\n");
	return discrete_heaviside;
}

//void kernel_rank_transform_bruteforce(
//		char *kernel_string,    // textual description of the kernel
//		char *heaviside_string, // textual description of Heaviside f.
//		float *v,               // output image data
//		float *u,               // input image data
//		int w,                  // image width
//		int h                   // image height
//		)
//{
//	int W, H; // kernel width, height
//	float *k = build_kernel_from_string(kernel_string, &W, &H);
//	float (*HH)(float) = get_heaviside_from_string(heaviside_string);
//
//	for (int i = 0; i < w*h; i++)
//		v[i] = 0;
//
//	for (int j = 0; j < h; j++)   // image line
//	for (int i = 0; i < w; i++)   // image column
//	for (int q = 0; q < H; q++)   // kernel line
//	for (int p = 0; p < W; p++)   // kernel column
//	if (p!=W/2 || q!=H/2)
//	{
//		float ux = pixel(u, w, h, i, j);
//		float uy = pixel(u, w, h, i + p - W/2, j + q - H/2);
//		float kxy = pixel(k, W, H, p, q);
//		v[j*w+i] += kxy * HH(ux - uy);
//		//v[j*w+i] += kxy * discrete_heaviside(ux - uy);
//		// TODO: verify that this is correct for asymmetric kernels
//	}
//
//	free(k);
//}

// krt                                                                      {{{1

void kernel_rank_transform(
		float *kappa,           // a 3xn array with kernel positions
		int n,                  //  and values
		float (*sigma)(float),  // comparison function
		float *v,               // output image data
		float *u,               // input image data
		int w,                  // image width
		int h                   // image height
		)
{
	for (int i = 0; i < w*h; i++)
		v[i] = 0;

	// Definition 2.
	for (int x1 = 0; x1 < w; x1++)   // image column
	for (int x2 = 0; x2 < h; x2++)   // image line
	for (int k = 0; k < n; k++)      // kernel pixel index
	{
		int y1 = x1 + kappa[3*k+0];
		int y2 = x2 + kappa[3*k+1];
		float kappa_yx = kappa[3*k+2];
		float ux = pixel(u, w, h, x1, x2);
		float uy = pixel(u, w, h, y1, y2);
		v[x1 + x2*w] += kappa_yx * sigma(ux - uy);
	}
}

void kernel_rank_transform_split(
		float *kappa,           // a 3xn array with kernel positions
		int n,                  //  and values
		float (*sigma)(float),  // comparison function
		float *v,               // output image data
		float *u,               // input image data
		int w,                  // image width
		int h,                  // image height
		int pd                  // pixel dimension
		)
{
	for (int i = 0; i < pd; i++)
		kernel_rank_transform(kappa, n, sigma,
				v + i*w*h, u + i*w*h, w, h);
}


// main                                                                     {{{1

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
"Usage: krt KERNEL [-h HEAVISIDE] in out\n"
"   or: krt KERNEL [-h HEAVISIDE] in > out\n"
"   or: cat in | krt KERNEL [-h HEAVISIDE] > out\n"
"\n"
"Kernels:\n"
// TODO: de-couple kernel type and kernel size in two separate optoins
" actualsquareN  discrete square of size N x N\n"
" diskR          discrete disk of radius R\n"
" squareX        square of size (2N+1) x (2N+1)\n"
" gaussX         gaussian kernel of sigma S\n"
" cauchyX        cauchy kernel of sigma S\n"
" laplaceX       laplace kernel of sigma S\n"
" landcX         truncated land kernel of sigma S\n"
" file.npy       read the kernel weights from an image file\n"
"\n"
"Heaviside:\n"
" h           Discrete Heaviside function\n"
" gapR        Heaviside function with a gap of radius R\n"
" logisticK   logistic function scaled by a factor K\n"
" arctanK     arc tangent function scaled by a factor K\n"
" erfK        error function scaled by a factor K\n"
"\n"
"Options:\n"
" --help      print help message\n"
" --version   print version number\n"
" -ox X       horizontal offset X pixels the center of the kernel\n"
" -oy X       vertical offset X pixels the center of the kernel\n"
" -nn         do not normalize the kernel weights\n"
" -kc         keep the center of the discrete kernel\n"
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

//int main_plot_heavisides(void)
//{
//	float p = 2;
//	gap_heaviside_parameter = p;
//	logistic_h_parameter = p;
//	arctan_h_parameter = p;
//	erf_h_parameter = p;
//	for (float x = -2; x <= 2; x += 0.01)
//	{
//		float y1 = discrete_heaviside(x);
//		float y2 = gap_heaviside(x);
//		float y3 = logistic_h(x);
//		float y4 = arctan_h(x);
//		float y5 = erf_h(x);
//		printf("%g %g %g %g %g %g\n", x, y1, y2, y3, y4, y5);
//	}
//}

int main(int c, char *v[])
{
	// process help options
	if (c == 2 && 0 == strcmp(v[1], "--help"   )) return 0*puts(help);
	if (c == 2 && 0 == strcmp(v[1], "--version")) return 0*puts(version);

	// extract_named_parameters
	char *heaviside_string = pick_option(&c, &v, "h", "h");
	int offset_x = atoi(pick_option(&c, &v, "ox", "0"));
	int offset_y = atoi(pick_option(&c, &v, "oy", "0"));
	bool normalize_kernel = !pick_option(&c, &v, "nn", NULL);
	bool remove_kernel_center = !pick_option(&c, &v, "kc", NULL);

	// process input arguments
	if (c != 2 && c != 3 && c != 4)
		return fprintf(stderr,
			"usage:\n\t%s KERNEL [-h HEAVISIDE] [in [out]]\n", *v);
			//          0 1                      2   3
	char *kernel_string = v[1];
	char *filename_in   = c > 2 ? v[2] : "-";
	char *filename_out  = c > 3 ? v[3] : "-";
	int n_kappa;
	float *kappa = build_kernel_from_string(kernel_string, &n_kappa);
	float (*sigma)(float) = get_heaviside_from_string(heaviside_string);

	// apply kernel offsets in-place
	for (int i = 0; i < n_kappa; i++)
	{
		kappa[3*i + 0] += offset_x;
		kappa[3*i + 1] += offset_y;
	}

	// remove kernel centers in-place
	if (remove_kernel_center)
	{
		for (int i = 0; i < n_kappa; i++)
			if (kappa[3*i+0]==0 && kappa[3*i+1]==0)
				kappa[3*i+2] = 0;
	}

	// normalize kernel in-place
	if (normalize_kernel)
	{
		float K = 0;
		for (int i = 0; i < n_kappa; i++)
			K += kappa[3*i+2];
		for (int i = 0; i < n_kappa; i++)
			kappa[3*i+2] /= K;
	}

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);

	// allocate space for output image
	float *y = malloc(w * h * pd * sizeof*y);

	// run the algorithm
	kernel_rank_transform_split(kappa, n_kappa, sigma, y, x, w, h, pd);

	// write result and exit
	setenv("IIO_REM", version, 0);
	iio_write_image_float_split(filename_out, y, w, h, pd);
}
#endif//KRT_MAIN

// vim:set foldmethod=marker:
