#pragma once

#include "image.hpp"
#include <fftw3.h>
#include <cstdlib>	// for debug
#include <cmath>

// macros for complex number
#define REAL 0
#define IMAG 1


struct ComplexPixel {
    std::complex<double> r, g, b;

    ComplexPixel() = default;

    ComplexPixel(uint8_t r_real, uint8_t g_real, uint8_t b_real): r(r_real), g(g_real), b(b_real) {}

    ComplexPixel(double real, double image): r(real, image), g(real, image), b(real, image) {}

    ComplexPixel(const std::complex<double> &r, const std::complex<double> &g, const std::complex<double> &b):
            r(r), g(g), b(b) {}

    [[nodiscard]] inline double real_sum() const {
        return r.real() + g.real() + b.real();
    }

    friend ComplexPixel operator + (const ComplexPixel &x, const ComplexPixel &y) {
        return ComplexPixel(x.r + y.r, x.g + y.g, x.b + y.b);
    }

    friend ComplexPixel operator - (const ComplexPixel &x, const ComplexPixel &y) {
        return ComplexPixel(x.r - y.r, x.g - y.g, x.b - y.b);
    }

    friend ComplexPixel operator * (const ComplexPixel &x, const ComplexPixel &y) {
        return ComplexPixel(x.r * y.r, x.g * y.g, x.b * y.b);
    }

    friend ComplexPixel operator / (const ComplexPixel &x, double k) {
        return ComplexPixel(x.r / k, x.g / k, x.b / k);
    }

    friend ComplexPixel operator * (const ComplexPixel &x, double k) {
        return ComplexPixel(x.r * k, x.g * k, x.b * k);
    }
};


/// Convert `Pixel` into `ComplexPixel`
[[nodiscard]] ComplexPixel to_complex_pixel(const Pixel &x) {
    return ComplexPixel(x.r, x.g, x.b);
}


/// Round length to 2^n
int dft_round(int x) {
    int len = 1;
    while (len < x) {
        len *= 2;
    }
    return len;
}


/// Get low-bit number
int dft_lowbit(int x) {
    return x & (-x);
}


/// Allocate DFT working space
void dft_alloc(const std::shared_ptr<Image> &image, int dft_w, int dft_h, ComplexPixel* &dft_space) {
    dft_space = static_cast<ComplexPixel*> (std::malloc(dft_w * dft_h * sizeof(ComplexPixel)));
    std::fill(dft_space, dft_space + dft_w * dft_h, ComplexPixel());
    for (int i = 0, index = 0; i < image->h; ++ i) {
        for (int j = 0; j < image->w; ++ j, ++ index) {
            dft_space[i * dft_w + j] = to_complex_pixel(image->data[index]);
        }
    }
}


/// Run DFT and IDFT
void dft(int dft_w, int dft_h, ComplexPixel* dft_space, bool inverse=false) {
    // Allocate space
    double coefficient = inverse ? -1 : 1;
    assert(dft_w > 0 and dft_h > 0 and dft_w == dft_lowbit(dft_w) and dft_h == dft_lowbit(dft_h));

    // Butterfly changes by w
    for (int row = 0; row < dft_h; ++ row) {
        ComplexPixel *base = dft_space + row * dft_w;
        for(int i = 0, j = 0; i < dft_w; ++ i){
            if (i > j) {
                std::swap(base[i], base[j]);
            }
            for (int t = dft_w / 2; (j ^= t) < t; t /= 2);
        }
    }

    // Butterfly changes by h
    for (int col = 0; col < dft_w; ++ col) {
        ComplexPixel *base = dft_space + col;
        for(int i = 0, j = 0; i < dft_h; ++ i){
            if (i > j) {
                std::swap(base[i * dft_w], base[j * dft_w]);
            }
            for (int t = dft_h / 2; (j ^= t) < t; t /= 2);
        }
    }

    // DFT by w
    for (int row = 0; row < dft_h; ++ row) {
        ComplexPixel *base = dft_space + row * dft_w;
        ComplexPixel wn, w, t, u;
        for(int m = 2; m <= dft_w; m *= 2){
            wn = ComplexPixel(cos(2.0 * M_PI / m), coefficient * sin(2.0 * M_PI / m));
            for (int i = 0 ; i < dft_w; i += m) {
                w = ComplexPixel(1, 0);
                for (int k = 0, p = m / 2; k < p; ++ k, w = w * wn) {
                    t = w * base[i + k + m / 2];
                    u = base[i + k];
                    base[i + k] = u + t;
                    base[i + k + m / 2] = u - t;
                }
            }
        }
    }

    // DFT by h
    for (int col = 0; col < dft_w; ++ col) {
        ComplexPixel *base = dft_space + col;
        ComplexPixel wn, w, t, u;
        for(int m = 2; m <= dft_h; m *= 2){
            wn = ComplexPixel(cos(2.0 * M_PI / m), coefficient * sin(2.0 * M_PI / m));
            for (int i = 0 ; i < dft_h; i += m) {
                w = ComplexPixel(1, 0);
                for (int k = 0, p = m / 2; k < p; ++ k, w = w * wn) {
                    t = w * base[(i + k + m / 2) * dft_w];
                    u = base[(i + k) * dft_w];
                    base[(i + k) * dft_w] = u + t;
                    base[(i + k + m / 2) * dft_w] = u - t;
                }
            }
        }
    }

    // Inverse
    if (inverse) {
        double inv = 1.0 / (dft_w * dft_h);
        for (int i = 0; i < dft_w * dft_h; ++ i) {
            dft_space[i] = dft_space[i] * inv;
        }
    }
}


/// Multiply DFT result 2 into result 1
void dft_multiply(int dft_w, int dft_h, ComplexPixel* dft_space1, ComplexPixel* dft_space2) {
    assert(dft_w > 0 and dft_h > 0 and dft_w == dft_lowbit(dft_w) and dft_h == dft_lowbit(dft_h));
    for (int i = 0; i < dft_w * dft_h; ++ i) {
        dft_space1[i] = dft_space1[i] * dft_space2[i];
    }
}


/// Free DFT result memory
void dft_free(ComplexPixel* dft_space) {
    assert(dft_space);
    std::free(dft_space);
}




// rewrite with fftw library
namespace DFT{
	// load in pixel
	void dft_load(const std::shared_ptr<Image> &image, int dft_w, int dft_h, 
			fftw_complex *matR, fftw_complex *matG, fftw_complex *matB){
		// matR/matG/matB has been initialized to zero before function
		// matR/matG/matB: (dft_w * dft_h) * 2, row major	

		for(int i = 0, index = 0; i < image -> h; ++i){
			for(int j = 0; j < image -> w; ++j, ++index){
				
				matR[i * dft_w + j][REAL] = (double)(image -> data[index].r);
				matG[i * dft_w + j][REAL] = (double)(image -> data[index].g);
				matB[i * dft_w + j][REAL] = (double)(image -> data[index].b);
				

				if(std::isnan((double)(image -> data[index].r))){ printf("NAN in data!!!! \n"); }

				// just to make sure no garbage in memory
				matR[i * dft_w + j][IMAG] = 0.0;
				matG[i * dft_w + j][IMAG] = 0.0;
				matB[i * dft_w + j][IMAG] = 0.0;
			}
		}
	}


	void dft(fftw_complex *in, fftw_complex *out, const int dft_w, const int dft_h){
		// initialize 2d fft
		fftw_plan p;
		p = fftw_plan_dft_2d(dft_h, dft_w, in, out, FFTW_FORWARD, FFTW_ESTIMATE); // if use FFTW_MEASURE, need to initialze input after set the plan
		fftw_execute(p);
		// clean up
		fftw_destroy_plan(p);
		fftw_cleanup();
	}


	void idft(fftw_complex *in, fftw_complex *out, const int dft_w, const int dft_h){
		// initialize 2d fft
		fftw_plan p;
		p = fftw_plan_dft_2d(dft_h, dft_w, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(p);
		// clean up
		fftw_destroy_plan(p);
		fftw_cleanup();
	}


	void dft_multiply(const int dft_w, const int dft_h, 
				fftw_complex *a, fftw_complex *b, fftw_complex *outFFT){
		// multiply complex dft result
		// TODO: apply omp here
		
		//printf("compled dft multiplication\n");
		for(int i = 0; i < dft_h; ++i){
			for(int j = 0; j < dft_w; ++j){
				// i for row, j for col
				/*							
				outFFT[i * dft_w + j][REAL] = a[i * dft_w + j][REAL]*b[i * dft_w + j][REAL] 
								- a[i * dft_w + j][IMAG]*b[i * dft_w + j][IMAG];
				outFFT[i * dft_w + j][IMAG] = a[i * dft_w + j][REAL]*b[i * dft_w + j][IMAG] 
								+ a[i * dft_w + j][IMAG]*b[i * dft_w + j][REAL];
				*/


					
				// compensate conjugate
				outFFT[i * dft_w + j][REAL] = a[i * dft_w + j][REAL]*b[i * dft_w + j][REAL]
								- a[i * dft_w + j][IMAG]*b[i * dft_w + j][IMAG];
				outFFT[i * dft_w + j][IMAG] = -1.0 * a[i * dft_w + j][REAL]*b[i * dft_w + j][IMAG]
								- a[i * dft_w + j][IMAG]*b[i * dft_w + j][REAL];
				
			}
		
		}
	}


} // end namespace FFT










