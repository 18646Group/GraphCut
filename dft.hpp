#pragma once

#include "image.hpp"
#include <fftw3.h>
#include <cstdlib> // for debug
#include <cmath>

// macros for complex number
#define REAL 0
#define IMAG 1

struct ComplexPixel {
    std::complex<double> r, g, b;

    ComplexPixel() = default;

    ComplexPixel(uint8_t r_real, uint8_t g_real, uint8_t b_real) : r(r_real), g(g_real), b(b_real) {}

    ComplexPixel(double real, double image) : r(real, image), g(real, image), b(real, image) {}

    ComplexPixel(const std::complex<double>& r, const std::complex<double>& g, const std::complex<double>& b) : r(r), g(g), b(b) {}

    [[nodiscard]] inline double real_sum() const {
        return r.real() + g.real() + b.real();
    }

    friend ComplexPixel operator+(const ComplexPixel& x, const ComplexPixel& y) {
        return ComplexPixel(x.r + y.r, x.g + y.g, x.b + y.b);
    }

    friend ComplexPixel operator-(const ComplexPixel& x, const ComplexPixel& y) {
        return ComplexPixel(x.r - y.r, x.g - y.g, x.b - y.b);
    }

    friend ComplexPixel operator*(const ComplexPixel& x, const ComplexPixel& y) {
        return ComplexPixel(x.r * y.r, x.g * y.g, x.b * y.b);
    }

    friend ComplexPixel operator/(const ComplexPixel& x, double k) {
        return ComplexPixel(x.r / k, x.g / k, x.b / k);
    }

    friend ComplexPixel operator*(const ComplexPixel& x, double k) {
        return ComplexPixel(x.r * k, x.g * k, x.b * k);
    }
};

/// Convert `Pixel` into `ComplexPixel`
[[nodiscard]] ComplexPixel to_complex_pixel(const Pixel& x) {
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
void dft_alloc(const std::shared_ptr<Image>& image, int dft_w, int dft_h, ComplexPixel*& dft_space) {
    dft_space = static_cast<ComplexPixel*>(std::malloc(dft_w * dft_h * sizeof(ComplexPixel)));
    std::fill(dft_space, dft_space + dft_w * dft_h, ComplexPixel());
    for (int i = 0, index = 0; i < image->h; ++i) {
        for (int j = 0; j < image->w; ++j, ++index) {
            dft_space[i * dft_w + j] = to_complex_pixel(image->data[index]);
        }
    }
}

/// Run DFT and IDFT
void dft(int dft_w, int dft_h, ComplexPixel* dft_space, bool inverse = false) {
    // Allocate space
    double coefficient = inverse ? -1 : 1;
    assert(dft_w > 0 and dft_h > 0 and dft_w == dft_lowbit(dft_w) and dft_h == dft_lowbit(dft_h));

    // Butterfly changes by w
    for (int row = 0; row < dft_h; ++row) {
        ComplexPixel* base = dft_space + row * dft_w;
        for (int i = 0, j = 0; i < dft_w; ++i) {
            if (i > j) {
                std::swap(base[i], base[j]);
            }
            for (int t = dft_w / 2; (j ^= t) < t; t /= 2)
                ;
        }
    }

    // Butterfly changes by h
    for (int col = 0; col < dft_w; ++col) {
        ComplexPixel* base = dft_space + col;
        for (int i = 0, j = 0; i < dft_h; ++i) {
            if (i > j) {
                std::swap(base[i * dft_w], base[j * dft_w]);
            }
            for (int t = dft_h / 2; (j ^= t) < t; t /= 2)
                ;
        }
    }

    // DFT by w
    for (int row = 0; row < dft_h; ++row) {
        ComplexPixel* base = dft_space + row * dft_w;
        ComplexPixel wn, w, t, u;
        for (int m = 2; m <= dft_w; m *= 2) {
            wn = ComplexPixel(cos(2.0 * M_PI / m), coefficient * sin(2.0 * M_PI / m));
            for (int i = 0; i < dft_w; i += m) {
                w = ComplexPixel(1, 0);
                for (int k = 0, p = m / 2; k < p; ++k, w = w * wn) {
                    t = w * base[i + k + m / 2];
                    u = base[i + k];
                    base[i + k] = u + t;
                    base[i + k + m / 2] = u - t;
                }
            }
        }
    }

    // DFT by h
    for (int col = 0; col < dft_w; ++col) {
        ComplexPixel* base = dft_space + col;
        ComplexPixel wn, w, t, u;
        for (int m = 2; m <= dft_h; m *= 2) {
            wn = ComplexPixel(cos(2.0 * M_PI / m), coefficient * sin(2.0 * M_PI / m));
            for (int i = 0; i < dft_h; i += m) {
                w = ComplexPixel(1, 0);
                for (int k = 0, p = m / 2; k < p; ++k, w = w * wn) {
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
        for (int i = 0; i < dft_w * dft_h; ++i) {
            dft_space[i] = dft_space[i] * inv;
        }
    }
}

/// Multiply DFT result 2 into result 1
void dft_multiply(int dft_w, int dft_h, ComplexPixel* dft_space1, ComplexPixel* dft_space2) {
    assert(dft_w > 0 and dft_h > 0 and dft_w == dft_lowbit(dft_w) and dft_h == dft_lowbit(dft_h));
    for (int i = 0; i < dft_w * dft_h; ++i) {
        dft_space1[i] = dft_space1[i] * dft_space2[i];
    }
}

/// Free DFT result memory
void dft_free(ComplexPixel* dft_space) {
    assert(dft_space);
    std::free(dft_space);
}

// rewrite with fftw library
namespace DFT {
    class Domain {
    private:
        int dft_w;
        int dft_h;
        int size;

    public:
        fftw_complex* R;
        fftw_complex* G;
        fftw_complex* B;

        Domain(int dft_w, int dft_h) : dft_w(dft_w), dft_h(dft_h) {
            size = dft_h * dft_w * sizeof(fftw_complex);
            R = (fftw_complex*)fftw_malloc(size);
            G = (fftw_complex*)fftw_malloc(size);
            B = (fftw_complex*)fftw_malloc(size);
        }

        ~Domain() {
            fftw_free(R);
            fftw_free(G);
            fftw_free(B);
        }

        void load(const std::shared_ptr<Image>& image) {
            memset(R, 0x00, size);
            memset(G, 0x00, size);
            memset(B, 0x00, size);
            for (int i = 0, index = 0; i < image->h; i++) {
                for (int j = 0; j < image->w; j++, index++) {
                    R[i * dft_w + j][REAL] = (double)(image->data[index].r);
                    G[i * dft_w + j][REAL] = (double)(image->data[index].g);
                    B[i * dft_w + j][REAL] = (double)(image->data[index].b);
                }
            }
        }

        void loadProductOf(Domain &d1, Domain &d2) {
            assert(dft_w == d1.dft_w && dft_w == d2.dft_w);
            assert(dft_w == d1.dft_h && dft_w == d2.dft_w);
            memset(R, 0x00, size);
            memset(G, 0x00, size);
            memset(B, 0x00, size);

            // TODO: apply omp/simd here
            for (int i = 0; i < dft_h; i++) {
                for (int j = 0; j < dft_w; j++) {
                    // R
                    R[i * dft_w + j][REAL] =
                        d1.R[i * dft_w + j][REAL] * d2.R[i * dft_w + j][REAL] -
                        d1.R[i * dft_w + j][IMAG] * d2.R[i * dft_w + j][IMAG];
                    R[i * dft_w + j][IMAG] =
                        d1.R[i * dft_w + j][REAL] * d2.R[i * dft_w + j][IMAG] +
                        d1.R[i * dft_w + j][IMAG] * d2.R[i * dft_w + j][REAL];
                    // G
                    G[i * dft_w + j][REAL] =
                        d1.G[i * dft_w + j][REAL] * d2.G[i * dft_w + j][REAL] -
                        d1.G[i * dft_w + j][IMAG] * d2.G[i * dft_w + j][IMAG];
                    G[i * dft_w + j][IMAG] =
                        d1.G[i * dft_w + j][REAL] * d2.G[i * dft_w + j][IMAG] +
                        d1.G[i * dft_w + j][IMAG] * d2.G[i * dft_w + j][REAL];
                    // B
                    B[i * dft_w + j][REAL] =
                        d1.B[i * dft_w + j][REAL] * d2.B[i * dft_w + j][REAL] -
                        d1.B[i * dft_w + j][IMAG] * d2.B[i * dft_w + j][IMAG];
                    B[i * dft_w + j][IMAG] =
                        d1.B[i * dft_w + j][REAL] * d2.B[i * dft_w + j][IMAG] +
                        d1.B[i * dft_w + j][IMAG] * d2.B[i * dft_w + j][REAL];
                }
            }
        }

        inline double realsum(int i, int j) {
            return (
                R[i * dft_w + j][REAL] + 
                G[i * dft_w + j][REAL] + 
                B[i * dft_w + j][REAL]
            ) / (
                dft_h * dft_w
            );
        }
    };

    static int planner_cnt = 0;

    class Planner {
        private:
        int dft_w;
        int dft_h;
        fftw_plan flippedPlanR;
        fftw_plan flippedPlanG;
        fftw_plan flippedPlanB;
        fftw_plan canvasPlanR;
        fftw_plan canvasPlanG;
        fftw_plan canvasPlanB;
        fftw_plan outputPlanR;
        fftw_plan outputPlanG;
        fftw_plan outputPlanB;

        public:
        Domain *flipped;
        Domain *flippedFFT;
        Domain *canvas;
        Domain *canvasFFT;
        Domain *output;
        Domain *outputFFT;

        Planner(int dft_w, int dft_h) : dft_w(dft_w), dft_h(dft_h) {
            // initialize domains
            flipped = new Domain(dft_w, dft_h);
            flippedFFT = new Domain(dft_w, dft_h);
            canvas = new Domain(dft_w, dft_h);
            canvasFFT = new Domain(dft_w, dft_h);
            output = new Domain(dft_w, dft_h);
            outputFFT = new Domain(dft_w, dft_h);
            // initialize plans
            flippedPlanR = fftw_plan_dft_2d(dft_h, dft_w, flipped->R, flippedFFT->R, FFTW_FORWARD, FFTW_ESTIMATE);
            flippedPlanG = fftw_plan_dft_2d(dft_h, dft_w, flipped->R, flippedFFT->R, FFTW_FORWARD, FFTW_ESTIMATE);
            flippedPlanB = fftw_plan_dft_2d(dft_h, dft_w, flipped->R, flippedFFT->R, FFTW_FORWARD, FFTW_ESTIMATE);
            canvasPlanR = fftw_plan_dft_2d(dft_h, dft_w, canvas->R, canvasFFT->R, FFTW_FORWARD, FFTW_ESTIMATE);
            canvasPlanG = fftw_plan_dft_2d(dft_h, dft_w, canvas->R, canvasFFT->R, FFTW_FORWARD, FFTW_ESTIMATE);
            canvasPlanB = fftw_plan_dft_2d(dft_h, dft_w, canvas->R, canvasFFT->R, FFTW_FORWARD, FFTW_ESTIMATE);
            outputPlanR = fftw_plan_dft_2d(dft_h, dft_w, outputFFT->R, output->R, FFTW_BACKWARD, FFTW_ESTIMATE);
            outputPlanG = fftw_plan_dft_2d(dft_h, dft_w, outputFFT->R, output->R, FFTW_BACKWARD, FFTW_ESTIMATE);
            outputPlanB = fftw_plan_dft_2d(dft_h, dft_w, outputFFT->R, output->R, FFTW_BACKWARD, FFTW_ESTIMATE);
            // add planner cnt
            planner_cnt ++;
        }

        ~Planner() {
            delete flipped;
            delete flippedFFT;
            delete canvas;
            delete canvasFFT;
            delete output;
            delete outputFFT;

            fftw_destroy_plan(flippedPlanR);
            fftw_destroy_plan(flippedPlanG);
            fftw_destroy_plan(flippedPlanB);

            fftw_destroy_plan(canvasPlanR);
            fftw_destroy_plan(canvasPlanG);
            fftw_destroy_plan(canvasPlanB);

            fftw_destroy_plan(outputPlanR);
            fftw_destroy_plan(outputPlanG);
            fftw_destroy_plan(outputPlanB);

            planner_cnt --;
            if (!planner_cnt) {
                fftw_cleanup();
            }
        }

        void main(const std::shared_ptr<Image>& img_flipped, const std::shared_ptr<Image>& img_canvas) {
            flipped->load(img_flipped);
            canvas->load(img_canvas);
            // TODO: OpenMP
            fftw_execute(flippedPlanR);
            fftw_execute(flippedPlanG);
            fftw_execute(flippedPlanB);
            fftw_execute(canvasPlanR);
            fftw_execute(canvasPlanG);
            fftw_execute(canvasPlanB);
            outputFFT->loadProductOf(*flippedFFT, *canvasFFT);
            // TODO: OpenMP
            fftw_execute(outputPlanR);
            fftw_execute(outputPlanG);
            fftw_execute(outputPlanB);
        }
    };

    static Planner* planner = NULL;
} // end namespace FFT
