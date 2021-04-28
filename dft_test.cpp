#include <iostream>

#include "dft.hpp"


int main() {
    // Set pixels
    auto image1 = std::make_shared<Image>(3, 2);
    auto image2 = std::make_shared<Image>(2, 5);
    for (int y = 0; y < 2; ++ y) {
        for (int x = 0; x < 3; ++ x) {
            image1->set(x, y, Pixel(0, 1, 2));
        }
    }
    for (int y = 0; y < 5; ++ y) {
        for (int x = 0; x < 2; ++ x) {
            image2->set(x, y, Pixel(0, 1, 2));
        }
    }

    // DFT
    int dft_w = dft_round(5), dft_h = dft_round(7);
    ComplexPixel *dft_space1, *dft_space2;
    dft_alloc(image1, dft_w, dft_h, dft_space1);
    dft_alloc(image2, dft_w, dft_h, dft_space2);
    dft(dft_w, dft_h, dft_space1);
    dft(dft_w, dft_h, dft_space2);
    dft_multiply(dft_w, dft_h, dft_space1, dft_space2);
    dft(dft_w, dft_h, dft_space1, true);
    
     

    fftw_complex  *fftw_spaceIn= (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex));
    fftw_complex  *fftw_spaceOut= (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex));

    // TODO: test multiplication result with added minus sign for conjugate

    /*
    double *fftw_real = fftw_alloc_real(dft_h * dft_w);
    fftw_complex *fftw_imag = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex));	    
    */

    for(int i = 0, index = 0; i < image1 -> h; ++i){
	for(int j = 0; j < image1 -> w; ++j, ++index){
				
	    fftw_spaceIn[i * dft_w + j][REAL] = (double)(image1 -> data[index].b);
	
	    fftw_real[i * dft_w + j] = (double)(image1 -> data[index].b);
	}
    }
    
    DFT::dft(fftw_spaceIn, fftw_spaceOut, dft_w, dft_h);
    DFT::dft_real(fftw_real, fftw_imag, dft_w, dft_h);

    // Show result
    std::cout << "dft_space1: " << std::endl;
    for (int y = 0, index = 0; y < dft_h; ++ y) {
        for (int x = 0; x < dft_w; ++ x, ++ index) {
            std::cout << dft_space1[index].b << " ";
        }
        std::cout << std::endl;
    }


    std::cout << "fftw_space：" << std::endl;
    for (int i = 0; i < dft_h; ++ i) {
        for (int j = 0; j < dft_w; ++ j) {
            std::cout << "(" << fftw_spaceOut[i * dft_w + j][REAL] << "," 
		    << fftw_spaceOut[i * dft_w + j][IMAG] << ") ";
        }
        std::cout << std::endl;
    }

    std::cout << "fftw_imag：" << std::endl;
    for (int i = 0; i < dft_h; ++ i) {
	for (int j = 0; j < dft_w; ++ j) {
	    std::cout << "(" << fftw_imag[i * dft_w + j][REAL] << ","
		    << fftw_imag[i * dft_w + j][IMAG] << ") ";
	}
	std::cout << std::endl;
    }


    // Free
    dft_free(dft_space1);
    dft_free(dft_space2);


    fftw_free(fftw_spaceIn);
    fftw_free(fftw_spaceOut);
    
    
    //fftw_free(fftw_imag);

}
