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

    std::cout << "dft_space1 before dft: " << std::endl;
    for (int y = 0, index = 0; y < dft_h; ++ y) {
        for (int x = 0; x < dft_w; ++ x, ++ index) {
            std::cout << dft_space1[index].b << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    dft_alloc(image2, dft_w, dft_h, dft_space2);
    dft(dft_w, dft_h, dft_space1);
    dft(dft_w, dft_h, dft_space2);
    
    
    std::cout << "dft_space1 after dft: " << std::endl;
    for (int y = 0, index = 0; y < dft_h; ++ y) {
        for (int x = 0; x < dft_w; ++ x, ++ index) {
            std::cout << dft_space1[index].b << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    
    
    
    std::cout << std::endl;
    dft_multiply(dft_w, dft_h, dft_space1, dft_space2);
    dft(dft_w, dft_h, dft_space1, true);
     
    std::cout << "dft_space1 after multiplication: " << std::endl;
    for (int y = 0, index = 0; y < dft_h; ++ y) {
        for (int x = 0; x < dft_w; ++ x, ++ index) {
            std::cout << dft_space1[index].b << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;





    fftw_complex *fftw_space1In= (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex));
    fftw_complex *fftw_space1Out= (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex));

    // TODO: test multiplication result with added minus sign for conjugate
    fftw_complex *fftw_space2In = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex));
    fftw_complex *fftw_space2Out = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex));

    fftw_complex *fftw_resultIn = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex));
    fftw_complex *fftw_resultOut = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex));




    for(int i = 0, index = 0; i < image1 -> h; ++i){
	for(int j = 0; j < image1 -> w; ++j, ++index){
				
	    fftw_space1In[i * dft_w + j][REAL] = (double)(image1 -> data[index].b);
	}
    }
   
    for(int i = 0, index = 0; i < image2 -> h; ++i){
        for(int j = 0; j < image2 -> w; ++j, ++index){
	    fftw_space2In[i * dft_w + j][REAL] = (double)(image2 -> data[index].b);
	}
    }

   


    std::cout << "fftw_space1In：" << std::endl;
    for (int i = 0; i < dft_h; ++ i) {
        for (int j = 0; j < dft_w; ++ j) {
            std::cout << "(" << fftw_space1In[i * dft_w + j][REAL] << ","
                    << fftw_space1In[i * dft_w + j][IMAG] << ") ";
        }
        std::cout << std::endl;
    }



    DFT::dft(fftw_space1In, fftw_space1Out, dft_w, dft_h);
    DFT::dft(fftw_space2In, fftw_space2Out, dft_w, dft_h);


    std::cout << "fftw_space1In：" << std::endl;
    for (int i = 0; i < dft_h; ++ i) {
        for (int j = 0; j < dft_w; ++ j) {
            std::cout << "(" << fftw_space1Out[i * dft_w + j][REAL] << ","
                    << fftw_space1Out[i * dft_w + j][IMAG] << ") ";
        }
        std::cout << std::endl;
    }



    DFT::dft_multiply(dft_w, dft_h, fftw_space1Out, fftw_space2Out, fftw_resultIn);
    
    DFT::dft(fftw_resultIn, fftw_resultOut, dft_w, dft_h);


    std::cout << std::endl;
    std::cout << "fftw_resultOut normalized：" << std::endl;
    for (int i = 0; i < dft_h; ++ i) {
        for (int j = 0; j < dft_w; ++ j) {
            std::cout << "(" << fftw_resultOut[i * dft_w + j][REAL] / (dft_h * dft_w) << ","
                    << fftw_resultOut[i * dft_w + j][IMAG] / (dft_h * dft_w) << ") ";
        }
        std::cout << std::endl;
    }




    // Show result

    /*
    std::cout << std::endl;
    std::cout << "dft_space2: " << std::endl;
    for (int y = 0, index = 0; y < dft_h; ++ y) {
        for (int x = 0; x < dft_w; ++ x, ++ index) {
            std::cout << dft_space2[index].b << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "fftw_space2Out：" << std::endl;
    for (int i = 0; i < dft_h; ++ i) {
        for (int j = 0; j < dft_w; ++ j) {
            std::cout << "(" << fftw_space2Out[i * dft_w + j][REAL] << ","
                    << fftw_space2Out[i * dft_w + j][IMAG] << ") ";
        }
        std::cout << std::endl;
    }
    */

    /*
    std::cout << "dft_space1: " << std::endl;
    for (int y = 0, index = 0; y < dft_h; ++ y) {
        for (int x = 0; x < dft_w; ++ x, ++ index) {
            std::cout << dft_space1[index].b << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "fftw_resultOut：" << std::endl;
    for (int i = 0; i < dft_h; ++ i) {
        for (int j = 0; j < dft_w; ++ j) {
            std::cout << "(" << fftw_resultOut[i * dft_w + j][REAL] << ","
                    << fftw_resultOut[i * dft_w + j][IMAG] << ") ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "fftw_resultOut normalized：" << std::endl;
    for (int i = 0; i < dft_h; ++ i) {
        for (int j = 0; j < dft_w; ++ j) {
            std::cout << "(" << fftw_resultOut[i * dft_w + j][REAL] / (dft_h * dft_w) << ","
                    << fftw_resultOut[i * dft_w + j][IMAG] / (dft_h * dft_w) << ") ";
        }
        std::cout << std::endl;
    }
    */

    // Free
    dft_free(dft_space1);
    dft_free(dft_space2);


    fftw_free(fftw_space1In);
    fftw_free(fftw_space1Out);
    fftw_free(fftw_space2In);
    fftw_free(fftw_space2Out);

    fftw_free(fftw_resultIn);
    fftw_free(fftw_resultOut);
    

}
