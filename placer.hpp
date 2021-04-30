#pragma once

#include "cherry.hpp"
#include "dft.hpp"
#include "image.hpp"
#include <fftw3.h>

#include <cmath>	// for debug

class Placer {
public:
    static void init(const std::shared_ptr<Canvas> &canvas, const std::shared_ptr<Image> &texture) {
        // apply patches on canvas
        auto random_y = Random(texture->h / 3, texture->h * 2 / 3);
        auto random_x = Random(texture->w / 3, texture->w * 2 / 3);
        for (int y = 0; y < canvas->h; y += random_y()) {
            for (int x = 0; x < canvas->w; x += random_x()) {
                auto patch = std::make_shared<Patch>(texture, x, y);
                canvas->apply(patch);
            }
        }
    }

    static void random(const std::shared_ptr<Canvas> &canvas, const std::shared_ptr<Image> &texture) {
        auto random_x = Random(0, texture->w - 1);
        auto random_y = Random(0, texture->h - 1);
        auto patch = std::make_shared<Patch>(texture, random_x(), random_y());
        canvas->apply(patch);
    }

    // Bigger means more randomness
    static constexpr double possibility_k = 0.3;

    static void entire_matching(const std::shared_ptr<Canvas> &canvas, const std::shared_ptr<Image> &texture, 
		    		bool random=false, int times=100) {
        
	//printf("entire matching start \n");
	std::shared_ptr<Patch> best_patch;

        if (random) {
            auto random_x = Random(0, canvas->w - 1);
            auto random_y = Random(0, canvas->h - 1);

            uint64_t best_ssd = UINT64_MAX;
            for (int i = 0; i < times; ++ i) {
                auto patch = std::make_shared<Patch>(texture, random_x(), random_y());
                uint64_t ssd = canvas->ssd(patch);
                if (ssd < best_ssd) {
                    best_ssd = ssd;
                    best_patch = patch;
                }
            }
        } else {
            // FFT-based acceleration
            // Prefix sum
            assert(canvas->none_empty());
            auto *texture_sum = static_cast<uint64_t*> (std::malloc(texture->w * texture->h * sizeof(uint64_t)));
            auto *canvas_sum = static_cast<uint64_t*> (std::malloc(canvas->w * canvas->h * sizeof(uint64_t)));
            auto query = [](const uint64_t *sum, int x, int y, int size_x, int size_y, int w, int h) {
                int last_x = x + size_x - 1, last_y = y + size_y - 1;
                uint64_t result = sum[last_y * w + last_x];
                result += (x > 0 and y > 0) ? sum[(y - 1) * w + x - 1] : 0;
                result -= x > 0 ? sum[last_y * w + x - 1] : 0;
                result -= y > 0 ? sum[(y - 1) * w + last_x] : 0;
                return result;
            };
            auto do_prefix_sum = [](int w, int h, Pixel *pixels, uint64_t *sum) {
                for (int y = 0, index = 0; y < h; ++ y) {
                    for (int x = 0; x < w; ++ x, ++ index) {
                        auto up = y > 0 ? sum[index - w] : 0;
                        auto left = x > 0 ? sum[index - 1] : 0;
                        auto left_up = (y > 0 and x > 0) ? sum[index - w - 1] : 0;
                        sum[index] = up + left + pixels[index].sqr_sum() - left_up;
                    }
                }
            };
            do_prefix_sum(texture->w, texture->h, texture->data, texture_sum);
            do_prefix_sum(canvas->w, canvas->h, canvas->data, canvas_sum);

            // FFT
	    auto flipped = texture->flip();	// why flip?
            int dft_w = dft_round(texture->w + canvas->w);	// pad length to 2^n
            int dft_h = dft_round(texture->h + canvas->h); 	

	    
	    //printf("allocate dft buffer \n");
	    fftw_complex  *flipped_RIn = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
		    	  *flipped_ROut = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *flipped_GIn = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *flipped_GOut = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *flipped_BIn = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *flipped_BOut = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex));
	
	    fftw_complex  *canvas_RIn = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
		    	  *canvas_ROut = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *canvas_GIn = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *canvas_GOut = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *canvas_BIn = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *canvas_BOut = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex));

	    fftw_complex  *outputFFT_R = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *output_R = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *outputFFT_G = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *output_G = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *outputFFT_B = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex)),
			  *output_B = (fftw_complex*)fftw_malloc(dft_h * dft_w * sizeof(fftw_complex));

        /*
	    // load in pixel values
	    DFT::dft_load(flipped, dft_w, dft_h, 
			    flipped_RIn, flipped_GIn, flipped_BIn);

	    DFT::dft_load(canvas, dft_w, dft_h, 
			    canvas_RIn, canvas_GIn, canvas_BIn);
        */

	   // TODO: there are nan value when loading

	    // perform fft for each channel
        /*
	    DFT::dft(flipped_RIn, flipped_ROut, dft_w, dft_h);
	    DFT::dft(flipped_GIn, flipped_GOut, dft_w, dft_h);
	    DFT::dft(flipped_BIn, flipped_BOut, dft_w, dft_h);

	    DFT::dft(canvas_RIn, canvas_ROut, dft_w, dft_h);
	    DFT::dft(canvas_GIn, canvas_GOut, dft_w, dft_h);
	    DFT::dft(canvas_BIn, canvas_BOut, dft_w, dft_h);
        */




        // dft
        // initiazte plan -> load data to memory -> execute plan -> clean up
        // TODO: reuse planner
        fftw_plan pR, pB, pG;
	// if use FFTW_MEASURE, need to initialze input after set the plan
	pR = fftw_plan_dft_2d(dft_h, dft_w, flipped_RIn, flipped_ROut, FFTW_FORWARD, FFTW_ESTIMATE);
	pG = fftw_plan_dft_2d(dft_h, dft_w, flipped_GIn, flipped_GOut, FFTW_FORWARD, FFTW_ESTIMATE);
        pB = fftw_plan_dft_2d(dft_h, dft_w, flipped_BIn, flipped_BOut, FFTW_FORWARD, FFTW_ESTIMATE);

        DFT::dft_load(flipped, dft_w, dft_h, 
			    flipped_RIn, flipped_GIn, flipped_BIn);

        fftw_execute(pR);
        fftw_execute(pG);
        fftw_execute(pB);
		
        pR = fftw_plan_dft_2d(dft_h, dft_w, canvas_RIn, canvas_ROut, FFTW_FORWARD, FFTW_ESTIMATE);
        pG = fftw_plan_dft_2d(dft_h, dft_w, canvas_GIn, canvas_GOut, FFTW_FORWARD, FFTW_ESTIMATE);
        pB = fftw_plan_dft_2d(dft_h, dft_w, canvas_BIn, canvas_BOut, FFTW_FORWARD, FFTW_ESTIMATE);
        
        DFT::dft_load(canvas, dft_w, dft_h, 
			    canvas_RIn, canvas_GIn, canvas_BIn);

        fftw_execute(pR);
        fftw_execute(pG);
        fftw_execute(pB);

        // clean up
		fftw_destroy_plan(pR);
        fftw_destroy_plan(pG);
        fftw_destroy_plan(pB);
		fftw_cleanup();

        // idft
        // initiazte plan -> load data to memory -> execute plan -> clean up
	/*
	fftw_plan p_inv;
        p_inv = fftw_plan_dft_2d(dft_h, dft_w, canvas_RIn, canvas_ROut, FFTW_BACKWARD, FFTW_ESTIMATE);
	DFT::dft_multiply(dft_w, dft_h, flipped_ROut, canvas_ROut, outputFFT_R);
        fftw_execute(p_inv);

        p_inv = fftw_plan_dft_2d(dft_h, dft_w, canvas_GIn, canvas_GOut, FFTW_BACKWARD, FFTW_ESTIMATE);
        DFT::dft_multiply(dft_w, dft_h, flipped_GOut, canvas_GOut, outputFFT_G);
        fftw_execute(p_inv);
     
        p_inv = fftw_plan_dft_2d(dft_h, dft_w, canvas_BIn, canvas_BOut, FFTW_BACKWARD, FFTW_ESTIMATE);
	// complex multiplication
	DFT::dft_multiply(dft_w, dft_h, flipped_BOut, canvas_BOut, outputFFT_B);
        fftw_execute(p_inv);

        // clean up
	fftw_destroy_plan(p_inv);
	fftw_cleanup();
	*/

		
	fftw_plan p_inv;
        p_inv = fftw_plan_dft_2d(dft_h, dft_w, canvas_RIn, canvas_ROut, FFTW_FORWARD, FFTW_ESTIMATE);
        DFT::dft_multiply(dft_w, dft_h, flipped_ROut, canvas_ROut, outputFFT_R);
        fftw_execute(p_inv);

        p_inv = fftw_plan_dft_2d(dft_h, dft_w, canvas_GIn, canvas_GOut, FFTW_FORWARD, FFTW_ESTIMATE);
        DFT::dft_multiply(dft_w, dft_h, flipped_GOut, canvas_GOut, outputFFT_G);
        fftw_execute(p_inv);
     
        p_inv = fftw_plan_dft_2d(dft_h, dft_w, canvas_BIn, canvas_BOut, FFTW_FORWARD, FFTW_ESTIMATE);
        // complex multiplication
        DFT::dft_multiply(dft_w, dft_h, flipped_BOut, canvas_BOut, outputFFT_B);
        fftw_execute(p_inv);

        // clean up
        fftw_destroy_plan(p_inv);
        fftw_cleanup();
		

	    // ifft
        /*
	    DFT::idft(outputFFT_R, output_R, dft_w, dft_h);
	    DFT::idft(outputFFT_G, output_G, dft_w, dft_h);
	    DFT::idft(outputFFT_B, output_B, dft_w, dft_h);
        */
	    /*
    	    DFT::dft(outputFFT_R, output_R, dft_w, dft_h);
            DFT::dft(outputFFT_G, output_G, dft_w, dft_h);
            DFT::dft(outputFFT_B, output_B, dft_w, dft_h);	    
	    */

        


	    // TODO: there are nan in output_R, output_G, and output_B

        std::cout << std::endl;
            std::cout << "flipped_RIn：" << std::endl;
            for (int i = 0; i < 5; ++ i) {
                for (int j = 0; j < 5; ++ j) {
                    std::cout << "(" << flipped_RIn[i * dft_w + j][REAL]<< ","
                        << flipped_RIn[i * dft_w + j][IMAG] << ") ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;


	    if(std::isnan(flipped_ROut[0][REAL])){
	    	printf("nan in flipped_ROut\n");
	    }
	    

	    std::cout << std::endl;
            std::cout << "flipped_ROut：" << std::endl;
            for (int i = 0; i < 5; ++ i) {
                for (int j = 0; j < 5; ++ j) {
                    std::cout << "(" << flipped_ROut[i * dft_w + j][REAL]<< ","
                        << flipped_ROut[i * dft_w + j][IMAG] << ") ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;

	    /*
	    std::cout << std::endl;
            std::cout << "canvas_ROut：" << std::endl;
            for (int i = 0; i < 5; ++ i) {
                for (int j = 0; j < 5; ++ j) {
                    std::cout << "(" << canvas_ROut[i * dft_w + j][REAL]<< ","
                        << canvas_ROut[i * dft_w + j][IMAG] << ") ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;



	    std::cout << std::endl;
            std::cout << "outputFFT_R：" << std::endl;
            for (int i = 0; i < 5; ++ i) {
                for (int j = 0; j < 5; ++ j) {
                    std::cout << "(" << outputFFT_R[i * dft_w + j][REAL]<< ","
                        << outputFFT_R[i * dft_w + j][IMAG] << ") ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;


	    std::cout << std::endl;
	    std::cout << "output_R normalized：" << std::endl;
    	    for (int i = 0; i < 5; ++ i) {
        	for (int j = 0; j < 5; ++ j) {
            	    std::cout << "(" << output_R[i * dft_w + j][REAL] / (dft_h * dft_w) << ","
                    	<< output_R[i * dft_w + j][IMAG] / (dft_h * dft_w) << ") ";
        	}
        	std::cout << std::endl;
    	    }
	    std::cout << std::endl;
	    */

	    //printf("dft multiply finished \n");
	    
        // Get results
        uint64_t variance = texture->variance();
        uint64_t ssd = 0;
	    auto *possibility = static_cast<double*> (std::malloc(canvas->h * canvas->w * sizeof(double)));
            
	   
	    double normalize = 1.0 / (dft_w * dft_h); 
	    for (int y = 0, index = 0; y < canvas->h; ++ y) {
                for (int x = 0; x < canvas->w; ++ x, ++ index) {
                    int overlapped_w = std::min(texture->w, canvas->w - x);
                    int overlapped_h = std::min(texture->h, canvas->h - y);
                    
		    ssd = 0;
                    ssd += texture_sum[(overlapped_h - 1) * texture->w + overlapped_w - 1];
                    ssd += query(canvas_sum, x, y, overlapped_w, overlapped_h, canvas->w, canvas->h);
                    
		    // dft output used here
		    ssd -= (uint64_t) std::floor(2.0 * normalize * (output_R[(texture->h + y - 1)* dft_w + texture->w + x - 1][REAL] 
					    		+ output_G[(texture->h + y - 1)* dft_w + texture->w + x - 1][REAL]
							+ output_B[(texture->h + y - 1)* dft_w + texture->w + x - 1][REAL]) );

		    ssd /= overlapped_w * overlapped_h;
                    possibility[index] = std::exp(-1.0 * ssd / (possibility_k * variance));
                }
            }

	    //printf("get ssd \n");

            double possibility_sum = 0;
            for (int i = 0; i < canvas->h * canvas->w; ++ i) {
                possibility_sum += possibility[i];
            }

	    //printf("possibility sum: %f\n", possibility_sum);

            double position = Random<double>(0, 1)(), up = 0;
            for (int y = 0, index = 0; y < canvas->h and not best_patch; ++ y) {
                for (int x = 0; x < canvas->w; ++ x, ++ index) {
                    possibility[index] /= possibility_sum;
                    if (up + possibility[index] >= position) {
                        best_patch = std::make_shared<Patch>(texture, x, y);
                        break;
                    }
                    up += possibility[index];
                }
            }
            assert(best_patch);
	    
	    //printf("find current best patch\n");
	    

            // Free resources
	    //printf("start free resources\n");
            std::free(possibility);
            std::free(texture_sum);
            std::free(canvas_sum);
	  
	    
	    fftw_free(flipped_RIn); fftw_free(flipped_ROut);
	    fftw_free(flipped_GIn); fftw_free(flipped_GOut);
	    fftw_free(flipped_BIn); fftw_free(flipped_BOut);

	    fftw_free(canvas_RIn); fftw_free(canvas_ROut);
	    fftw_free(canvas_GIn); fftw_free(canvas_GOut);
	    fftw_free(canvas_BIn); fftw_free(canvas_BOut);

	    fftw_free(outputFFT_R); fftw_free(output_R);
	    fftw_free(outputFFT_G); fftw_free(output_G); 
	    fftw_free(outputFFT_B); fftw_free(output_B); 
	    
	    //printf("free all resources");
        }
	
	/*
	if((best_patch->x <= 0) || (best_patch->y <= 0) 
		|| (best_patch->x_end() > canvas->w) || (best_patch->y_end() > canvas->h)){
		printf("patch bigger than canvas!!! \n");
	}
	*/

	canvas->apply(best_patch);
	
	//printf("entire matching finished \n");
    }

    static void sub_patch_matching(const std::shared_ptr<Canvas> &canvas, const std::shared_ptr<Image> &texture, int times=100) {
        int sub_patch_w = texture->w / 3, sub_patch_h = texture->h / 3;
        auto random_canvas_x = Random(0, canvas->w - sub_patch_w), random_canvas_y = Random(0, canvas->h - sub_patch_h);
        int canvas_x = random_canvas_x(), canvas_y = random_canvas_y();

        auto random_x = Random(0, texture->w - sub_patch_w);
        auto random_y = Random(0, texture->h - sub_patch_h);
        uint64_t best_ssd = UINT64_MAX;
        std::shared_ptr<Patch> best_patch;
        for (int i = 0; i < times; ++ i) {
            int x = random_x(), y = random_y();
            auto patch = std::make_shared<Patch>(texture, canvas_x - x, canvas_y - y);
            uint64_t ssd = canvas->ssd(patch, canvas_x, canvas_y, sub_patch_w, sub_patch_h); // SSD only calculates the sub-patch region
            if (ssd < best_ssd) {
                best_ssd = ssd;
                best_patch = patch;
            }
        }
        canvas->apply(best_patch);
    }
};
