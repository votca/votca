// 
// File:   crosscorrelate.cc
// Author: ruehle
//
// Created on May 22, 2007, 5:42 PM
//

#include <fftw3.h>
#include "crosscorrelate.h"

/**
    \todo clean implementation!!!
*/
void CrossCorrelate::AutoCorrelate(DataCollection<double>::selection *data, bool average)
{
    size_t N = (*data)[0].size();
    _corrfunc.resize(N);

    fftw_complex *tmp;
    fftw_plan fft, ifft;
    
    tmp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));

    fft = fftw_plan_dft_r2c_1d(N, &(*data)[0][0], tmp,
                                    FFTW_ESTIMATE);
    ifft = fftw_plan_dft_c2r_1d(N, tmp, &_corrfunc[0],
                                    FFTW_ESTIMATE);
    fftw_execute(fft);
    
    tmp[0][0] = tmp[0][1] = 0;
    for(int i=1; i<N/2+1; i++) {
        tmp[i][0] = tmp[i][0]*tmp[i][0] + tmp[i][1]*tmp[i][1];
        tmp[i][1] = 0;       
    }
    fftw_execute(ifft);
    
    /*double m=0;
    for(int i=0; i<N; i++) {
        _corrfunc[i] = 0;
        m+=(*data)[0][i];
    }
    m=m/(double)N;
    for(int i=0;i<N; i++)
        for(int j=0; j<N-i-1; j++)
            _corrfunc[i]+=((*data)[0][j]-m)*((*data)[0][(i+j)]-m);
    */
    double d = _corrfunc[0];
    for(int i=0; i<N; i++)
        _corrfunc[i] = _corrfunc[i]/d;
    //cout << *data << endl;
    fftw_destroy_plan(fft);
    fftw_destroy_plan(ifft);
    fftw_free(tmp);
}

void CrossCorrelate::AutoFourier(vector <double>& ivec){
    size_t N = ivec.size();
    _corrfunc.resize(N);

    fftw_complex *tmp;
    fftw_plan fft;
    
    tmp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));

    fft = fftw_plan_dft_r2c_1d(N, &ivec[0], tmp, FFTW_ESTIMATE);
    fftw_execute(fft);
    
    tmp[0][0] = tmp[0][1] = 0;
    for(int i=1; i<N/2+1; i++) {
        tmp[i][0] = tmp[i][0]*tmp[i][0] + tmp[i][1]*tmp[i][1];
        tmp[i][1] = 0;       
    }
    
    // copy the real component of temp to the _corrfunc vector
    for(int i=0; i<N; i++){
        _corrfunc[i] = tmp[i][0];
    }
    
    fftw_destroy_plan(fft);
    fftw_free(tmp);
}

void CrossCorrelate::AutoCorr(vector <double>& ivec){
    size_t N = ivec.size();
    _corrfunc.resize(N);

    fftw_complex *tmp;
    fftw_plan fft, ifft;
    
    tmp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));

    fft = fftw_plan_dft_r2c_1d(N, &ivec[0], tmp, FFTW_ESTIMATE);
    ifft = fftw_plan_dft_c2r_1d(N, tmp, &_corrfunc[0], FFTW_ESTIMATE);
    fftw_execute(fft);
    
    tmp[0][0] = tmp[0][1] = 0;
    for(int i=1; i<N/2+1; i++) {
        tmp[i][0] = tmp[i][0]*tmp[i][0] + tmp[i][1]*tmp[i][1];
        tmp[i][1] = 0;       
    }
    
    fftw_execute(ifft);
    
    double d = _corrfunc[0];
    for(int i=0; i<N; i++)
        _corrfunc[i] = _corrfunc[i]/d;
    
    fftw_destroy_plan(fft);
    fftw_destroy_plan(ifft);
    fftw_free(tmp);
}