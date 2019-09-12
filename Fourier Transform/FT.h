#pragma once
#include <iostream>
#include <algorithm> 
class Complex {
public:
	double rl, im;
	Complex();
	//Complex(double r, double i);

};

class FT
{
private:

public:
	FT();
	void DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void DFT(double** pFreqReal, double** pFreqImag, int** InputImage, int h, int w, int u, int v);

	void InverseDiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	void InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y);

	void FastFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	//void FFT(Complex* pVec, int len);
	void FFT(double *xreal, double* ximag, int n);
	void InverseFastFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
	//void InverseFFT(Complex* pVec, int len);
	void InverseFFT(double *xreal, double* ximag, int n);

	void LowpassFilter(double** Real, double** Img, int** output, double d, int h, int w);
	void HighpassFilter(double** Real, double** Img, int** output, double d, int h, int w);

};

int get_computation_layers(int num);
void swap(double &a, double &b);
void bitrp(double* real, double*imag, int n);