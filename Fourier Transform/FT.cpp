#include "FT.h"
#define PI 3.1415926
#include <algorithm>
#include<cmath>
using namespace std;

//#define DEBUG

FT::FT()
{
}

void FT::DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt < M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // 傅立葉頻率陣列
	}
	for (int forzero_i = 0; forzero_i < M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j < N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			DFT(FreqReal, FreqImag, InputImage, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
		}
	}
	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}

void FT::DFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int y = 0; y < M; y++)
	{
		for (int x = 0; x < N; x++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleDFT = (-1.0f * 2.0f * 3.14159 * (double)(u*x + v * y) / (double)M);
			double c = cos(angleDFT);
			double s = sin(angleDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			pFreqReal[u][v] += (double)InputImage[y][x] * c;
			pFreqImag[u][v] -= (double)InputImage[y][x] * s;
		}
	}

	pFreqReal[u][v] = pFreqReal[u][v] / (double)(M);
	pFreqImag[u][v] = pFreqImag[u][v] / (double)(M);
}

void FT::InverseDiscreteFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** InverseReal = new double*[M];
	double** InverseImag = new double*[M];
	double** pFreq = new double*[M];

	for (int i = 0; i < M; i++)
	{
		InverseReal[i] = new double[N];
		InverseImag[i] = new double[N];
		pFreq[i] = new double[N]; // 傅立葉頻率陣列
	}

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			InverseReal[i][j] = 0.0f;
			InverseImag[i][j] = 0.0f;
			pFreq[i][j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			InverseDFT(InverseReal, InverseImag, FreqReal, FreqImag, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(InverseReal[i][j], (double) 2.0) + pow(InverseImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
			//存下反傅立葉實數與虛數部分
			FreqReal[i][j] = InverseReal[i][j];
			FreqImag[i][j] = InverseImag[i][j];
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		delete[] pFreq[i];
		delete[] InverseReal[i];
		delete[] InverseImag[i];

	}
	delete[] pFreq;
	delete[] InverseReal;
	delete[] InverseImag;

}

void FT::InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int v = 0; v < M; v++)
	{
		for (int u = 0; u < N; u++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleIDFT = (2.0f * 3.14159 * (double)(u*x + v * y) / (double)M);
			double c = cos(angleIDFT);
			double s = sin(angleIDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			InverseReal[x][y] += (pFreqReal[v][u] * c - pFreqImag[v][u] * s);
			InverseImag[x][y] += (pFreqReal[v][u] * s + pFreqImag[v][u] * c);
		}
	}
	InverseReal[x][y] = InverseReal[x][y] / (float)M;
	InverseImag[x][y] = InverseImag[x][y] / (float)M;
}

void FT::FastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;
	Complex *pVec;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt < M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // 傅立葉頻率陣列
	}
	for (int forzero_i = 0; forzero_i < M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j < N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	pVec = new Complex[N];
	double* xreal = new double[N];
	double* ximag = new double[N];
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++) {
#ifdef DEBUG
			pVec[j].rl = InputImage[j][i];
			pVec[j].im = 0;
#endif // DEBUG

#ifndef DEBUG
			xreal[j] = InputImage[j][i];
			ximag[j] = 0;
#endif // !DEBUG
		}
		//FFT(pVec, N);

		FFT(xreal, ximag, M);
		for (int j = 0; j < N; j++) {
#ifdef DEBUG
			FreqReal[j][i] = pVec[j].rl / (double)(M);
			FreqImag[j][i] = pVec[j].im / (double)(M);
#endif // DEBUG
#ifndef DEBUG
			FreqReal[j][i] = xreal[j] / (double)(M);
			FreqImag[j][i] = ximag[j] / (double)(M);
#endif // !DEBUG

		}
	}
	delete[] pVec;
	pVec = new Complex[M];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++) {
#ifdef DEBUG
			pVec[j].rl = FreqReal[i][j];
			pVec[j].im = FreqImag[i][j];
#endif // DEBUG

#ifndef DEBUG
			xreal[j] = FreqReal[i][j];
			ximag[j] = FreqImag[i][j];
#endif // !DEBUG
		}
		//(pVec, M);
		FFT(xreal, ximag, N);
		for (int j = 0; j < M; j++) {
#ifdef DEBUG
			FreqReal[i][j] = pVec[j].rl;
			FreqImag[i][j] = pVec[j].im;
#endif // DEBUG
#ifndef DEBUG
			FreqReal[i][j] = xreal[j];
			FreqImag[i][j] = ximag[j];
#endif // !DEBUG
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = std::sqrt(std::pow(FreqReal[i][j], (double) 2.0) + std::pow(FreqImag[i][j], (double) 2.0));

			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
		}
	}

	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}
//Complex *pVec, int len
//double *xreal, double* ximag, int n
void FT::FFT(double *xreal, double* ximag, int n)
{
	double *wreal = new double[n / 2];
	double* wimag = new double[n / 2];
	double  treal, timag, ureal, uimag, arg;
	int m, k, j, t, index1, index2;

	bitrp(xreal, ximag, n);

	// 计算 1 的前 n / 2 个 n 次方根的共轭复数 W'j = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
	arg = -2 * PI / n;
	treal = cos(arg);
	timag = sin(arg);
	wreal[0] = 1.0;
	wimag[0] = 0.0;
	for (j = 1; j < n / 2; j++)
	{
		wreal[j] = wreal[j - 1] * treal - wimag[j - 1] * timag;
		wimag[j] = wreal[j - 1] * timag + wimag[j - 1] * treal;
	}

	for (m = 2; m <= n; m *= 2)
	{
		for (k = 0; k < n; k += m)
		{
			for (j = 0; j < m / 2; j++)
			{
				index1 = k + j;
				index2 = index1 + m / 2;
				t = n * j / m;
				treal = wreal[t] * xreal[index2] - wimag[t] * ximag[index2];
				timag = wreal[t] * ximag[index2] + wimag[t] * xreal[index2];
				ureal = xreal[index1];
				uimag = ximag[index1];
				xreal[index1] = ureal + treal;
				ximag[index1] = uimag + timag;
				xreal[index2] = ureal - treal;
				ximag[index2] = uimag - timag;
			}
		}
	}

	/*
		Complex *Weights = new Complex[len];
		Complex *X = new Complex[len];
		int *pnInvBits = new int[len];
		double fixed_factor = (-2 * PI) / len;
		for (int i = 0; i < len / 2; i++) {
			double angle = i * fixed_factor;
			Weights[i].rl = cos(angle);
			Weights[i].im = sin(angle);
		}
		for (int i = len / 2; i < len; i++) {
			Weights[i].rl = -(Weights[i - len / 2].rl);
			Weights[i].im = -(Weights[i - len / 2].im);
		}

		int r = get_computation_layers(len);

		// 计算倒序位码
		int index = 0;
		for (int i = 0; i < len; i++) {
			index = 0;
			for (int m = r - 1; m >= 0; m--) {
				index += (1 && (i & (1 << m))) << (r - m - 1);
			}
			pnInvBits[i] = index;
			X[i].rl = pVec[pnInvBits[i]].rl;
			X[i].im = pVec[pnInvBits[i]].im;
		}

		// 计算快速傅里叶变换
		for (int L = 1; L <= r; L++) {
			int distance = 1 << (L - 1);
			int W = 1 << (r - L);

			int B = len >> L;
			int N = len / B;

			for (int b = 0; b < B; b++) {
				int mid = b * N;
				for (int n = 0; n < N / 2; n++) {
					int index = n + mid;
					int dist = index + distance;
					pVec[index].rl = X[index].rl + (Weights[n*W].rl * X[dist].rl - Weights[n*W].im * X[dist].im);                      // Fe + W*Fo
					pVec[index].im = X[index].im + Weights[n*W].im * X[dist].rl + Weights[n*W].rl * X[dist].im;
				}
				for (int n = N / 2; n < N; n++) {
					int index = n + mid;
					pVec[index].rl = X[index - distance].rl + Weights[n*W].rl * X[index].rl - Weights[n*W].im * X[index].im;        // Fe - W*Fo
					pVec[index].im = X[index - distance].im + Weights[n*W].im * X[index].rl + Weights[n*W].rl * X[index].im;
				}
			}

			memcpy(X, pVec, len * sizeof(Complex));
		}

		if (Weights)      delete[] Weights;
		if (X)                 delete[] X;
		if (pnInvBits)    delete[] pnInvBits;*/
}

void FT::InverseFastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;
	Complex *pVec;

	double** pFreq = new double*[M];
	for (int newcnt = 0; newcnt < M; newcnt++)
	{
		pFreq[newcnt] = new double[N]; // 傅立葉頻率陣列
	}
	for (int forzero_i = 0; forzero_i < M; forzero_i++)
	{
		for (int forzero_j = 0; forzero_j < N; forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	pVec = new Complex[N];
	double* xreal = new double[N];
	double* ximag = new double[N];
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++) {
#ifdef DEBUG
			pVec[j].rl = FreqReal[i][j];
			pVec[j].im = FreqImag[i][j];
#endif // DEBUG
#ifndef DEBUG
			xreal[j] = FreqReal[i][j];
			ximag[j] = FreqImag[i][j];
#endif // !DEBUG

		}
		//InverseFFT(pVec, M);
		InverseFFT(xreal, ximag, N);
		for (int j = 0; j < N; j++) {
#ifdef DEBUG
			FreqReal[i][j] = pVec[j].rl;
			FreqImag[i][j] = pVec[j].im;
#endif // DEBUG
#ifndef DEBUG
			FreqReal[i][j] = xreal[j];
			FreqImag[i][j] = ximag[j];
#endif // !DEBUG
		}
	}
	delete[] pVec;
	pVec = new Complex[M];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++) {
#ifdef DEBUG
			pVec[j].rl = FreqReal[j][i];
			pVec[j].im = FreqImag[j][i];
#endif // DEBUG
#ifndef DEBUG
			xreal[j] = FreqReal[j][i];
			ximag[j] = FreqImag[j][i];
#endif // !DEBUG
		}
		//InverseFFT(pVec, N);
		InverseFFT(xreal, ximag, M);
		for (int j = 0; j < M; j++) {
#ifdef DEBUG
			FreqReal[j][i] = pVec[j].rl;
			FreqImag[j][i] = pVec[j].im;
#endif // DEBUG

#ifndef DEBUG
			FreqReal[j][i] = xreal[j];
			FreqImag[j][i] = ximag[j];
#endif // !DEBUG

		}
	}

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = std::sqrt(std::pow(FreqReal[i][j], (double) 2.0) + std::pow(FreqImag[i][j], (double) 2.0));

			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j] / h;
		}
	}

	//-------------------------------------------
	for (int delcnt = 0; delcnt < M; delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[] pFreq;
}
//Complex* pVec, int len
//
void FT::InverseFFT(double *xreal, double* ximag, int n)
{
	double* wreal = new double[n / 2];
	double* wimag = new double[n / 2];
	double treal, timag, ureal, uimag, arg;
	int m, k, j, t, index1, index2;

	bitrp(xreal, ximag, n);

	// 计算 1 的前 n / 2 个 n 次方根 Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
	arg = 2 * PI / n;
	treal = cos(arg);
	timag = sin(arg);
	wreal[0] = 1.0;
	wimag[0] = 0.0;
	for (j = 1; j < n / 2; j++)
	{
		wreal[j] = wreal[j - 1] * treal - wimag[j - 1] * timag;
		wimag[j] = wreal[j - 1] * timag + wimag[j - 1] * treal;
	}

	for (m = 2; m <= n; m *= 2)
	{
		for (k = 0; k < n; k += m)
		{
			for (j = 0; j < m / 2; j++)
			{
				index1 = k + j;
				index2 = index1 + m / 2;
				t = n * j / m;
				treal = wreal[t] * xreal[index2] - wimag[t] * ximag[index2];
				timag = wreal[t] * ximag[index2] + wimag[t] * xreal[index2];
				ureal = xreal[index1];
				uimag = ximag[index1];
				xreal[index1] = ureal + treal;
				ximag[index1] = uimag + timag;
				xreal[index2] = ureal - treal;
				ximag[index2] = uimag - timag;
			}
		}
	}

	/*
	double         *W_rl = new double[len];
	double         *W_im = new double[len];
	double         *X_rl = new double[len];
	double         *X_im = new double[len];
	double         *X2_rl = new double[len];
	double         *X2_im = new double[len];

	double fixed_factor = (-2 * PI) / len;
	for (int i = 0; i < len / 2; i++) {
		double angle = i * fixed_factor;
		W_rl[i] = cos(angle);
		W_im[i] = sin(angle);
	}
	for (int i = len / 2; i < len; i++) {
		W_rl[i] = -(W_rl[i - len / 2]);
		W_im[i] = -(W_im[i - len / 2]);
	}

	// 时域数据写入X1
	for (int i = 0; i < len; i++) {
		X_rl[i] = pVec[i].rl;
		X_im[i] = pVec[i].im;
	}
	memset(X2_rl, 0, sizeof(double)*len);
	memset(X2_im, 0, sizeof(double)*len);

	int r = get_computation_layers(len);
	for (int L = r; L >= 1; L--) {
		int distance = 1 << (L - 1);
		int W = 1 << (r - L);

		int B = len >> L;
		int N = len / B;

		for (int b = 0; b < B; b++) {
			for (int n = 0; n < N / 2; n++) {
				int index = n + b * N;
				X2_rl[index] = (X_rl[index] + X_rl[index + distance]) / 2;
				X2_im[index] = (X_im[index] + X_im[index + distance]) / 2;
			}
			for (int n = N / 2; n < N; n++) {
				int index = n + b * N;
				X2_rl[index] = (X_rl[index] - X_rl[index - distance]) / 2;           // 需要再除以W_rl[n*W]
				X2_im[index] = (X_im[index] - X_im[index - distance]) / 2;
				double square = W_rl[n*W] * W_rl[n*W] + W_im[n*W] * W_im[n*W];          // c^2+d^2
				double part1 = X2_rl[index] * W_rl[n*W] + X2_im[index] * W_im[n*W];         // a*c+b*d
				double part2 = X2_im[index] * W_rl[n*W] - X2_rl[index] * W_im[n*W];          // b*c-a*d
				if (square > 0) {
					X2_rl[index] = part1 / square;
					X2_im[index] = part2 / square;
				}
			}
		}
		memcpy(X_rl, X2_rl, sizeof(double)*len);
		memcpy(X_im, X2_im, sizeof(double)*len);
	}

	// 位码倒序
	int index = 0;
	for (int i = 0; i < len; i++) {
		index = 0;
		for (int m = r - 1; m >= 0; m--) {
			index += (1 && (i & (1 << m))) << (r - m - 1);
		}
		pVec[i].rl = X_rl[index];
		pVec[i].im = X_im[index];
		//sprintf(msg + 6, "X_rl[i]: %lf, %lf,  index: %d", out_rl[i], out_im[i], index);
		//OutputDebugStringA(msg);
	}

	delete[] W_rl;
	delete[] W_im;
	delete[] X_rl;
	delete[] X_im;
	delete[] X2_rl;
	delete[] X2_im;*/
}


void FT::LowpassFilter(double** Real, double** Img, int** output, double d, int h, int w)
{
	int m = h / 2;
	int n = w / 2;
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			if ((i - m)*(i - m) + (j - n)*(j - n) > d)
			{
				Real[i][j] = 0;
				Img[i][j] = 0;
			}
			output[i][j] = std::sqrt(std::pow(Real[i][j], (double) 2.0) + std::pow(Img[i][j], (double) 2.0));
		}
	}
}

void FT::HighpassFilter(double** Real, double** Img, int** output, double d, int h, int w)
{
	int m = h / 2;
	int n = w / 2;
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			if ((i - m)*(i - m) + (j - n)*(j - n) < d)
			{
				Real[i][j] = 0;
				Img[i][j] = 0;
			}
			output[i][j] = std::sqrt(std::pow(Real[i][j], (double) 2.0) + std::pow(Img[i][j], (double) 2.0));
		}
	}

}


int get_computation_layers(int num) {
	int nLayers = 0;
	int len = num;
	if (len == 2)
		return 1;
	while (true) {
		len = len / 2;
		nLayers++;
		if (len == 2)
			return nLayers + 1;
		if (len < 1)
			return -1;
	}
}

Complex::Complex() {
	rl = 0;
	im = 0;
}

void swap(double &a, double &b) {
	double tmp = a;
	a = b;
	b = tmp;
}

void bitrp(double* real, double*imag, int n) {
	int i, j, a, b, p;
	for (i = 1, p = 0; i < n; i *= 2) {
		p++;
	}
	for (i = 0; i < n; i++) {
		a = i;
		b = 0;
		for (j = 0; j < p; j++) {
			;			b = (b << 1) + (a & 1);
			a >>= 1;
		}
		if (b > i) {
			swap(real[i], real[b]);
			swap(imag[i], imag[b]);
		}
	}
}