//Computes the sum of a series.Which series is it,though...
//tip:don't set both compile time flags,though apparently ChatGPT thinks
//i could use a "mixed partial derivative" or smth...
//derivatives not tested..yet.
double FUNC_NAME(double coeff[], double th, double ph, int level) {
	double cos_th = cos(th), nsin_th = -sin(th);
	double rt2 = 0, zero = 0;
#ifdef DERIVATIVE_PH
	double cosph[MAX_SH_LEVEL + 1] = { 0 };
	double sinph[MAX_SH_LEVEL + 1] = { 0 };
#else
	double cosph[MAX_SH_LEVEL + 1] = { 1 };
	double sinph[MAX_SH_LEVEL + 1] = { 0 };
#endif
	double expct[MAX_SH_LEVEL + 1] = { 1 };
	double expox[MAX_SH_LEVEL + 1] = { 1 };
#ifdef DERIVATIVE_TH
	double dexpct[MAX_SH_LEVEL + 1] = { 0 };
	double dexpox[MAX_SH_LEVEL + 1] = { 0 };
#endif
	double preve = 1, prevx = 1;
	for (int m = 1; m <= level; ++m) {
#ifdef DERIVATIVE_PH
		cosph[m] = -m * sin(m * ph);
		sinph[m] = m * cos(m * ph);
#else
		cosph[m] = cos(m * ph);
		sinph[m] = sin(m * ph);
#endif
#ifdef DERIVATIVE_TH
		dexpct[m] = preve * nsin_th * m;
		dexpox[m] = -prevx * cos_th * m;
#endif
		preve = expct[m] = preve * cos_th;
		prevx = expox[m] = prevx * nsin_th;
	}
	double* lpc = lpoly_coeff, * lpn = lpoly_norm, * tc = coeff,rest;
	for (int l = 0; l <= level; l++) {

#ifdef DERIVATIVE_PH
		//zero is a constant,constants have 0 derivative
		//still have to advance the iterators tho
		lpc += l + 1;
		lpn++;
#else
		rest = 0;
		for (int i = 0; i < l + 1; i++) {
#ifdef DERIVATIVE_TH
			rest += *(lpc++) * dexpct[i];
#else
			rest += *(lpc++) * expct[i];
#endif
		}
		zero += *(lpn++) * *tc * rest;
#endif
		for (int abs_m = 1; abs_m <= l; abs_m++) {
			rest = 0;
			double E = expox[abs_m];
#ifdef DERIVATIVE_TH
			double dE = dexpox[abs_m];
#endif
			for (int i = 0; i < l + 1; i++) {
#ifdef DERIVATIVE_TH
				rest += *(lpc++) * (expct[i] * dE + dexpct[i] * E);//product rule.annoying here...
#else
				rest += *(lpc++) * expct[i] * E;
#endif
			}
			double nplm = rest * *(lpn++);
			rt2 += nplm * (tc[abs_m] * cosph[abs_m] + tc[-abs_m] * sinph[abs_m]);
		}
		tc += 2 * l + 2;
	}
	return numbers::sqrt2 * rt2 + zero;
}