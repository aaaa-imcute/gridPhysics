//on default detail,points,and terrain level settings,
//nearly a million calls to SH are made when a planet is rendered
//(128*64*11**2)
//this makes it a critical performance bottleneck
#include "ofMain.h"

glm::dvec2 sphericalCoordinates(glm::dvec3 v) {
	glm::dvec3 t = glm::normalize(v);
	return glm::dvec2(atan2(t.z, t.x), acos(t.y));
}
double fac_t[16] = {
	1,
	1,
	2,
	6,
	24,
	120,
	720,
	5040,
	40320,
	362880,
	3628800,
	39916800,
	479001600,
	6227020800,
	87178291200,
	1307674368000
};
double dfac_t[26] = {
	1,
	1,
	2,
	3,
	8,
	15,
	48,
	105,
	384,
	945,
	3840,
	10395,
	46080,
	135135,
	645120,
	2027025,
	10321920,
	34459425,
	185794560,
	654729075,
	3715891200,
	13749310575,
	81749606400,
	316234143225,
	1961990553600,
	7905853580625,
};
//okay i thought i would use these tables but apparently not
double fac_ratio(int x, int y) {
	if (x == y) return 1;
	if (x < 16 && y < 16)return fac_t[x] / fac_t[y];
	//yeah yeah i could use a double lookup table but whatever
	double ret=1;
	for (int i = y; i < x; i++) {
		ret *= i + 1;
	}
	return ret;
}
double qpow(double x, int p) {
	
	double r = 1;
	for (int i = 0; i < p; i++) {
		r *= x;
	}
	return r;
	
	/*
	if (p == 0)return 1;
	double s = qpow(x, p / 2);
	if (p & 1)return s * s * x;//x&1 means x is odd
	return s * s;
	*/
	//do not use recursive power for small exponents
	//iterative solution gives ~1 more FPS(27 vs 28)
}
double double_factorial(int n) {
	//sometimes i don't need integer precision and so not those tables
	if (n == 0 || n == -1) return 1.0;
	if (n & 1) {
		int k = (n - 1) / 2;
		return tgamma(n + 1) / (pow(2, k) * tgamma(k + 1));
	}
	return pow(2, n / 2) * tgamma(n / 2 + 1);
}
vector<vector<vector<double>>> lpoly_coeff = {{{1}},{{0,1},{1,0}}};
vector<vector<double>> lpoly_norm = {
	{sqrt(1 / (4.0 * PI))},
	{sqrt(3 / (4.0 * PI)),sqrt(3 / (8.0 * PI))}
};//unfortunately we are not calculating these in the compute function so i have to
void compute_legendre_coeff(int s) {
	int l = lpoly_coeff.size();
	while (l <= s) {
		vector<vector<double>> temp;
		for (int m = 0; m < l-1; m++) {
			vector<double> t;
			vector<double>& lft = lpoly_coeff[l - 1][m], rgt = lpoly_coeff[l - 2][m];
			for (int d = 0; d <= l; d++) {
				double left = (d == 0) ? 0 : (2 * l - 1) * lft[d - 1];
				double right = (d >= l - 1) ? 0 : (l + m - 1) * rgt[d];
				t.emplace_back((left - right) / (l - m));
			}
			temp.emplace_back(t);
		}
		vector<double> last0(l + 1, 0), last1(l + 1, 0);
		last0[0]=double_factorial(2 * l - 1);
		last1[1] = last0[0];
		temp.emplace_back(last1);
		temp.emplace_back(last0);
		lpoly_coeff.emplace_back(temp);
		vector<double> tn;
		for (int m = 0; m <= l; m++) {
			tn.emplace_back(sqrt((2.0 * l + 1.0) / (4.0 * PI) / fac_ratio(l + m, l - m)));
		}
		lpoly_norm.emplace_back(tn);
		l++;
	}
}
double assoclegendre(int l, int m, double x) {
	//compute_legendre_coeff(l);
	double p = qpow(sqrt(1 - x * x), m);
	double rest = 0,xx=1;
	const vector<double>& coeff = lpoly_coeff[l][m];
	for (double c:coeff) {
		rest += c*xx;
		xx *= x;
	}
	return p * rest;
}
const double SQRT_2 = 1.4142135623730951;
double sphericalHarmonics(int l, int m, double th, double ph) {
	int abs_m = abs(m);
	double plm = ((abs_m & 1) * -2 + 1) * assoclegendre(l, abs_m, cos(th));
	plm *= lpoly_norm[l][abs_m];
	if (m == 0) return plm;
	if (m > 0) return SQRT_2 * plm * cos(m * ph);
	return -SQRT_2 * plm * sin(m * ph);
}