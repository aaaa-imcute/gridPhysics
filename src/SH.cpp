//on default detail,points,and terrain level settings,
//nearly a million calls to SH are made when a planet is rendered
//(128*64*11**2)
//this makes it a critical performance bottleneck
#include "ofMain.h"
const int MAX_SH_LEVEL = 11;
#define LC_FIND_INDEX(l,m) ((l)*((l)+1)*(2*(l)+1)/6+(l+1)*(m))
#define LN_FIND_INDEX(l) ((l)*((l)+1)/2)
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
double double_factorial(int n) {
	//sometimes i don't need integer precision and so not those tables
	if (n == 0 || n == -1) return 1.0;
	if (n & 1) {
		int k = (n - 1) / 2;
		return tgamma(n + 1) / (pow(2, k) * tgamma(k + 1));
	}
	return pow(2, n / 2) * tgamma(n / 2 + 1);
}
double SQRT_2 = sqrt(2.0);
double lpoly_coeff[LC_FIND_INDEX(MAX_SH_LEVEL + 1, 0)] = { 1,0,1,1,0 };
double lpoly_norm[LN_FIND_INDEX(MAX_SH_LEVEL + 1)] = {
	sqrt(1 / (4.0 * PI)),
	sqrt(3 / (4.0 * PI)),
	sqrt(3 / (8.0 * PI))
};//unfortunately we are not calculating these in the compute function so i have to
void compute_legendre_coeff() {
	int l = 2, index = 5, nndex = 3;
	while (l <= MAX_SH_LEVEL) {
		for (int m = 0; m < l-1; m++) {
			double *lft = &lpoly_coeff[LC_FIND_INDEX(l - 1, m)], *rgt = &lpoly_coeff[LC_FIND_INDEX(l - 2, m)];
			for (int d = 0; d <= l; d++) {
				double left = (d == 0) ? 0 : (2 * l - 1) * lft[d-1];
				double right = (d >= l - 1) ? 0 : (l + m - 1) * rgt[d];
				lpoly_coeff[index++]=(left - right) / (l - m);
			}
		}
		double c = double_factorial(2 * l - 1);
		lpoly_coeff[index++] = 0;
		lpoly_coeff[index++] = c;
		for (int i = 0; i < l - 1; i++)lpoly_coeff[index++] = 0;
		lpoly_coeff[index++] = c;
		for (int i = 0; i < l; i++)lpoly_coeff[index++] = 0;
		for (int m = 0; m <= l; m++) {
			lpoly_norm[nndex++]=sqrt((2.0 * l + 1.0) / (4.0 * PI) / fac_ratio(l + m, l - m));
		}
		l++;
	}
}
//using macros we define three variations of the same function
#define FUNC_NAME get_terrain_height
#include "del_sh.inl"
#undef FUNC_NAME
#define FUNC_NAME get_terrain_height_dth
#define DERIVATIVE_TH
#include "del_sh.inl"
#undef DERIVATIVE_TH
#undef FUNC_NAME
#define FUNC_NAME get_terrain_height_dph
#define DERIVATIVE_PH
#include "del_sh.inl"
#undef DERIVATIVE_PH
#undef FUNC_NAME
double dRaycastingSH(double coeff[], glm::dvec3 origin, glm::dvec3 normal, double t, int level) {
	//TODO testing
	double EP = 1e-8;
	glm::dvec3 ray = origin + normal * t;
	double r = glm::length(ray), th = acos(ray.y / r), ph = atan2(ray.z, ray.x);
	double dr = glm::dot(normal, ray) / r;
	double dth = -(r * normal.y - ray.y * dr) / (r * r * sqrt(1 - (ray.y / r) * (ray.y / r)));
	double dph = ray.x * ray.x + ray.z * ray.z;
	if (abs(dph) < 1e-8) dph = 0;
	else dph = (ray.x * normal.z - ray.z * normal.x) / dph;
	return dr - get_terrain_height_dth(coeff, th, ph, level) * dth - get_terrain_height_dph(coeff, th, ph, level) * dph;
}
bool raycast_SH(glm::dvec3& result,double& length,double coeff[], glm::dvec3 origin, glm::dvec3 normal, int level) {
	//TODO:Fix
	double EP = 1e-8;
	double t = 0;//hope that we start outside the sphere ig
	double dt = 0;//uninitialized memory warning blah blah
	const int MAX_ITER = 20;
	int i = 0;
	do {
		glm::dvec3 ray = origin + normal * t;
		double rad = glm::length(ray), th = acos(ray.y / rad), ph = atan2(ray.z, ray.x);
		double f = rad - get_terrain_height(coeff, th, ph, level)-1;
		double df = dRaycastingSH(coeff, origin, normal, t, level);
		dt = f / df;
		t -= dt;
		if (isnan(dt) || isinf(dt) || abs(dt) > 1.0 / EP || t<0 || ++i>MAX_ITER)return false;
	} while (abs(dt) > EP);
	result = origin + normal * t;
	length = t;
	return true;
}