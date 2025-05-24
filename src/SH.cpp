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
double sphericalHarmonics(int l, int m, double th, double ph) {
	int abs_m = abs(m);
	double norm = sqrt((2.0 * l + 1.0) / (4.0 * PI) / fac_ratio(l + abs_m, l - abs_m));
	double plm = ((abs_m % 2 == 0) ? 1 : -1) * assoc_legendre(l, abs_m, cos(th));
	if (m == 0) return norm * plm;
	if (m > 0) return sqrt(2.0) * norm * plm * cos(m * ph);
	else       return sqrt(2.0) * norm * plm * sin(-m * ph);
}