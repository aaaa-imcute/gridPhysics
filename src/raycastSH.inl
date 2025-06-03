//raymarching,but different returns.
//no longer supports gradient calculation
//as sphere marching is derivative-free
bool FUNC_NAME(
#ifdef RET_CONTACT
	glm::dvec3& result, 
#endif
#ifdef RET_LENGTH
	double& length, 
#endif
	double coeff[], 
	double lipschitz,
	double maxHeight,
	glm::dvec3 origin, 
	glm::dvec3 normal
) {
	double EP = 1e-8;
	double t = 0, dt = 0;
	double cor = sqrt(1 + lipschitz * lipschitz);
	//"correction factor",corrects the naive guess of SDF(v)=v-r(v)
	//derived from geometric construction of a circle and a chord(distance from chord to center,
	//chord represents worse case according to lipschitz)
	double r, th, ph;
	glm::dvec3 ray;
	do {
		ray = origin + normal * t;
		r = glm::length(ray);
		if (r > 1 + maxHeight && glm::dot(ray, normal) > 0)return false;//out of range
		th = acos(ray.y / r);
		ph = atan2(ray.z, ray.x);
		dt = (r - 1 - get_terrain_height(coeff, th, ph, MAX_SH_LEVEL)) / cor;
		t += dt;
	} while (dt > EP);//don't need the abs() here,ray marching doesn't go back like NR does.
	//if (abs(r - 1 - get_terrain_height(coeff, th, ph, MAX_SH_LEVEL))>1e-4) throw "oh no";
#ifdef RET_CONTACT
	result = origin + normal * t;
#endif
#ifdef RET_LENGTH
	length = t;
#endif
	return true;
}