//the same NR method but different return
bool FUNC_NAME(
#ifdef RET_CONTACT
	glm::dvec3& result, 
#endif
#ifdef RET_LENGTH
	double& length, 
#endif
#ifdef RET_GRADIENT
	glm::dvec2& gradient,
#endif
	double coeff[], 
	glm::dvec3 origin, 
	glm::dvec3 normal
) {
	double EP = 1e-8;
	double t = 0;//hope that we start outside the sphere ig
	double dt = 0;//uninitialized memory warning blah blah
	const int MAX_ITER = 20;
	int i = 0;
	double r, th, ph, dr, dth, dph, f, dsdth, dsdph, df;
	do {
		glm::dvec3 ray = origin + normal * t;
		r = glm::length(ray), th = acos(ray.y / r), ph = atan2(ray.z, ray.x);
		dr = glm::dot(normal, ray) / r;
		dth = -(r * normal.y - ray.y * dr) / (r * r * sqrt(1 - (ray.y / r) * (ray.y / r)));
		dph = ray.x * ray.x + ray.z * ray.z;
		if (abs(dph) < 1e-8) dph = 0;
		else dph = (ray.x * normal.z - ray.z * normal.x) / dph;
		//TODO:if raycasting is the performance bottleneck,modify del_sh.inl to allow for
		//multiple calculations at the same time(just like this function)
		//nothing stopping dth mixing with anything else rn,and ph only needs new tables
		//ofc there would be new variables to store the results
		f = r - get_terrain_height(coeff, th, ph, MAX_SH_LEVEL) - 1;
		dsdth = get_terrain_height_dth(coeff, th, ph, MAX_SH_LEVEL);
		dsdph = get_terrain_height_dph(coeff, th, ph, MAX_SH_LEVEL);
		df = dr - dsdth * dth - dsdph * dph;
		dt = f / df;
		t -= dt;
		if (isnan(dt) || isinf(dt) || abs(dt) > 1.0 / EP || t<0 || ++i>MAX_ITER)return false;
	} while (abs(dt) > EP);
#ifdef RET_CONTACT
	result = origin + normal * t;
#endif
#ifdef RET_LENGTH
	length = t;
#endif
#ifdef RET_GRADIENT
	gradient = glm::dvec2(dsdph,dsdth);
#endif
	return true;
}