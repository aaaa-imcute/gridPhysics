#pragma once

#include "ofMain.h"
#include "SH.cpp"
#include <variant>
double ffmod(double x, double y) {
	return fmod(fmod(x, y) + y, y);
}
unordered_map<int, bool> keys;
unordered_map<int, bool> mouse;
glm::dvec2 untransform2D(glm::dvec2 vec) {
	ofRectangle viewport = ofGetCurrentViewport();
	glm::vec3 projected = glm::project(
		glm::vec3(0, 0, 0), // or wherever you draw your object
		ofGetCurrentMatrix(OF_MATRIX_MODELVIEW),
		ofGetCurrentMatrix(OF_MATRIX_PROJECTION),
		glm::vec4(viewport.x, viewport.y, viewport.width, viewport.height)
	);
	glm::vec3 windowCoords(vec.x, viewport.height - (vec.y - viewport.y), projected.z);
	glm::dvec3 worldPos = glm::unProject(
		windowCoords,
		ofGetCurrentMatrix(OF_MATRIX_MODELVIEW),
		ofGetCurrentMatrix(OF_MATRIX_PROJECTION),
		glm::vec4(viewport.x, viewport.y, viewport.width, viewport.height)
	);
	return glm::dvec2(worldPos.x, worldPos.y);
}
glm::dvec2 mousePos, pmousePos;
bool mousePressedOnLastFrame = false;
bool dontDragMap = false;
int keyPressedOnLastFrame = 0;
bool dontDetectKeys = false;
ofTrueTypeFont font;
class easyCam : public ofCamera {
public:
	glm::dvec3 angle;//pitch yaw radius
	glm::dvec3 pos;
	easyCam() {
		angle = { 0,0,1000 };
		updateOrientation();
	}
	void update() {
		if (!mouse[2])return;
		glm::dvec2 rel = mousePos - pmousePos;
		angle.x += -rel.y * PI / ofGetHeight();
		angle.y += -rel.x * PI * 2 / ofGetWidth();
		angle.x = ofClamp(angle.x, -PI / 2, PI / 2);
		updateOrientation();
	}
	void mouseScrolled(double sx, double sy) {
		angle.z += -sy * 50;
		updateOrientation();
	}
	void updateOrientation() {
		glm::dquat q;
		q = glm::rotate(q, angle.y, { 0,1,0 });
		q = glm::rotate(q, angle.x, { 1,0,0 });
		glm::dvec3 p = { 0,0,angle.z };
		p = glm::rotate(q, p);
		setGlobalOrientation(q);
		setGlobalPosition(pos + p);
	}
};
glm::dvec2 mapPos;//sceneDisplayed==0||sceneDisplayed==2
double mapScale = 16;
int dm2_layer = 0;
easyCam camera1,camera3;//sceneDisplayed==1,sceneDisplayed==3
shared_ptr<ofLight> sunLight;//pointer,otherwise error
//TODO:move this into the any planet without an orbit(i.e the sun)
namespace std {
	template <>
	struct hash<glm::dvec3> {
		auto operator()(const glm::dvec3& x) const -> size_t {
			return hash<glm::i32>{}(x.x) ^ (hash<glm::i32>{}(x.y) << 1) ^ (hash<glm::i32>{}(x.z) << 2);
		}
	};
}//make it so that this can be the key of a unordered map(used in class PhysicsGrid)
glm::dvec3 orbitPos;//sceneDisplayed==1
double orbitScale = 1.0 / 600;
vector<string> faceNames = { "right","left","top","bottom","front","back" };
unordered_map<string, shared_ptr<ofImage>> imageCache;
shared_ptr<ofImage> loadImage(string img) {
	auto result = imageCache.find(img);
	shared_ptr<ofImage> ret;
	if (result == imageCache.end()) {
		ret = make_shared<ofImage>();
		if (!ret->load(img))throw "image " + img + " not loaded correctly";
		ret->getTexture().setTextureMinMagFilter(GL_NEAREST, GL_NEAREST);
		imageCache[img] = ret;
	}
	else {
		ret = result->second;
	}
	return ret;
}
void drawImage(string img, float x, float y) {
	shared_ptr<ofImage> temp = loadImage(img);
	ofSetColor(255, 255, 255);
	temp->draw(x, y);
}
ofFbo atlas;
unordered_map<string, glm::dvec2> atlasMap;
bool atlasReady = false;
ofFbo planetAtlas;//map of planet colors per altitude
bool planetAtlasReady = false;
constexpr double DM3_SCALE = 256;
constexpr int ATLAS_SIZE = 4096;//needs to be divisible by 16
//dk why it was 256 before,mb i thought of it as the number of textures on a row of the atlas
//instead of how it was actually used,the number of pixels.
void walkTextures(string path, int& x, int& y) {
	ofDirectory dir(path);
	if (!dir.isDirectory()) {
		shared_ptr<ofImage> temp = loadImage(path);//fixed dumb getOriginalDirectory call and string manipulation
		atlas.begin();
		temp->draw(x, y);
		atlas.end();
		filesystem::path relPath = filesystem::relative(path, "textures");
		relPath.replace_extension();
		atlasMap[relPath.string()] = { (float)x, (float)y };
		x += 16;
		if (x == ATLAS_SIZE) {
			x = 0;
			y += 16;
			if (y == ATLAS_SIZE)throw "Too many textures";
		}
		return;
	}
	dir.listDir();
	for (int i = 0; i < dir.size(); i++) {
		walkTextures(dir.getPath(i), x, y);
	}
}
void createAtlas() {
	if (atlasReady)return;
	atlasReady = true;
	ofFbo::Settings fboSettings;
	fboSettings.width = ATLAS_SIZE;
	fboSettings.height = ATLAS_SIZE;
	fboSettings.internalformat = GL_RGBA;
	fboSettings.numSamples = 0;
	atlas.allocate(fboSettings);
	atlas.getTexture().setTextureMinMagFilter(GL_NEAREST, GL_NEAREST);
	atlas.begin();
	ofClear(0);
	atlas.end();
	int x = 0, y = 0;
	walkTextures(".\\textures\\", x, y);
}
void drawLine(glm::dvec2 p1, glm::dvec2 p2, float width) {//replacement method for line width
	//expects color already set
	ofPushMatrix();
	ofTranslate(p1);
	ofRotateRad(atan2f(p2.y - p1.y, p2.x - p1.x));
	ofDrawRectangle(-glm::dvec2(0, width), glm::length(p1 - p2), width * 2);
	ofPopMatrix();
}
void drawFaceNormal2D(glm::dvec3 v, glm::dvec2 pos) {//sceneDisplay==2
	//expects color already set
	glm::dvec2 d(-v.z, -v.x);
	if (abs(v.y) < 0.01) {
		drawLine(pos, pos + d * 16, 1);
	}
	else if (v.y > 0) {
		ofDrawCircle(pos, 2);
	}
	else if (v.y < 0) {
		drawLine(pos + glm::dvec2(-4, -4), pos + glm::dvec2(4, 4), 1);
		drawLine(pos + glm::dvec2(-4, 4), pos + glm::dvec2(4, -4), 1);
	}
}
constexpr double FONT_SCALE = 0.25;
void drawText(string t, double x, double y) {
	ofPushMatrix();
	ofTranslate(x, y);
	ofScale(FONT_SCALE, FONT_SCALE);
	font.drawStringAsShapes(t, 0, 0);
	ofPopMatrix();
}
void drawText(string t, glm::dvec2 pos) {
	drawText(t, pos.x, pos.y);
}
double textWidth(string t) {
	return font.stringWidth(t) * FONT_SCALE;
}
double textHeight(string t) {
	return font.stringHeight(t) * FONT_SCALE;
}
void getTriangleMesh(ofMesh& m, glm::dvec3 a, glm::dvec3 b, glm::dvec3 c) {
	double scale = DM3_SCALE;
	int baseIndex = (int)m.getNumVertices();//compiler unhappy about "size_t to int is downgrading conversion" or smth
	glm::dvec3 normal = glm::cross(c - b, a - b);
	m.addVertex(a * scale);
	m.addVertex(b * scale);
	m.addVertex(c * scale);
	for (int i = 0; i < 3; i++)m.addNormal(normal);
	m.addTriangle(baseIndex, baseIndex + 1, baseIndex + 2);
}
void getTriangleMesh(ofMesh& m, glm::dvec3 a, glm::dvec3 b, glm::dvec3 c, glm::dvec2 aa, glm::dvec2 ab, glm::dvec2 ac) {
	getTriangleMesh(m, a, b, c);
	m.addTexCoord(aa);
	m.addTexCoord(ab);
	m.addTexCoord(ac);
}
void getRectangleMesh(ofMesh& m, glm::dvec3 pos, glm::dvec3 width, glm::dvec3 height, glm::dvec2 apos, glm::dvec2 aw, glm::dvec2 ah) {
	getTriangleMesh(m, pos, pos + height, pos + width, apos, apos + ah, apos + aw);
	getTriangleMesh(m, pos + width + height, pos + width, pos + height, apos + aw + ah, apos + aw, apos + ah);
}
void drawPlaneWithNormal(glm::dvec3 normal, float distance, float size) {
	glm::vec3 unitNormal = glm::normalize(normal);
	glm::vec3 center = -unitNormal * distance;

	// Create a transformation that rotates Z-axis to the desired normal
	glm::vec3 defaultNormal = glm::vec3(0, 0, 1);
	glm::quat rotation = glm::rotation(defaultNormal, unitNormal);

	// Setup the ofNode for transformation
	ofNode node;
	node.setPosition(center);
	node.setOrientation(rotation);

	// Create the plane primitive
	ofPlanePrimitive plane;
	plane.set(size, size, 2, 2);
	plane.setPosition(0, 0, 0); // Will be drawn relative to node

	// Apply transformation
	ofPushMatrix();
	ofMultMatrix(node.getGlobalTransformMatrix());

	ofSetColor(150, 200, 255, 100);
	plane.draw();
	ofPopMatrix();
	// Draw normal vector from center
	ofSetColor(255, 0, 0);
	ofDrawLine(center, center + unitNormal * (size / 2));
}
void make_orthonormal_basis(const glm::dvec3& N, glm::dvec3& T, glm::dvec3& B) {
	if (fabs(N.x) > fabs(N.z)) {
		T = glm::normalize(glm::dvec3(-N.y, N.x, 0.0));
	}
	else {
		T = glm::normalize(glm::dvec3(0.0, -N.z, N.y));
	}
	B = glm::cross(N, T);
}
class Terrain {
public:
	double coeff[(MAX_SH_LEVEL + 1) * (MAX_SH_LEVEL + 1)] = {};
	bool generated = false;
	double maxHeight = 0, lipschitz = 0;
	void generate(unsigned int seed, double A, double a) {
		//if a planet isn't generated,it will render as a sphere(gas giant),so don't if it actually is!
		//seed fineness amplitude(ratio of mountains to planet radius) smoothness(more than 2)
		mt19937 rng(seed);
		int index = 0;
		maxHeight = 0;
		lipschitz = 0;
		for (int l = 0; l <= MAX_SH_LEVEL; l++) {
			double var = A / pow(1 + l, a);
			normal_distribution<double> d(0, var);
			for (int m = -l; m <= l; m++) {
				double height = d(rng);
				coeff[index++] = height;
				maxHeight += height * height;
				lipschitz += l * l * height * height;
			}
		}
		maxHeight *= 2 * sqrt(log(MAX_SH_LEVEL + 1));//blah blah blah extreme value of gaussian
		lipschitz *= 2 * sqrt(log(MAX_SH_LEVEL + 1));//this has a 1/((L+1)^2*sqrt(8pilog(L+1))) chance to fail
		//current settings it is about 1/1138
		generated = true;
	}
	void generateTest() {
		//dummy function for generating very round planets(for testing)
		int index = 0;
		maxHeight = 0;
		lipschitz = 0;
		for (int l = 0; l <= MAX_SH_LEVEL; l++) {
			for (int m = -l; m <= l; m++) {
				coeff[index++] = 0;
			}
		}
		generated = true;
	}
	double get(double th, double ph, int level = -1) {
		//returns relative altitude not distance to planet center,need to +1 then *radius
		if (level == -1)level = MAX_SH_LEVEL;
		return get_terrain_height(coeff, th, ph, MAX_SH_LEVEL);
	}
	double get(glm::dvec3 v, int level = -1) {
		glm::dvec2 c = sphericalCoordinates(v);
		return get(c.y, c.x, level);
	}
	void mesh(ofMesh& m, glm::dvec3 ref, int planetI, bool centered = false) {
		//ref is relative to the planet,scaled down by the radius
		//(so is the resulting mesh)
		//TODO:centered makes the mesh centered around the camera
		double ter = 1 + get(ref);
		int points = 64;
		double r = glm::length(ref);
		glm::dvec3 norm = ref / r;
		if (!generated || r > 16) {
			ofSpherePrimitive sphere;
			sphere.setRadius(256);
			sphere.setResolution(points);
			vector<glm::vec2> texCoords(m.getNumVertices(), { ATLAS_SIZE / 2,planetI });
			m.clearTexCoords();
			m.addTexCoords(texCoords);
			m = sphere.getMesh();
			return;
		}
		double limit = PI / 2 + 1e-5;//so apparently some calculations might be off by one if
		//limit is divisible by detail
		double detail = PI / 256;//now it is relative to the camera not the planet
		int length = ceil(limit / detail) + 1;
		vector<glm::dvec3> vertices(points * length);
		vector<int> texCoords(points * length);
		int i = 0, k = 0;
		glm::dvec3 T, B;
		make_orthonormal_basis(-norm, T, B);
		for (double ph = 0;; ph += detail) {
			if (ph > limit)ph = limit;
			//yep the angle convention thing is backwards here wrt the
			//SH math,oops!(th came first and represented angles on a 2d system hence the name)
			double h = r * cos(ph) - sqrt(ter * ter - r * r * sin(ph) * sin(ph));
			double dist = h * cos(ph);
			double radius = h * sin(ph);
			double end = ((i++) & 1) * PI / points;
			double th = end;
			for (int j = 0; j < points; j++) {
				glm::dvec3 V = norm * (r - dist);
				V += (cos(th) * T + sin(th) * B) * radius;
				V /= ter;
				if (abs(glm::length(V) - 1) > 0.001)V = -norm;//remove problem dots
				//mb i could optimize it if i had more brain capacity for all the geometry
				//and do this calculation in the outer loop
				//but whatever
				double h = get(V);
				V *= h + 1;//max(h + 1,1.0);//temporary;replace with water rendering with actual option in planet
				if (centered)V -= ref;
				vertices[k] = V;//faster than push_back,but i think the SH is the actual performance issue
				//if (h > maxHeight * 1.1)throw "way off";
				texCoords[k] = (ofClamp(h / maxHeight, -1, 1) + 1) / 2 * ATLAS_SIZE;
				k++;
				th += 2 * PI / points;
			}
			if (ph == limit)break;
		}
		for (int i = 0; i < length - 1; i++) {
			for (int j = 0; j < points; j++) {
				int a, b, c;
				a = i * points + j;
				b = i * points + ffmod(j + 1, points);
				c = (i + 1) * points + ffmod(j + i % 2, points);
				//if the current layer is even then +1 otherwise don't +1
				//this is for the layer shifting thing
				if (i != 0)getTriangleMesh(
					m,
					vertices[a], vertices[b], vertices[c],
					{ texCoords[a], planetI },
					{ texCoords[b], planetI },
					{ texCoords[c], planetI }
				);
				a = i * points + j;
				b = (i + 1) * points + ffmod(j + i % 2, points);
				c = (i + 1) * points + ffmod(j + i % 2 - 1, points);
				getTriangleMesh(
					m,
					vertices[a], vertices[b], vertices[c],
					{ texCoords[a], planetI },
					{ texCoords[b], planetI },
					{ texCoords[c], planetI }
				);
			}
		}
	}
	glm::dvec3 getSurfaceNormal(glm::dvec3 pos,double& H) {
		double EP = 1e-8;
		double th, ph;
		glm::dvec3 rb, thb, phb, sn;
		double alt, rho;
		alt = glm::length(pos);
		th = acos(pos.y / alt);
		ph = atan2(pos.z, pos.x);
		if (ph < 0.0f)
			ph += glm::two_pi<float>();
		rho = sqrt(pos.x * pos.x + pos.z * pos.z);
		rb = glm::dvec3(
			sin(th) * cos(ph),
			cos(th),
			sin(th) * sin(ph)
		);
		thb = glm::dvec3(
			cos(th) * cos(ph),
			-sin(th),
			cos(th) * sin(ph)
		);
		phb = glm::dvec3(
			-sin(ph),
			0.0,
			cos(ph)
		);
		H = 1 + get_terrain_height(coeff, th, ph, MAX_SH_LEVEL);
		sn = rb -
			thb * get_terrain_height_dth(coeff, th, ph, MAX_SH_LEVEL) / H -
			phb * get_terrain_height_dph(coeff, th, ph, MAX_SH_LEVEL) / H / sin(th);
		sn = glm::normalize(sn);
		return sn;
	}
	glm::dvec3 getSurfaceNormal(glm::dvec3 pos) {
		double H;
		return getSurfaceNormal(pos, H);
	}
};
class Planet;
class OrbitalElements;
constexpr double ATMO_TOL = 0.001;
class Planet {
public:
	shared_ptr<OrbitalElements> o;
	double gravity;//m^3/s^2
	double radius;
	pair<double, double> atmos;//sea level,scale height(unit m,"how much height does it take for pressure to drop e times")
	function<ofColor(double)> color;
	string name;
	Terrain terrain;
	of3dPrimitive brush;
	ofMesh mesh;
	Planet(string n, double g, double r, double p0, double ph);
	Planet(string n, double g, double r, double p0, double ph, shared_ptr<OrbitalElements> o1);
	glm::dvec3 pos(double t);
	glm::dvec3 vel(double t);
	glm::dvec3 apos(double t);
	glm::dvec3 avel(double t);
	double SOI();
	double atmoPressure(double h);
	double atmoHeight();
	void displayMode1(double t, int planetI);
	void displayMode3(glm::dvec3 shipPos, int planetI);
};
class OrbitalElements {
public:
	double a, e, i, o, w, v, periapsis;
	//o is LAN(capital omega)
	//v is mean anomaly at t=0 NOT TRUE ANOMALY,unless if the orbit is parabolic
	//then a is infinity,e is 1,periapsis is used,mean motion is 1,
	//and v is how much time have passed since periapsis at t=0
	//(so that the equation for mean anomaly still makes sense)
	shared_ptr<Planet> p;
	OrbitalElements(shared_ptr<Planet> p1, glm::dvec3 r, glm::dvec3 V, double t) {
		p = p1;
		set(r, V, t);
	}
	OrbitalElements(
		shared_ptr<Planet> p1,
		double a1, double e1,
		double i1, double o1, double w1,
		double v1
	)
		:p(p1), a(a1), e(e1), i(i1), o(o1), w(w1), v(v1), periapsis(a* (1 - e)) {
	};
	OrbitalElements(
		shared_ptr<Planet> p1,
		double periapsis1, double e1,
		double i1, double o1, double w1,
		double v1, bool parabolic
	)
		:p(p1), a(numeric_limits<double>::infinity()), e(e1), i(i1), o(o1), w(w1), v(v1), periapsis(periapsis1) {
	};
	double period() {
		//yep this and meanMotion multiply to make 2PI
		//apart from the sign part which I basically invented to make hyperbolic
		//orbits have negative periods just like how they have negative semi major axes
		return 2 * PI * ofSign(a) * sqrt(abs(a * a * a) / p->gravity);
	}
	double meanMotion() {
		if (e == 1)return 1;
		return sqrt(p->gravity / abs(a * a * a));
	}
	double trueAnomaly(double t) {
		double M = v + meanMotion() * t;
		if (e < 1) {
			M = ffmod(M, 2 * PI);
			if (e < 1e-3) return M;//that is a big epsilon because
			//the error is epsilon squared here
			double E = M + e * sin(M) + 0.5 * e * e * sin(2 * M);
			double dE = 10000;
			while (abs(dE) > 1e-12) {
				dE = (E - e * sin(E) - M) / (1 - e * cos(E));
				E -= dE;
				E = ffmod(E, 2 * PI);
			}
			return 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));
		}
		if (e > 1) {
			double H = asinh(M / e);//ofSign(M) * log(2 * abs(M) / e + 1.8);
			double dH = 10000;
			while (abs(dH) > 1e-12) {
				dH = (e * sinh(H) - H - M) / (e * cosh(H) - 1);
				H -= dH;
			}
			return 2 * atan2(sqrt(e + 1) * sinh(H / 2), sqrt(e - 1) * cosh(H / 2));
		}
		double B = -1.5 * M * sqrt(p->gravity / (2 * periapsis * periapsis * periapsis));
		double rtdelta = sqrt(B * B - 1);
		double D = cbrt(B + rtdelta) + cbrt(B - rtdelta);
		return 2 * atan(D);
	}
	glm::dmat4 mat() {
		glm::dquat q =
			glm::angleAxis(w, glm::dvec3(0, 1, 0)) *
			glm::angleAxis(i, glm::dvec3(1, 0, 0)) *
			glm::angleAxis(o, glm::dvec3(0, 1, 0));
		return glm::mat4_cast(q);
		/*
		glm::dquat q_o = glm::angleAxis(o, glm::dvec3(0, 1, 0));
		glm::dquat q_i = glm::angleAxis(i, glm::rotate(q_o, glm::dvec3(1, 0, 0)));
		glm::dquat q_w = glm::angleAxis(w, glm::rotate(q_i * q_o, glm::dvec3(0, 1, 0)));

		glm::dquat q = q_w * q_i * q_o;
		return glm::mat4_cast(q);
		*/
	}
	glm::dvec3 tpos(double n) {
		//input:true anomaly
		double r;
		if (e == 1) {
			r = 2 * periapsis / (1 + cos(n));
		}
		else {
			r = a * (1 - e * e) / (1 + e * cos(n));
		}
		glm::dvec4 ret(r * cos(n), 0.0, r * sin(n), 1.0);
		ret = mat() * ret;
		return glm::dvec3(ret);
	}
	glm::dvec3 tvel(double n) {
		//input:true anomaly
		double r, th, coeff, vcos = 1 + e * cos(n);//a bit off from what "vercosin" is,but whatever
		glm::dvec4 ret;
		if (e == 1) {
			coeff = sqrt(p->gravity / (a * (1 - e * e)));
			r = 2 * periapsis / vcos;
			th = sqrt(2 * p->gravity / r);
			ret = { -sin(n), 0.0, vcos, 1.0 };
		}
		else {
			coeff = sqrt(p->gravity / (2 * periapsis));
			r = coeff * e * sin(n);
			th = coeff * vcos;
			ret = { r * cos(n) - th * sin(n), 0.0, r * sin(n) + th * cos(n), 1.0 };
		}
		ret = mat() * ret;
		//TODO:testing.this piece of code has not been used anywhere.
		if (glm::length(glm::dvec3(ret) - (tpos(n + 0.0000000001) - tpos(n)) / 0.0000000001) > 0.01)throw "far off";
		return glm::dvec3(ret);
	}
	glm::dvec3 pos(double t) {
		//input:time elapsed since t=0
		double n = trueAnomaly(t);
		return tpos(n);
	}
	glm::dvec3 vel(double t) {
		double n = trueAnomaly(t);
		return tvel(n);
		//double ep = 1.0/60;
		//return (pos(t+ep)-pos(t))/ep;
		//yeah chat gpt keeps using the vis viva equation so here's my (practical) answer
	}
	glm::dvec3 apos(double t) {
		return pos(t) + p->apos(t);
	};
	glm::dvec3 avel(double t) {
		return vel(t) + p->avel(t);
	};
	void set(glm::dvec3 r, glm::dvec3 V, double t) {
		//relative position and velocity
		glm::dvec3 h = glm::cross(r, V);
		glm::dvec3 N = glm::cross(glm::dvec3(0, 1, 0), h);
		glm::dvec3 E = glm::cross(V, h) / p->gravity - glm::normalize(r);
		double ep = glm::length2(V) / 2 - p->gravity / glm::length(r);
		a = -p->gravity / (2 * ep);
		e = glm::length(E);
		i = ffmod(acos(glm::clamp(h.y / glm::length(h), -1.0, 1.0)), 2 * PI);
		if (isnan(i))i = 0;
		o = ffmod(acos(glm::clamp(N.x / glm::length(N), -1.0, 1.0)), 2 * PI);
		if (N.z < 0)o = 2 * PI - o;
		w = ffmod(acos(glm::clamp(glm::dot(N, E) / (glm::length(N) * glm::length(E)), -1.0, 1.0)), 2 * PI);
		if (E.y < 0)w = 2 * PI - w;
		if (isnan(o)) {
			o = 0;
			w = ffmod(atan2(E.z, E.x), 2 * PI);
			if (h.y < 0)w = 2 * PI - w;
		}
		if (isnan(w))w = 0;
		double n = -ffmod(acos(glm::clamp(glm::dot(E, r) / (glm::length(E) * glm::length(r)), -1.0, 1.0)), 2 * PI);
		//its negative for some weird reason
		if (glm::dot(r, V) < 0) n = 2 * PI - n;
		if (isnan(n)) {
			n = ffmod(acos(glm::clamp(glm::dot(N, r) / (glm::length(N) * glm::length(r)), -1.0, 1.0)), 2 * PI);
			if (r.y < 0) n = 2 * PI - n;
			if (isnan(n)) {
				n = ffmod(acos(glm::clamp(r.x / glm::length(r), -1.0, 1.0)), 2 * PI);
				if (V.x > 0) n = 2 * PI - n;
			}
		}
		/*
		double n = atan2(
			glm::dot(r, glm::cross(h, E)) / (glm::length(h) * e * glm::length(r)),
			glm::dot(E, r) / (glm::length(h) * glm::length(r))
		);
		n = ffmod(n, 2 * PI);
		*/
		double M;
		if (e < 1) {
			//M = atan2(sqrt(1 - e * e) * sin(n), e + cos(n)) - e * sqrt(1 - e * e) * sin(n) / (1 + e * cos(n));
			double E = atan2(sqrt(1 - e * e) * sin(n), e + cos(n));
			M = E - e * sin(E);
		}
		else if (e > 1) {
			double F = 2 * atanh(tan(n / 2) * sqrt((e - 1) / (e + 1)));
			M = e * sinh(F) - F;
		}
		else {
			//parabolic orbits are annoying because a*(1-e) is NaN so we will have to store that
			periapsis = glm::length2(h) / (2 * p->gravity);
			double D = tan(n / 2);
			M = sqrt(2 * periapsis * periapsis * periapsis / p->gravity) * (D + D * D * D / 3);
		}
		v = M - meanMotion() * t;
		if (e < 1)v = ffmod(v, 2 * PI);
	}
	void displayMode1(double t) {
		double div = 2 * PI / 256;
		ofNoFill();
		ofSetColor(0, 255, 0);
		ofBeginShape();
		double th = PI - (e > 1) * PI / 2;
		for (double i = -th; i <= th; i += div) {
			ofVertex(glm::vec3((tpos(i) + p->apos(t) - orbitPos) * orbitScale));
		}
		ofEndShape(e < 1);
	}
};
vector<shared_ptr<Planet>> planets = {
	make_shared<Planet>(Planet("Kerbin",3.5316e12,600000,101325,70000/log(101325/ATMO_TOL)))
};
void createPlanetAtlas() {
	if (planetAtlasReady)return;
	planetAtlasReady = true;
	ofFbo::Settings fboSettings;
	fboSettings.width = ATLAS_SIZE;
	fboSettings.height = ATLAS_SIZE;
	fboSettings.internalformat = GL_RGBA;
	fboSettings.numSamples = 0;
	planetAtlas.allocate(fboSettings);
	planetAtlas.getTexture().setTextureMinMagFilter(GL_NEAREST, GL_NEAREST);
	planetAtlas.begin();
	ofClear(0);
	for (int y = 0; y < planets.size(); y++) {
		for (int x = 0; x < ATLAS_SIZE; x++) {
			ofFill();
			ofSetColor(planets[y]->color((double(x) / ATLAS_SIZE) * 2 - 1));
			ofDrawRectangle(x, y, 1, 1);
		}
	}
	planetAtlas.end();
}
Planet::Planet(string n, double g, double r, double p0, double ph) {
	name = n;
	gravity = g;
	radius = r;
	atmos = make_pair(p0, ph);
}
Planet::Planet(string n, double g, double r, double p0, double ph, shared_ptr<OrbitalElements> o1) {
	name = n;
	gravity = g;
	radius = r;
	atmos = make_pair(p0, ph);
	o = o1;
};
glm::dvec3 Planet::pos(double t) {
	if (o == nullptr)return glm::dvec3(0, 0, 0);
	return o->pos(t);
}
glm::dvec3 Planet::vel(double t) {
	if (o == nullptr)return glm::dvec3(0, 0, 0);
	return o->vel(t);
}
glm::dvec3 Planet::apos(double t) {
	if (o == nullptr)return glm::dvec3(0, 0, 0);
	return o->apos(t);
}
glm::dvec3 Planet::avel(double t) {
	if (o == nullptr)return glm::dvec3(0, 0, 0);
	return o->avel(t);
}
double Planet::SOI() {
	if (o == nullptr)return std::numeric_limits<double>::infinity();
	return o->a * pow(gravity / o->p->gravity, 0.4);
};
double Planet::atmoPressure(double h) {
	//with a linear temperature gradient(assumed) atmospheric pressure is always an exponential of this form.
	if (h <= 0)return atmos.first;//prevent pressure from being a lot
	return atmos.first * exp(-h / atmos.second);
}
double Planet::atmoHeight() {
	//remember this is relative to sea level not the center of the planet
	return -atmos.second * log(ATMO_TOL / atmos.first);
}
void Planet::displayMode1(double t, int planetI) {
	glm::dvec3 p = (apos(t) - orbitPos) * orbitScale;
	ofPushMatrix();
	ofTranslate(glm::vec3(p));
	ofFill();
	//ofDrawSphere(radius * orbitScale);
	glm::dvec3 ref = (glm::dvec3)camera1.getGlobalPosition() - p;
	double r = glm::length(ref);
	double ter = 1 + terrain.get(ref);
	if (r <= radius * orbitScale * ter) {
		orbitScale = r / (radius * ter);
	}
	ref /= radius * orbitScale;
	//ref = { 0,1.5,0 };
	ofSetColor(255, 255, 255);
	mesh.clear();//probably superfluous
	terrain.mesh(mesh, ref, planetI);
	brush.getMesh() = std::move(mesh);
	brush.setScale(radius * orbitScale / DM3_SCALE);
	createPlanetAtlas();
	ofTexture tex = planetAtlas.getTexture();
	tex.setTextureMinMagFilter(GL_NEAREST, GL_NEAREST);
	tex.bind();
	brush.draw();
	tex.unbind();
	/*
	glm::dvec3 rayDir = glm::normalize(-ref), rayHit;
	if (raycastSH_c(rayHit, terrain.coeff, terrain.lipschitz, terrain.maxHeight, ref, rayDir)) {
		ofSetColor(255, 0, 0);
		ofDrawSphere(glm::vec3(radius * orbitScale * rayHit), 16);
	}
	*/
	ofTranslate(-glm::vec3(p));//superfluous but i'll keep it here because why not
	ofPopMatrix();
	if (o != nullptr)o->displayMode1(t);
}
void Planet::displayMode3(glm::dvec3 shipPos, int planetI) {
	ofPushMatrix();//superfluous but i'll keep it here because why not
	ofFill();
	glm::dvec3 cameraPos = glm::dvec3(0, 0, 0);//camera.getGlobalPosition()/DM3_SCALE;
	glm::dvec3 ref = cameraPos + shipPos;
	double r = glm::length(ref);
	double ter = 1 + terrain.get(ref);
	if (r <= radius * ter) {
		//TODO:Avoid problems by moving camera away from terrain
	}
	ref /= radius;
	//ref = { 0,1.001,0 };
	ofSetColor(255, 255, 255);//TODO
	mesh.clear();//probably superfluous
	terrain.mesh(mesh, ref, planetI, true);
	brush.getMesh() = std::move(mesh);
	brush.setScale(radius);
	createPlanetAtlas();
	ofTexture t = planetAtlas.getTexture();
	t.setTextureMinMagFilter(GL_NEAREST, GL_NEAREST);
	t.bind();
	brush.draw();
	t.unbind();
	/*
	glm::dvec3 rayDir = glm::normalize(-cameraPos), rayHit;
	if (raycastSH_c(rayHit, terrain.coeff, terrain.lipschitz, terrain.maxHeight, ref, rayDir)) {
		ofSetColor(255, 0, 0);
		glm::dvec3 actualHit = radius * rayHit - shipPos;
		ofDrawSphere(glm::vec3(actualHit), 16000);
		//remember that the radius is in real coordinates rn
		//a sphere of 16 is not visible lol
	}
	*/
	ofPopMatrix();
}
class MenuElement;
class GridElement;
class PhysicsGrid;
shared_ptr<GridElement> selectedPart;//sceneDisplayed==2
glm::dvec3 selectedPos;
glm::dvec2 menuPos;
constexpr double MENU_MARGIN = 4;
const ofColor MENU_COLOR1 = ofColor(0, 0, 0);//background
const ofColor MENU_COLOR2 = ofColor(255, 255, 255);//buttons etc
const ofColor MENU_COLOR3 = ofColor(0, 255, 0);//text
const ofColor MENU_COLOR4 = ofColor(0, 0, 255);//text on white background
shared_ptr<GridElement> resourceTransferOrigin;
string selectedElementType = "fire-tank";
const unordered_map<string, unordered_map<string, double>> elementBuildCosts = {
	{"fire-tank",{{"metal",1000}}},
	{"metal-tank",{{"metal",1000}}},
	{"solid-rocket-engine",{{"metal",1000}}},
	{"command-core",{{"spirits",1000}}},
	{"reaction-wheel",{{"metal",500},{"spirits",500}}}
};
enum class ShipSituation {
	LANDED,
	SlIDING,
	FLYING
};
string printSituation(ShipSituation situ) {
	switch (situ) {
	case ShipSituation::LANDED:
		return "Landed";
	case ShipSituation::SlIDING:
		return "Sliding";
	case ShipSituation::FLYING:
		return "Flying";
	default:
		return "idk";
	}
}
glm::dvec3 getAlignedAxis(glm::dvec3 v, vector<pair<glm::dvec3, double>>& c) {
	//returns normal of array of vectors aligned with the normal v
	//if the array is colinear,returns {nan,nan,nan}

	//proof that selecting any two points from any array with atleast three points
	//non-colinear,there is atleast one point that the three points are not colinear:
	//suppose the three non-colinear points are a,b,c
	//if the two points are colinear to a and b,choose c
	//if one is colinear to a and b,choose another point colinear to 
	//a and b that is not the point.the three points are not colinear
	//because the other point does not belong in the line of the point and the chosen point.
	//if neither are colinear to a and b,choose a if a is not colinear,otherwise b
	//if the resultant three points are colinear then a b and the line the two points
	//are on are actually the same line,which is wrong
	//"simple" trick to save time complexity lol
	if (c.size() < 3)return glm::dvec3(0) / 0.0;
	glm::dvec3 x = c[0].first, y = c[1].first, r;
	for (int i = 2; i < c.size(); i++) {
		r = glm::cross(y - x, c[i].first - x);
		if (r != glm::dvec3(0, 0, 0))break;
	}
	r = glm::normalize(r);
	if (glm::dot(r, v) < 0)return r;
	return -r;
}
const double SF_BONUS = 1.25;
glm::dvec3 corners_of_unit_cube[8] = {
	{-0.5,-0.5,-0.5},
	{-0.5,-0.5,0.5},
	{-0.5,0.5,-0.5},
	{-0.5,0.5,0.5},
	{0.5,-0.5,-0.5},
	{0.5,-0.5,0.5},
	{0.5,0.5,-0.5},
	{0.5,0.5,0.5}
};
struct MonoPropSpec {
	double Esp = 0;//chemical energy per kg
	int f = 0;//DOF of exhaust,if there are multiple types of exhaust give molar average
	double Ae = 0;//how big the bell is
	double Ar = 0;//if you look into the bell,how likely are you to see into the hole?
	double maxR = 0;//maximum injection rate
	string R;
};
struct BiPropSpec {
	int a = 0;
	int b = 0;
	int c = 0;
	double H = 0;//aA+bB->cC+H
	double Ma = 0;//molar masses
	double Mb = 0;
	double Mc = 0;
	int Fa = 0;//DOFs
	int Fb = 0;
	int Fc = 0;
	double Ae = 0;
	double Ar = 0;
	double maxX = 0;//maximum injection rates
	double maxY = 0;
	string A;
	string B;
};
constexpr double IDEAL_GAS_COEFF = 8.314;
double solveForMachNumber2(double r2, double f) {
	//(after manipulating the equation so it doesn't include x,only x^2,and moving a factor to the left side:)
	//r2x^2/(f/(f+1))^(f+1)=(1+x^2/f)^(f+1)
	//this equation has atmost two roots,take the one where x^2 is larger.
	//that solution represents the case where flow is choked and supersonic.
	//also return x^2,nobody needs x itsself
	double a = r2 / pow(f / (f + 1), f + 1), b = 1 / f, c = f + 1;
	double x = 100, dx = 0;//mach 10.apparently irl the highest-scoring ones here score ~8->x=64
	//NR method works great on this kind of equation where the first derivative is monotonous
	//apparently if x is ridiculously large then each step cuts it down by 1/(1+1/c)
	if (a * x - pow(1 + b * x, c) >= 0)throw "unreasonable mach number(>10)";
	double EP = 1e-8;
	do {
		dx = (a * x - pow(1 + b * x, c)) / (a - b * c * pow(1 + b * x, c - 1));
		x -= dx;
	} while (abs(dx) > EP);//dx shouldn't actually be less than 0 because of the monotonicity thing
	if (x < 1)throw "still unreasonable mach number(<1)";
	return x;
}
double calculateEjectionVelocityInner(double f, double Esp, double r, double Ar, double Ae, double Pa) {
	//calculates effective ejection velocity
	double y = 1 + 2 / f, fa = f / (f + 1), fb = (f + 1) / 2, fc = (f + 2) / 2, fd = 1 / f;
	double Me2 = solveForMachNumber2(Ar * Ar, f);//proportional to something about ar
	double Tr = 1 / (1 + Me2 * fd);//proportion of chemical energy turned into thermal energy(the rest is kinetic energy)
	double Ve = pow(Tr, fc) / (sqrt(y * fc / Esp) * pow(fa, fb));//pressure at exit=Ve*Ar/Ae*r
	//proportional to sqrt(Esp) and something about ar
	//"velocity the pressure would accelerate the mass to at throat if it was at exit pressure" yeah that's weird
	return max(0.0,sqrt(2 * Esp * (1 - Tr)) + Ve * Ar - Pa * Ae / r);//clamp isp so it doesn't go negative(engine no longer works at all)
	//first term proportional to sqrt(Esp) and something about ar
	//second term proportional to sqrt(Esp) and Ar/something about ar
	//(doesn't seem right but remember real rocket engines don't have a fixed mass flow rate and get choked)
	//third term proportional to atmospheric pressure
	//third term proportional to exit area and inversely proportional to mass flow rate
	//so:
	//atmospheric pressure linearly decreases isp
	//mass flow rate inversely decreases isp
	//exit area linearly decreases isp
	//fuel energy per unit mass ^1/2ly increases isp
	//area ratio either linearly or logarithmically increases isp depending on how incorrectly expanded the engine is
}
double calculateEjectionVelocity(MonoPropSpec spec, double Pa, double r) {
	return calculateEjectionVelocityInner(spec.f, spec.Esp, r * spec.maxR, spec.Ar, spec.Ae, Pa);
}
double calculateEjectionVelocity(BiPropSpec spec, double Pa, double x, double y) {
	double s = min(x / (spec.a * spec.Ma), y / (spec.b * spec.Mb));
	//reaction stochiometricness=rc*spec.Mc/(x+y) (?)
	double ra = x * spec.maxX / spec.Ma - spec.a * s;
	double rb = y * spec.maxY / spec.Mb - spec.b * s;
	double rc = spec.c * s;
	double mols = ra + rb + rc;//linear to the reagent depleted.
	//if the mols of product is greater than the OTHER reagent
	//then it increases with the mass flow rate of the depleted reagent
	//otherwise it decreases
	double f = (ra * spec.Fa + rb * spec.Fb + rc * spec.Fc) / mols;
	//complex relation between non-depleted reagent and f
	//something like (r/M+something)/(r/M+something else)
	double M = (ra * spec.Ma + rb * spec.Mb + rc * spec.Mc) / mols;
	double m = mols / s;
	//m*M increases with the ratio of the non-depleted reagent to the depleted reagent
	//effective energy density is proportional to the stochiometric energy density
	//and inversely proportional to m*M
	return calculateEjectionVelocityInner(f, spec.H / m / M, x * spec.maxX + y * spec.maxY, spec.Ar, spec.Ae, Pa);
}
class MenuElement {
public:
	virtual ~MenuElement() = default;

	virtual double width() const = 0;
	virtual double height() const = 0;
	virtual void display() = 0;
	virtual void update(shared_ptr<GridElement> parent, glm::dvec2 mp) = 0;
	ofRectangle boundingBox() const {
		return ofRectangle(0, 0, width(), height());
	}
};
class LabelMenuElement : public MenuElement {
public:
	string label;
	function<string(shared_ptr<GridElement>)> content;
	string value;//so updates only happen when update() is called
	LabelMenuElement(string l, function<string(shared_ptr<GridElement>)> c) : label(l),content(c) {}
	void update(shared_ptr<GridElement> parent, glm::dvec2 mp) {
		value = content(parent);
	}
	string text() const {
		return label + value;
	}
	void display() {
		ofSetColor(MENU_COLOR3);
		drawText(text(), 0, height());
	}
	double width() const {
		return textWidth(text());
	}
	double height() const {
		return textHeight(text());
	}
};
class ButtonMenuElement : public MenuElement, public enable_shared_from_this<ButtonMenuElement> {
public:
	string label;
	function<void(shared_ptr<GridElement>, shared_ptr<ButtonMenuElement>)> callback;
	double buttonWidth;
	ButtonMenuElement(string l, function<void(shared_ptr<GridElement>, shared_ptr<ButtonMenuElement>)> c, double w) : label(l), callback(c), buttonWidth(w) {}
	void update(shared_ptr<GridElement> parent, glm::dvec2 mp) {
		if (mousePressedOnLastFrame && mouse[0] && buttonBoundingBox().inside(mp)) {
			mousePressedOnLastFrame = false;
			callback(parent, shared_from_this());
		}
		else if(buttonBoundingBox().inside(mp)){
		}
	}
	ofRectangle buttonBoundingBox() const {
		return ofRectangle(textWidth(label) + MENU_MARGIN, 0, buttonWidth, height());
	}
	void display() {
		ofSetColor(MENU_COLOR3);
		drawText(label, 0, height());
		ofSetColor(MENU_COLOR2);
		ofDrawRectRounded(buttonBoundingBox(), height() / 2);
	}
	double width() const {
		return textWidth(label) + MENU_MARGIN + buttonWidth;
	}
	double height() const {
		return textHeight(label);
	}
};
class InputMenuElement : public MenuElement, public enable_shared_from_this<InputMenuElement> {
public:
	string label;
	function<void(shared_ptr<GridElement>, shared_ptr<InputMenuElement>,string)> callback;
	string input;
	bool focused = false;
	int pointer = 0;
	InputMenuElement(string l, function<void(shared_ptr<GridElement>, shared_ptr<InputMenuElement>,string)> c) : label(l), callback(c) {}
	void update(shared_ptr<GridElement> parent, glm::dvec2 mp) {
		if (mousePressedOnLastFrame && mouse[0]) {
			if (inputBoundingBox().inside(mp)) {
				dontDetectKeys = true;
				mousePressedOnLastFrame = false;
				focused = true;
				pointer = 0;
			}
			else if (focused) {
				dontDetectKeys = false;
				focused = false;
			}
		}
		if (!focused) return;
		if (keyPressedOnLastFrame == 0) return;
		if (keyPressedOnLastFrame == OF_KEY_BACKSPACE) {
			if (input.size() != 0 && pointer != 0) {
				input.erase(--pointer, 1);
			}
		}
		else if (keyPressedOnLastFrame == OF_KEY_LEFT) {
			pointer--;
		}
		else if (keyPressedOnLastFrame == OF_KEY_RIGHT) {
			pointer++;
		}
		else if (isprint(keyPressedOnLastFrame)) {
			input.insert(pointer++, 1, keyPressedOnLastFrame);
		}
		pointer = max(pointer, 0);
		pointer = min(pointer, (int)input.size());
		callback(parent, shared_from_this(), input);
	}
	ofRectangle inputBoundingBox() const {
		return ofRectangle(textWidth(label) + MENU_MARGIN, 0, max(textWidth(input), 2 * MENU_MARGIN), height());
	}
	void display() {
		ofSetColor(MENU_COLOR3);
		drawText(label, 0, height());
		ofSetColor(MENU_COLOR2);
		ofDrawRectangle(inputBoundingBox());
		ofSetColor(MENU_COLOR4);
		drawText(input, textWidth(label) + MENU_MARGIN, height());
		if (focused) {
			ofSetColor(0, 0, 0);
			ofDrawRectangle(textWidth(label) + MENU_MARGIN + textWidth("A") * pointer, 0, 1, height());
		}
	}
	double width() const {
		return textWidth(label) + MENU_MARGIN + textWidth(input);
	}
	double height() const {
		return max(textHeight(label), textHeight(input));
	}
};
class ProgressMenuElement : public MenuElement {
public:
	string label;
	function<pair<double, double>(shared_ptr<GridElement>)> content;
	pair<double, double> value;//so updates only happen when update() is called
	ProgressMenuElement(string l, function<pair<double, double>(shared_ptr<GridElement>)> c) : label(l), content(c) {}
	void update(shared_ptr<GridElement> parent, glm::dvec2 mp) {
		value = content(parent);
	}
	ofRectangle progressBoundingBox() const {
		return ofRectangle(textWidth(label) + MENU_MARGIN, 0, 256, height());
	}
	void display() {
		ofSetColor(MENU_COLOR3);
		drawText(label, 0, height());
		ofSetColor(MENU_COLOR2);
		ofDrawRectRounded(progressBoundingBox(), height() / 2);
		ofSetColor(MENU_COLOR3);
		ofDrawRectRounded(textWidth(label) + MENU_MARGIN, 0, 256.0 * value.first / value.second, height(), height() / 2);
		ofSetColor(MENU_COLOR4);
		string t = to_string(value.first) + "/" + to_string(value.second);
		drawText(t, textWidth(label) + MENU_MARGIN + 128 - textWidth(t) / 2, height());
	}
	double width() const {
		return textWidth(label) + MENU_MARGIN + 256;
	}
	double height() const {
		return textHeight(label);
	}
};
class SliderMenuElement : public MenuElement, public enable_shared_from_this<SliderMenuElement> {
public:
	string label;
	function<pair<double, double>(shared_ptr<GridElement>, shared_ptr<SliderMenuElement>, double v)> callback;
	pair<double, double> value;//so updates only happen when update() is called
	bool focused = false;
	SliderMenuElement(string l, function<pair<double, double>(shared_ptr<GridElement>, shared_ptr<SliderMenuElement>, double v)> c) : label(l), callback(c) {}
	void update(shared_ptr<GridElement> parent, glm::dvec2 mp) {
		if (progressBoundingBox().inside(mp)&&mouse[0]) {
			focused = true;
		}
		if (focused) {
			if (!mouse[0]) {
				focused = false;
			}
			else {
				mousePressedOnLastFrame = false;
				double pos = (mp.x - progressBoundingBox().getX()) / 256.0;
				pos = glm::clamp(pos, 0.0, 1.0);
				double factor = pow(2, glm::round(log(mapScale)) + 1);
				pos = glm::round(pos * factor) / factor;
				value.first = value.second * pos;
			}
		}
		value = callback(parent, shared_from_this(), value.first);
	}
	ofRectangle progressBoundingBox() const {
		return ofRectangle(textWidth(label) + MENU_MARGIN, 0, 256, height());
	}
	void display() {
		ofSetColor(MENU_COLOR3);
		drawText(label, 0, height());
		ofSetColor(MENU_COLOR2);
		ofDrawRectRounded(progressBoundingBox(), height() / 2);
		ofSetColor(MENU_COLOR3);
		ofDrawRectRounded(textWidth(label) + MENU_MARGIN, 0, 256.0 * value.first / value.second, height(), height() / 2);
		ofSetColor(MENU_COLOR4);
		string t = to_string(value.first) + "/" + to_string(value.second);
		drawText(t, textWidth(label) + MENU_MARGIN + 128 - textWidth(t) / 2, height());
	}
	double width() const {
		return textWidth(label) + MENU_MARGIN + 256;
	}
	double height() const {
		return textHeight(label);
	}
};
class GridElement : public enable_shared_from_this<GridElement> {
public:
	GridElement(string t);
	string type;
	double dryMass, COR, COF;//COF is kinetic friction,static friction is approximated by COF*1.25
	int frontFace, topFace, rightFace;
	glm::dvec3 faceNormal(int face);
	glm::dvec3 front();
	glm::dvec3 top();
	glm::dvec3 right();
	vector<int> faces();
	glm::dvec3 pointingFace(int p);
	void rotateFrontFace();
	void rotateTopFace();
	void rotateRightFace();
	vector<shared_ptr<MenuElement>> menu;
	void displayMenu();
	void displayMode2();
	unordered_map<string, pair<double, double>> fluids;
	unordered_map<string, double> fluidChanges;
	double tankTransferRate();
	MonoPropSpec monopropEngineData();//std::variant is too troublesome for this
	BiPropSpec bipropEngineData();
	glm::dvec3 reactionWheelTorque();
	double throttleA = 0;
	double throttleB = 0;
	double throttleC = 0;
	glm::dvec3 thrust;
	glm::dvec3 torque;
	//TODO:possible off-center thrust?
	double mass();
	void update(glm::dvec3 pos, double t, double dt);
	void integrate(double t, double dt);
	shared_ptr<PhysicsGrid> parent;
};
class PhysicsGrid : public enable_shared_from_this<PhysicsGrid> {
	//stupid c++ requiring me to declare such a large class
public:
	shared_ptr<Planet> soi;
	glm::dvec3 position;//position of COM now
	glm::dvec3 velocity;
	glm::dvec3 accel;//verlet internal
	glm::dquat angle;
	glm::dvec3 avel;//along axis of rotation,magnitude is amount and direction of rotation
	//in local space apparently
	OrbitalElements orbit;
	ShipSituation situ = ShipSituation::FLYING;
	vector<pair<glm::dvec3, double>> contacts;//local position/friction coefficient
	//remember after each state switch justify that the subsequent code does handle the state again
	PhysicsGrid(glm::dvec3 p, glm::dvec3 v, shared_ptr<Planet> planet, double t);
	shared_ptr<GridElement> getItem(int x, int y, int z);
	shared_ptr<GridElement> getItem(glm::dvec3 index);
	shared_ptr<GridElement> setItem(shared_ptr<GridElement> e, int x, int y, int z);
	shared_ptr<GridElement> setItem(shared_ptr<GridElement> e, glm::dvec3 index);
	void removeItem(int x, int y, int z);
	void removeItem(glm::dvec3 v);
	void displayMode1(double t);
	void mousePressed(int x, int y, int button);
	void displayMode2(int y);
	void displayMode3();
	void updateGrid();
	glm::dvec3 calculateContacts(glm::dvec3 sn);
	glm::dvec3 computeNFDistribution(glm::dvec3 kh, glm::dvec3 force, glm::dvec3 torque, glm::dvec3& ih, glm::dvec3& jh);
	void updatePhysics(double t, double dt);
	void updateContactPhysics(double t, double dt);
	void updatePhysicsLanded(double t, double dt);
	void updatePhysicsSliding(double t, double dt);
	void updatePhysicsFlying(double t, double dt);
	void integrate(double t, double dt, glm::dvec3 force, glm::dvec3 torque);
	pair<glm::dvec3, glm::dvec3> calculateLoads(double t, double dt);
	double checkCollision(double maxDT, double H, glm::dvec3 sn, glm::dvec3& ret_position, glm::dvec3& contact_normal, shared_ptr<GridElement>& part);
	double totalEnergy();
	void updateInternal(double t, double dt);
private:
	unordered_map<glm::dvec3, shared_ptr<GridElement>> contents;
	//y points down,z points right,x points forwards in that order
	double mass, radius;//radius of bounding sphere
	glm::dvec3 COM;//ship coordinates
	glm::dmat3 inertialTensor;
	of3dPrimitive brush;
	ofMesh mesh;
};
GridElement::GridElement(string t) {
	type = t;
	dryMass = -1;
	COR = 0.7;
	COF = 0.15;
	frontFace = 4;
	topFace = 2;
	rightFace = 0;
	menu = {};
	if (type == "fire-tank") {
		dryMass = 1000;
		fluids = { {"fire",{0,1000}} };
	}
	else if (type == "metal-tank") {
		dryMass = 1000;
		fluids = { {"metal",{0,1000}} };
		throttleA = 1;
		/*menu.push_back(make_shared<SliderMenuElement>(SliderMenuElement("Injection rate", [&](auto g, auto s, double v) {
			g->throttleA = v;
			return make_pair(v, 1);
			}
		)));*///having flow rate control here is unrealistic
	}
	else if (type == "solid-rocket-engine") {
		dryMass = 1000;
		//TODO:thrust
		fluids = { {"metal",{0,1000}} };
		menu.push_back(make_shared<SliderMenuElement>(SliderMenuElement("Injection rate", [&](auto g, auto s, double v) {
			g->throttleA = v;
			return make_pair(v, 1);
			}
		)));
		menu.push_back(make_shared<LabelMenuElement>(LabelMenuElement("Isp:", [&](auto g) {
			double h = glm::length(g->parent->position) - g->parent->soi->radius;
			double Pa = (h > g->parent->soi->atmoHeight()) ? 0 : g->parent->soi->atmoPressure(h);
			return to_string(calculateEjectionVelocity(g->monopropEngineData(), Pa, g->throttleA));
			}
		)));
	}
	else if (type == "reaction-wheel") {
		dryMass = 1000;
		menu.push_back(make_shared<SliderMenuElement>(SliderMenuElement("Pitch authority", [&](auto g, auto s, double v) {
			g->throttleA = v;
			return make_pair(v, 1);
			}
		)));
		menu.push_back(make_shared<SliderMenuElement>(SliderMenuElement("Yaw authority", [&](auto g, auto s, double v) {
			g->throttleB = v;
			return make_pair(v, 1);
			}
		)));
		menu.push_back(make_shared<SliderMenuElement>(SliderMenuElement("Roll authority", [&](auto g, auto s, double v) {
			g->throttleC = v;
			return make_pair(v, 1);
			}
		)));
	}
	else if (type == "command-core") {
		dryMass = 1000;
		fluids = {
			{ "metal",{0,2000} },
			{ "spirits",{0,1000} }
		};
		menu.push_back(make_shared<SliderMenuElement>(SliderMenuElement("Pitch authority", [&](auto g, auto s, double v) {
			g->throttleA = v;
			return make_pair(v, 1);
			}
		)));
		menu.push_back(make_shared<SliderMenuElement>(SliderMenuElement("Yaw authority", [&](auto g, auto s, double v) {
			g->throttleB = v;
			return make_pair(v, 1);
			}
		)));
		menu.push_back(make_shared<SliderMenuElement>(SliderMenuElement("Roll authority", [&](auto g, auto s, double v) {
			g->throttleC = v;
			return make_pair(v, 1);
			}
		)));
		for (auto& pair : elementBuildCosts) {
			string key = pair.first;
			menu.push_back(make_shared<ButtonMenuElement>(ButtonMenuElement("Create " + key, [key](auto g, auto s) {
				/*if (resourceTransferOrigin == nullptr)*/ resourceTransferOrigin = g;
				selectedElementType = key;
				}
			, textWidth("aa"))));
		}
		menu.push_back(make_shared<ButtonMenuElement>(ButtonMenuElement("Recycle element", [](auto g, auto s) {
			/*if (resourceTransferOrigin == nullptr)*/resourceTransferOrigin = g;
			selectedElementType = "";
			}
		, textWidth("aa"))));
	}
	if (fluids.size() != 0) {
		menu.push_back(make_shared<ButtonMenuElement>(ButtonMenuElement("Set as transfer origin", [&](auto g, auto s) {
			resourceTransferOrigin = g;
			}
		, textWidth("aa"))));
		menu.push_back(make_shared<ButtonMenuElement>(ButtonMenuElement("Reset transfer origin", [&](auto g, auto s) {
			resourceTransferOrigin = nullptr;
			}
		, textWidth("aa"))));
	}
	for (auto& pair : fluids) {
		string key = pair.first;
		menu.push_back(make_shared<ProgressMenuElement>(ProgressMenuElement(key, [key](auto g) {return g->fluids[key]; })));
		menu.push_back(make_shared<ButtonMenuElement>(ButtonMenuElement("Transfer resource "+key, [key](auto g, auto s) {
			if (
				resourceTransferOrigin == nullptr ||
				resourceTransferOrigin == g ||
				resourceTransferOrigin->fluids.find(key) == resourceTransferOrigin->fluids.end()
				)return;
			double delta = min(g->fluids[key].second - g->fluids[key].first, resourceTransferOrigin->fluids[key].first);
			resourceTransferOrigin->fluids[key].first -= delta;
			g->fluids[key].first += delta;
			}
		, textWidth("aa"))));
	}
	if (dryMass == -1)throw"fix your code";
};
glm::dvec3 GridElement::faceNormal(int face) {//returns ship coordinates
	switch (face) {
	case 0:
		return { 1,0,0 };//right
	case 1:
		return { -1,0,0 };//left
	case 2:
		return { 0,1,0 };//up
	case 3:
		return { 0,-1,0 };//down
	case 4:
		return { 0,0,1 };//forwards
	case 5:
		return { 0,0,-1 };//backwards
	}
	throw "screw you compiler warning";
	return { 0,0,0 };
}
glm::dvec3 GridElement::front() {//returns normal vector of front face(in ship coordinates)
	return faceNormal(frontFace);
}
glm::dvec3 GridElement::top() {
	//alright,now we need to determine where the top face is.
	//it cannot be the front face or the back face.
	int i = 0, j = topFace;
	for (; i <= j; i++) {
		if (i == frontFace || i == (frontFace ^ 1))j++;
	}
	return faceNormal(j);
}
glm::dvec3 GridElement::right() {
	//the right face cannot be any of the other four faces.
	int j = topFace;
	for (int i = 0; i <= j; i++) {
		if (i == frontFace || i == (frontFace ^ 1))j++;
	}
	int k = rightFace;
	for (int i = 0; i <= k; i++) {
		if (i == frontFace || i == (frontFace ^ 1))k++;
		if (i == j || i == (j ^ 1))k++;
	}
	return faceNormal(k);
}
vector<int> GridElement::faces() {
	vector<int> result(6);
	result[frontFace] = 4;
	result[frontFace ^ 1] = 5;
	int j = topFace;
	for (int i = 0; i <= j; i++) {
		if (i == frontFace || i == (frontFace ^ 1))j++;
	}
	result[j] = 2;
	result[j ^ 1] = 3;
	int k = rightFace;
	for (int i = 0; i <= k; i++) {
		if (i == frontFace || i == (frontFace ^ 1))k++;
		if (i == j || i == (j ^ 1))k++;
	}
	result[k] = 0;
	result[k ^ 1] = 1;
	return result;
}
glm::dvec3 GridElement::pointingFace(int p) {
	//what face "up" points to on a face
	//front if not front or back,
	//front points to top,back points to bottom
	//input is in craft coordinates NOT part coordinates
	vector<int> f = faces();
	if (f[p] < 4) return front();
	if (f[p] == 4)return top();
	return -top();
}
void GridElement::rotateFrontFace() {
	frontFace++;
	frontFace %= 6;
}
void GridElement::rotateTopFace() {
	topFace++;
	topFace %= 4;
}
void GridElement::rotateRightFace() {
	rightFace++;
	rightFace %= 2;
}
void GridElement::displayMenu() {
	double w = 0, h = MENU_MARGIN;
	for (auto& item : menu) {
		w = max(w, item->width()) + MENU_MARGIN * 2;
		h += item->height() + MENU_MARGIN;
	}
	ofSetColor(MENU_COLOR1);
	ofDrawRectangle(0, 0, w, h);
	double x = MENU_MARGIN, y = MENU_MARGIN;
	glm::dvec2 mp = untransform2D(mousePos);
	for (auto& item : menu) {
		ofPushMatrix();
		ofTranslate(x, y);
		glm::dvec2 lmp = mp - glm::dvec2(x, y);
		y += item->height() + MENU_MARGIN;
		item->update(shared_from_this(),lmp);
		item->display();
		ofPopMatrix();
	}
}
void GridElement::displayMode2() {
	vector<int> f = faces();
	glm::dvec3 p = pointingFace(2);
	glm::dvec2 n(-p.z, -p.x);
	ofPushMatrix();
	ofTranslate(8, 8);
	ofRotateRad(atan2f(n.y, n.x) + PI / 2);
	if (glm::cross(front(), top()) == right())ofScale(-1, 1);
	//problem:does not consider float imprecision,though in theory it doesn't need to
	ofTranslate(-8, -8);
	drawImage(".\\textures\\parts\\" + type + "_" + faceNames[f[2]] + ".png", 0, 0);
	ofPopMatrix();
}
double GridElement::tankTransferRate() {
	if (type == "fire-tank")return 1000;
	if (type == "metal-tank")return 1000;
	return 0;
}
MonoPropSpec GridElement::monopropEngineData() {
	//v=sqrt(2fcRT/M*(1-(pa/p0)^(1/fc)))+(Pe-Pa)/Ae/r
	//p0=r/(An*sqrt(yM/(RT))*fa^fb)
	//MeAe/An=(kfa)^fb
	//k=1+fdMe^2
	//Pe=P0*k^-fc
	//y=1+2/f
	//fa=f/(f+1)
	//fb=(f+1)/2
	//fc=(f+2)/2
	//fd=1/f
	//T=H/(Rfcm)->RT/M=H/m*(fcM)

	//for monopropellant engines:
	//m=mol of exhaust per reaction
	//f=3*number of atoms in exhaust-number of DOF not activated
	//(look up wavenumber then mT=hcv/kb where v is the wavenumber.
	//a DOF is not activated if its equivalent temperature is more than the actual temperature.)
	//r=kg of fuel per second
	//needs:m,H,M,f,Ae,An,maxR
	//needs r and Pa from elsewhere

	//this is more like a sugar/KNO3 rocket but i'm still going to stick to the elemental thing
	//and call it metal.also the nozzle/throat sizes are dubious and come from chatgpt's janky math.
	//although i made it have about the right amount of expansion at atmospheric pressure
	if (type == "solid-rocket-engine")return MonoPropSpec(3.5e6, 10, 0.15, 15, 30, "metal");
	return MonoPropSpec();
}
BiPropSpec GridElement::bipropEngineData() {
	//v=sqrt(2fcRT/M*(1-(pa/p0)^(1/fc)))+(Pe-Pa)/Ae/r
	//p0=r/(An*sqrt(yM/(RT))*fa^fb)
	//MeAe/An=(kfa)^fb
	//k=1+fdMe^2
	//Pe=P0*k^-fc
	//y=1+2/f
	//fa=f/(f+1)
	//fb=(f+1)/2
	//fc=(f+2)/2
	//fd=1/f
	//T=H/(Rfcm)->RT/M=H/m*(fcM)

	//for bipropellant engines:
	//if one mol of reaction needs a mol of A and b mol of B
	//producing c mol of exhaust,and actually
	//x kg and y kg of reactant is supplied every second,
	//then s=min(x/(aMa),y/(bMb)) sets of reagents have reacted,
	//ra=x/Ma-a*s
	//rb=y/Mb-b*s
	//rc=c*s
	//f=(ra*Fa+rb*Fb+rc*Fc)/(ra+rb+rc)
	//m=(ra+rb+rc)/s
	//r=x+y
	//M=(ra*Ma+rb*Mb+rc*Mc)/(ra+rb+rc)
	//needs:a,b,c,H,Ma,Mb,Mc,Fa,Fb,Fc,Ae,An,maxX,maxY
	//needs x,y,and Pa from somewhere else(how do we let the user control the fuel ratio?)
	return BiPropSpec();
}
glm::dvec3 GridElement::reactionWheelTorque() {
	if (type == "reaction-wheel")return glm::dvec3(10000, 10000, 10000);
	if (type == "command-core")return glm::dvec3(5000, 5000, 5000);
	return glm::dvec3(0);
}
double GridElement::mass() {
	double m = dryMass;
	for (auto& pair : fluids) {
		m += pair.second.first;
	}
	return m;
}
void GridElement::update(glm::dvec3 pos, double t, double dt) {
	double rate = tankTransferRate();
	if (rate != 0) {
		shared_ptr<GridElement> next = parent->getItem(pos + front());
		if (next == nullptr || type == "metal-tank" &&
			(front() != next->front() || next->type != "metal-tank" && next->type != "solid-rocket-engine")) return;
		//TODO:test solid fuel transfer stuff
		for (auto& pair : fluids) {
			auto target = next->fluids.find(pair.first);
			if (target == next->fluids.end())continue;
			double delta = min(rate * throttleA * dt, pair.second.first);
			delta = min(delta, target->second.second - target->second.first);
			delta = max(0.0, delta);//just in case
			next->fluidChanges[pair.first] += delta;
			fluidChanges[pair.first] -= delta;
			//TODO:fix obvious problem of multiple fuel tanks pointing to a single one
			//being able to push too much fuel into said tank
			//(not a major issue:problem is self limiting as the tank cannot get fuller
			//after it is more than full on the next tick)
		}
		return;
	}
	double h = glm::length(parent->position) - parent->soi->radius;
	double Pa = (h > parent->soi->atmoHeight()) ? 0 : parent->soi->atmoPressure(h);
	auto specA = monopropEngineData();
	auto specB = bipropEngineData();
	if (specA.maxR) {
		shared_ptr<GridElement> next = parent->getItem(pos + front());
		if (next != nullptr) return;//oops,obstructed
		//TODO:better obstruction check
		//(check for all blocks in a cone?and the angle is related to how much the atmosphere harms the isp?)
		double r = throttleA;
		pair<double, double> pair = fluids[specA.R];
		double fuel = pair.first;
		double delta = min(fuel / dt, specA.maxR * r);
		thrust = -front() * delta * calculateEjectionVelocity(specA, Pa, r);
		fluidChanges[specA.R] -= delta * dt;
	}
	else if (specB.maxX) {
		shared_ptr<GridElement> next = parent->getItem(pos + front());
		if (next != nullptr) return;//oops,obstructed
		//TODO:better obstruction check
		//(check for all blocks in a cone?and the angle is related to how much the atmosphere harms the isp?)
		double x = specB.maxX * throttleA;
		double y = specB.maxY * throttleB;
		pair<double, double> pairA = fluids[specB.A];
		pair<double, double> pairB = fluids[specB.B];
		double A = pairA.first;
		double B = pairB.first;
		double deltaA = min(A / dt, x);
		double deltaB = min(B / dt, y);
		double delta = deltaA + deltaB;
		thrust = -front() * delta * calculateEjectionVelocity(specB, Pa, x, y);
		fluidChanges[specB.A] -= deltaA * dt;
		fluidChanges[specB.B] -= deltaB * dt;
	}
	else {
		thrust = glm::dvec3(0, 0, 0);
	}
	glm::dvec3 maxT = reactionWheelTorque();
	if (maxT != glm::dvec3(0)) {
		glm::dmat3 ori(right(), top(), front());//very convenient :D
		glm::dvec3 T = (ori * maxT) * glm::dvec3(throttleA, throttleB, throttleC);
		torque = T * glm::dvec3(keys['w'] - keys['s'], keys['d'] - keys['a'], keys['q'] - keys['e']);
	}
	else {
		torque = glm::dvec3(0);
	}
}
void GridElement::integrate(double t, double dt) {
	//actually does changes
	//this is to ensure that the order of processing does not matter
	for (auto& pair : fluids) {
		pair.second.first += fluidChanges[pair.first];
	}
	fluidChanges.clear();
}
PhysicsGrid::PhysicsGrid(glm::dvec3 p, glm::dvec3 v, shared_ptr<Planet> planet, double t)
	:orbit(planet, p, v, t) {
	//soi = planets[0];//assuming this is the one holding the global frame of reference
	//if (soi->o != nullptr)throw "hey wrong order of planets";
	soi = planet;
	position = p;
	velocity = v;
	accel = { 0,0,0 };
	avel = { 0,0,0 };
	angle = { 1,0,0,0 };
	brush = { mesh };
	updateGrid();
}
shared_ptr<GridElement> PhysicsGrid::getItem(int x, int y, int z) {
	return getItem(glm::dvec3(x, y, z));
}
shared_ptr<GridElement> PhysicsGrid::getItem(glm::dvec3 index) {
	if (contents.find(index) == contents.end())return nullptr;
	return contents[index];
}
shared_ptr<GridElement> PhysicsGrid::setItem(shared_ptr<GridElement> e, int x, int y, int z) {
	return setItem(e, glm::dvec3(x, y, z));
}
shared_ptr<GridElement> PhysicsGrid::setItem(shared_ptr<GridElement> e, glm::dvec3 index) {
	e->parent = shared_from_this();
	contents[index] = e;
	updateGrid();
	return e;
}
void PhysicsGrid::removeItem(int x, int y, int z) {
	removeItem({ (float)x,(float)y,(float)z });
}
void PhysicsGrid::removeItem(glm::dvec3 v) {
	getItem(v)->parent = nullptr;
	contents.erase(v);
	updateGrid();
}
void PhysicsGrid::displayMode1(double t) {
	glm::dvec3 pos = (orbit.apos(t) - orbitPos) * orbitScale;
	ofPushMatrix();
	ofTranslate(glm::vec3(pos));//draw ship(TODO)
	ofFill();
	ofSetColor(127, 127, 127);
	ofDrawSphere(16);
	ofTranslate(-glm::vec3(pos));
	ofPopMatrix();
	orbit.displayMode1(t);
}
void PhysicsGrid::displayMode2(int y) {
	if (mousePressedOnLastFrame && mouse[2]) selectedPart = nullptr;
	for (auto& ptr : contents) {
		glm::dvec3 pos = ptr.first;
		if (round(pos.y) != y) continue;
		//in this scene,x+ goes up while z+ goes left
		glm::dvec3 sPos(-pos.z * 16, -pos.x * 16, 0);
		glm::dvec2 m = untransform2D(mousePos);
		ofPushMatrix();
		ofTranslate(glm::vec3(sPos));
		ptr.second->displayMode2();
		if (mousePressedOnLastFrame && mouse[2]) {
			if (m.x > sPos.x && m.x < sPos.x + 16 && m.y > sPos.y && m.y < sPos.y + 16) {
				selectedPart = ptr.second;
				selectedPos = pos;
				menuPos = m - glm::dvec2(sPos);
			}
		}
		ofPopMatrix();
	}
	if (selectedPart != nullptr && round(selectedPos.y) == y) {
		glm::dvec2 sPos(-selectedPos.z * 16, -selectedPos.x * 16);
		ofFill();
		ofSetColor(0, 0, 255);
		drawFaceNormal2D(selectedPart->front(), sPos + glm::dvec2(8, 8));
		ofSetColor(0, 255, 0);
		drawFaceNormal2D(selectedPart->top(), sPos + glm::dvec2(8, 8));
		ofSetColor(255, 0, 0);
		drawFaceNormal2D(selectedPart->right(), sPos + glm::dvec2(8, 8));
		ofPushMatrix();
		ofTranslate(sPos + menuPos);
		ofScale(1.0 / 16);
		selectedPart->displayMenu();
		ofPopMatrix();
	}
	if (mousePressedOnLastFrame && mouse[0] && resourceTransferOrigin != nullptr) {
		glm::dvec2 m = untransform2D(mousePos);
		glm::dvec2 p = glm::floor(m / 16);
		glm::dvec3 pos(-p.y, y, -p.x);
		shared_ptr<GridElement> item = getItem(pos);
		if (item == nullptr && elementBuildCosts.find(selectedElementType) != elementBuildCosts.end()) {
			unordered_map<string,double> cost = elementBuildCosts.at(selectedElementType);
			bool afford = true;
			for (auto& pair : cost) {
				if (
					resourceTransferOrigin->fluids.find(pair.first) == resourceTransferOrigin->fluids.end() ||
					resourceTransferOrigin->fluids[pair.first].first < pair.second
					) {
					afford = false;
					break;
				}
			}
			if (afford) {
				for (auto& pair : cost) {
					resourceTransferOrigin->fluids[pair.first].first -= pair.second;
				}
				setItem(make_shared<GridElement>(GridElement(selectedElementType)), pos);
			}
		}
		else if (
			item != nullptr &&
			item != resourceTransferOrigin &&
			elementBuildCosts.find(selectedElementType) == elementBuildCosts.end() &&
			elementBuildCosts.find(item->type) != elementBuildCosts.end()
			) {
			auto cost = elementBuildCosts.at(item->type);
			for (auto& pair : item->fluids) {
				cost[pair.first] += pair.second.first;
			}
			bool afford = true;
			for (auto& pair : cost) {
				auto it = resourceTransferOrigin->fluids.find(pair.first);
				if (
					it == resourceTransferOrigin->fluids.end() ||
					it->second.second - it->second.first < pair.second
					) {
					afford = false;
					break;
				}
			}
			if (afford) {
				for (auto& pair : cost) {
					resourceTransferOrigin->fluids[pair.first].first += pair.second;
				}
				removeItem(pos);
				if (selectedPos == pos)selectedPart = nullptr;//this should always be true btw
			}
		}
	}
}
void PhysicsGrid::displayMode3() {
	brush.setGlobalOrientation(glm::dquat(1, 0, 0, 0));
	brush.setGlobalPosition(0, 0, 0);
	ofSetColor(255, 255, 255);
	ofPushMatrix();
	glm::mat4 matr = glm::toMat4(angle);
	ofMultMatrix(matr);
	ofTranslate(glm::vec3(-COM * DM3_SCALE));
	createAtlas();
	ofTexture t = atlas.getTexture();
	t.bind();
	brush.draw();
	t.unbind();
	//brush.drawNormals(100, true);
	ofPopMatrix();
}
void PhysicsGrid::updateGrid() {
	bool emptyShip = contents.size() == 0;
	mass = 0;
	if (emptyShip) {
		COM = glm::dvec3(0, 0, 0);
	}
	else {
		position -= glm::rotate(angle, COM);//so that the positions of the old blocks do not change after editing
		//need testing tho
		COM = glm::dvec3(0, 0, 0);
		for (auto& ptr : contents) {
			if (ptr.second == nullptr)continue;
			mass += ptr.second->mass();
			COM += (ptr.first + glm::dvec3(0.5, 0.5, 0.5)) * ptr.second->mass();
		}
		COM /= mass;
		position += glm::rotate(angle, COM);
	}
	inertialTensor = glm::dmat3(0.0);
	for (auto& ptr : contents) {
		if (ptr.second == nullptr)continue;
		glm::dvec3 pos = ptr.first;
		glm::dmat3 c(1.0 / 6);
		//should also have side length squared but that is 1.mass is accounted for later
		glm::dvec3 p = pos + glm::dvec3(0.5, 0.5, 0.5) - COM;
		double x = p.x, y = p.y, z = p.z;
		glm::dmat3 s(
			glm::dvec3(y * y + z * z, -x * y, -x * z),  // First row
			glm::dvec3(-y * x, x * x + z * z, -y * z),  // Second row
			glm::dvec3(-z * x, -z * y, x * x + y * y)   // Third row
		);
		inertialTensor += (c + s) * ptr.second->mass();
	}
	createAtlas();
	mesh.clear();
	for (auto& ptr : contents) {
		if (ptr.second == nullptr)continue;
		glm::dvec3 pos = ptr.first;
		bool flipped = glm::cross(ptr.second->front(), ptr.second->top()) == ptr.second->right();
		string path = "parts\\" + ptr.second->type + "_";
		glm::dvec3 normals[6] = {};
		for (int i = 0; i < 6; i++) {
			normals[i] = ptr.second->faceNormal(i);
		}
		for (int i = 0; i < 6; i++) {
			glm::dvec3 normal = normals[i];
			if (getItem(pos + normal) != nullptr)continue;
			glm::dvec3 w = normals[(i / 2 * 2 + 4) % 6],
				h = normals[(i / 2 * 2 + 2) % 6];
			int bf = ptr.second->faces()[i];
			glm::dvec3 pf = ptr.second->pointingFace(i);
			glm::dvec3 rp, rw, rh;
			if (i % 2 == 0) {
				rp = pos + normal;
				rw = w;
				rh = h;
			}
			else {
				rp = pos;
				rw = h;
				rh = w;
			}
			//this is a very long function
			glm::dvec2 ap, aw, ah, t1, t2, t3;
			t1 = atlasMap[path + faceNames[bf]];
			t2 = { 16,0 };
			t3 = { 0,16 };
			if (!flipped || !(bf > 1)) {
				t1 += t2;
				t2 *= -1;
			}
			if (pf == -rh) {
				ap = t1;
				aw = t2;
				ah = t3;
			}
			else if (pf == rw) {
				ap = t1 + t3;
				aw = -t3;
				ah = t2;
			}
			else if (pf == rh) {
				ap = t1 + t2 + t3;
				aw = -t2;
				ah = -t3;
			}
			else if (pf == -rw) {
				ap = t1 + t2;
				aw = t3;
				ah = -t2;
			}
			getRectangleMesh(mesh, rp, rw, rh, ap, aw, ah);
		}
	}
	brush.getMesh() = std::move(mesh);
	radius = 0;
	for (auto& ptr : contents) {
		if (ptr.second == nullptr)continue;
		glm::dvec3 pos = ptr.first + glm::dvec3(0.5, 0.5, 0.5);
		pos.x += glm::sign(pos.x) / 2;
		pos.y += glm::sign(pos.y) / 2;
		pos.z += glm::sign(pos.z) / 2;
		if (pos.x == 0)pos.x = 0.5;//edge case,sign()==0
		if (pos.y == 0)pos.y = 0.5;
		if (pos.z == 0)pos.z = 0.5;
		if (glm::length(pos) > radius)radius = glm::length(pos);
	}
}
glm::dvec3 PhysicsGrid::calculateContacts(glm::dvec3 sn) {
	contacts.clear();
	double alt = glm::length(position) - soi->radius * (1 + soi->terrain.get(position));
	double EP_CONTACT = 0.01;
	for (auto& ptr : contents) {
		for (auto& cc : corners_of_unit_cube) {
			glm::dvec3 localpos = ptr.first + glm::dvec3(0.5, 0.5, 0.5) - COM + cc;
			//note that here localpos is actually local,as opposed to the other copies of this loop
			//also note that because we are seeking all answers,
			//not just the one closest to the plane,we need to have an inner loop
			if (glm::dot(-sn, angle * localpos) > alt - EP_CONTACT) {
				contacts.push_back({ localpos,ptr.second->COF });
			}
		}
	}
	glm::dvec3 axis = getAlignedAxis(glm::inverse(angle) * sn, contacts);
	bool isSliding = contacts.size() >= 3 && !glm::any(glm::isnan(axis));
	if (!isSliding && situ != ShipSituation::FLYING) {
		situ = ShipSituation::FLYING;
	}
	return axis;//TODO:vary friction coefficient based on shape of contact(three points,two lines,one face)
}
vector<pair<glm::dvec3, glm::dvec3>> arrows;
glm::dvec3 PhysicsGrid::computeNFDistribution(glm::dvec3 kh, glm::dvec3 force, glm::dvec3 torque, glm::dvec3& ih, glm::dvec3& jh) {
	//calculates (affine) normal force distribution
	//also trims contacts,returns nan vector if invalid
	//force torque and normal(kh) must be local,i.e probably glm::inverse(angle)*force
	//output:{a,b,c}
	//also kh points into the planet i think?i.e probably -glm::inverse(angle)*sn

	if (glm::any(glm::isnan(getAlignedAxis(-kh, contacts))))return glm::dvec3(0) / 0.0;
	make_orthonormal_basis(kh, ih, jh);
	double z = glm::dot(contacts[0].first, kh);//all points have the same z coordinate
	double d = z * glm::dot(force, ih) - glm::dot(torque, jh);
	double e = glm::dot(torque, ih) + z * glm::dot(force, jh);
	double f = glm::dot(force, kh);
	double sum_x = 0, sum_y = 0;
	double sum_xx = 0, sum_yy = 0, sum_xy = 0;
	double n = static_cast<double>(contacts.size());
	for (auto& pair : contacts) {
		double x = glm::dot(pair.first, ih);
		double y = glm::dot(pair.first, jh);
		sum_x += x;
		sum_y += y;
		sum_xx += x * x;
		sum_yy += y * y;
		sum_xy += x * y;
	}
	glm::dmat3 J = {
		{sum_xx, sum_xy, sum_x},
		{sum_xy, sum_yy, sum_y},
		{sum_x, sum_y, n}
	};
	glm::dvec3 rhs = { d,e,f };
	glm::dvec3 abc = glm::inverse(J) * rhs;
	bool ok = true;
	for (auto& pair : contacts) {
		if (glm::dot(pair.first, ih) * abc.x + glm::dot(pair.first, jh) * abc.y + abc.z <= 0){
			ok = false;
			contacts.erase(remove(contacts.begin(), contacts.end(), pair), contacts.end());
		}
	}
	if (ok)return abc;
	return computeNFDistribution(kh, force, torque, ih, jh);
}
void PhysicsGrid::updatePhysics(double t, double dt) {
	//TODO:Many parts assume that only contact forces can be off-center
	glm::dvec3 gp = position + soi->apos(t);
	shared_ptr<Planet> p;
	double dist = std::numeric_limits<double>::infinity();
	for (auto& ptr : planets) {
		double nd = glm::length(ptr->apos(t) - gp);
		if (nd < dist) {
			dist = nd;
			p = ptr;
		}
	}
	while (dist > p->SOI()) {
		p = p->o->p;
		dist = glm::length(p->apos(t) - gp);
	}
	if (soi != p) {
		position = gp - p->apos(t);
		velocity += soi->avel(t) - p->avel(t);
		soi = p;
		orbit.p = p;
		orbit.set(position, velocity, t);
	}
	updateContactPhysics(t, dt);
	orbit.set(position, velocity, t);
}
void PhysicsGrid::updateContactPhysics(double t, double dt) {
	calculateContacts(soi->terrain.getSurfaceNormal(position));
	if (situ == ShipSituation::LANDED) {
		updatePhysicsLanded(t, dt);
	}
	else if (situ == ShipSituation::SlIDING) {
		updatePhysicsSliding(t, dt);
	}
	else {
		updatePhysicsFlying(t, dt);
	}
}
void PhysicsGrid::updatePhysicsLanded(double t, double dt) {
	glm::dvec3 sn = soi->terrain.getSurfaceNormal(position);//wf
	glm::dvec3 lsn = glm::inverse(angle) * sn;//lf
	auto [force, torque] = calculateLoads(t, dt);//[wf,lf]
	glm::dvec3 ih, jh;
	glm::dvec3 fd = computeNFDistribution(-lsn, glm::inverse(angle) * force, torque, ih, jh);//<ih,jh,kh>
	if (glm::any(glm::isnan(fd))) {
		situ = ShipSituation::FLYING;
		updatePhysicsFlying(t, dt);
		return;
	}
	glm::dvec3 nf = sn * glm::dot(sn, force);//wf
	glm::dvec3 tf = force - nf;//wf
	glm::dvec3 friction;
	for (auto& pair : contacts) {
		glm::dvec3 test = angle * torque, test2 = angle * (glm::cross(pair.first, torque) / glm::length2(pair.first));
		glm::dvec3 vel = force + angle * (glm::cross(pair.first, torque) / glm::length2(pair.first));//wf
		double ef = glm::length(vel - sn * glm::dot(sn, vel)) / contacts.size();//wf
		//temporary patch so friction doesn't overshoot
		//the friction limit for not overshooting should probably not be divided equally across contacts
		double f = (glm::dot(pair.first, ih) * fd.x + glm::dot(pair.first, jh) * fd.y + fd.z) * pair.second;//kh
		if (f * SF_BONUS < ef) situ = ShipSituation::SlIDING;
	}
	if (situ == ShipSituation::SlIDING) {
		updatePhysicsSliding(t, dt);
		return;
	}
	glm::dvec3 axis = calculateContacts(sn);
	velocity = glm::dvec3(0);
	avel = glm::dvec3(0);
	angle = angle * glm::rotation(-axis, lsn);
	double dist = soi->radius * (1 + soi->terrain.get(position)) + glm::dot(contacts[0].first, axis);
	position = dist * glm::normalize(position);
}
void PhysicsGrid::updatePhysicsSliding(double t, double dt) {
	glm::dvec3 sn = soi->terrain.getSurfaceNormal(position);//wf
	glm::dvec3 lsn = glm::inverse(angle) * sn;//lf
	auto [force, torque] = calculateLoads(t, dt);//[wf,lf]
	glm::dvec3 ih, jh;
	glm::dvec3 fd = computeNFDistribution(-lsn, glm::inverse(angle) * force, torque, ih, jh);//<ih,jh,kh>
	if (glm::any(glm::isnan(fd))) {
		situ = ShipSituation::FLYING;
		updatePhysicsFlying(t, dt);
		return;
	}
	glm::dvec3 nf = sn * glm::dot(sn, force);//wf
	glm::dvec3 tf = force - nf;//wf
	glm::dvec3 friction, ft;
	if (glm::length(avel) * radius < 0.01)situ = ShipSituation::LANDED;
	for (auto& pair : contacts) {
		glm::dvec3 vel = force + angle * (glm::cross(pair.first, torque) / glm::length2(pair.first)) +
			(velocity + angle * glm::cross(avel, pair.first)) * mass / dt;//wf
		double ef = glm::length(vel - sn * glm::dot(sn, vel)) / contacts.size();//wf
		//temporary patch so friction doesn't overshoot
		//the friction limit for not overshooting should probably not be divided equally across contacts
		double f = (glm::dot(pair.first, ih) * fd.x + glm::dot(pair.first, jh) * fd.y + fd.z) * pair.second;//kh
		if (f * SF_BONUS < ef) situ = ShipSituation::SlIDING;
		glm::dvec3 fr = -f * glm::normalize(vel);//wf
		fr = fr - sn * glm::dot(sn, fr);//wf
		glm::dvec3 gc = angle * pair.first;//wf
		glm::dvec3 sgc = gc - sn * glm::dot(sn, gc);//wf
		//glm::dvec3 nfr = sgc * glm::dot(sgc, fr) / glm::length2(sgc);
		//arrows.push_back({ gc,gc + nfr });
		//arrows.push_back({ gc,gc + fr - nfr });
		arrows.push_back({ gc,gc + fr });
		friction += fr;//wf
		ft += glm::cross(gc, fr);//wf
		//sometimes the ship tilts in some odd direction
		//sometimes it makes its position nan for some reason
	}
	if (situ == ShipSituation::SlIDING) {
		integrate(t, dt, tf + friction, lsn * glm::dot(lsn, torque) + glm::inverse(angle) * ft);
	}
	else {
		velocity = glm::dvec3(0);
		avel = glm::dvec3(0);
	}
	glm::dvec3 axis = calculateContacts(sn);
	velocity = velocity - sn * glm::dot(sn, velocity);
	avel = lsn * glm::dot(lsn, avel);
	angle = angle * glm::rotation(-axis, lsn);
	double dist = soi->radius * (1 + soi->terrain.get(position)) + glm::dot(contacts[0].first, axis);
	position = dist * glm::normalize(position);
}
void PhysicsGrid::updatePhysicsFlying(double t, double dt) {
	double H;
	glm::dvec3 sn = soi->terrain.getSurfaceNormal(position,H);//wf
	glm::dvec3 lsn = glm::inverse(angle) * sn;//lf
	double elapsed = 0;
	while (dt - elapsed > 1e-6) {
		//check for collision
		glm::dvec3 lp;//wf
		shared_ptr<GridElement> cpart;
		glm::dvec3 temp;
		double step = min(dt - elapsed, checkCollision(dt - elapsed, H, sn, lp, temp, cpart));
		auto [force, torque] = calculateLoads(t + elapsed, step);//[wf,lf]
		integrate(t + elapsed, step, force, torque);
		if (step == dt - elapsed) break;
		elapsed += step;
		glm::dvec3 axis = calculateContacts(sn);//lf
		glm::dvec3 vel = velocity + glm::cross(glm::rotate(angle, avel), lp);//wf
		glm::dvec3 ih, jh;
		glm::dvec3 fd = computeNFDistribution(-lsn, glm::inverse(angle) * force, torque, ih, jh);//<ih,jh,kh>
		if (
			glm::dot(sn, velocity) <= 0 &&
			glm::dot(sn, velocity) > -1e-2 &&
			!glm::any(glm::isnan(fd))
			) {
			situ = ShipSituation::SlIDING;
			updatePhysicsSliding(t + elapsed, dt - elapsed);
			break;
		}
		glm::dmat3 R = glm::toMat3(angle);
		glm::dmat3 it = R * inertialTensor * glm::transpose(R);
		temp = glm::inverse(it) * glm::cross(lp, sn);
		glm::dvec3 nf = (sn * (glm::dot(vel, sn) * -(1 + cpart->COR))) /
			(1 / mass + glm::dot(sn, glm::cross(temp, lp)));//wf
		velocity += nf / mass;
		avel += glm::inverse(inertialTensor) * (glm::inverse(angle) * glm::cross(lp, nf));
	}
}
void PhysicsGrid::integrate(double t, double dt, glm::dvec3 force, glm::dvec3 torque) {
	velocity += force * dt / mass;
	position += velocity * dt - force * (dt * dt / mass / 2);
	avel += glm::inverse(inertialTensor) * torque * dt;
	glm::dvec3 v = avel - (glm::inverse(inertialTensor) * torque) * (dt / 2.0);
	double l = glm::length(v);
	if (l != 0) {
		angle = angle * glm::angleAxis(l * dt, v / l);
		angle = glm::normalize(angle);
	}
}
pair<glm::dvec3, glm::dvec3> PhysicsGrid::calculateLoads(double t, double dt) {
	glm::dvec3 force, torque;
	force += position * (-mass * soi->gravity / pow(glm::length(position), 3));//gravity
	for (auto& ptr : contents) {
		force += angle * ptr.second->thrust;
		torque += glm::cross(ptr.first + glm::dvec3(0.5, 0.5, 0.5) - COM, ptr.second->thrust);
		torque += ptr.second->torque;
	}
	return { force,torque };
}
double PhysicsGrid::checkCollision(double maxDT, double H, glm::dvec3 sn, glm::dvec3& ret_position, glm::dvec3& contact_normal, shared_ptr<GridElement>& part) {
	//return +inf for no result,otherwise return time to collision
	//btw not to be called when situ is sliding or landed for obvious reasons
	double EP = 1e-8;
	glm::dvec3 pos, localpos;
	glm::dquat ang;
	double l_of_avel = glm::length(avel), speed = glm::length(velocity);
	glm::dvec3 cc;
	double alt, dts;
	double step = maxDT / 10;//artificially introduce tunnelling :(
	//maxDT is expected to be PHYSICS_DT
	//design the tunnelling around that
	//at the same time the divisor(here 10) represents how many steps the linear searching does
	//in total in the worst case
	//actually step should be tb/10 i think
	bool last = false;
	double ta = 0, tb = maxDT, time = tb;
	bool bisection = false;
	//bisect the time and do actual specific calculations to compute exact time
	//since we are thinking the planet as a plane locally anyways
	//wouldnt hurt to cut some costly get_terrain_height() calls
	pos = position + velocity * time;
	alt = glm::length(pos);
	do {
		//initialization
		pos = position + velocity * time;
		if (l_of_avel != 0) {
			ang = angle * glm::angleAxis(l_of_avel * time, avel / l_of_avel);
		}
		else {
			ang = angle;
		}
		//calculate planet plane normal and distance to COM
		//(not assuming the planet is locally a plane is impractical)
		alt = glm::length(pos);
		dts = alt - H * soi->radius;
		//closest corner to sn.effectively this is figuring out which octant sn points against
		//and finding the corner there
		cc = -0.5 * glm::sign(glm::inverse(ang) * sn);
		//loop over parts and calculate intersection
		for (auto& ptr : contents) {
			localpos = glm::rotate(ang, ptr.first + glm::dvec3(0.5, 0.5, 0.5) - COM + cc);
			if (glm::dot(-sn, localpos) < dts) continue;
			//a collision has happened
			//for the linear search it means it is okay to do bisection now if tb is time
			bisection = true;
			tb = time;
			ret_position = localpos;
			contact_normal = sn;
			part = ptr.second;
			break;
		}
		if (bisection) {
			if (tb != time)ta = time;
			time = (ta + tb) / 2;
		}
		else time -= step;
	} while (bisection ? (tb - ta > EP) : (time - step > ta));
	if (!bisection)return numeric_limits<double>::infinity();
	//check if collision is at or before 0
	if (ta != 0)return ta;
	//pull out of surface and try again
	do {
		position += sn * step;
		alt = glm::length(position);
		dts = alt - H * soi->radius;
		cc = -0.5 * glm::sign(glm::inverse(angle) * sn);
		bisection = false;
		for (auto& ptr : contents) {
			localpos = glm::rotate(angle, ptr.first + glm::dvec3(0.5, 0.5, 0.5) - COM + cc);
			if (glm::dot(-sn, localpos) < dts) continue;
			bisection = true;
			break;
		}
	} while (bisection);
	return checkCollision(maxDT, H, sn, ret_position, contact_normal, part);
}
double PhysicsGrid::totalEnergy() {
	//uses dry mass for calculation
	//idk if using dry mass for rotational kinetic energy makes any sense
	double m = 0;
	glm::dmat3 i;
	glm::dvec3 com;
	for (auto& ptr : contents) {
		if (ptr.second == nullptr)continue;
		m += ptr.second->dryMass;
		com += (ptr.first + glm::dvec3(0.5, 0.5, 0.5)) * ptr.second->dryMass;
	}
	com /= m;
	for (auto& ptr : contents) {
		if (ptr.second == nullptr)continue;
		glm::dvec3 pos = ptr.first;
		glm::dmat3 c(1.0 / 6);
		glm::dvec3 p = pos + glm::dvec3(0.5, 0.5, 0.5) - com;
		double x = p.x, y = p.y, z = p.z;
		glm::dmat3 s(
			glm::dvec3(y * y + z * z, -x * y, -x * z),
			glm::dvec3(-y * x, x * x + z * z, -y * z),
			glm::dvec3(-z * x, -z * y, x * x + y * y)
		);
		i += (c + s) * ptr.second->dryMass;
	}
	double KE = 0.5 * m * glm::dot(velocity, velocity);
	double rotKE = 0.5 * glm::dot(avel, i * avel);
	double GE = -soi->gravity * m / glm::length(position);
	return KE + rotKE + GE;
}
void PhysicsGrid::updateInternal(double t, double dt) {
	for (auto& ptr : contents) {
		if (ptr.second == nullptr)continue;
		ptr.second->update(ptr.first, t, dt);
	}
	for (auto& ptr : contents) {
		if (ptr.second == nullptr)continue;
		ptr.second->integrate(t, dt);
	}
	updateGrid();
}
int sceneDisplayed = 2;//0=instruments,1=map,2=craft part details,3=camera
void dragMap() {
	if (sceneDisplayed != 0 && sceneDisplayed != 2)return;
	if (!mouse[1] || dontDragMap)return;
	mapPos += (mousePos - pmousePos) / mapScale;
}
class ofApp : public ofBaseApp {

public:
	void setup();
	void update();
	void draw();

	void keyPressed(int key);
	void keyReleased(int key);
	void mouseMoved(int x, int y);
	void mouseDragged(int x, int y, int button);
	void mousePressed(int x, int y, int button);
	void mouseReleased(int x, int y, int button);
	void mouseEntered(int x, int y);
	void mouseExited(int x, int y);
	void mouseScrolled(int x, int y, float sx, float sy);
	void windowResized(int w, int h);
	void dragEvent(ofDragInfo dragInfo);
	void gotMessage(ofMessage msg);

};