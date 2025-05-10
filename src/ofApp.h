#pragma once

#include "ofMain.h"
#include "easyCam.cpp"
unordered_map<int, bool> keys;
unordered_map<int, bool> mouse;
glm::dvec2 untransform2D(glm::dvec2 vec) {
	glm::dvec2 temp = vec-glm::dvec2(ofGetWidth() / 2, ofGetHeight() / 2);
	glm::dmat4 modelView = ofGetCurrentMatrix(OF_MATRIX_MODELVIEW);
	glm::dmat4 inverse = glm::inverse(modelView);
	return glm::dvec4(temp.x, temp.y, 0, 1) * inverse;
}
glm::dvec2 mousePos, pmousePos;
glm::dvec2 mapPos;//sceneDisplayed==0||sceneDisplayed==2
double mapScale = 16;
easyCam camera;//sceneDisplayed==3
glm::dvec3 orbitPos;//sceneDisplayed==1
double orbitScale=1.0/6000000;
vector<string> faceNames = { "right","left","top","bottom","front","back" };
unordered_map<string, ofImage> imageCache;
ofImage loadImage(string img) {
	auto result = imageCache.find(img);
	ofImage temp;
	if (result == imageCache.end()) {
		if (!temp.load(img))throw "image " + img + " not loaded correctly";
		temp.getTexture().setTextureMinMagFilter(GL_NEAREST, GL_NEAREST);
		imageCache[img] = temp;
	}
	else {
		temp = result->second;
	}
	return temp;
}
void drawImage(string img, float x, float y) {
	ofImage temp = loadImage(img);
	ofSetColor(255,255,255);
	temp.draw(x, y);
}
ofFbo atlas;
unordered_map<string, glm::dvec2> atlasMap;
void walkTextures(string path, int& x, int& y) {
	ofDirectory dir(path);
	if (!dir.isDirectory()) {
		string n = dir.getOriginalDirectory();
		ofImage temp = loadImage(n.substr(0,n.size()-1));
		atlas.begin();
		temp.draw(x, y);
		atlas.end();
		//length of .\\textures\\ is 11 and length of .png is 4,also account for extra \\ at end
		atlasMap[n.substr(11, n.size() - 16)] = { (float)x, (float)y };
		x += 16;
		if (x == 4096) {
			x = 0;
			y += 16;
			if (y == 4096)throw "Too many textures";
		}
		return;
	}
	dir.listDir();
	for (int i = 0; i < dir.size(); i++) {
		walkTextures(dir.getPath(i), x, y);
	}
}
void createAtlas() {
	int SIZE = 1024;//so max.4096 textures
	atlas.allocate(SIZE, SIZE);
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
	ofDrawRectangle(- glm::dvec2(0, width), glm::length(p1-p2), width * 2);
	ofPopMatrix();
}
void drawFaceNormal2D(glm::dvec3 v, glm::dvec2 pos) {//sceneDisplay==2
	//expects color already set
	glm::dvec2 d(-v.z, -v.x);
	if (abs(v.y) < 0.01) {
		drawLine(pos, pos + d * 16,1);
	}
	else if (v.y > 0) {
		ofDrawCircle(pos, 2);
	}
	else if (v.y < 0) {
		drawLine(pos + glm::dvec2(-4, -4), pos + glm::dvec2(4, 4),1);
		drawLine(pos + glm::dvec2(-4, 4), pos + glm::dvec2(4, -4),1);
	}
}
static constexpr double DM3_SCALE = 256;
ofMesh getRectangleMesh(glm::dvec3 pos, glm::dvec3 width, glm::dvec3 height, glm::dvec2 apos, glm::dvec2 aw, glm::dvec2 ah) {
	double scale = DM3_SCALE;
	glm::dvec3 normal = glm::cross(height, width);
	ofMesh m;
	m.addVertex(pos * scale);
	m.addTexCoord(apos);
	m.addVertex((pos + width) * scale);
	m.addTexCoord(apos + aw);
	m.addVertex((pos + height) * scale);
	m.addTexCoord(apos + ah);
	m.addVertex((pos + width + height) * scale);
	m.addTexCoord(apos + aw + ah);
	m.addVertex((pos + width) * scale);
	m.addTexCoord(apos + aw);
	m.addVertex((pos + height) * scale);
	m.addTexCoord(apos + ah);
	for (int i = 0; i < 6; i++)m.addNormal(normal);
	m.addTriangle(0, 1, 2);
	m.addTriangle(3, 4, 5);
	return m;
}
class Planet;
class OrbitalElements;
class Planet {
public:
	shared_ptr<OrbitalElements> o;
	double gravity;//m^3/s^2
	double radius;
	Planet(double g, double r);
	Planet(double g, double r, shared_ptr<Planet> p1, double a1, double e1, double i1, double o1, double w1, double v1);
	glm::dvec3 pos(double t);
	glm::dvec3 vel(double t);
	glm::dvec3 apos(double t);
	glm::dvec3 avel(double t);
	double SOI();
	void displayMode1(double t);
};
class OrbitalElements {
public:
	double a, e, i, o, w, v;//v is mean anomaly at t=0 NOT TRUE ANOMALY,o is LAN(capital omega)
	shared_ptr<Planet> p;
	OrbitalElements(shared_ptr<Planet> p1,glm::dvec3 r,glm::dvec3 V,double t) {
		p = p1;
		set(r, V, t);
	}
	OrbitalElements(shared_ptr<Planet> p1, double a1, double e1, double i1, double o1, double w1, double v1)
		:p(p1), a(a1), e(e1), i(i1), o(o1), w(w1), v(v1) {};
	double period() {
		//yep this and meanMotion multiply to make 2PI
		//apart from the sign part which I basically invented to make hyperbolic
		//orbits have negative periods just like how they have negative semi major axes
		return 2 * PI * ofSign(a) * sqrt(abs(a * a * a) / p->gravity);
	}
	double meanMotion() {
		return sqrt(p->gravity / abs(a * a * a));
	}
	double trueAnomaly(double t) {
		double M = v + meanMotion() * t;
		if (e < 1) {
			double E = M + e * sin(M) + 0.5 * e * e * sin(2 * M);
			double dE = 10000;
			while (abs(dE) > 1e-16) {
				dE = (E - e * sin(E) - M) / (e * cos(E) - 1);
				E += dE;
			}
			return 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));
		}
		if (e > 1) {
			double H = ofSign(M) * log(2 * abs(M) / e + 1.8);
			double dH = 10000;
			while (abs(dH) > 1e-16) {
				dH = (e * sinh(H) - H - M) / (1 - e * cosh(H));
				H += dH;
			}
			return 2 * atan2(sqrt(e + 1) * sinh(H / 2), sqrt(e - 1) * cosh(H / 2));
		}
		throw "hey I'm not programming that case";
	}
	glm::dmat4 mat() {
		glm::dmat4 ri = glm::rotate(glm::dmat4(1.0f), glm::radians(i), glm::dvec3(1.0, 0.0, 0.0));
		glm::dmat4 ro = glm::rotate(glm::dmat4(1.0f), glm::radians(o), glm::dvec3(0.0f, 0.0f, 1.0f));
		glm::dmat4 rw = glm::rotate(glm::dmat4(1.0f), glm::radians(w), glm::dvec3(0.0f, 0.0f, 1.0f));
		return ro * ri * rw;
	}
	glm::dvec3 tpos(double n) {
		//input:true anomaly
		double r = a * (1 - e * e) / (1 + e * cos(n));
		glm::dvec4 ret(r * sin(n), r * cos(n), 0.0, 1.0);
		ret = mat() * ret;
		return glm::dvec3(ret);
	}
	glm::dvec3 pos(double t) {
		//input:time elapsed since t=0
		double n = trueAnomaly(t);
		return tpos(n);
	}
	glm::dvec3 vel(double t) {
		double ep = 1.0/60;
		return (pos(t+ep)-pos(t))/ep;
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
		glm::dvec3 h = glm::cross(r, V), N(-h.y, h.x, 0);
		glm::dvec3 E = glm::cross(V, h) / p->gravity - glm::normalize(r);
		double ep = glm::length2(V) / 2 - p->gravity / glm::length(r);
		a = -p->gravity / (2 * ep);
		e = glm::length(E);
		i = acos(h.z / glm::length(h));
		o = acos(N.x / glm::length(N));
		if (N.y < 0)o = 2 * PI - o;
		w = acos(glm::dot(N, E) / glm::length(glm::dot(N, E)));
		if (E.z < 0)w = 2 * PI - w;
		double n = acos(glm::dot(E, r) / glm::length(glm::dot(E, r)));
		if (glm::dot(r, V) < 0)n = 2 * PI - n;
		double M;
		if (e < 1) {
			M = atan2(sqrt(1 - e * e) * sin(n), e + cos(n)) - e * sqrt(1 - e * e) * sin(n) / (1 + e * cos(n));
			/*double cosE = (e + cos(n)) / (1 + e * cos(n));
			double sinE = sqrt(1 - e * e) * sin(n) / (1 + e * cos(n));
			double E_anom = atan2(sinE, cosE);
			M = E_anom - e * sin(E_anom);*/
		}
		else {
			double F = 2 * atanh(tan(n / 2) * sqrt((e - 1) / (e + 1)));
			M = e * sinh(F) - F;
			/*double coshF = (e + cos(n)) / (1 + e * cos(n));
			double sinhF = sqrt(e * e - 1) * sin(n) / (1 + e * cos(n));
			double F = asinh(sinhF);
			M = e * sinh(F) - F;*/
		}
		v = M;// - t * meanMotion();//TODO:Fix mean anomaly at epoch calculation
	}
	void displayMode1(double t) {
		double div = 2 * PI / 256;
		ofNoFill();
		ofSetColor(0, 255, 0);
		ofBeginShape();
		for (double i = 0; i < 2 * PI; i += div) {
			ofVertex(glm::vec3((tpos(i) + p->apos(t) - orbitPos) * orbitScale));
		}
		ofEndShape(true);
	}
};
vector<shared_ptr<Planet>> planets = { make_shared<Planet>(Planet(3.5316e12,600000)) };
Planet::Planet(double g, double r) {
	gravity = g;
	radius = r;
}
Planet::Planet(double g, double r, shared_ptr<Planet> p1, double a1, double e1, double i1, double o1, double w1, double v1) {
	gravity = g;
	radius = r;
	o = make_shared<OrbitalElements>(OrbitalElements(p1, a1, e1, i1, o1, w1, v1));
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
void Planet::displayMode1(double t) {
	glm::dvec3 p = (apos(t) - orbitPos) * orbitScale;
	ofPushMatrix();
	ofTranslate(glm::vec3(p));
	ofSetColor(127, 127, 127);//TODO
	ofDrawSphere(radius * orbitScale);
	ofTranslate(-glm::vec3(p));
	ofPopMatrix();
	if(o!=nullptr)o->displayMode1(t);
}
class GridElement {
public:
	GridElement(string t, double m) {
		type = t;
		mass = m;
		frontFace = 4;
		topFace = 2;
		rightFace = 0;
	}
	string type;
	double mass;
	int frontFace, topFace, rightFace;
	glm::dvec3 faceNormal(int face) {//returns ship coordinates
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
	glm::dvec3 front() {//returns normal vector of front face(in ship coordinates)
		return faceNormal(frontFace);
	}
	glm::dvec3 top() {
		//alright,now we need to determine where the top face is.
		//it cannot be the front face or the back face.
		int i = 0, j=topFace;
		for (; i <= j; i++) {
			if (i == frontFace || i == (frontFace ^ 1))j++;
		}
		return faceNormal(j);
	}
	glm::dvec3 right() {
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
	vector<int> faces() {
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
	glm::dvec3 pointingFace(int p) {
		//what face "up" points to on a face
		//front if not front or back,
		//front points to top,back points to bottom
		//input is in craft coordinates NOT part coordinates
		vector<int> f = faces();
		if (f[p] < 4) return front();
		if (f[p] == 4)return top();
		return -top();
	}
	void rotateFrontFace() {
		frontFace++;
		frontFace %= 6;
	}
	void rotateTopFace() {
		topFace++;
		topFace %= 4;
	}
	void rotateRightFace() {
		rightFace++;
		rightFace %= 2;
	}
	void displayMode2() {
		vector<int> f = faces();
		glm::dvec3 p = pointingFace(2);
		glm::dvec2 n(-p.z, -p.x);
		ofPushMatrix();
		ofTranslate(8, 8);
		ofRotateRad(atan2f(n.y, n.x) + PI / 2);
		if (glm::cross(front(),top()) == right() && f[2] > 1)ofScale(-1, 1);
		//problem:does not consider float imprecision,though in theory it doesn't need to
		ofTranslate(-8, -8);
		drawImage(".\\textures\\parts\\" + type + "_" + faceNames[f[2]] + ".png", 0, 0);
		ofPopMatrix();
	}
	glm::dvec3 getThrust(double dt) {
		//returns impulse from center of this part in dt
		//only call once per time interval(because it removes fuel among other things)
		return { 0,0,0 };
	}
};
shared_ptr<GridElement> selectedPart;//sceneDisplayed==2
glm::dvec3 selectedPos;
class PhysicsGrid {
public:
	shared_ptr<Planet> soi;
	glm::dvec3 position;//need high-precision
	glm::dvec3 velocity;//doesn't need high-precision but unfortunately glm::dvec3 are floats
	//bottleneck on fixedpoint_prec is the handling of small accelerations
	glm::dvec3 accel;//verlet internal
	glm::dquat angle;
	glm::dvec3 avel;//along axis of rotation,magnitude is amount and direction of rotation
	OrbitalElements orbit;
	PhysicsGrid(shared_ptr<GridElement> root, glm::dvec3 p, glm::dvec3 v,shared_ptr<Planet> planet,double t)
	:orbit(planet,p,v,t){
		//soi = planets[0];//assuming this is the one holding the global frame of reference
		//if (soi->o != nullptr)throw "hey wrong order of planets";
		soi = planet;
		position = p;
		velocity = v;
		accel = { 0,0,0 };
		avel = { 0,0,0 };
		angle = { 1,0,0,0 };
		setItem(root, 0, 0, 0);
		brush = { mesh };
		updateGrid();
	}
	shared_ptr<GridElement> getItem(int x, int y, int z) {
		return getItem(glm::dvec3(x,y,z));
	}
	shared_ptr<GridElement> getItem(glm::dvec3 index) {
		string i = keyString(index);
		if (contents.find(i) == contents.end())return nullptr;
		return contents[keyString(index)];
	}
	shared_ptr<GridElement> setItem(shared_ptr<GridElement> e, int x, int y, int z) {
		return setItem(e, glm::dvec3(x, y, z));
	}
	shared_ptr<GridElement> setItem(shared_ptr<GridElement> e, glm::dvec3 index) {
		contents[keyString(index)] = e;
		updateGrid();
		return e;
	}
	void removeItem(int x, int y, int z) {
		removeItem({ (float)x,(float)y,(float)z });
	}
	void removeItem(glm::dvec3 v) {
		contents.erase(keyString(v));
		updateGrid();
	}
	void displayMode1(double t) {
		glm::dvec3 pos = (orbit.apos(t) - orbitPos) * orbitScale;
		ofPushMatrix();
		ofTranslate(glm::vec3(pos));//draw ship
		ofSetColor(127, 127, 127);
		ofDrawSphere(16);
		ofTranslate(-glm::vec3(pos));
		ofPopMatrix();
		orbit.displayMode1(t);
	}
	void displayMode2(int y) {
		for (auto& ptr : contents) {
			glm::dvec3 pos = keyPos(ptr.first);
			if (round(pos.y) != y) continue;
			//in this scene,x+ goes up while z+ goes left
			glm::dvec3 sPos(-pos.z * 16, -pos.x * 16, 0);
			ofPushMatrix();
			ofTranslate(glm::vec3(sPos));
			ptr.second->displayMode2();
			ofPopMatrix();
			glm::dvec2 m = untransform2D(mousePos);
			if (mouse[2] && m.x > sPos.x && m.x<sPos.x + 16 && m.y>sPos.y && m.y < sPos.y + 16) {
				selectedPart = ptr.second;
				selectedPos = pos;
			}
		}
		if (selectedPart != nullptr) {
			glm::dvec2 sPos(-selectedPos.z * 16, -selectedPos.x * 16);
			ofSetColor(0, 0, 255);
			drawFaceNormal2D(selectedPart->front(), sPos + glm::dvec2(8, 8));
			ofSetColor(0, 255, 0);
			drawFaceNormal2D(selectedPart->top(), sPos + glm::dvec2(8, 8));
			ofSetColor(255, 0, 0);
			drawFaceNormal2D(selectedPart->right(), sPos + glm::dvec2(8, 8));
		}
	}
	void displayMode3() {
		brush.setGlobalOrientation(glm::dquat(1, 0, 0, 0));
		brush.setGlobalPosition(0, 0, 0);
		ofSetColor(255, 255, 255);
		ofPushMatrix();
		ofScale(1, 1, -1);
		glm::mat4 matr = glm::toMat4(angle);
		ofMultMatrix(matr);
		ofTranslate(glm::vec3(-COM)*DM3_SCALE);
		ofTexture t = atlas.getTexture();
		t.bind();
		brush.draw();
		t.unbind();
		ofPopMatrix();
	}
	void updateGrid() {
		if (contents.size() == 0)throw "empty ship";
		mass = 0;
		COM = glm::dvec3(0, 0, 0);//TODO:Fix COM calculation
		for (auto& ptr : contents) {
			if (ptr.second == nullptr)continue;
			mass += ptr.second->mass;
			COM += (keyPos(ptr.first)+glm::dvec3(0.5,0.5,0.5)) * ptr.second->mass;
		}
		COM /= mass;
		inertialTensor = glm::dmat3(0.0);
		for (auto& ptr : contents) {
			if (ptr.second == nullptr)continue;
			glm::dvec3 pos = keyPos(ptr.first);
			glm::dmat3 c(1.0 / 6);
			//should also have side length squared but that is 1.mass is accounted for later
			glm::dvec3 p = pos + glm::dvec3(0.5, 0.5, 0.5)-COM;
			double x = p.x, y = p.y, z = p.z;
			glm::dmat3 s(
				glm::vec3(y * y + z * z, -x * y, -x * z),  // First row
				glm::vec3(-y * x, x * x + z * z, -y * z),  // Second row
				glm::vec3(-z * x, -z * y, x * x + y * y)   // Third row
			);
			inertialTensor += (c + s) * ptr.second->mass;
		}
		mesh.clear();
		for (auto& ptr : contents) {
			if (ptr.second == nullptr)continue;
			glm::dvec3 pos = keyPos(ptr.first);
			bool flipped = glm::cross(ptr.second->front(), ptr.second->top()) == ptr.second->right();
			for (int i = 0; i < 6; i++) {
				glm::dvec3 normal = ptr.second->faceNormal(i);
				if (getItem(pos + normal) != nullptr)continue;
				glm::dvec3 w = ptr.second->faceNormal((i / 2 * 2 + 4) % 6),
					h = ptr.second->faceNormal((i / 2 * 2 + 2) % 6);
				int bf = ptr.second->faces()[i];
				string fn = faceNames[bf];
				string an = "parts\\" + ptr.second->type + "_" + fn;
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
				glm::dvec2 ap, aw, ah,t1,t2,t3;
				t1 = atlasMap[an];
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
				mesh.append(getRectangleMesh(rp, rw, rh, ap, aw, ah));
			}
		}
		brush = { mesh };

	}
	void updatePhysics(double t, double dt) {
		/*
		set<shared_ptr<Planet>> p;
		for (auto& ptr : planets) {
			glm::dvec3 pos = ptr->pos(t);
			if (glm::length(position - pos) < ptr->SOI()) {
				p.insert(ptr);
			}
		}
		*/
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
		glm::dvec3 torque;
		//torque = { 0.1,0,0 };
		//TODO
		avel += glm::inverse(inertialTensor) * torque * dt;
		double l = glm::length(avel) * (dt / 2);
		if (l != 0) {
			glm::dvec3 n = glm::normalize(avel);
			glm::dquat d(0, n.x, n.y, n.z);
			d *= sin(l);
			d.w = cos(l);
			angle = d * angle;
			angle = glm::normalize(angle);
		}
		glm::dvec3 acc;
		position = position + velocity * dt + glm::dvec3(acc * (dt * dt / 2));
		velocity = velocity + accel * dt;//estimate(for drag and other speed related forces)
		//TODO
		acc += position * (-soi->gravity / pow(glm::length(position),3));//gravity
		velocity = velocity + (acc - accel) * (dt / 2.0);//factor in actual acceleration
		accel = acc;
		orbit.set(position, velocity, t);
	}
private:
	unordered_map<string,shared_ptr<GridElement>> contents;
	//y points down,z points right,x points forwards in that order
	double mass;
	glm::dvec3 COM;//ship coordinates
	glm::dmat3 inertialTensor;
	of3dPrimitive brush;
	ofMesh mesh;
	string keyString(glm::dvec3 index) {
		return keyString((int)round(index.x), (int)round(index.y), (int)round(index.z));
	}
	string keyString(int x, int y, int z) {
		return to_string(x) + "," + to_string(y) + "," + to_string(z);
	}
	glm::dvec3 keyPos(string str) {
		stringstream s(str);
		string temp;
		glm::dvec3 res;
		getline(s, temp, ',');
		res.x = stod(temp);
		getline(s, temp, ',');
		res.y = stod(temp);
		getline(s, temp, ',');
		res.z = stod(temp);
		return res;
	}
};
int sceneDisplayed = 2;//0=instruments,1=map,2=craft part details,3=camera(future)
void dragMap() {
	if (sceneDisplayed != 0 && sceneDisplayed != 2)return;
	if (!mouse[0])return;
	mapPos += (mousePos - pmousePos)/mapScale;
}
class ofApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
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
