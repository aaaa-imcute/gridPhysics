#pragma once

#include "ofMain.h"

unordered_map<int, bool> keys;
unordered_map<int, bool> mouse;
ofVec2f untransform2D(ofVec2f vec) {
	ofVec2f temp = vec-ofVec2f(ofGetWidth() / 2, ofGetHeight() / 2);
	ofMatrix4x4 modelView = ofGetCurrentMatrix(OF_MATRIX_MODELVIEW);
	ofMatrix4x4 inverse = modelView.getInverse();
	return ofVec4f(temp.x, temp.y, 0, 1) * inverse;
}
ofVec2f mousePos, pmousePos;
ofVec2f mapPos;//sceneDisplayed==0||sceneDisplayed==2
double mapScale = 16;
ofEasyCam camera;//sceneDisplayed==1||sceneDisplayed==3
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
unordered_map<string, ofVec2f> atlasMap;
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
void drawLine(ofVec2f p1, ofVec2f p2, float width) {//replacement method for line width
	ofPushMatrix();
	ofTranslate(p1);
	ofRotateRad(atan2f(p2.y - p1.y, p2.x - p1.x));
	ofDrawRectangle(- ofVec2f(0, width), p1.distance(p2), width * 2);
	ofPopMatrix();
}
void drawFaceNormal2D(ofVec3f v, ofVec2f pos) {//sceneDisplay==2
	ofVec2f d(-v.z, -v.x);
	if (abs(v.y) < 0.01) {
		drawLine(pos, pos + d * 16,1);
	}
	else if (v.y > 0) {
		ofDrawCircle(pos, 2);
	}
	else if (v.y < 0) {
		drawLine(pos + ofVec2f(-4, -4), pos + ofVec2f(4, 4),1);
		drawLine(pos + ofVec2f(-4, 4), pos + ofVec2f(4, -4),1);
	}
}
ofMesh getRectangleMesh(ofVec3f pos, ofVec3f width, ofVec3f height, ofVec2f apos, ofVec2f aw, ofVec2f ah) {
	float scale = 256;
	ofVec3f normal = height.getCrossed(width);
	ofMesh m;
	m.addVertex(ofPoint(pos * scale));
	m.addTexCoord(apos);
	m.addVertex(ofPoint((pos + width) * scale));
	m.addTexCoord(apos + aw);
	m.addVertex(ofPoint((pos + height) * scale));
	m.addTexCoord(apos + ah);
	m.addVertex(ofPoint((pos + width + height) * scale));
	m.addTexCoord(apos + aw + ah);
	m.addVertex(ofPoint((pos + width) * scale));
	m.addTexCoord(apos + aw);
	m.addVertex(ofPoint((pos + height) * scale));
	m.addTexCoord(apos + ah);
	for (int i = 0; i < 6; i++)m.addNormal(normal);
	m.addTriangle(0, 1, 2);
	m.addTriangle(3, 4, 5);
	return m;
}
static constexpr int FIXEDPOINT_PREC = 1000;
class Vec3 {
	//ofVec3f but much more precision
	//todo:angle stuff
private:
	long long _x, _y, _z;
	class fixedPoint {
		long long& internal;
	public:
		fixedPoint(long long& ref) :internal(ref) {}
		operator double() const {
			return double(internal) / FIXEDPOINT_PREC;
		}
		fixedPoint& operator=(double o) {
			internal = o * FIXEDPOINT_PREC;
			return *this;
		}
	};
	Vec3(long long x1, long long y1, long long z1, bool raw) : _x(0), _y(0), _z(0), x(_x), y(_y), z(_z) {
		_x = x1;
		_y = y1;
		_z = z1;
	}
public:
	fixedPoint x, y, z;
	Vec3() :_x(0), _y(0), _z(0), x(_x), y(_y), z(_z) {}
	Vec3(double x1, double y1, double z1) :_x(0), _y(0), _z(0), x(_x), y(_y), z(_z) {
		x = x1;
		y = y1;
		z = z1;
	}
	Vec3(ofVec3f v) :_x(0), _y(0), _z(0), x(_x), y(_y), z(_z) {
		x = v.x;
		y = v.y;
		z = v.z;
	}
	Vec3 operator=(const Vec3& o) {
		_x = o._x;
		_y = o._y;
		_z = o._z;
		return *this;
	}
	Vec3 operator+(Vec3 o) {
		return Vec3(_x + o._x, _y + o._y, _z + o._z, 1);
	}
	Vec3 operator-(Vec3 o) {
		return Vec3(_x - o._x, _y - o._y, _z - o._z, 1);
	}
	Vec3 operator*(double o) {
		return Vec3(_x * o, _y * o, _z * o, 1);
	}
	Vec3 operator/(double o) {
		return Vec3(_x / o, _y / o, _z / o, 1);
	}
	double sqMag() {
		return (*this) * (*this);
	}
	double mag() {
		return sqrt(sqMag());
	}
	double operator*(Vec3 o) {
		return x * o.x + y * o.y + z * o.z;
	}
	Vec3 cross(ofVec3f o) {//used to determine orbital angular momentum
		return Vec3(y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x);
	}
	ofVec3f downgrade() {//only downgrade if the result does not represent positions relative to the world origin
		return ofVec3f(x, y, z);
	}
};
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
	ofVec3f faceNormal(int face) {//returns ship coordinates
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
	}
	ofVec3f front() {//returns normal vector of front face(in ship coordinates)
		return faceNormal(frontFace);
	}
	ofVec3f top() {
		//alright,now we need to determine where the top face is.
		//it cannot be the front face or the back face.
		int i = 0, j=topFace;
		for (; i <= j; i++) {
			if (i == frontFace || i == (frontFace ^ 1))j++;
		}
		return faceNormal(j);
	}
	ofVec3f right() {
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
	ofVec3f pointingFace(int p) {
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
		ofVec3f p = pointingFace(2);
		ofVec2f n(-p.z, -p.x);
		ofPushMatrix();
		ofTranslate(8, 8);
		ofRotateRad(atan2f(n.y, n.x) + PI / 2);
		if (front().getCrossed(top()) == right() && f[2] > 1)ofScale(-1, 1);
		//problem:does not consider float imprecision,though in theory it doesn't need to
		ofTranslate(-8, -8);
		drawImage(".\\textures\\parts\\" + type + "_" + faceNames[f[2]] + ".png", 0, 0);
		ofPopMatrix();
	}
	ofVec3f getThrust(double dt) {
		//returns impulse from center of this part in dt
		//only call once per time interval(because it removes fuel among other things)
		return { 0,0,0 };
	}
};
shared_ptr<GridElement> selectedPart;//sceneDisplayed==2
ofVec3f selectedPos;
class PhysicsGrid {
public:
	Vec3 position;//need high-precision
	Vec3 velocity;//doesn't need high-precision but unfortunately ofvec3f are floats
	ofVec3f accel;//verlet internal
	ofQuaternion angle;
	ofVec3f avel;//along axis of rotation,magnitude is amount and direction of rotation
	PhysicsGrid(shared_ptr<GridElement> root,Vec3 p, Vec3 v) {
		position = p;
		velocity = v;
		setItem(root, 0, 0, 0);
		brush = { mesh };
		updateGrid();
	}
	shared_ptr<GridElement> getItem(int x, int y, int z) {
		return getItem(ofVec3f(x,y,z));
	}
	shared_ptr<GridElement> getItem(ofVec3f index) {
		string i = keyString(index);
		if (contents.find(i) == contents.end())return nullptr;
		return contents[keyString(index)];
	}
	shared_ptr<GridElement> setItem(shared_ptr<GridElement> e, int x, int y, int z) {
		return setItem(e, ofVec3f(x, y, z));
	}
	shared_ptr<GridElement> setItem(shared_ptr<GridElement> e, ofVec3f index) {
		contents[keyString(index)] = e;
		updateGrid();
		return e;
	}
	void removeItem(int x, int y, int z) {
		removeItem({ (float)x,(float)y,(float)z });
	}
	void removeItem(ofVec3f v) {
		contents.erase(keyString(v));
		updateGrid();
	}
	void displayMode2(int y) {
		for (auto& ptr : contents) {
			ofVec3f pos = keyPos(ptr.first);
			if (round(pos.y) != y) continue;
			//in this scene,x+ goes up while z+ goes left
			ofPoint sPos(-pos.z*16, -pos.x*16);
			ofPushMatrix();
			ofTranslate(sPos);
			ptr.second->displayMode2();
			ofPopMatrix();
			ofVec2f m = untransform2D(mousePos);
			if (mouse[2] && m.x > sPos.x && m.x<sPos.x + 16 && m.y>sPos.y && m.y < sPos.y + 16) {
				selectedPart = ptr.second;
				selectedPos = pos;
			}
		}
		if (selectedPart != nullptr) {
			ofPoint sPos(-selectedPos.z * 16, -selectedPos.x * 16);
			ofSetColor(0, 0, 255);
			drawFaceNormal2D(selectedPart->front(), sPos + ofVec2f(8, 8));
			ofSetColor(0, 255, 0);
			drawFaceNormal2D(selectedPart->top(), sPos + ofVec2f(8, 8));
			ofSetColor(255, 0, 0);
			drawFaceNormal2D(selectedPart->right(), sPos + ofVec2f(8, 8));
		}
	}
	void displayMode3() {
		brush.setGlobalOrientation(angle);
		brush.setGlobalPosition(0, 0, 0);
		ofSetColor(255, 255, 255);
		ofPushMatrix();
		ofScale(1, 1, -1);
		ofTexture t = atlas.getTexture();
		t.bind();
		brush.draw();
		t.unbind();
		ofPopMatrix();
	}
	void updateGrid() {
		mass = 0;
		COM.set(0, 0, 0);
		for (auto& ptr : contents) {
			if (ptr.second == nullptr)continue;
			mass += ptr.second->mass;
			COM += keyPos(ptr.first) * ptr.second->mass;
		}
		//TODO
		//blah blah blah inertial tensor blah blah
		mesh.clear();
		for (auto& ptr : contents) {
			if (ptr.second == nullptr)continue;
			ofVec3f pos = keyPos(ptr.first);
			bool flipped = ptr.second->front().getCrossed(ptr.second->top()) == ptr.second->right();
			for (int i = 0; i < 6; i++) {
				ofVec3f normal = ptr.second->faceNormal(i);
				if (getItem(pos + normal) != nullptr)continue;
				ofVec3f w = ptr.second->faceNormal((i / 2 * 2 + 4) % 6),
					h = ptr.second->faceNormal((i / 2 * 2 + 2) % 6);
				int bf = ptr.second->faces()[i];
				string fn = faceNames[bf];
				string an = "parts\\" + ptr.second->type + "_" + fn;
				ofVec3f pf = ptr.second->pointingFace(i);
				ofVec3f rp, rw, rh;
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
				ofVec2f ap, aw, ah,t1,t2,t3;
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
	void updatePhysics(double dt) {
		ofVec3f torque;
		//TODO
		avel += inertialTensor.getInverse()*ofVec4f(torque)*dt;
		angle += ofQuaternion(avel.x,avel.y,avel.z,0)*angle*(dt/2);
		ofVec3f acc;
		position = position + velocity * dt + Vec3(acc * (dt * dt / 2));
		velocity = velocity + Vec3(accel);//estimate(for drag and other speed related forces)
		//TODO
		velocity = velocity + Vec3((acc - accel) / 2);//factor in actual acceleration
		accel = acc;
	}
private:
	unordered_map<string,shared_ptr<GridElement>> contents;
	//y points down,z points right,x points forwards in that order
	double mass;
	ofVec3f COM;//ship coordinates
	ofMatrix4x4 inertialTensor;
	of3dPrimitive brush;
	ofMesh mesh;
	string keyString(ofVec3f index) {
		return keyString((int)round(index.x), (int)round(index.y), (int)round(index.z));
	}
	string keyString(int x, int y, int z) {
		return to_string(x) + "," + to_string(y) + "," + to_string(z);
	}
	ofVec3f keyPos(string str) {
		stringstream s(str);
		string temp;
		ofVec3f res;
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
