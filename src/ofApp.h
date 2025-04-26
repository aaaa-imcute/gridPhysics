#pragma once

#include "ofMain.h"

unordered_map<int, bool> keys;
unordered_map<int, bool> mouse;
ofVec2f mousePos, pmousePos;
ofVec2f mapPos;//sceneDisplayed==0||sceneDisplayed==2
double mapScale;
class Vec3 {
	//ofVec3f but double precision
	//todo:angle stuff
public:
	double x, y, z;
	Vec3() {
		x = 0;
		y = 0;
		z = 0;
	}
	Vec3(double x1, double y1, double z1) {
		x = x1;
		y = y1;
		z = z1;
	}
	Vec3 operator+(Vec3 o) {
		return Vec3(x + o.x, y + o.y, z + o.z);
	}
	Vec3 operator-(Vec3 o) {
		return Vec3(x - o.x, y - o.y, z - o.z);
	}
	Vec3 operator*(double o) {
		return Vec3(x * o, y * o, z * o);
	}
	Vec3 operator/(double o) {
		return Vec3(x / o, y / o, z / o);
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
	Vec3 cross(Vec3 o) {
		return Vec3(y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x);
	}
};
class GridElement {
public:
	GridElement(string t, double m) {
		type = t;
		mass = m;
	}
	string type;
	double mass;
	void displayMode2() {
		//ofDrawBitmapString(type, 0, 0);
		ofDrawRectangle(0, 0, 1, 1);
	}
};
class PhysicsGrid {
public:
	Vec3 pos;//world position of origin
	Vec3 vel;
	ofQuaternion angle = { 1,0,0,0 };
	ofQuaternion avel = { 1,0,0,0 };
	PhysicsGrid(shared_ptr<GridElement> root,Vec3 p, Vec3 v) {
		pos = p;
		vel = v;
		setItem(root, 0, 0, 0);
		updateGrid();
	}
	shared_ptr<GridElement> getItem(int x, int y, int z) {
		return getItem(ofVec3f(x,y,z));
	}
	shared_ptr<GridElement> getItem(ofVec3f index) {
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
	void displayMode2(int y) {
		for (auto& ptr : contents) {
			ofVec3f pos = keyPos(ptr.first);
			if (round(pos.y) != y) continue;
			ofPoint sPos(pos.x*16, pos.z*16);
			ofPushMatrix();
			ofTranslate(sPos);
			ofScale(16);
			ptr.second->displayMode2();
			ofPopMatrix();
		}
	}
private:
	unordered_map<string,shared_ptr<GridElement>> contents;
	//y points down,z points right,x points forwards in that order
	double mass;
	ofVec3f COM;//logical coordinates
	string keyString(ofVec3f index) {
		return keyString((int)index.x, (int)index.y, (int)index.z);
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
		res.y = stod(temp);
		return res;
	}
	void updateGrid() {
		mass = 0;
		COM.set(0, 0, 0);
		for(auto& ptr:contents){
			mass += ptr.second->mass;
			COM += keyPos(ptr.first) * ptr.second->mass;
		}
	}
};
int sceneDisplayed = 2;//0=instruments,1=map,2=craft part details,3=camera(future)
void dragMap() {
	if (sceneDisplayed != 0 && sceneDisplayed != 2)return;
	mapPos += mousePos - pmousePos;
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
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
		
};
