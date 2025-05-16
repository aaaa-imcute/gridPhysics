#include "ofApp.h"

GridElement root("test", 1);
GridElement t("test2", 1);
PhysicsGrid p(make_shared<GridElement>(root), {680000,0,0}, { 0,0,4000/*2278.9316*/ },planets[0],0);
//more precise orbital velocity than 2279 so I don't confuse actual errors with
//p.orbit.a being about 680040
void ofApp::setup(){
	//ofSetFrameRate(60);
	//ofSetVerticalSync(true);
	ofDisableAntiAliasing();
	ofDisableBlendMode();
	createAtlas();
	planets[0]->terrain.generate(0, 5, 0.5, 1.5);
	//planets[0]->mesh = planets[0]->terrain.mesh({0,0,1},2);
	t.rotateRightFace();
	p.setItem(make_shared<GridElement>(t), 0, 1, 0);
}
static constexpr double DT = 1.0/600;
double remainingSimulation = 0.0;
double totalTime = 0.0;
void ofApp::update(){
	dragMap();
	if (sceneDisplayed != 2) {
		remainingSimulation += ofGetLastFrameTime();
		for (; remainingSimulation > DT; remainingSimulation -= DT) {
			p.updatePhysics(totalTime,DT);
			totalTime += DT;
		}
	}
}
void ofApp::draw(){
	pmousePos=mousePos;
	mousePos=glm::vec2(ofGetMouseX(), ofGetMouseY());
	if(sceneDisplayed == 1 || sceneDisplayed == 3)camera.update();
	if (keys['1'])sceneDisplayed = 1;
	if (keys['2'])sceneDisplayed = 2;
	if (keys['3'])sceneDisplayed = 3;
	switch(sceneDisplayed){
	case 1:
		camera.begin();
		ofEnableDepthTest();
		p.displayMode1(totalTime);
		for (auto& ptr : planets) {
			ptr->displayMode1(totalTime);
		}
		//ofDrawAxis(256);
		ofDisableDepthTest();
		camera.end();
		break;
	case 2:
		ofPushMatrix();
		ofTranslate(ofGetWidth() / 2, ofGetHeight() / 2);
		ofScale(mapScale);
		ofTranslate(mapPos);
		if (mouse[2])selectedPart = nullptr;
		p.displayMode2(0);
		//ofSetColor(255, 255, 255);
		//ofDrawCircle(untransform2D(mousePos), 4.0 / mapScale);
		ofPopMatrix();
		break;
	case 3:
		camera.begin();
		ofEnableDepthTest();
		p.displayMode3();
		ofDrawAxis(256);
		//ofDrawArrow({0,0,0}, p.avel*100.0, 20.0);
		ofDisableDepthTest();
		camera.end();
		break;
	}
}

void ofApp::keyPressed(int key){
	keys[key] = true;
	if (selectedPart != nullptr) {
		if (key == 'r')selectedPart->rotateFrontFace();
		if (key == 't')selectedPart->rotateTopFace();
		if (key == 'y')selectedPart->rotateRightFace();
		p.updateGrid();
	}
}

void ofApp::keyReleased(int key){
	keys[key] = false;
}


void ofApp::mouseMoved(int x, int y ){
}

void ofApp::mouseDragged(int x, int y, int button){

}

void ofApp::mousePressed(int x, int y, int button){
	mouse[button] = true;
}

void ofApp::mouseReleased(int x, int y, int button){
	mouse[button] = false;
}

void ofApp::mouseEntered(int x, int y){

}

void ofApp::mouseExited(int x, int y){

}

void ofApp::mouseScrolled(int x, int y, float sx, float sy) {
	if (sceneDisplayed == 0 || sceneDisplayed == 2)mapScale *= pow(2, (sy / 4));
	if (sceneDisplayed == 1)orbitScale *= pow(2, (sy / 4));
	if (sceneDisplayed == 3)camera.mouseScrolled(sx, sy);
}

void ofApp::windowResized(int w, int h){

}

void ofApp::gotMessage(ofMessage msg){

}

void ofApp::dragEvent(ofDragInfo dragInfo){ 

}

int main() {

	//Use ofGLFWWindowSettings for more options like multi-monitor fullscreen
	ofGLFWWindowSettings settings;
	settings.setSize(1920, 1080);
	settings.windowMode = OF_WINDOW;
	settings.decorated = false;
	auto window = ofCreateWindow(settings);
	ofRunApp(window, make_shared<ofApp>());
	ofRunMainLoop();

}
