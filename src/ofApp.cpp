#include "ofApp.h"

GridElement root("test", 1);
PhysicsGrid p(make_shared<GridElement>(root), { 0,0,0 }, { 0,0,0 });
void ofApp::setup(){
	//ofSetFrameRate(60);
	//ofSetVerticalSync(true);
	p.getItem(0, 0, 0);
}

void ofApp::update(){
	dragMap();
}

void ofApp::draw(){
	pmousePos.set(mousePos);
	mousePos.set(mouseX, mouseY);
	ofDrawBitmapString(to_string(mapPos.x) + "," + to_string(mapPos.y) + "," + to_string(mapScale), 0, 100);
	ofPushMatrix();
	ofTranslate(ofGetWidth() / 2, ofGetHeight() / 2);
	ofScale(mapScale);
	ofTranslate(mapPos);
	p.displayMode2(0);
	ofPopMatrix();
}

void ofApp::keyPressed(int key){
	keys[key] = true;
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
	mapScale *= pow(2,(sy / 4));
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
