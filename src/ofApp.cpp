#include "ofApp.h"

GridElement root("test", 1);
PhysicsGrid p(make_shared<GridElement>(root), { 0,0,0 }, { 0,0,0 });
void ofApp::setup(){
	ofSetFrameRate(60);
	p.getItem(0, 0, 0);
}

void ofApp::update(){
	dragMap();
}

void ofApp::draw(){
	pmousePos.set(mousePos);
	mousePos.set(mouseX, mouseY);
	ofPushMatrix();
	ofTranslate(mapPos);
	ofScale(mapScale);
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

void ofApp::windowResized(int w, int h){

}

void ofApp::gotMessage(ofMessage msg){

}

void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
