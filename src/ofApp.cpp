#include "ofApp.h"
shared_ptr<PhysicsGrid> p;
double initialEnergy;
//orbital speed at 680000=2278.9316
void ofApp::setup(){
	compute_legendre_coeff();
	ofDisableAntiAliasing();
	ofDisableBlendMode();
	if (!font.load("cour.ttf", 12, true, true, true, 0.0f, 96))throw "bad font";
	sunLight = make_shared<ofLight>(ofLight());
	sunLight->setPointLight();
	sunLight->setSpecularColor(ofColor::white);
	sunLight->setPosition(0, 10000, 0);
	camera3.setFarClip(16777215.0);
	planets[0]->color = [&](double h)->ofColor {
		if (h <= 0)return ofColor(0, 0, 255);
		ofColor beach = ofColor(240, 220, 130);     // sandy yellow
		ofColor plains = ofColor(70, 180, 90);      // rich green
		ofColor mountains = ofColor(255, 255, 255); // snow white
		if (h < 0.5) {
			return beach.getLerped(plains, h * 2);
		}
		return plains.getLerped(mountains, h * 2 - 1);
		};
	planets[0]->terrain.generate(436, 0.1, 2);
	planetAtlasReady = false;//do this after every change made to the terrain or color scheme of a planet
	createPlanetAtlas();//for some reason this doesnt work in the display methods of the planet class
	//it only works here
	//TODO:find out why
	GridElement t1("metal-tank", 1, 0.7, 0.15);
	GridElement t2("metal-tank", 1, 0.7, 0.15);
	GridElement t3("solid-rocket-engine", 1, 0.7, 0.15);
	PhysicsGrid p1(make_shared<GridElement>(t1), { 574134.9,0,0 }, { 0,0,10.0 }, planets[0], 0);
	p = make_shared<PhysicsGrid>(p1);
	p->setItem(make_shared<GridElement>(t2), 0, 1, 0);
	p->setItem(make_shared<GridElement>(t3), 0, -1, 0);
	p->getItem(0, -1, 0)->rotateFrontFace();
	p->getItem(0, -1, 0)->rotateFrontFace();
	p->getItem(0, -1, 0)->rotateFrontFace();
	p->getItem(0, -1, 0)->rotateFrontFace();
	p->getItem(0, -1, 0)->rotateFrontFace();
	p->getItem(0, 0, 0)->rotateFrontFace();
	p->getItem(0, 0, 0)->rotateFrontFace();
	p->getItem(0, 0, 0)->rotateFrontFace();
	p->getItem(0, 0, 0)->rotateFrontFace();
	p->getItem(0, 0, 0)->rotateFrontFace();
	p->getItem(0, 1, 0)->rotateFrontFace();
	p->getItem(0, 1, 0)->rotateFrontFace();
	p->getItem(0, 1, 0)->rotateFrontFace();
	p->getItem(0, 1, 0)->rotateFrontFace();
	p->getItem(0, 1, 0)->rotateFrontFace();
	p->updateGrid();
	initialEnergy = p->totalEnergy();
	/*auto spec = p->getItem(0, -1, 0)->engineData();
	visit(
		[&](const auto& spec) {
			using T = std::decay_t<decltype(spec)>;
			if constexpr (std::is_same_v<T, MonoPropSpec>) {
				double test = 0.0;
				test = calculateEjectionVelocity(spec, 101325, spec.maxR);
				test *= 1.0;
			}
		}, *spec
	);*/
}
constexpr double PHYSICS_DT = 1.0/600;
double remainingSimulation = 0.0;
double totalTime = 0.0;
void ofApp::update(){
	dragMap();
	if (sceneDisplayed != 2/* && keys[' ']*/) {
		remainingSimulation += ofGetLastFrameTime();
		for (; remainingSimulation > PHYSICS_DT; remainingSimulation -= PHYSICS_DT) {
			p->updatePhysics(totalTime,PHYSICS_DT);
			p->updateInternal(totalTime, PHYSICS_DT);
			totalTime += PHYSICS_DT;
		}
	}
	/*if (keys['r']) {
		p->angle = glm::inverse(camera.getGlobalOrientation());
	}*/
	
}
void ofApp::draw(){
	pmousePos=mousePos;
	mousePos=glm::vec2(ofGetMouseX(), ofGetMouseY());
	if(sceneDisplayed == 1)camera1.update();
	if(sceneDisplayed == 3)camera3.update();
	if (keys['1'])sceneDisplayed = 1;
	if (keys['2'])sceneDisplayed = 2;
	if (keys['3'])sceneDisplayed = 3;
	switch(sceneDisplayed){
	case 1:
		ofEnableLighting();
		sunLight->enable();
		camera1.begin();
		ofEnableDepthTest();
		for (int i = 0; i < planets.size();i++) {
			planets[i]->displayMode1(totalTime,i);
		}
		ofDisableLighting();
		sunLight->disable();
		p->displayMode1(totalTime);
		//ofDrawAxis(256);
		ofDisableDepthTest();
		camera1.end();
		break;
	case 2:
		ofPushMatrix();
		ofTranslate(ofGetWidth() / 2, ofGetHeight() / 2);
		ofScale(mapScale);
		ofTranslate(mapPos);
		//atlas.draw(0, 0);
		p->displayMode2(dm2_layer);
		//ofSetColor(255, 255, 255);
		//ofDrawCircle(untransform2D(mousePos), 4.0 / mapScale);
		ofPopMatrix();
		break;
	case 3:
		//ofEnableLighting();
		//sunLight->enable();
		ofPushMatrix();
		camera3.begin();
		ofEnableDepthTest();
		ofScale(1, 1, -1);
		for (int i = 0; i < planets.size(); i++) {
			if (p->soi != planets[i])continue;
			planets[i]->displayMode3(p->position, i);
		}
		p->displayMode3();
		ofDrawAxis(256);
		glm::dvec3 test = planets[0]->terrain.getSurfaceNormal(p->position);
		double d = glm::length(p->position) - planets[0]->radius * (1 + planets[0]->terrain.get(p->position));
		//drawPlaneWithNormal(test, d*DM3_SCALE, 1000);
		for (auto& pair : p->contacts) {
			ofSetColor(255, 255, 255);
			ofDrawSphere(DM3_SCALE * (p->angle * pair.first), 16);
		}
		for (auto& pair : arrows) {
			ofSetColor(255, 255, 0);
			ofDrawArrow(pair.first * DM3_SCALE, pair.second * DM3_SCALE);
		}
		arrows.clear();
		ofDisableDepthTest();
		camera3.end();
		ofPopMatrix();
		//ofDisableLighting();
		//sunLight->disable();
		ofPushStyle();
		ofSetColor(255, 255, 0);
		font.drawString(
			to_string(ofGetFrameRate()) + "\n" +
			printSituation(p->situ) + "\n" +
			to_string(glm::length(p->position)) + "\n" +
			to_string(glm::length(p->velocity)) + "\n" +
			to_string(p->totalEnergy()-initialEnergy) + "\n" +
			to_string(planets[0]->radius * (planets[0]->terrain.get(p->position) + 1))
			, 100, 100);
		ofPopStyle();
		break;
	}
}

void ofApp::keyPressed(int key){
	keys[key] = true;
	if (sceneDisplayed == 2) {
		if (selectedPart != nullptr && round(selectedPos.y) == dm2_layer) {
			if (key == 'r')selectedPart->rotateFrontFace();
			if (key == 't')selectedPart->rotateTopFace();
			if (key == 'y')selectedPart->rotateRightFace();
			p->updateGrid();
		}
		if (key == 'w')dm2_layer++;
		if (key == 's')dm2_layer--;
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
	if (sceneDisplayed == 2)p->mousePressed(x, y, button);
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
	if (sceneDisplayed == 3)camera3.mouseScrolled(sx, sy);
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
