#pragma once

#include "ofMain.h"
static const unsigned long doubleclickTime = 200;
class easyCam : public ofCamera {//can of worms.
public:
	bool disableScroll = false;
	easyCam() {
		reset();
		sensitivityTranslate = { 1,1,1 };
		sensitivityRot = { 1,1,1 };

		addInteraction(TRANSFORM_TRANSLATE_XY, OF_MOUSE_BUTTON_LEFT, doTranslationKey);
		addInteraction(TRANSFORM_ROTATE, OF_MOUSE_BUTTON_LEFT);
		addInteraction(TRANSFORM_TRANSLATE_Z, OF_MOUSE_BUTTON_RIGHT);
		addInteraction(TRANSFORM_TRANSLATE_XY, OF_MOUSE_BUTTON_MIDDLE);

	}
	virtual void begin(const ofRectangle& _viewport) {
		if (!bEventsSet) {
			setEvents(ofEvents());
		}
		viewport = _viewport;
		ofCamera::begin(viewport);
	}
	virtual void begin() {
		begin(getViewport());
	}
	void reset() {
		target.resetTransform();

		target.setPosition(0, 0, 0);
		lookAt(target);

		resetTransform();
		setPosition(0, 0, lastDistance);

		rot = { 0,0,0 };
		translate = { 0,0,0 };

		if (bAutoDistance) {
			bDistanceSet = false;
		}
		currentTransformType = TRANSFORM_NONE;
	}
	void setTargetsetTarget(const glm::vec3& targetPoint) {
		target.setPosition(targetPoint);
		lookAt(target);
	}
	void setTargetsetTarget(ofNode& targetNode) {
		target = targetNode;
		lookAt(target);
	}
	const ofNode& getTarget() const {
		return target;
	}
	void setDistance(float distance) {
		setDistance(distance, true);
	}
	float getDistance() const {
		return glm::distance(target.getPosition(), getPosition());
	}
	void setAutoDistance(bool bAutoDistance) {
		this->bAutoDistance = bAutoDistance;
		if (bAutoDistance) {
			bDistanceSet = false;
		}
	}

	void setEvents(ofCoreEvents& _events) {
		// If en/disableMouseInput were called within ofApp::setup(),
		// bMouseInputEnabled will tell us about whether the camera
		// mouse input needs to be initialised as enabled or disabled.
		// we will still set `events`, so that subsequent enabling
		// and disabling can work.

		// we need a temporary copy of bMouseInputEnabled, since it will 
		// get changed by disableMouseInput as a side-effect.
		bool wasMouseInputEnabled = bMouseInputEnabled;// || !events;
		disableMouseInput();
		events = &_events;
		if (wasMouseInputEnabled) {
			// note: this will set bMouseInputEnabled to true as a side-effect.
			enableMouseInput();
		}
		bEventsSet = true;
	}
	void setRotationSensitivity(const glm::vec3& sensitivity) {
		sensitivityRot = sensitivity;
	}
	void setRotationSensitivity(float x, float y, float z) {
		setRotationSensitivity({ x,y,z });
	}
	void setTranslationSensitivity(const glm::vec3& sensitivity) {
		sensitivityTranslate = sensitivity;
	}
	void setTranslationSensitivity(float x, float y, float z) {
		sensitivityTranslate = { x,y,z };
	}
	void setTranslationKey(char key) {
		doTranslationKey = key;
	}
	char getTranslationKey() const {
		return doTranslationKey;
	}
	void enableMouseInput() {
		if (!bMouseInputEnabled && events) {
			listeners.push(events->update.newListener(this, &easyCam::update));
			listeners.push(events->mousePressed.newListener(this, &easyCam::mousePressed));
			listeners.push(events->mouseReleased.newListener(this, &easyCam::mouseReleased));
			listeners.push(events->mouseScrolled.newListener(this, &easyCam::mouseScrolled));
		}
		// if enableMouseInput was called within ofApp::setup()
		// `events` will still carry a null pointer, and bad things
		// will happen. Therefore we only update the flag. 
		bMouseInputEnabled = true;
		// setEvents() is called upon first load, and will make sure 
		// to enable the mouse input once the camera is fully loaded.
	}
	void disableMouseInput() {
		if (bMouseInputEnabled && events) {
			listeners.unsubscribeAll();
		}
		// if disableMouseInput was called within ofApp::setup()
		// `events` will still carry a null pointer, and bad things
		// will happen. Therefore we only update the flag. 
		bMouseInputEnabled = false;
		// setEvents() is called upon first load, and will make sure 
		// to enable the mouse input once the camera is fully loaded.
	}
	bool getMouseInputEnabled() const {
		return bMouseInputEnabled;
	}
	void enableMouseMiddleButton() {
		bEnableMouseMiddleButton = true;
	}
	void disableMouseMiddleButton() {
		bEnableMouseMiddleButton = false;
	}
	bool getMouseMiddleButtonEnabled() const {
		return bEnableMouseMiddleButton;
	}
	void setRelativeYAxis(bool relative) {
		bRelativeYAxis = relative;
	}
	bool getRelativeYAxis() const {
		return bRelativeYAxis;
	}
	void setUpAxis(const glm::vec3& _up) {
		upAxis = _up;
	}
	const glm::vec3& getUpAxis() const {
		return upAxis;
	}
	enum TransformType {
		TRANSFORM_NONE,
		TRANSFORM_ROTATE,
		TRANSFORM_TRANSLATE_XY,
		TRANSFORM_TRANSLATE_Z,
		TRANSFORM_SCALE
	};
	void addInteraction(TransformType type, int mouseButton, int key=-1) {
		if (!hasInteraction(mouseButton, key)) {
			interactions.push_back(interaction(type, mouseButton, key));
		}
		else {
			ofLogNotice("easyCam") << "Can not add interaction. It already exists";
		}
	}
	void removeInteraction(TransformType type, int mouseButton, int key=-1) {
		ofRemove(interactions, [&](interaction& i) { return i.transformType == type && i.mouseButton == mouseButton && i.key == key; });
	}
	bool hasInteraction(int mouseButton, int key=-1) {
		for (const auto& i : interactions) {
			if (i.mouseButton == mouseButton && i.key == key) {
				return true;
			}
		}
		return false;
	}
	bool hasInteraction(TransformType type, int mouseButton, int key) {
		for (const auto& i : interactions) {
			if (i.transformType == type && i.mouseButton == mouseButton && i.key == key) {
				return true;
			}
		}
		return false;
	}
	void removeAllInteractions() {
		interactions.clear();
	}
protected:
	virtual void onPositionChanged() {
		if (!bDistanceSet && bAutoDistance) {
			bDistanceSet = true;
		}
	}

private:
	void setDistance(float distance, bool save) {//should this be the distance from the camera to the target?
		if (distance > 0.0f) {
			if (save) {
				this->lastDistance = distance;
			}
			setPosition(target.getPosition() + (distance * getZAxis()));
			bDistanceSet = true;
		}
	}

	ofNode target;

	bool bEnableMouseMiddleButton = true;
	bool bMouseInputEnabled = true;
	bool bDistanceSet = false;
	bool bAutoDistance = true;
	bool bEventsSet = false;
	bool bIsScrolling = false;
	float lastDistance = 0.f;
	glm::vec3 rot;
	glm::vec3 translate;
	glm::vec3 sensitivityTranslate;
	glm::vec3 sensitivityRot;
	float     sensitivityScroll = 1.0f;
	glm::vec2 prevMouse;
	glm::vec2 mouseVel;
	void updateRotation() {
		curRot = glm::angleAxis(rot.z, lastPressAxisZ) * glm::angleAxis(rot.y, up()) * glm::angleAxis(rot.x, lastPressAxisX);
		setOrientation(curRot * lastPressOrientation);
		setPosition(curRot * (lastPressPosition - target.getGlobalPosition()) + target.getGlobalPosition());
	}
	void updateTranslation() {
		if (currentTransformType == TRANSFORM_TRANSLATE_XY ||
			currentTransformType == TRANSFORM_TRANSLATE_Z ||
			currentTransformType == TRANSFORM_SCALE) {
			if (getOrtho()) {
				//In ortho mode moving along the z axis has no effect besides clipping.
				// Instead, scale is applied to achieve the effect of getting near or far from the target.
				glm::vec3 mousePre;
				bool bDoScale = (currentTransformType == TRANSFORM_SCALE || currentTransformType == TRANSFORM_TRANSLATE_Z);
				if (bDoScale) {
					mousePre = screenToWorld(glm::vec3((bIsScrolling ? mouseAtScroll : lastPressMouse), 0));
				}
				move(glm::vec3(lastPressAxisX * translate.x) + (lastPressAxisY * translate.y));
				if (bDoScale) {
					setScale(getScale() + translate.z);
					// this move call is to keep the scaling centered below the mouse.
					move(mousePre - screenToWorld(glm::vec3((bIsScrolling ? mouseAtScroll : lastPressMouse), 0)));
				}
			}
			else {
				move(glm::vec3(lastPressAxisX * translate.x) + (lastPressAxisY * translate.y) + (lastPressAxisZ * translate.z));
			}
		}
		if (bIsScrolling) {
			//this it to avoid the transformation to keep on after scrolling ended.
			currentTransformType = TRANSFORM_NONE;
			bIsScrolling = false;
		}
	}
	void update(ofEventArgs& args) {
		if (this->viewport.isZero()) {
			viewport = getViewport();
		}
		if (!bDistanceSet && bAutoDistance) {
			setDistance(getImagePlaneDistance(viewport), true);
		}
		if (bMouseInputEnabled && events) {
			if (events->getMousePressed()) {
				updateMouse(glm::vec2(events->getMouseX(), events->getMouseY()));
			}
			if (currentTransformType == TRANSFORM_ROTATE) {
				updateRotation();
			}
			else if (currentTransformType == TRANSFORM_TRANSLATE_XY ||
				currentTransformType == TRANSFORM_TRANSLATE_Z ||
				currentTransformType == TRANSFORM_SCALE) {
				updateTranslation();
			}
		}
	}
	void mousePressed(ofMouseEventArgs& mouse) {
		lastPressMouse = mouse;
		prevMouse = mouse;
		lastPressAxisX = getXAxis();
		lastPressAxisY = getYAxis();
		lastPressAxisZ = getZAxis();
		lastPressPosition = ofCamera::getGlobalPosition();
		lastPressOrientation = ofCamera::getGlobalOrientation();

		currentTransformType = TRANSFORM_NONE;
		if (events) {
			for (const auto& i : interactions) {
				if (i.mouseButton == mouse.button && ((i.key == -1) || events->getKeyPressed(i.key))) {
					currentTransformType = i.transformType;
					break;
				}
			}
		}
	}
	void mouseReleased(ofMouseEventArgs& mouse) {
		// Check if it's double click
		unsigned long curTap = ofGetElapsedTimeMillis();
		if (lastTap != 0 && curTap - lastTap < doubleclickTime) {
			reset();
			return;
		}
		lastTap = curTap;
		currentTransformType = TRANSFORM_NONE;
		rot = { 0,0,0 };
		translate = { 0,0,0 };
	}
	void mouseScrolled(ofMouseEventArgs& mouse) {
		if (disableScroll)return;
		mouseVel = mouse - prevMouse;
		prevMouse = mouse;
		lastPressPosition = ofCamera::getGlobalPosition();
		lastPressAxisZ = getZAxis();
		if (getOrtho()) {
			translate.z = sensitivityScroll * mouse.scrollY / viewport.height;
			mouseAtScroll = mouse;
		}
		else {
			translate.z = mouse.scrollY * 30 * sensitivityTranslate.z * (getDistance() + std::numeric_limits<float>::epsilon()) / ofGetHeight();
		}
		currentTransformType = TRANSFORM_SCALE;
		bIsScrolling = true;
	}
	void updateMouse(const glm::vec2& mouse) {
		int vFlip = (isVFlipped() ? -1 : 1);

		rot = { 0,0,0 };
		translate = { 0,0,0 };
		switch (currentTransformType) {
		case TRANSFORM_ROTATE:
			mouseVel = mouse - lastPressMouse;
			rot.x = vFlip * -mouseVel.y * sensitivityRot.x * glm::pi<float>() / std::min(ofGetWidth(), ofGetHeight());
			rot.y = -mouseVel.x * sensitivityRot.y * glm::pi<float>() / std::min(ofGetWidth(), ofGetHeight());
			break;
		case TRANSFORM_TRANSLATE_XY:
			mouseVel = mouse - prevMouse;
			if (getOrtho()) {
				translate.x = -mouseVel.x * getScale().z;
				translate.y = vFlip * mouseVel.y * getScale().z;
			}
			else {
				translate.x = -mouseVel.x * sensitivityTranslate.x * 0.5f * (getDistance() + std::numeric_limits<float>::epsilon()) / ofGetWidth();
				translate.y = vFlip * mouseVel.y * sensitivityTranslate.y * 0.5f * (getDistance() + std::numeric_limits<float>::epsilon()) / ofGetHeight();
			}
			break;
		case TRANSFORM_TRANSLATE_Z:
			mouseVel = mouse - prevMouse;
			if (getOrtho()) {
				translate.z = mouseVel.y * sensitivityScroll / ofGetHeight();
			}
			else {
				translate.z = mouseVel.y * (sensitivityTranslate.z * 0.7f) * (getDistance() + std::numeric_limits<float>::epsilon()) / ofGetHeight();
			}
			break;
		default:
			break;
		}
		prevMouse = mouse;
	}
	glm::vec3 up() const {
		if (bRelativeYAxis) {
			return lastPressAxisY;
		}
		else {
			return upAxis;
		}
	}
	char doTranslationKey = 'm';
	unsigned long lastTap = 0;
	glm::quat curRot;
	glm::vec3 lastPressAxisX;
	glm::vec3 lastPressAxisY;
	glm::vec3 lastPressAxisZ;
	glm::vec3 lastPressPosition;
	glm::quat lastPressOrientation;
	glm::vec2 lastPressMouse;
	ofRectangle viewport;
	ofRectangle controlArea;
	ofEventListeners listeners;
	ofCoreEvents* events = nullptr;
	bool bRelativeYAxis = false;
	bool doInertia = false;
	glm::vec3 upAxis{ 0,1,0 };
	glm::vec2 mouseAtScroll;
	TransformType currentTransformType;
	struct interaction {
		interaction() :mouseButton(0), key(-1), transformType(TRANSFORM_NONE) {}
		interaction(TransformType type, int _mouseButton, int _key = -1) :mouseButton(_mouseButton), key(_key), transformType(type) {}
		int mouseButton;
		int key;
		TransformType transformType;
	};
	std::vector<interaction> interactions;
};
