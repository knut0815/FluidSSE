#pragma once

#include "ofMain.h"
#include "Simulator.h"
#include "ofxSSAO.h"

class testApp : public ofBaseApp{
    Simulator *simulator;
    ofEasyCam cam;
    ofImage image;
    
    ofFbo pointsFbo;
    ofxSSAO ssao;
	public:
		void setup();
		void update();
		void draw();

		void keyPressed  (int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
        void enableFog(float near, float far);
        void disableFog();
};
