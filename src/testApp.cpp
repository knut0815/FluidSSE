#include "testApp.h"

//--------------------------------------------------------------
void testApp::setup(){
    simulator.InitializeGrid(100, 100, 100);
    simulator.AddParticles();
    
    glEnable(GL_DEPTH_TEST);
    cam.setDistance(150);
}

//--------------------------------------------------------------
void testApp::update(){
    simulator.Update();
}

//--------------------------------------------------------------
void testApp::draw(){
    Particle *particles = simulator.particles;
    cam.begin();
    ofPushMatrix();
    ofTranslate(-50, -50, -50);
    glBegin(GL_POINTS);
    for (int i = 0; i < simulator.nParticles; i++) {
        Particle &p = particles[i];
        glVertex3f(p.x[0], p.x[1], p.x[2]);
    }
    glEnd();
    ofPopMatrix();
    cam.end();
    
    ofDrawBitmapString(ofToString(ofGetFrameRate()), 10, 10);
}

//--------------------------------------------------------------
void testApp::keyPressed(int key){

}

//--------------------------------------------------------------
void testApp::keyReleased(int key){

}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){ 

}