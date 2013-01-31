#include "testApp.h"

//--------------------------------------------------------------
void testApp::setup(){
    simulator = new Simulator();
    simulator->InitializeGrid(128, 96, 40);
    simulator->AddParticles();
    
    glEnable(GL_DEPTH_TEST);
    cam.setDistance(96);
    
    enableFog(100, 250);
    ofBackground(0);
}

//--------------------------------------------------------------
void testApp::update(){\
}

//--------------------------------------------------------------
void testApp::draw(){
    simulator->startThread();
    Particle *particles = simulator->particles;
    cam.begin();
    ofPushMatrix();
    ofTranslate(-64, -48, -20);
    glBegin(GL_POINTS);
    for (int i = 0; i < simulator->nParticles; i++) {
        Particle &p = particles[i];
        float v = fminf(sqrtf(p.u[0]*p.u[0]+p.u[1]*p.u[1]+p.u[2]*p.u[2]), 1);
        glColor3f(v, .6+.4*v, 1);
        glVertex3f(p.x[0], p.x[1], p.x[2]);
    }
    glEnd();
    ofPopMatrix();
    cam.end();
    
    ofDrawBitmapString(ofToString(ofGetFrameRate()), 10, 10);
    simulator->waitForThread();
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

void testApp::enableFog(float near, float far) {
    glEnable(GL_FOG);
    glFogi(GL_FOG_MODE, GL_LINEAR);
    GLfloat fogColor[4]= {0, 0, 0, 1};
    glFogfv(GL_FOG_COLOR, fogColor);
    glHint(GL_FOG_HINT, GL_FASTEST);
    glFogf(GL_FOG_START, near);
    glFogf(GL_FOG_END, far);
}