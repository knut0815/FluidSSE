#include "testApp.h"

//--------------------------------------------------------------
void testApp::setup(){
    simulator = new Simulator();
    simulator->InitializeGrid(256, 144, 120);
    //simulator->AddParticles();
    
    glEnable(GL_DEPTH_TEST);
    cam.setDistance(180);
    
    
    ofBackground(0);
    
    ofFbo::Settings fboSettings;
    fboSettings.width = ofGetWidth();
    fboSettings.height = ofGetHeight();
    fboSettings.internalformat = GL_RGBA;
    fboSettings.useDepth = true;
    fboSettings.depthStencilAsTexture = true;
    fboSettings.depthStencilInternalFormat = GL_DEPTH_COMPONENT32;
    pointsFbo.allocate(fboSettings);
    
    ssao.setup();
    
    cam.setNearClip(100);
    cam.setFarClip(280);
}

//--------------------------------------------------------------
void testApp::update(){
}

//--------------------------------------------------------------
void testApp::draw(){
    simulator->startThread();
    Particle *particles = simulator->particles;
    enableFog(120, 280);
    pointsFbo.begin();
    ofClear(0);
    cam.begin();
    ofPushMatrix();
    ofTranslate(-128, -72, -60);
    glBegin(GL_POINTS);
    for (int i = 0; i < simulator->nParticles; i++) {
        Particle &p = particles[i];
        float v = fminf(.5*sqrtf(p.u[0]*p.u[0]+p.u[1]*p.u[1]+p.u[2]*p.u[2]), 1);
        int mod = i/100%4;
        if (mod == 0) {
            glColor3f(v, .5+.5*v, 1);
        } else if (mod == 1) {
            glColor3f(v, 1, .5+.5*v);
        } else if (mod == 2) {
            glColor3f(1, .5+.5*v, v);
        } else {
            glColor3f(1, v, .5+.5*v);
        }
        glVertex3f(p.x[0], p.x[1], p.x[2]);
    }
    glEnd();
    ofPopMatrix();
    cam.end();
    pointsFbo.end();
    disableFog();
    pointsFbo.getDepthTexture().draw(0, 0);
    //pointsFbo.draw(0, 0);
    ofDrawBitmapString(ofToString(ofGetFrameRate()), 10, 10);
    //image.grabScreen(0, 0, ofGetWidth(), ofGetHeight());
    //image.saveImage("frame"+ofToString(ofGetFrameNum())+".jpg");
    simulator->waitForThread();
    if (ofGetFrameNum()>2000) {
        exit();
    }
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
void testApp::disableFog() {
    glDisable(GL_FOG);
}