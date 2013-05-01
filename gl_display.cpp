#include <vector>
#include <deque>
#include <string>
#include <sstream>
#include <iomanip>
#include "freeglut.h"
#pragma comment( lib, "freeglut.lib" )

#include "particles.h"
#include "bitmapUtil.h"

using namespace particles;
const int Windowsize = 600;
const int PointNum = 50;
particleClass pc;
std::deque<std::vector<std::pair<double,double> > > locus;
bool showLocus = true;
char windowPixels[3*Windowsize*Windowsize];
int tickCount = 0;

void locusAdd(){
	locus.push_front(pc.vertices());
	if(locus.size() > PointNum){
		locus.pop_back();
	}
}

void timertick(int value)
{
	const int repeat = 100;
	const int sleepms = 30;
	const double speedup = 100;
	for(int i = 0; i < repeat; i++){
		pc.move(sleepms / 1000.0 / static_cast<double>(repeat) * speedup);
		locusAdd();
	}
    glutPostRedisplay();
    glutTimerFunc(sleepms, timertick, 0);

	//glReadPixels(0, 0, Windowsize, Windowsize, GL_BGR_EXT, GL_UNSIGNED_BYTE, windowPixels);
	//std::stringstream ss;
	//ss << "bmp/case3/case_" << std::setw(10) << std::setfill('0') << tickCount << ".bmp";
	//saveBMP(ss.str(), Windowsize, Windowsize, windowPixels);

	tickCount++;
	/*
	if(tickCount > 7000){
		exit(0);
	}
	*/
}


void keydown(unsigned char key, int x, int y){
	switch(key){
	case '1':
	case '2':
	case '3':
	case '4':
		pc.setIntegrate(key - '0');
		break;
	case 'i':
		pc.initialize();
		locus.clear();
		break;
	case 'l':
		showLocus = !showLocus;
		break;
	case 'q':
		exit(0);
	}
}

void display_string(float x, float y, std::string str){
	glDisable(GL_TEXTURE_2D);
	glColor4f(0.7f,0.8f,0.8f,1.0f);
	glRasterPos3f(x, y, -1.0f);
	const char* p = str.c_str();
	while (*p != '\0') glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *p++);
 	glEnable(GL_TEXTURE_2D);
}

void display_params(){
	std::stringstream ss;

	display_string(-0.9f, 0.9f, "<Keys> 1/2/3:Methods, i:init, l:locus q:quit");

	ss.setf(std::ios::showpoint);
	ss.str("");
	ss  << "Current Method: " << pc.getIntegrateName();
	display_string(-0.9f, 0.82f, ss.str());
	ss.str("");
	ss << "energy: "
		<< std::setprecision(10) << std::setw(6) << pc.energyDiff()*1.0e2;
	display_string(-0.9f, 0.74f, ss.str());
	ss.str("");
	ss << "Moment: "
		<< std::setprecision(10) << std::setw(6) << pc.mx()*1.0e16 << "  " << pc.my()*1.0e16;
	display_string(-0.9f, 0.66f, ss.str());
	ss.str("");
	ss << "Count: "
		<< tickCount;
	display_string(-0.9f, 0.58f, ss.str());
}

void display(void)
{
	const double dR = 1.0;
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	float colors[3][3] = {{0.275f, 0.510f, 0.706f},{0.196f, 0.804f, 0.196f},{0.933f, 0.902f, 0.522f}};

	if(showLocus){
		glPointSize(1.0f);
		glLineWidth(1.0f);
		for(unsigned int i = 0; i < locus.front().size(); i++){
			glBegin(GL_LINE_STRIP);
			int pCount = 0;
			for(auto it = locus.rbegin(); it != locus.rend(); ++it){
				float colorMul = static_cast<float>(static_cast<double>(pCount)/static_cast<double>(locus.size()));
				glColor3f(colors[i%3][0]*colorMul, colors[i%3][1]*colorMul, colors[i%3][2]*colorMul);
				glVertex2d(it->at(i).first*dR, it->at(i).second*dR);
				pCount++;
			}
			glEnd();
		}
	}
	glPointSize(5.0f);
	glColor3f(1.000f, 1.000f, 1.000f);
	glBegin(GL_POINTS);
	for(unsigned int i = 0; i < locus.front().size(); i++){
		glVertex2d(locus.front().at(i).first*dR, locus.front().at(i).second*dR);
	}
	glEnd();
	
	display_params();

	glutSwapBuffers();
}
int main(int argc, char** argv){
	std::vector<std::pair<double, partState> > particles;
	particles.push_back(std::pair<double, partState>(5000.0, partState(0.0, 0.0, 0.0, 0.0)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.2, 1.3, 0.4)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.24, 1.2, 0.4)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.28, 1.1, 0.4)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.32, 1.0, 0.4)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.36, 0.9, 0.4)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.4, 0.8, 0.4)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.44, 0.7, 0.4)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.48, 0.7, 0.4)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.52, 0.7, 0.4)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.56, 0.7, 0.4)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.6, 0.7, 0.4)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.64, 0.7, 0.4)));
	particles.push_back(std::pair<double, partState>(1.0, partState(0.0, -0.68, 0.7, 0.4)));
	pc = particleClass(particles);
	pc.setIntegrate(4);
	locusAdd();

    glutInitWindowPosition(100, 100);
	glutInitWindowSize(Windowsize, Windowsize);
 
    glutInit(&argc, argv);
 
    glutCreateWindow("Triple Pendulums");
 
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
 
    glutDisplayFunc(display);
    glutTimerFunc(16, timertick, 0);
	glutKeyboardFunc(keydown);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    glutMainLoop();

	return 0;
}
