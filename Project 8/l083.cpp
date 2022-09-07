#include <stdio.h>
#include <opencv2/opencv.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/core/matx.hpp"
#include <iostream>
#include <string>
#include <cmath>
#include <math.h>

using namespace std;
using namespace cv;

float phi = 0.5*(1+sqrt(5));
ofstream myfile;
bool firstLine = false;

class HomogeneousPoint
{
private:
    double x, y, z;
    int w;
public:
    HomogeneousPoint() : x(0.0), y(0.0), z(0.0), w(0) {};
    HomogeneousPoint(Vec3f v)
    {
        x = v[0];
        y = v[1];
        z = v[2];
    }
    HomogeneousPoint(double num1, double num2, double num3, int weight)
	{
		x = num1;
		y = num2;
        z = num3;
        w = weight;
	}
    double getX()
	{
		return x;
	}
	double getY()
	{
		return y;
	}
    double getZ()
	{
		return z;
	}
    double getW()
	{
		return w;
	}
	void setX(double newX)
	{
		x = newX;
	}
	void setY(double newY)
	{
		y = newY;
	}
    void setZ(double newZ)
	{
		z = newZ;
	}
    HomogeneousPoint add(HomogeneousPoint& other)
    {
        return HomogeneousPoint(getX()+other.getX(), getY()+other.getY(), getZ()+other.getZ(), 1);
    }
    HomogeneousPoint subtract(HomogeneousPoint& other)
    {
        return HomogeneousPoint(getX()-other.getX(), getY()-other.getY(), getZ()-other.getZ(), 1);
    }
    double dot (HomogeneousPoint& other)
    {
        return getX()*other.getX() + getY()*other.getY() + getZ()*other.getZ();
    }
    bool operator == (const HomogeneousPoint& p) const
    {
        return (x/w==p.x/p.w and y/w==p.y/p.w and z/w==p.z/p.w);
    }
    double sqrd_distance(HomogeneousPoint other)
    {
        return pow(getX()-other.getX(), 2) + pow(getY()-other.getY(), 2) + pow(getZ()-other.getZ(), 2);
    }
    Point to_orthographic()
    {
        return Point(x,y);
    }
    Vec3f to_vector()
    {
        return Vec3f(x, y, z);
    }
	string toString()
	{
		return "(" + to_string(int(x)) + "," + to_string(int(y)) + "," + to_string(int(z)) + ")";
	}
};

vector<double> create(string method, double theta, int d)
{
    vector<double> ultimate = {};
    return ultimate;
}

HomogeneousPoint translation(HomogeneousPoint initial, double tX, double tY, double tZ)
{
    double p0[1][4] = {{initial.getX(), initial.getY(), initial.getZ(), 1}};
    double matrix[4][4] = {{1, 0, 0, tX},
                           {0, 1, 0, tY},
                           {0, 0, 1, tZ},
                           {0, 0, 0, 1}};
    Mat result = Mat(4, 4, CV_64FC1, matrix) * Mat(1, 4, CV_64FC1, p0).t();
    result = result.t();
    HomogeneousPoint pf(result.at<double>(0,0), result.at<double>(0,1), result.at<double>(0,2), 1);
    return pf;
}

HomogeneousPoint scale(HomogeneousPoint initial, double sX, double sY, double sZ)
{
    double p0[1][4] = {{initial.getX(), initial.getY(), initial.getZ(), 1}};
    double matrix[4][4] = {{sX, 0, 0, 0},
                           {0, sY, 0, 0},
                           {0, 0, sZ, 0},
                           {0, 0, 0, 1}};
    Mat result = Mat(4, 4, CV_64FC1, matrix) * Mat(1, 4, CV_64FC1, p0).t();
    result = result.t();
    HomogeneousPoint pf(result.at<double>(0,0), result.at<double>(0,1), result.at<double>(0,2), 1);
    return pf;
}

HomogeneousPoint rotation(HomogeneousPoint initial, double X, double Y, double Z)
{
    double aX = X * M_PI/180;
    double aY = Y * M_PI/180;
    double aZ = Z * M_PI/180;
    double p0[1][4] = {{initial.getX(), initial.getY(), initial.getZ(), 1}};
    double matrixX[4][4] = {{1, 0, 0, 0},
                            {0, cos(aX), -sin(aX), 0},
                            {0, sin(aX), cos(aX), 0},
                            {0, 0, 0, 1}};
    double matrixY[4][4] = {{cos(aY), 0, sin(aY), 0},
                            {0, 1, 0, 0},
                            {-sin(aY), 0, cos(aY), 0},
                            {0, 0, 0, 1}};
    double matrixZ[4][4] = {{cos(aZ), -sin(aZ), 0, 0},
                            {sin(aZ), cos(aZ), 0, 0},
                            {0, 0, 1, 0},
                            {0, 0, 0, 1}};
    Mat R = Mat(4, 4, CV_64FC1, matrixX) * Mat(4, 4, CV_64FC1, matrixY) * Mat(4, 4, CV_64FC1, matrixZ);
    Mat result = R * Mat(1, 4, CV_64FC1, p0).t();
    result = result.t();
    HomogeneousPoint pf(result.at<double>(0,0), result.at<double>(0,1), result.at<double>(0,2), 1);
    return pf;
}

double operation(Vec3f x, Vec3f a, Vec3f n)
{
    return (x-a).dot(n);
}

HomogeneousPoint perspective(HomogeneousPoint initial, int value, int eyeX, int eyeY, int eyeZ)
{
    double t = (value-eyeX)/(initial.getX() - eyeX);
    double u = t*initial.getY() + (1-t)*eyeY;
    double v = t*initial.getZ() + (1-t)*eyeZ;
    return HomogeneousPoint(value, u, v, 1);
}

HomogeneousPoint projection(Vec3f v, Vec3f a, Vec3f n, Vec3f eye)
{
    double t = operation(a, eye, n)/operation(v, eye, n);
    return HomogeneousPoint(t*(v-eye) + eye);
}

Point proj_perspective(Vec3f v, Vec3f p0, Vec3f w1, Vec3f w2, double tx, double ty)
{
    return Point(((v-p0).dot(w1)/w1.dot(w1)) + tx, ((v-p0).dot(w2)/w2.dot(w2)) + ty);
}

void write_to_video(VideoWriter video, vector<HomogeneousPoint> points, int scaleFactor, string type, Vec3f a, Vec3f n, Vec3f eye)
{
    Vec3f p0(projection(Vec3f(0,0,0), a, n, eye).to_vector());
    for (int i = 0; i<points.size(); i++) {
        points[i] = scale(points[i], scaleFactor, scaleFactor, scaleFactor);
    } 
    int side = pow(scaleFactor, 2);
    if (type == "square" or type == "icosahedron") {
        side *= 4;
    } else if (type == "octahedron") {
        side *= 2; 
    } else if (type == "dodecahedron") {
        side *= (6 - 2*sqrt(5));
    }
    int i = 0;
    while (i != 360)
    {
        vector<Vec3f> threeD;
        vector<Point> twoD;
        
        Mat frame(800, 600, CV_8UC3, cv::Scalar(1, 1, 1));   
        
//         vector<Vec3f> p;
//         for (int j=0; j<points.size(); j++) {       
//             if (p.size() < 3) {
//                 for (int k=j+1; k < points.size(); k++) {
//                     if (abs(points[j].sqrd_distance(points[k])-side) < 1 ) {
//                         if (p.size() == 0) {
//                             p.push_back(-points[k].to_vector());
//                         } else {
//                             if (abs(points[k].sqrd_distance(HomogeneousPoint(p[0]))-side) < 1) {
//                                 p.push_back(points[k].to_vector());
//                             }
//                         }
//                     }
//                 }
//             } else {
//                 break;
//             }
//         }
        
//         Vec3f pA;
//         pA = p[0] - p[1];
//         Vec3f pB;
//         pB = p[2] - p[1];
//         Vec3f pC;
//         pC = pB - pA*(pB.dot(pA)/pA.dot(pA));
//         Vec3f w1;
//         w1 = pA * pow(pow(pA[0],2)+pow(pA[1],2)+pow(pA[2],2), -0.5);
//         Vec3f w2;
//         w2 = pC * pow(pow(pC[0],2)+pow(pC[1],2)+pow(pC[2],2), -0.5);
        
        Vec3f w1(0, 1, 0);
        Vec3f w2(0, 0, -1);
        for (int j=0; j<points.size(); j++) {               
            Vec3f temp = translation(rotation(points[j], i, i, 0), 0, 300, 400).to_vector();
            Point u(proj_perspective(projection(rotation(points[j], i, i, 0).to_vector(), a, n, eye).to_vector(), p0, w1, w2, 300, 400));
            for (int k=j+1; k < points.size(); k++) {
                if (abs(points[j].sqrd_distance(points[k])-side) < 1 ) {
                    Vec3f temp2 = translation(rotation(points[k], i, i, 0), 0, 300, 400).to_vector();
                    Point v(proj_perspective(projection(rotation(points[k], i, i, 0).to_vector(), a, n, eye).to_vector(), p0, w1, w2, 300, 400));
                    line(frame, u, v, Scalar(255, 255, 255), 1, 8);
                    threeD.push_back(temp);
                    threeD.push_back(temp2);
                    twoD.push_back(u);
                    twoD.push_back(v);
                }
            }
        }
        if (i < 4 and type == "square") {
            if (i==0) {
                myfile << "The coordinates in the plane x = p0 + u*w1 + v*w2 is:" << endl;
                myfile << "\t(where p0 is the origin, preferraby the projection of the center of the cube in first frame, w1 and w2 are 2 perpendicular vertices in the plane)" << endl;
                myfile << "\tp0 = (" <<p0[0]<<", "<<p0[1]<<", "<<p0[2]<<")" << endl;
                myfile << "\tw1 = (" <<w1[0]<<", "<<w1[1]<<", "<<w1[2]<<")" << endl;
                myfile << "\tw2 = (" <<w2[0]<<", "<<w2[1]<<", "<<w2[2]<<")" << endl;
            }
            myfile << "The frame"<<i+1<<" in 3d has the following edges:"<<endl;
            for (int j=0; j<threeD.size(); j++) {
                if (j%2==0) {
                    myfile << "\t("<<threeD[j][0]<<", "<<threeD[j][1]<<", "<<threeD[j][2]<<")";
                } else {
                    myfile << ", ("<<threeD[j][0]<<", "<<threeD[j][1]<<", "<<threeD[j][2]<<")"<<endl;
                }
            }
            myfile << "The frame"<<i+1<<" in 2d has the following edges:"<<endl;
            for (int j=0; j<twoD.size(); j++) {
                if (j%2==0) {
                    myfile << "\t("<<twoD[j].x<<", "<<twoD[j].y<<")";
                } else {
                    myfile << ", ("<<twoD[j].x<<", "<<twoD[j].y<<")"<<endl;
                }
            }
        }
        video.write(frame);
        i += 1;
    }
}

void part3 () 
{
    myfile.open("log.txt", ios::out);
    Vec3f a(500, 300, 200);
    Vec3f n(1, 2, 3);
    Vec3f eye(800, 50, 123);
    myfile << "The plane defined by (x-a)*n =0 is:" << endl;
    myfile << "\t"<<"a: ("<<a[0]<<", "<<a[1]<<", "<<a[2]<<")"<<endl;
    myfile << "\t"<<"n: ("<<n[0]<<", "<<n[1]<<", "<<n[2]<<")"<<endl;
    myfile << "The eye e is:"<<endl;
    myfile << "\t"<<"e: ("<<eye[0]<<", "<<eye[1]<<", "<<eye[2]<<")"<<endl;
    VideoWriter video("rotation.avi", cv::VideoWriter::fourcc('M','J','P','G'), 35, Size(600, 800)); 
    vector<HomogeneousPoint> points = {HomogeneousPoint(-1,-1,-1,1), HomogeneousPoint( 1, 1, 1,1),
                                       HomogeneousPoint( 1,-1,-1,1), HomogeneousPoint( 1, 1,-1,1),
                                       HomogeneousPoint(-1, 1,-1,1), HomogeneousPoint( 1,-1, 1,1),
                                       HomogeneousPoint(-1,-1, 1,1), HomogeneousPoint(-1, 1, 1,1)};
    write_to_video(video, points, 50, "square", a, n, eye);
    points = {HomogeneousPoint(pow(3,-0.5),0,0,1), HomogeneousPoint(0,0,2*pow(6,-0.5),1),
              HomogeneousPoint(-sqrt(3)/6,0.5,0,1), HomogeneousPoint(-sqrt(3)/6,-0.5,0,1)};
    write_to_video(video, points, 150, "tetrahedron", a, n, eye);
    points = {HomogeneousPoint(1,0,0,1), HomogeneousPoint(-1,0,0,1),
              HomogeneousPoint(0,1,0,1), HomogeneousPoint(0,-1,0,1),
              HomogeneousPoint(0,0,1,1), HomogeneousPoint(0,0,-1,1)};
    write_to_video(video, points, 75, "octahedron", a, n, eye);
    points = {HomogeneousPoint(0,phi,1,1), HomogeneousPoint(0,phi,-1,1),
              HomogeneousPoint(0,-phi,1,1), HomogeneousPoint(0,-phi,-1,1),
              HomogeneousPoint(phi,1,0,1), HomogeneousPoint(phi,-1,0,1),
              HomogeneousPoint(-phi,1,0,1), HomogeneousPoint(-phi,-1,0,1),
              HomogeneousPoint(1,0,phi,1), HomogeneousPoint(-1,0,phi,1),
              HomogeneousPoint(1,0,-phi,1), HomogeneousPoint(-1,0,-phi,1)};
    write_to_video(video, points, 50, "icosahedron", a, n, eye);
    points = {HomogeneousPoint(0,1/phi,phi,1), HomogeneousPoint(0,1/phi,-phi,1),
              HomogeneousPoint(0,-1/phi,phi,1), HomogeneousPoint(0,-1/phi,-phi,1),
              HomogeneousPoint(1/phi,phi,0,1), HomogeneousPoint(1/phi,-phi,0,1),
              HomogeneousPoint(-1/phi,phi,0,1), HomogeneousPoint(-1/phi,-phi,0,1),
              HomogeneousPoint(phi,0,1/phi,1), HomogeneousPoint(phi,0,-1/phi,1),
              HomogeneousPoint(-phi,0,1/phi,1), HomogeneousPoint(-phi,0,-1/phi,1),
              HomogeneousPoint(-1,-1,-1,1), HomogeneousPoint( 1, 1, 1,1),
              HomogeneousPoint( 1,-1,-1,1), HomogeneousPoint( 1, 1,-1,1),
              HomogeneousPoint(-1, 1,-1,1), HomogeneousPoint( 1,-1, 1,1),
              HomogeneousPoint(-1,-1, 1,1), HomogeneousPoint(-1, 1, 1,1)};
    write_to_video(video, points, 50, "dodecahedron", a, n, eye);
    video.release();
    myfile.close();
}

int main(int argc, char** argv )
{
    part3();
    waitKey(0);
    return 0;
}