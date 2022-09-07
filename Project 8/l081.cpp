#include <stdio.h>
#include <opencv2/opencv.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <string>
#include <cmath>
#include <math.h>

using namespace std;
using namespace cv;

float phi = 0.5*(1+sqrt(5));
ofstream myfile;

class HomogeneousPoint
{
private:
    double x, y, z;
    int w;
public:
    HomogeneousPoint() : x(0.0), y(0.0), z(0.0), w(0) {};
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
                            {0, sin(aX), cos(aZ), 0},
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

void write_to_video(VideoWriter video, vector<HomogeneousPoint> points, int scaleFactor, string type)
{
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
    while (i != 105)
    {
        Mat frame(800, 600, CV_8UC3, cv::Scalar(1, 1, 1)); // 600, 800
        for (int j=0; j<points.size(); j++) {
            HomogeneousPoint temp = translation(rotation(points[j], i*10, i*10, 0), 300, 400, 0);
            for (int k=j+1; k < points.size(); k++) {
                if (abs(points[j].sqrd_distance(points[k])-side) < 1 ) {
                    line(frame, temp.to_orthographic(), translation(rotation(points[k], i*10, i*10, 0), 300, 400, 0).to_orthographic(), Scalar(255, 255, 255), 1, 8);
                }
            }
            if (i < 4 and type == "square") {
                if (j != points.size()-1) 
                     myfile << temp.toString() << ", ";
                else
                     myfile << temp.toString() << endl;
            }
        }
        video.write(frame);
        i += 1;
    }
}

void part1 () 
{
    myfile.open("coordinates.txt", ios::out);
    VideoWriter video("rotation.avi", cv::VideoWriter::fourcc('M','J','P','G'), 15, Size(600, 800)); // should be (800,600)
    vector<HomogeneousPoint> points = {HomogeneousPoint(-1,-1,-1,1), HomogeneousPoint( 1, 1, 1,1),
                                       HomogeneousPoint( 1,-1,-1,1), HomogeneousPoint( 1, 1,-1,1),
                                       HomogeneousPoint(-1, 1,-1,1), HomogeneousPoint( 1,-1, 1,1),
                                       HomogeneousPoint(-1,-1, 1,1), HomogeneousPoint(-1, 1, 1,1)};
    write_to_video(video, points, 30, "square");
    points = {HomogeneousPoint(pow(3,-0.5),0,0,1), HomogeneousPoint(0,0,2*pow(6,-0.5),1),
              HomogeneousPoint(-sqrt(3)/6,0.5,0,1), HomogeneousPoint(-sqrt(3)/6,-0.5,0,1)};
    write_to_video(video, points, 80, "tetrahedron");
    points = {HomogeneousPoint(1,0,0,1), HomogeneousPoint(-1,0,0,1),
              HomogeneousPoint(0,1,0,1), HomogeneousPoint(0,-1,0,1),
              HomogeneousPoint(0,0,1,1), HomogeneousPoint(0,0,-1,1)};
    write_to_video(video, points, 45, "octahedron");
    points = {HomogeneousPoint(0,phi,1,1), HomogeneousPoint(0,phi,-1,1),
              HomogeneousPoint(0,-phi,1,1), HomogeneousPoint(0,-phi,-1,1),
              HomogeneousPoint(phi,1,0,1), HomogeneousPoint(phi,-1,0,1),
              HomogeneousPoint(-phi,1,0,1), HomogeneousPoint(-phi,-1,0,1),
              HomogeneousPoint(1,0,phi,1), HomogeneousPoint(-1,0,phi,1),
              HomogeneousPoint(1,0,-phi,1), HomogeneousPoint(-1,0,-phi,1)};
    write_to_video(video, points, 30, "icosahedron");
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
    write_to_video(video, points, 30, "dodecahedron");
    video.release();
    myfile.close();
}

int main(int argc, char** argv )
{
    part1();
    waitKey(0);
    return 0;
}