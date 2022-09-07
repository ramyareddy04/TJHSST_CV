// // Ramya Reddy
// // Period 4
// // 05/29/2022

#include <stdio.h>
#include <opencv2/opencv.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/core/matx.hpp"
#include <iostream>
#include <string>
#include <cmath>
#include <math.h>
#include <opencv2/core/types_c.h>

using namespace std;
using namespace cv;

void part1()
{
    vector<Point3f> objp = {Point3f(1, 1, 0), Point3f(-1, 1, 0), Point3f(1, -1, 0), Point3f(-1, -1, 0)};
    vector<Point3f> cube = {Point3f(1, 1, -2), Point3f(1, -1, -2), Point3f(-1, 1, -2), Point3f(-1, -1, -2), Point3f(1, 1, 0), Point3f(1, -1, 0), Point3f(-1, 1, 0), Point3f(-1, -1, 0)};
    vector<Point2f> projectionPoints, out, corners, cornersOF;
    vector<vector<Point3f>> obj, objSet;
    vector<vector<Point2f>> img, imgSet;
    
    Mat prev;
    VideoCapture cap("withChessBoard.MOV");
    VideoWriter video("vr.avi", cv::VideoWriter::fourcc('M','J','P','G'), 25, Size(int(cap.get(3)), int(cap.get(4)))); 
    int idx = 0;
    while(cap.isOpened()) {
        Mat frame;
        bool isSuccess = cap.read(frame);
        if (isSuccess == true) {
            vector<uchar> status;
            vector<float> err;
            flip(frame, frame, 0);
            flip(frame, frame, 1);
            cvtColor(frame, frame, COLOR_BGR2GRAY);
            out.clear();
            bool found = findChessboardCorners(frame, Size(7, 7), out, CALIB_CB_ADAPTIVE_THRESH + CALIB_CB_NORMALIZE_IMAGE + CALIB_CB_FAST_CHECK);
            if (found == true) {
                cornerSubPix(frame, out, Size(11, 11), Size(-1, -1), TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1));
                corners.clear();
                corners.push_back(out[16]);
                corners.push_back(out[18]);
                corners.push_back(out[30]);
                corners.push_back(out[32]);                
            } else {
                calcOpticalFlowPyrLK(prev, frame, corners, cornersOF, status, err, Size(15, 15), 2, TermCriteria((TermCriteria::COUNT) + (TermCriteria::EPS), 10, 0.03));
                corners = cornersOF;
            }
            obj.push_back(objp);
            img.push_back(corners);
            if (idx%10 == 0) {
                objSet.push_back(objp);
                imgSet.push_back(corners);
            }
            idx += 1;
            frame.copyTo(prev);
        }
        else
            break;
    }
    cap.release();
    
    Mat camera, distC, R, T;
    calibrateCamera(objSet, imgSet, prev.size(), camera, distC, R, T);
    vector<Point2i> edges;
    for (int i=0; i<cube.size()-1; i++) {
        for (int j=i+1; j<cube.size(); j++) {
            if (sqrt(pow(cube[i].x - cube[j].x, 2) + pow(cube[i].y - cube[j].y, 2) + pow(cube[i].z - cube[j].z, 2)) == 2) {
                edges.push_back(Point2i(i, j));
            }
        }
    }
    VideoCapture cap2("withChessBoard.MOV");
    idx = 0;
    while(cap2.isOpened()) {
        Mat frame;
        bool isSuccess = cap2.read(frame);
        if (isSuccess == true) {
            flip(frame, frame, 0);
            flip(frame, frame, 1);
            solvePnP(obj[idx], img[idx], camera, distC, R, T);
            projectPoints(cube, R, T, camera, distC, projectionPoints);
            for (Point2i e: edges) {
                line(frame, projectionPoints[e.x], projectionPoints[e.y], Scalar(0, 0, 255), 2);
            }
            idx+=1;
            video.write(frame);
        }
        else
            break;
    }
    cap2.release();
    video.release();
    
}

int main(int argc, char** argv )
{
    part1();
    waitKey(0);
    return 0;
}