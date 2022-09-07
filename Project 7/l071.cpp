#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <iostream>
#include <math.h>

using namespace std;
using namespace cv;

Mat create_grayscale( string src ) 
{
    Mat gray;
    gray = imread( src, IMREAD_GRAYSCALE );
    imwrite("./imageg.jpg", gray);    
    return gray;
}

int l = 0, h = 250, r = 3, ks = 3, ms = 100, minR = 70, maxR = 150, tc = 100, minQ = 103;
string filename = "image.jpg";
Mat canny_edge (Mat src)
{
    Mat canny; 
    blur( src, canny, Size(3,3) );
    Canny( canny, canny, l, h*r, ks );
    imwrite("./imagef.jpg", canny);
    return canny;
}

vector<Vec3f> hough_transform ( Mat gray )
{
    medianBlur( gray, gray, 5);
    vector<Vec3f> circles;
    if (maxR > 120) {
        HoughCircles(gray, circles, HOUGH_GRADIENT, 1, ms, h, tc , minR, 120);
        if (maxR >= 160) {
            vector<Vec3f> largeCircles;
            HoughCircles(gray, largeCircles, HOUGH_GRADIENT, 1, ms, h, max(tc-30, 0) , 160, maxR);
            circles.insert( circles.end(), largeCircles.begin(), largeCircles.end() );
        }
    } else {
        HoughCircles(gray, circles, HOUGH_GRADIENT, 1, ms, h, tc , minR, maxR);
    }
    
    return circles;
}

vector<int> averageRGB ( Mat src, int radius , Point center )
{
    vector<int> rgb = { 0 , 0 , 0};
    for (int i=0; i<3; i++) 
        rgb[i] += src.at<cv::Vec3b>( center.y , center.x )[2-i];
    bool right = center.x + radius < src.rows, left = center.x - radius >= 0, top = center.y + radius < src.cols, bottom = center.y - radius >= 0;
    if ( top and right) {
        for (int i=0; i<3; i++) 
            rgb[i] += src.at<cv::Vec3b>( center.y + radius , center.x + radius )[2-i];
    }
    if ( bottom and left ) {
        for (int i=0; i<3; i++) 
            rgb[i] += src.at<cv::Vec3b>( center.y - radius , center.x - radius )[2-i];
    }
    if ( top and left) {
        for (int i=0; i<3; i++) 
            rgb[i] += src.at<cv::Vec3b>( center.y + radius , center.x - radius )[2-i];
    }
    if ( bottom and right ) {
        for (int i=0; i<3; i++) 
            rgb[i] += src.at<cv::Vec3b>( center.y - radius , center.x + radius )[2-i];
    }
    return rgb;
}

void create_results(string filename, int count[5])
{
    ofstream myfile;
    myfile.open(filename, ios::out);
    double sum = 0.0;
    string type[5] = {"Pennies", "Nickels", "Dimes", "Quarters", "Silver Dollars" }; 
    double val[5] = {0.01, 0.05, 0.10, 0.25, 0.50 }; 
    for (int i=0; i<5; i++) {
        myfile << count[i] << " " << type[i] << endl;
        cout << count[i] << " " << type[i] << endl;
        sum += count[i]*val[i];
    }
    myfile << "Total: $" << sum << endl;
    cout << "Total: $" << sum << endl;
    myfile.close();
}

void classify ( Mat src, vector<Vec3f> identified )
{
    int type;
    int count[5] = { 0 , 0 , 0 , 0 };
    Scalar colors[5] = { Scalar(0,0,255) , Scalar(128,0,128) , Scalar(255,0,0) , Scalar(0,255,0) , Scalar(0,255,255) };
    for (int i=0; i<identified.size(); i++)
    {
        int radius = cvRound(identified[i][2]);
        Point center(cvRound(identified[i][0]), cvRound(identified[i][1]));
        vector<int> avg = averageRGB( src , radius/4 , center );
        int avgR = avg[0], avgB = avg[2];
        if ( radius < 89) {
            if ( avgR > 1.1*avgB ) {
                type = 0;
            } else {
                type = 2;
            }
        } else if ( radius < minQ ) {
            if ( avgR > 1.3*avgB ) {
                type = 0;
            } else {
                type = 1;
            }
        } else if ( radius < 121 ) {
            type = 3;
        } else if (radius > 160) {
            type = 4;
        }
        count[type] += 1;
        circle( src , center , 3 , colors[type] , -1 , 8 , 0 );
        for (int j = 0; j < 4; j++ ) 
            circle( src , center , radius+j , colors[type] , 3 , 8 , 0 );
    }
    create_results( "./results.txt" , count);
    imwrite( "./coins.jpg", src );
}

void part1 () 
{
    Mat img;
    img = imread( filename, 1 );

    if ( !img.data )
    {
        printf("No image data \n");
        return;
    }
    
    Mat gray = create_grayscale( filename );
    Mat canny = canny_edge( gray );
    vector<Vec3f> coins = hough_transform( gray );
    classify( img, coins );
}

int main(int argc, char** argv )
{
    if (argc > 1) {
        for (int i=1; i<argc-1;) {
            if (strcmp(argv[i],"-H") == 0) {
                h = stoi(argv[i+1]);
            } else if (strcmp(argv[i],"-L") == 0) {
                l = stoi(argv[i+1]);
            } else if (strcmp(argv[i],"-F") == 0) {
                filename = argv[i+1];
            } else if (strcmp(argv[i],"-KS") == 0) {
                ks = stod(argv[i+1]);
            } else if (strcmp(argv[i],"-MS") == 0) {
                ms = stod(argv[i+1]);
            } else if (strcmp(argv[i],"-R") == 0) {
                r = stoi(argv[i+1]);
            } else if (strcmp(argv[i],"-MinR") == 0) {
                minR = stoi(argv[i+1]);
            } else if (strcmp(argv[i],"-MaxR") == 0) {
                maxR = stoi(argv[i+1]);
            } else if (strcmp(argv[i],"-TC") == 0) {
                tc = stoi(argv[i+1]);
            } else if (strcmp(argv[i],"-MinQ") == 0) {
                minQ = stoi(argv[i+1]);
            } 
            i+=2;
        }
    }
    
    part1();
    waitKey(0);

    return 0;
}