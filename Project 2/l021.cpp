// Ramya Reddy
// Period 4
// 9/9/2021

#include<iostream>
#include <iomanip> 
#include <fstream>
#include<math.h>
#include <vector>
using namespace std;

// method that calculates the area of a triangle 
double calculateArea(vector<pair<double, double>> points)
{
    double a = sqrt((pow(points[0].first - points[1].first, 2) + pow(points[0].second - points[1].second, 2)));
    double b = sqrt((pow(points[1].first - points[2].first, 2) + pow(points[1].second - points[2].second, 2)));
    double c = sqrt((pow(points[0].first - points[2].first, 2) + pow(points[0].second - points[2].second, 2)));
    double s = (a + b + c) / 2;
    return sqrt(s*(s-a)*(s-b)*(s-c));
}
// method that determines whether point is within triangle or not using the area method discussed in class
bool isInTriangle(vector<pair<double, double>> points, pair<double, double> target)
{
    double tot = calculateArea(points); 
    double a = calculateArea({points[0], points[1], target});
    double b = calculateArea({points[0], points[2], target});
    double c = calculateArea({points[2], points[1], target});
    double eps = pow(2, -52);
    return (a+b+c>=tot-eps) and (a+b+c<=tot+eps);
}
// method that goes through all possible triangles within 4 points to ensure a point is not within another triangle
bool isInIntersection(vector<pair<double, double>> points)
{
    return (isInTriangle({points[0], points[1], points[2]}, points[3]) 
            or isInTriangle({points[3], points[1], points[2]}, points[0]) 
            or isInTriangle({points[0], points[3], points[2]}, points[1]) 
            or isInTriangle({points[0], points[1], points[3]}, points[2]));
}
// method that writes points to .txt file
void create_file(vector<pair<double,double>> points)
{
    ofstream myfile;
    myfile.open("points.txt", ios::out);
    for (int i=0; i<points.size() - 1; i++) {
        myfile << setprecision(17) << "(" << points[i].first << ", " << points[i].second << ") ," << fixed << " ";
    }
    myfile << setprecision(17) << "(" << points.back().first << ", " << points.back().second << ")" << fixed;
    myfile.close();
}
// method that generates four random points
vector<pair<double,double>> part1()
{
    srand(time(0));
    vector<pair<double,double>> points;
    for (int i = 0; i < 4; i++) {
        points.push_back(make_pair(((double)rand() / (RAND_MAX)), ((double)rand() / (RAND_MAX))));
    }
    bool inIntersection = isInIntersection(points);
    while(inIntersection) {
        points[3].first = ((double)rand() / (RAND_MAX));
        points[3].second = ((double)rand() / (RAND_MAX));
        inIntersection = isInIntersection(points);
    }
    create_file(points);
    return points;
}

int main()
{   
    part1();
}