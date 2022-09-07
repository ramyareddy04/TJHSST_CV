// Ramya Reddy
// Period 4
// 10/18/2021

#include <bits/stdc++.h>
#include <chrono>
#include<iostream>
#include <iomanip> 
#include <fstream>
#include<math.h>
#include <list>
#include <string>
#include <vector>
using namespace std;
using namespace std::chrono;

vector<duration<float>> times;

class Point
{
private:
	double x;
	double y;
public:
	Point() : x(0.0), y(0.0) {};
	Point(double num1, double num2)
	{
		x = num1;
		y = num2;
	}
	// gets X-coordinate
	double getX()
	{
		return x;
	}
	// gets Y-coordinate
	double getY()
	{
		return y;
	}
	// sets X-coordinate value
	void setX(double newX)
	{
		x = newX;
	}
	// sets Y-coordinate value
	void setY(double newY)
	{
		y = newY;
	}
	// determines length, or distance, between this and given point
	double find_length(Point other)
	{
		return sqrt(pow(getX() - other.getX(), 2) + pow(getY() - other.getY(), 2));
	}
	// determines length, or distance, between this and given point
	double find_length_sqrd(Point other)
	{
		return pow(getX() - other.getX(), 2) + pow(getY() - other.getY(), 2);
	}
	// prints Point as string to console
	string toString()
	{
		return "(" + to_string(x) + "," + to_string(y) + ")";
	}
};

class TwoPoints
{
private:
	double distance;
	Point one;
	Point two;
	int mid_pt;
public:
	TwoPoints();
	TwoPoints(Point p1, Point p2)
	{
		distance = p1.find_length_sqrd(p2);
		one = p1;
		two = p2;
	}
	double get_distance()
	{
		return get_first_point().find_length(get_second_point());
	}
	double get_sqrd_distance()
	{
		return distance;
	}
	Point get_first_point()
	{
		return one;
	}
	Point get_second_point()
	{
		return two;
	}
	string toString()
	{
		return to_string(get_distance()) + " " + get_first_point().toString() + " " + get_second_point().toString();
	}
};

class Line
{
private:
	Point one;
	Point two;
	vector<double> standardForm;
	vector<double> interceptForm;
public:
	Line();
	Line(Point p1, Point p2)
	{
		standardForm = find_equation_of_line({ p1, p2 });
		if (p1.getX() == p2.getX()) {
			one = Point(p1.getX(), 0.0);
			two = Point(p1.getX(), 1.0);
			interceptForm = {};
		}
		else {
			one = Point(0.0, standardForm[2] / standardForm[1]);
			two = Point(1.0, (standardForm[2] - standardForm[0]) / standardForm[1]);
			interceptForm = { two.getY() - one.getY(), one.getY() };
		}
	}
	// determines coefficients of equation of line in standard form
	vector<double> find_equation_of_line(vector<Point> points)
	{
		vector<double> coefficients;
		coefficients.push_back(points[1].getY() - points[0].getY());
		coefficients.push_back(points[0].getX() - points[1].getX());
		coefficients.push_back((coefficients[0] * points[0].getX() + coefficients[1] * points[0].getY()));
		return coefficients;
	}
	// gets coefficents for equation of line in standard form
	vector<double> get_standard_form()
	{
		return standardForm;
	}
	// gets coefficients for equation of line slope-intercept form
	vector<double> get_intercept_form()
	{
		return interceptForm;
	}
	// find intersection of this and other line
	Point find_intersection(Line other)
	{
		vector<double> otherIntercept = other.get_intercept_form();
		double X = (otherIntercept[1] - interceptForm[1]) / (interceptForm[0] - otherIntercept[0]);
		double Y = interceptForm[0] * X + interceptForm[1];
		return Point(X, Y);
	}
	// gets y-intercept of line
	Point get_first_point()
	{
		return one;
	}
	// gets intersection of this and x=800
	Point get_second_point()
	{
		return two;
	}
};

class Graph
{
private:
	char(*arr)[800];
	list<Point> coordinates;

public:
	const int size = 800;
	// helper method that instantiates graphs with white pixels
	void fill_graph()
	{
		arr = new char[800][800];
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				arr[i][j] = 'W';
			}
		}
	}
	// custom rounding method
	int round(double val)
	{
		return int(val) + int(abs(int(val) - val) >= 0.5);
	}
	// helper method that scales and rounds given points
	list<int> scale(list<Point> points)
	{
		list<int> scaled_points;
		list<Point>::iterator it;
		for (it = points.begin(); it != points.end(); ++it) {
			scaled_points.push_back(round(it->getX() * size));
			scaled_points.push_back(round(it->getY() * size));
		}
		return scaled_points;
	}
	list<Point> convert(vector<Point> points)
	{
		list<Point> scaled_points;
		vector<Point>::iterator it;
		for (it = points.begin(); it != points.end(); ++it) {
			scaled_points.push_back(*it);
		}
		return scaled_points;
	}
	// helper method that changes color of pixel
	void illuminate(int x, int y, char color)
	{
		if ((x >= 0 and x <= size) and (y >= 0 and y <= size)) {
			arr[x][y] = color;
		}
	}
	// helper method that creates line with negative slope and driving x axis
	void negative_x(int x1, int y1, int x2, int y2, int deltaX, int deltaY, int error, char color)
	{
		int j = y1;
		for (int i = x1; i >= x2; i--) {
			illuminate(i, j, color);
			if (error >= 0) {
				if (deltaY < 0) { // if true move down else move up
					j -= 1;
				}
				else if (deltaY > 0) {
					j += 1;
				} error -= abs(deltaX);
			} error += abs(deltaY);
		}
	}
	// helper method that creates line with positive slope and driving x axis
	void positive_x(int x1, int y1, int x2, int y2, int deltaX, int deltaY, int error, char color)
	{
		int j = y1;
		for (int i = x1; i <= x2; i++) {
			illuminate(i, j, color);
			if (error >= 0) {
				if (deltaY < 0) { // if true move down else move up
					j -= 1;
				}
				else if (deltaY > 0) {
					j += 1;
				} error -= abs(deltaX);
			} error += abs(deltaY);
		}
	}
	// helper method that creates line with positive slope and driving y axis
	void positive_y(int x1, int y1, int x2, int y2, int deltaX, int deltaY, int error, char color)
	{
		int i = x1;
		for (int j = y1; j <= y2; j++) { // if true move down else move up
			illuminate(i, j, color);
			if (error >= 0) {
				if (deltaX < 0) { // if true move left else move right
					i -= 1;
				}
				else if (deltaX > 0) {
					i += 1;
				} error -= abs(deltaY);
			} error += abs(deltaX);
		}
	}
	// helper method that creates line with negative slope and driving y axis
	void negative_y(int x1, int y1, int x2, int y2, int deltaX, int deltaY, int error, char color)
	{
		int i = x1;
		for (int j = y1; j >= y2; j--) { // if true move down else move up
			illuminate(i, j, color);
			if (error >= 0) {
				if (deltaX < 0) { // if true move left else move right
					i -= 1;
				}
				else if (deltaX > 0) {
					i += 1;
				} error -= abs(deltaY);
			} error += abs(deltaX);
		}
	}
	// modified Bresenham's algorithm that draws line for all twelve cases
	void create_line(int x1, int y1, int x2, int y2, char color)
	{
		int deltaX = x2 - x1;
		int deltaY = y2 - y1;
		int error = deltaY - deltaX;
		if (abs(deltaX) > abs(deltaY)) { // if true driving axis is x-axis else driving axis is y-axis
			if (deltaX < 0) { // if true move right else move left
				negative_x(x1, y1, x2, y2, deltaX, deltaY, error, color);
			}
			else {
				positive_x(x1, y1, x2, y2, deltaX, deltaY, error, color);
			}
		}
		else {
			if (deltaY < 0) {
				negative_y(x1, y1, x2, y2, deltaX, deltaY, error, color);
			}
			else {
				positive_y(x1, y1, x2, y2, deltaX, deltaY, error, color);
			}
		}
	}
	// rasterized algorithm that draws circle
	void create_circle(int r, int x0, int y0, char color) {
		int x, y, xmax, y2, y2_new, ty;
		xmax = (int)(r / sqrt(2));
		y = r;
		y2 = y * y;
		ty = (2 * y) - 1;
		y2_new = y2;
		for (x = 0; x < xmax + 2; x++) {
			if ((y2 - y2_new) >= ty) {
				y2 -= ty;
				y -= 1;
				ty -= 2;
			}
			illuminate(x + x0, y + y0, color);
			illuminate(x + x0, -y + y0, color);
			illuminate(-x + x0, y + y0, color);
			illuminate(-x + x0, -y + y0, color);
			illuminate(y + x0, x + y0, color);
			illuminate(-y + x0, x + y0, color);
			illuminate(y + x0, -x + y0, color);
			illuminate(-y + x0, -x + y0, color);
			y2_new -= (2 * x) - 3;
		}
	}
	// draws original 4 points onto to graph
	void draw_figures()
	{
		list<Point>::iterator it;
		for (it = coordinates.begin(); it != prev(coordinates.end(), 2); ++it) {
			create_circle(3, round((it->getX()) * 800), round((it->getY()) * 800), 'B');
		}
		for (it = prev(coordinates.end(), 2); it != coordinates.end(); ++it) {
			create_circle(3, round((it->getX()) * 800), round((it->getY()) * 800), 'R');
			create_circle(2, round((it->getX()) * 800), round((it->getY()) * 800), 'R');
		}
	}
	// creates ppm file with origin on the bottom left
	void create_file()
	{
		ofstream myfile;
		myfile.open("output.ppm", ios::out);
		myfile << "P3 " << size << " " << size << " 1\n";
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (arr[j][i] == 'W') {
					myfile << 1 << " " << 1 << " " << 1 << " ";
				}
				else if (arr[j][i] == 'B') {
					myfile << 0 << " " << 0 << " " << 0 << " ";
				}
				else if (arr[j][i] == 'R') {
					myfile << 1 << " " << 0 << " " << 0 << " ";
				}
			}
			myfile << endl;
		}
		myfile.close();
	}
	// gets coordinates vector
	list<Point> get_coords()
	{
		return coordinates;
	}
	// sets coordinates vector
	void set_coords(list<Point> p)
	{
		coordinates = p;
	}
	void set_coords(vector<Point> p)
	{
		coordinates = convert(p);
	}
	// destructor for Graph class
	~Graph()
	{
		if (arr != NULL) {
			delete[] arr;
		}
	}
};

// method that writes points to .txt file
void create_pt_file(list<Point> points)
{
	ofstream myfile;
	myfile.open("points.txt", ios::out);
	list<Point>::iterator it;
	for (it = points.begin(); it != points.end(); ++it) {
		myfile << setprecision(23) << it->getX() << "  " << it->getY() << fixed << endl;
	}
	myfile.close();
}

// brute force algorithim to get points of min distance from each other
vector<Point> brute_force(list<Point> points)
{
	double curr_min = 640000;
	Point pt1;
	Point pt2;
	list<Point>::iterator it;
	list<Point>::iterator it2;
	double temp;
	for (it = points.begin(); it != prev(points.end()); ++it) {
		for (it2 = next(it); it2 != points.end(); ++it2) {
			temp = (*it).find_length_sqrd(*it2);
			if (temp < curr_min) {
				curr_min = temp;
				pt1 = *it;
				pt2 = *it2;
			}
		}
	}
	return { pt1, pt2 };
}

// method that generates random points and does brute force approach
TwoPoints part1()
{
	srand(time(0));
	list<Point> points;
	for (int i = 0; i < 1000; i++) {
		points.push_back(Point((double)rand() / (RAND_MAX), (double)rand() / (RAND_MAX)));
	}
	create_pt_file(points);
	auto start = high_resolution_clock::now();
	vector<Point> finPts = brute_force(points);
	auto stop = high_resolution_clock::now();
	times.push_back(duration_cast<microseconds>(stop - start));
	points.push_back(finPts[0]);
	points.push_back(finPts[1]);
	Graph g;
	g.fill_graph();
	g.set_coords(points);
	g.draw_figures();
	g.create_file();
	return TwoPoints(finPts[0], finPts[1]);
}

// reads points.txt file
vector<Point> read_pt_file()
{
	string line;
	vector<Point> coords;
	ifstream myfile("points.txt");
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			int i = 0;
			while (i<int(line.size())) {
				if (line[i] == ' ') {
					coords.push_back(Point(stod(line.substr(0, i - 1)), stod(line.substr(i + 2, line.length()))));
					i = (int)line.size();
				}
				i += 1;
			}
		}
	}
	myfile.close();
	return coords;
}

// method that sorts vector of Points
bool compare_X_coord(Point p1, Point p2)
{
	return (p1.getX() < p2.getX());
}

// method that returns minimum of two given points
TwoPoints find_min(TwoPoints num1, TwoPoints num2)
{
	if (num1.get_sqrd_distance() < num2.get_sqrd_distance()) {
		return num1;
	}
	else {
		return num2;
	}
}

// implementation of brute force approach for strip in part 2
TwoPoints brute_force2(Point mid, int midPt, vector<Point> v, int pts, int pts2)
{
	double curr_min = 640000;
	Point pt1;
	Point pt2;
	double temp;
	vector<Point>::iterator it;
	vector<Point>::iterator it2;
	if (midPt != -1) {
		for (it = v.begin() + pts; it != v.begin() + midPt; ++it) {
			for (it2 = v.begin() + midPt + 1; it2 != v.begin() + pts2; ++it2) {
				temp = (*it).find_length_sqrd(*it2);
				if (temp < curr_min) {
					curr_min = temp;
					pt1 = *it;
					pt2 = *it2;
				}
			}
		}
	}
	else {
		for (it = v.begin() + pts; it != prev(v.begin() + pts2); ++it) {
			for (it2 = next(it); it2 != v.begin() + pts2; ++it2) {
				temp = (*it).find_length_sqrd(*it2);
				if (temp < curr_min) {
					curr_min = temp;
					pt1 = *it;
					pt2 = *it2;
				}
			}
		}
	}
	return TwoPoints(pt1, pt2);
}

// creates vectors to use for brute force in strip
TwoPoints val_check(TwoPoints smallest, vector<Point> pts, int i, int k)
{
	int mid = (i + k) / 2;
	double midX = pts[mid].getX();
	double posDif = (midX + smallest.get_distance());
	double negDif = (midX - smallest.get_distance());

	int idxL = mid - 1;
	while (pts[idxL].getX() > negDif and idxL >= i) {
		idxL -= 1;
	}
	int idxR = mid + 1;
	while (pts[idxR].getX() < posDif and idxR <= k) {
		idxR += 1;
	}

	if (idxL != mid - 1 and idxR != mid + 1) {
		smallest = find_min(smallest, brute_force2(pts[mid], mid, pts, idxL, idxR));
	}
	else {
		smallest = find_min(smallest, brute_force2(pts[mid], -1, pts, idxL, idxR));
	}
	return smallest;
}

// recursive merge sort method
TwoPoints recur(vector<Point> v, int i, int k)
{
	int s = k - i + 1;
	if (s == 2) {
		return TwoPoints(v[i], v[k]);
	}
	else if (s == 3) {
		return find_min(find_min(TwoPoints(v[i], v[k]), TwoPoints(v[i], v[i + 1])), TwoPoints(v[i + 1], v[k]));
	}
	else {
		double j = (i + k) / 2;
		TwoPoints l = recur(v, i, j);
		TwoPoints r = recur(v, j, k);
		return val_check(find_min(l, r), v, i, k);
	}
}

// method that reads points from file and implements intermediate recursive approach
TwoPoints part2()
{
	vector<Point> pts = read_pt_file();
	auto start2 = high_resolution_clock::now();
	sort(pts.begin(), pts.end(), compare_X_coord);
	TwoPoints smallest = recur(pts, 0, pts.size());
	auto stop2 = high_resolution_clock::now();
	times.push_back(duration_cast<microseconds>(stop2 - start2));
	return smallest;
}

// writes results to terminal and output file
void create_results_file(vector<TwoPoints> points)
{
	ofstream myfile;
	myfile.open("results.txt", ios::out);
	for (int i = 0; i < (int)points.size(); i++) {
		myfile << "Results for Part " << (i + 1) << ":" << endl;
		cout << "Results for Part " << (i + 1) << ":" << endl;
		myfile << "Points: " << points[i].get_first_point().toString() << " " << points[i].get_second_point().toString() << " Distance: " << points[i].get_distance() << " Time (seconds): " << times[i].count() << endl;
		cout << "Points: " << points[i].get_first_point().toString() << " " << points[i].get_second_point().toString() << " Distance: " << points[i].get_distance() << " Time (seconds): " << times[i].count() << endl;
	}
	myfile.close();
}

int main()
{
	TwoPoints res1 = part1();
	TwoPoints res2 = part2();
	create_results_file({ res1, res2 });
}