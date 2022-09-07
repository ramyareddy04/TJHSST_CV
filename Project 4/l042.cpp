// Ramya Reddy
// Period 4
// 12/6/2021

# include <algorithm>
#include <chrono>
#include<iostream>
#include <iomanip> 
#include <fstream>
#include<math.h>
#include <list>
#include <string>
#include <vector>
#include <stack>
using namespace std;
using namespace std::chrono;

int extraPoints;

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
	bool operator == (const Point& p) const
	{
		return (x == p.x and y == p.y);
	}
	// prints Point as string to console
	string toString()
	{
		return "(" + to_string(x) + "," + to_string(y) + ")";
	}
};

class Line
{
private:
	Point one;
	Point two;
	Point orgOne;
	Point orgTwo;
	vector<double> standardForm;
	vector<double> interceptForm;
public:
	Line();
	Line(Point p1, Point p2)
	{
		orgOne = p1;
		orgTwo = p2;
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
	Point get_first_pt()
	{
		return orgOne;
	}
	// gets intersection of this and x=800
	Point get_second_point()
	{
		return two;
	}
	Point get_second_pt()
	{
		return orgTwo;
	}
};

class Graph
{
private:
	char(*arr)[400];
	list<Point> coordinates;

public:
	const int size = 400;
	// helper method that instantiates graphs with white pixels
	void fill_graph()
	{
		arr = new char[400][400];
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
	// draws original points onto to graph
	void draw_figures()
	{
		list<Point>::iterator it;
		for (it = coordinates.begin(); it != prev(coordinates.end(), extraPoints); ++it) {
			create_circle(3, round((it->getX()) * 400), round((it->getY()) * 400), 'R');
		}
		for (it = prev(coordinates.end(), extraPoints); it != prev(coordinates.end()); ++it) {
			create_line(round((it->getX()) * 400), round((it->getY()) * 400), round((next(it)->getX()) * 400), round((next(it)->getY()) * 400), 'B');
		}
	}
	// creates ppm file with origin on the bottom left
	void create_file(string filename)
	{
		ofstream myfile;
		myfile.open(filename, ios::out);
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
void create_pt_file(vector<Point> points, string fileName)
{
	ofstream myfile;
	myfile.open(fileName, ios::out);
	for (int i = 0; i < (int)points.size(); i++) {
		myfile << setprecision(23) << points[i].getX() << " " << points[i].getY() << fixed << endl;
	}
	myfile.close();
}

// reads points.txt file
vector<Point> read_pt_file(string fileName)
{
	string line;
	vector<Point> coords;
	ifstream myfile(fileName);
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

bool compare_X_coord(Point& p1, Point& p2)
{
	return (p1.getX() < p2.getX());
}

// splits points into those left or right of a line
void divide_points(Line l, vector<Point>& v, vector<Point>& direction)
{
	vector<double> coeffs = l.get_intercept_form();
	double temp = 0.0;
	bool isRight = l.get_first_pt().getX() < l.get_second_pt().getX();
	for (int i = 0; i < (int)v.size(); i++) {
		temp = (coeffs[0] * v[i].getX()) + coeffs[1];
		if (isRight == true) {
			if (v[i].getY() > temp) {
				direction.push_back(v[i]);
			}
		}
		else {
			if (v[i].getY() < temp) {
				direction.push_back(v[i]);
			}
		}
	}
}

// returns point that is furthest away from a line
int determine_farthest(vector<Point>& v, vector<Point> hull, Line l)
{
	int p;
	double temp;
	double temp2;
	vector<double> coeffs = l.get_standard_form();
	for (int i = 0; i < (int)v.size(); i++) {
		temp2 = abs(coeffs[0] * v[i].getX() + coeffs[1] * v[i].getY() - coeffs[2]) / (pow(coeffs[0], 2) + pow(coeffs[1], 2));
		if (temp2 >= temp) {
			temp = temp2;
			p = i;
		}
	}
	return p;
}

// adds a point to the convex hull
void add_to_hull(vector<Point>& v, vector<Point>& hull, Point& p, Point& q, Point& temp)
{
	if (p.getX() < temp.getX() and q.getX() < temp.getX()) {
		hull.insert(find(hull.begin(), hull.end(), p), temp);
	}
	else if (p.getX() > temp.getX() and q.getX() > temp.getX()) {
		hull.insert(find(hull.begin(), hull.end(), q) + 1, temp);
	}
	else {
		hull.insert(find(hull.begin(), hull.end(), p) + 1, temp);
	}
}

// recursive method to find the convex hull
void convex_hull(vector<Point>& v, vector<Point>& hull, Point p, Point q)
{
	int size = (int)v.size();
	if (size == 0) {
		return;
	}
	else if (size == 1) {
		add_to_hull(v, hull, p, q, v[0]);
	}
	else {
		int farthest = determine_farthest(v, hull, Line(p, q));
		Point temp = v[farthest];
		v.erase(v.begin() + farthest);
		add_to_hull(v, hull, p, q, temp);

		vector<Point> left, right;
		divide_points(Line(temp, q), v, right);
		divide_points(Line(p, temp), v, left);

		convex_hull(right, hull, temp, q);
		convex_hull(left, hull, p, temp);
	}
}

// sets up the convex hull for quick hull method
void part1()
{
	Graph g;
	g.fill_graph();
	srand(time(0));
	vector<Point> pts;
	for (int i = 0; i < 60; i++) {
		pts.push_back(Point((double)rand() / (RAND_MAX), (double)rand() / (RAND_MAX)));
	}
	create_pt_file(pts, "points.txt");

	vector<Point> sorted = pts;
	sort(sorted.begin(), sorted.end(), compare_X_coord);
	Point a = sorted[0], b = sorted[sorted.size() - 1];
	vector<Point> hull = { a, b, a };

	vector<Point> left, right;
	sorted.pop_back();
	sorted.erase(sorted.begin());
	divide_points(Line(b, a), sorted, left);
	divide_points(Line(a, b), sorted, right);

	convex_hull(left, hull, b, a);
	convex_hull(right, hull, a, b);

	extraPoints = hull.size();
	pts.insert(pts.end(), hull.begin(), hull.end());
	g.set_coords(pts);
	g.draw_figures();
	g.create_file("quickhull.ppm");
}

int find_smallest_y(vector<Point>& v)
{
	int smallest = 0;
	double smallestY = v[0].getY(), temp;
	for (int i = 1; i < (int)v.size(); i++) {
		temp = v[i].getY();
		if (temp < smallestY) {
			smallestY = temp;
			smallest = i;
		}
		else if (temp == smallestY) {
			if (v[i].getX() < v[smallest].getX()) {
				smallest = i;
			}
		}
	}
	return smallest;
}

Point p0;

double get_polar_angle(Point& p1, Point& p2)
{
	return atan2(p2.getY() - p1.getY(), p2.getX() - p1.getX());
}

bool compare_polar_angle(Point& p1, Point& p2)
{
	return (get_polar_angle(p1, p0) < get_polar_angle(p2, p0));
}

bool isRightTurn(Point& prev, Point& top, Point& current)
{
	if ((top.getX() - prev.getX()) * (current.getY() - top.getY()) - (top.getY() - prev.getY()) * (current.getX() - top.getX()) < 0.0) {
		return true;
	}
	else {
		return false;
	}
}

void part2()
{
	Graph g;
	g.fill_graph();

	vector<Point> orgPts = read_pt_file("points.txt");
	vector<Point> pts = orgPts;
	int idx = find_smallest_y(pts);
	p0 = pts[idx];
	vector<Point> hull = { p0 };

	pts.erase(pts.begin() + idx);
	sort(pts.begin(), pts.end(), compare_polar_angle);
	stack<Point> s;
	s.push(p0);
	s.push(pts[0]);
	pts.erase(pts.begin());

	Point top, prev;
	for (Point& point : pts) {
		while (s.size() >= 2) {
			top = s.top();
			s.pop();
			prev = s.top();
			if (isRightTurn(prev, top, point) == false) {
				s.push(top);
				break;
			}
		}
		s.push(point);
	}

	while (s.size() > 0) {
		hull.push_back(s.top());
		s.pop();
	}

	extraPoints = hull.size();
	orgPts.insert(orgPts.end(), hull.begin(), hull.end());
	g.set_coords(orgPts);
	g.draw_figures();
	g.create_file("grahamscan.ppm");
}

int main()
{
	// 	part1();
	part2();
}