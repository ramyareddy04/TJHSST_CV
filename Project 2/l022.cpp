// Ramya Reddy
// Period 4
// 9/18/2021

#include<iostream>
#include <iomanip> 
#include <fstream>
#include<math.h>
#include <vector>
using namespace std;

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
	vector<double> standardForm;
	vector<double> interceptForm;
public:
	Line();
	Line(Point num1, Point num2)
	{
		standardForm = find_equation_of_line({ num1, num2 });
		one = Point(0.0, standardForm[2] / standardForm[1]);
		two = Point(1.0, (standardForm[2] - standardForm[0]) / standardForm[1]);
		interceptForm = { two.getY() - one.getY(), one.getY() };
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
	// determine perpendicular line from given point
	Line find_perp(Point p)
	{
		double slope = -1 / interceptForm[0];
		return Line(p, Point(p.getX() + 0.01, p.getY() + slope / 100));
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
	int(*arr)[800];
	vector<Point> coordinates;

public:
	const int size = 800;
	// helper method that instantiates graphs with white pixels
	void fill_graph()
	{
		arr = new int[800][800];
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				arr[i][j] = 1;
			}
		}
	}
	// custom rounding method
	int round(double val)
	{
		return int(val) + int(abs(int(val) - val) >= 0.5);
	}
	// helper method that scales and rounds given points
	vector<int> scale(vector<Point> points)
	{
		vector<int> scaled_points;
		for (int i = 0; i < int(points.size()); i++) {
			scaled_points.push_back(round(points[i].getX() * size));
			scaled_points.push_back(round(points[i].getY() * size));
		}
		return scaled_points;
	}
	// helper method that changes color of pixel
	void illuminate(int x, int y)
	{
		if ((x >= 0 and x <= size) and (y >= 0 and y <= size)) {
			arr[x][y] = 0;
		}
	}
	// helper method that creates line with negative slope and driving x axis
	void negative_x(int x1, int y1, int x2, int y2, int deltaX, int deltaY, int error)
	{
		int j = y1;
		for (int i = x1; i >= x2; i--) {
			illuminate(i, j);
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
	void positive_x(int x1, int y1, int x2, int y2, int deltaX, int deltaY, int error)
	{
		int j = y1;
		for (int i = x1; i <= x2; i++) {
			illuminate(i, j);
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
	void positive_y(int x1, int y1, int x2, int y2, int deltaX, int deltaY, int error)
	{
		int i = x1;
		for (int j = y1; j <= y2; j++) { // if true move down else move up
			illuminate(i, j);
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
	void negative_y(int x1, int y1, int x2, int y2, int deltaX, int deltaY, int error)
	{
		int i = x1;
		for (int j = y1; j >= y2; j--) { // if true move down else move up
			illuminate(i, j);
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
	void create_line(int x1, int y1, int x2, int y2)
	{
		int deltaX = x2 - x1;
		int deltaY = y2 - y1;
		int error = deltaY - deltaX;
		if (abs(deltaX) > abs(deltaY)) { // if true driving axis is x-axis else driving axis is y-axis
			if (deltaX < 0) { // if true move right else move left
				negative_x(x1, y1, x2, y2, deltaX, deltaY, error);
			}
			else {
				positive_x(x1, y1, x2, y2, deltaX, deltaY, error);
			}
		}
		else {
			if (deltaY < 0) {
				negative_y(x1, y1, x2, y2, deltaX, deltaY, error);
			}
			else {
				positive_y(x1, y1, x2, y2, deltaX, deltaY, error);
			}
		}
	}
	// rasterized algorithm that draws circle
	void create_circle(int r, int x0, int y0) {
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
			illuminate(x + x0, y + y0);
			illuminate(x + x0, -y + y0);
			illuminate(-x + x0, y + y0);
			illuminate(-x + x0, -y + y0);
			illuminate(y + x0, x + y0);
			illuminate(-y + x0, x + y0);
			illuminate(y + x0, -x + y0);
			illuminate(-y + x0, -x + y0);
			y2_new -= (2 * x) - 3;
		}
	}
	// method that finds 'Point E' to create square
	Point find_e(Line curr, Point org, Point other, double length)
	{
		double slope = curr.get_intercept_form()[0];
		double magnitude = sqrt(pow(slope, 2) + 1);
		double dt = length / magnitude;
		double dX = dt;
		double dY = slope * dt;
		return Point(org.getX() + dX, org.getY() + dY);
	}
	// draws original 4 points onto to graph
	void draw_figures()
	{
		vector<int> temp = scale(coordinates);
		for (int i = 0; i<int(coordinates.size()); i++) {
			create_circle(2, temp[2 * i], temp[2 * i + 1]);
		}
	}
	// creates vectors of line, that when intersect, form a square
	vector<Line> create_square(Point a, Point b, Point c, Point d)
	{
		Line step1(a, b);
		Line step2 = step1.find_perp(d);
		Point e = find_e(step2, d, c, a.find_length(b));
		Line step4(c, e);
		Line step5a = step4.find_perp(a);
		Line step5b = step4.find_perp(b);
		Line step6 = step5a.find_perp(d);
		vector<Line> lines = { step4, step5a, step5b, step6 };
		return lines;
	}
	// iterates through all possible distinct squares that can be made
	vector<vector<Line>> create_squares(vector<Point> p)
	{
		vector<vector<int>> ways = { {0,1,2,3},{0,1,3,2},{0,2,1,3},{0,2,3,1},{0,3,1,2},{0,3,2,1},{1,0,2,3},{1,0,3,2},{1,2,3,0},{1,2,0,3},{1,3,2,0},{1,3,0,2},{2,0,3,1},{2,3,0,1},{2,0,1,3},{2,1,0,3},{2,1,3,0},{2,3,1,0},{3,0,1,2},{3,0,2,1},{3,1,2,0},{3,1,0,2},{3,2,0,1},{3,2,1,0} };
		vector<vector<Line>> lines;
		for (int i = 0; i < (int)ways.size(); i++) {
			lines.push_back(create_square(p[ways[i][0]], p[ways[i][1]], p[ways[i][2]], p[ways[i][3]]));
		}
		return lines;
	}
	// draws the smallest square onto graph
	void create_smallest_square(vector<Line> lines)
	{
		vector<int> temp;
		for (int i = 0; i < (int)lines.size(); i++) {
			temp = scale({ lines[i].get_first_point(), lines[i].get_second_point() });
			create_line(temp[0], temp[1], temp[2], temp[3]);
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
				myfile << arr[j][i] << " " << arr[j][i] << " " << arr[j][i] << " ";
			}
			myfile << endl;
		}
		myfile.close();
	}
	// gets coordinates vector
	vector<Point> get_coords()
	{
		return coordinates;
	}
	// sets coordinates vector
	void set_coords(vector<Point> p)
	{
		coordinates = p;
	}
	// destructor for Graph class
	~Graph()
	{
		if (arr != NULL) {
			delete[] arr;
		}
	}
};

double calculateArea(vector<Point> points)
{
	double a = sqrt((pow(points[0].getX() - points[1].getX(), 2) + pow(points[0].getY() - points[1].getY(), 2)));
	double b = sqrt((pow(points[1].getX() - points[2].getX(), 2) + pow(points[1].getY() - points[2].getY(), 2)));
	double c = sqrt((pow(points[0].getX() - points[2].getX(), 2) + pow(points[0].getY() - points[2].getY(), 2)));
	double s = (a + b + c) / 2;
	return sqrt(s * (s - a) * (s - b) * (s - c));
}
// method that determines whether point is within triangle or not using the area method discussed in class
bool isInTriangle(vector<Point> points, Point target)
{
	double tot = calculateArea(points);
	double a = calculateArea({ points[0], points[1], target });
	double b = calculateArea({ points[0], points[2], target });
	double c = calculateArea({ points[2], points[1], target });
	double eps = pow(2, -52);
	return (a + b + c >= tot - eps) and (a + b + c <= tot + eps);
}
// method that goes through all possible triangles within 4 points to ensure a point is not within another triangle
bool isInIntersection(vector<Point> points)
{
	return (isInTriangle({ points[0], points[1], points[2] }, points[3])
		or isInTriangle({ points[3], points[1], points[2] }, points[0])
		or isInTriangle({ points[0], points[3], points[2] }, points[1])
		or isInTriangle({ points[0], points[1], points[3] }, points[2]));
}
// method that writes points to .txt file
void create_pt_file(vector<Point> points)
{
	ofstream myfile;
	myfile.open("points.txt", ios::out);
	for (int i = 0; i<int(points.size()) - 1; i++) {
		myfile << setprecision(17) << "(" << points[i].getX() << "," << points[i].getY() << ") ," << fixed << " ";
	}
	myfile << setprecision(17) << "(" << points.back().getX() << "," << points.back().getY() << ")" << fixed;
	myfile.close();
}
// method that generates four random points
vector<Point> part1()
{
	srand(time(0));
	vector<Point> points;
	for (int i = 0; i < 4; i++) {
		points.push_back(Point((double)rand() / (RAND_MAX), (double)rand() / (RAND_MAX)));
	}
	bool inIntersection = isInIntersection(points);
	while (inIntersection) {
		points[3].setX((double)rand() / (RAND_MAX));
		points[3].setY((double)rand() / (RAND_MAX));
		inIntersection = isInIntersection(points);
	}
	create_pt_file(points);
	return points;
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
			int beginIdx = 0;
			int endIdx = 0;
			int split = 0;
			for (int i = 0; i<int(line.size()); i++) {
				if (line[i] == '(') {
					beginIdx = i;
				}
				if (line[i] == ')') {
					endIdx = i;
					coords.push_back(Point(stod(line.substr(beginIdx + 1, split - beginIdx - 1)), stod(line.substr(split + 1, endIdx - split - 1))));
				}
				if (line[i] == ',') {
					split = i;
				}
			}
		}
	}
	myfile.close();
	return coords;
}
// determines intersections of given lines and returns points in a vector
vector<vector<Point>> get_intersections(vector<vector<Line>> lines)
{
	vector<vector<Point>> p;
	vector<Line> line;
	for (int i = 0; i < (int)lines.size(); i++) {
		line = lines[i];
		p.push_back({ line[0].find_intersection(line[1]), line[0].find_intersection(line[2]), line[1].find_intersection(line[3]), line[2].find_intersection(line[3]) });
	}
	return p;
}
// finds the area of a square by squaring the length of one side
double get_square_area(vector<Point> p)
{
	return pow(p[0].find_length(p[1]), 2);
}
// checks if area is distinct
bool in_vector(vector<double> v, double d)
{
	double eps = pow(10, -5);
	for (int i = 0; i < (int)v.size(); i++) {
		if (abs(v[i] - d) <= eps) {
			return true;
		}
	}
	return false;
}
// determines distinct squares as well as square with smallest area
vector<vector<double>> compareArea(vector<vector<Point>> p)
{
	vector<vector<double>> sorted = { {},{} };
	double curr_min = 640000;
	double temp;
	int in;
	for (int i = 0; i < (int)p.size(); i++) {
		temp = get_square_area(p[i]);
		if (in_vector(sorted[0], temp) == false) {
			if (temp < curr_min) {
				curr_min = temp;
				in = i;
			}
			sorted[0].push_back(temp);
			sorted[1].push_back(i);
		}
		//         }
	}
	sorted[1].push_back(double(in));
	return sorted;
}
// method that writes coordinates and vertices to .txt file
void create_final_file(vector<Point> coords, vector<vector<Point>> points, vector<double> areas)
{
	ofstream myfile;
	myfile.open("output.txt", ios::out);
	for (int j = 0; j<int(coords.size()) - 1; j++) {
		myfile << setprecision(17) << "(" << coords[j].getX() << "," << coords[j].getY() << ")," << fixed << " ";
	}
	myfile << setprecision(17) << "(" << coords.back().getX() << "," << coords.back().getY() << ")" << fixed << endl;
	for (int i = 0; i<int(areas.size()); i++) {
		for (int j = 0; j<int(points[i].size()) - 1; j++) {
			myfile << setprecision(17) << "(" << points[i][j].getX() << "," << points[i][j].getY() << ") ," << fixed << " ";
		}
		myfile << setprecision(17) << "(" << points[i].back().getX() << "," << points[i].back().getY() << ") Area = " << areas[i] << fixed << endl;
	}
	myfile.close();
}
// method that generates all possible squares from four coordinates made in part1() method and draws smallest square
void part2()
{
	Graph g;
	g.fill_graph();
	g.set_coords(read_pt_file());
	g.draw_figures();
	vector<vector<Line>> allLines = g.create_squares(g.get_coords());
	vector<vector<Point>> intersections = get_intersections(allLines);
	vector<vector<double>> areas = compareArea(intersections);
	g.create_smallest_square(allLines[int(areas[1].back())]);
	g.create_file();
	create_final_file(g.get_coords(), intersections, areas[0]);
}

int main()
{
	//     auto result = part1();
	part2();
}