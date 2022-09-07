// Ramya Reddy
// Period 4
// 11/16/2021

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
	double find_length_sqrd_spec(Point* other)
	{
		return pow(getX() - other->getX(), 2) + pow(getY() - other->getY(), 2);
	}
	// prints Point as string to console
	string toString()
	{
		return "(" + to_string(x) + "," + to_string(y) + ")";
	}
};

class hash_pair {
public:
	size_t operator()(const pair<int, int> p) const
	{
		return p.first + p.second;
	}
};

class TwoPoints
{
private:
	double distance;
	Point one;
	Point two;
public:
	TwoPoints();
	TwoPoints(Point p1, Point p2)
	{
		distance = p1.find_length_sqrd(p2);
		one = p1;
		two = p2;
	}
	TwoPoints(double dist, Point p1, Point p2)
	{
		distance = dist;
		one = p1;
		two = p2;
	}
	TwoPoints(Point* p1, Point* p2, double dist)
	{
		distance = dist;
		one = *p1;
		two = *p2;
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
		return to_string(pow(get_sqrd_distance(), 2)) + " " + get_first_point().toString() + " " + get_second_point().toString();
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
void part1()
{
	srand(time(0));
	list<Point> points;
	for (int i = 0; i < 200000; i++) {
		points.push_back(Point((double)rand() / (RAND_MAX), (double)rand() / (RAND_MAX)));
	}
	create_pt_file(points);
	// 	auto start = high_resolution_clock::now();
	// 	vector<Point> finPts = brute_force(points);
	// 	auto stop = high_resolution_clock::now();
	// 	times.push_back(duration_cast<microseconds>(stop - start));
	// 	points.push_back(finPts[0]);
	// 	points.push_back(finPts[1]);
	// 	return TwoPoints(finPts[0], finPts[1]);
}

// reads points.txt file
vector<Point> read_pt_file()
{
	string line;
	vector<Point> coords;
	ifstream myfile("points10k.txt");
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
bool compare_X_coord(Point& p1, Point& p2)
{
	return (p1.getX() < p2.getX());
}


// method that sorts vector of Points
bool compare_Y_coord(Point& p1, Point& p2)
{
	return (p1.getY() < p2.getY());
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
TwoPoints brute_force2(int midPt, vector<Point>& v, int pts, int pts2)
{
	double curr_min = 1000;
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
	return TwoPoints(curr_min, pt1, pt2);
}

// creates vectors to use for brute force in strip
TwoPoints val_check(TwoPoints smallest, vector<Point>& pts, int i, int k)
{
	int mid = (i + k) / 2;
	double midX = pts[mid].getX();
	double dist = pow(smallest.get_sqrd_distance(), 0.5);
	double posDif = (midX + dist);
	double negDif = (midX - dist);

	int idxL = mid - 1;
	while (pts[idxL].getX() > negDif and idxL >= i) {
		idxL -= 1;
	}
	int idxR = mid + 1;
	while (pts[idxR].getX() < posDif and idxR <= k) {
		idxR += 1;
	}

	if (idxL != mid - 1 and idxR != mid + 1) {
		smallest = find_min(smallest, brute_force2(mid, pts, idxL, idxR));
	}
	else {
		smallest = find_min(smallest, brute_force2(-1, pts, idxL, idxR));
	}
	return smallest;
}

// recursive merge sort method
TwoPoints recur(vector<Point>& v, int i, int k)
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

// implementation of brute force approach for strip in part 3
void brute_force2_opt(double distance, double margin, TwoPoints& smallest, vector<Point>& v)
{
	double curr_min = distance;
	Point pt1;
	Point pt2;
	double temp;
	vector<Point>::iterator it;
	vector<Point>::iterator it2;
	vector<Point>::iterator it3;
	bool isChanged = false;
	for (it = v.begin(); it != prev(v.end()); ++it) {
		if (it + 15 < v.end()) {
			it3 = next(it, 15);
		}
		else {
			it3 = v.end();
		}
		for (it2 = next(it); it2 != it3; ++it2) {
			if (abs((*it).getY() - (*it2).getY()) < 2 * margin) {
				temp = (*it).find_length_sqrd(*it2);
				if (temp < curr_min) {
					curr_min = temp;
					pt1 = *it;
					pt2 = *it2;
					isChanged = true;
				}
			}
			else {
				break;
			}
		}
	}
	if (isChanged == true) {
		smallest = TwoPoints(curr_min, pt1, pt2);
	}
}

// creates vectors to use for brute force in strip for part 3
TwoPoints val_check_opt(TwoPoints smallest, vector<Point>& pts, int i, int k)
{
	int mid = (i + k) / 2;
	double midX = pts[mid].getX();
	double dist = pow(smallest.get_sqrd_distance(), 0.5);
	double margin = dist;
	double posDifX = (midX + margin);
	double negDifX = (midX - margin);

	int idxL = mid - 1;
	while (pts[idxL].getX() > negDifX and idxL >= i) {
		idxL -= 1;
	}
	int idxR = mid + 1;
	while (pts[idxR].getX() < posDifX and idxR <= k) {
		idxR += 1;
	}
	vector<Point> pts2(pts.begin() + idxL, pts.begin() + idxR);
	sort(pts2.begin(), pts2.end(), compare_Y_coord);
	brute_force2_opt(smallest.get_sqrd_distance(), margin, smallest, pts2);
	return smallest;
}

// recursive merge sort method for part 3
TwoPoints recur_opt(vector<Point>& v, int i, int k)
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
		TwoPoints l = recur_opt(v, i, j);
		TwoPoints r = recur_opt(v, j, k);
		return val_check_opt(find_min(l, r), v, i, k);
	}
}

// method that reads points from file and implements recursive approach
TwoPoints part3()
{
	vector<Point> pts = read_pt_file();
	auto start3 = high_resolution_clock::now();
	sort(pts.begin(), pts.end(), compare_X_coord);
	TwoPoints smallest = recur_opt(pts, 0, pts.size());
	auto stop3 = high_resolution_clock::now();
	times.push_back(duration_cast<microseconds>(stop3 - start3));
	return smallest;
}

// implementation of Knuth Shuffle for part 4
void knuth_shuffle(vector<Point>& v)
{
	int n = (int)v.size();
	int randIdx;
	Point temp;
	for (int i = 0; i < n - 1; i++) {
		randIdx = rand() % (n - i) + i;
		temp = v[randIdx];
		v[randIdx] = v[i];
		v[i] = temp;
	}
}

// assigning a row and column to each subsquare
pair<unsigned long long int, unsigned long long int> assign_subsquare(Point* p, double delta)
{
	double val = (p->getX() / delta);
	if (val != 0 and val - int(val) == 0) {
		val += 1;
	}
	unsigned long long int x = val;
	val = (p->getY() / delta);
	if (val != 0 and val - int(val) == 0) {
		val += 1;
	}
	unsigned long long int y = val;
	return make_pair(x, y);
}

// changes row and col of each subsquare
void reassign_coords(unordered_map<pair<unsigned long long int, unsigned long long int>, Point*, hash_pair>& map, vector<Point>& v, double subsqrLen, int idx)
{
	Point* p;
	for (int i = 0; i <= idx; i++) {
		p = &v[i];
		map[assign_subsquare(p, subsqrLen)] = &v[i];
	}
}

// implementation of Khuller Algorithim for part 4
TwoPoints khuller_algorithm(vector<Point>& v)
{
	Point* p1 = &v[0];
	Point* p2 = &v[1];
	double dist = p1->find_length_sqrd_spec(p2);
	double delta = pow(dist, 0.5);
	unsigned long long int maxSqr = 1 / (delta / 2);

	unordered_map<pair<unsigned long long int, unsigned long long int>, Point*, hash_pair> map;
	map[assign_subsquare(p1, delta / 2)] = &v[0];
	map[assign_subsquare(p2, delta / 2)] = &v[1];

	double temp;
	pair<unsigned long long int, unsigned long long int> coord, tempCoord;
	unsigned long long int upperX = 0, upperY = 0, lowerX = 0, lowerY = 0;
	bool isChanged;
	for (int i = 2; i < (int)v.size(); i++) {
		coord = assign_subsquare(&v[i], delta / 2);
		lowerX = max(coord.first - 2, (unsigned long long int)0);
		lowerY = max(coord.second - 2, (unsigned long long int)0);
		upperX = min(coord.first + 2, maxSqr);
		upperY = min(coord.second + 2, maxSqr);
		for (unsigned long long int x = lowerX; x <= upperX; x++) {
			for (unsigned long long int y = lowerY; y <= upperY; y++) {
				tempCoord = make_pair(x, y);
				if (map.count(tempCoord) > 0) {
					temp = map[tempCoord]->find_length_sqrd(v[i]);
					if (temp < dist) {
						dist = temp;
						p2 = map[tempCoord];
						isChanged = true;
					}
				}
			}
		}
		if (isChanged == true) {
			map.clear();
			p1 = &v[i];
			delta = pow(dist, 0.5);
			maxSqr = 1 / (delta / 2);
			reassign_coords(map, v, delta / 2, i);
			isChanged = false;
		}
		else {
			map[coord] = &v[i];
		}
	}
	return TwoPoints(p1, p2, delta);
}

TwoPoints part4()
{
	vector<Point> pts = read_pt_file();
	auto start4 = high_resolution_clock::now();
	knuth_shuffle(pts);
	TwoPoints smallest = khuller_algorithm(pts);
	auto stop4 = high_resolution_clock::now();
	times.push_back(duration_cast<microseconds>(stop4 - start4));
	return smallest;
}

// writes results to terminal and output file
void create_results_file(vector<TwoPoints> points)
{
	ofstream myfile;
	myfile.open("results.txt", ios::out);
	for (int i = 0; i < (int)points.size(); i++) {
		myfile << "Results for Part " << (i + 3) << ":" << endl;
		cout << "Results for Part " << (i + 3) << ":" << endl;
		if (i == 0) {
			myfile << setprecision(23) << "Points: " << "(" << points[i].get_first_point().getX() << "," << points[i].get_first_point().getY() << ") (" << points[i].get_second_point().getX() << "," << points[i].get_second_point().getY() << ")" << " Distance: " << points[i].get_distance() << " Time (seconds): " << times[i].count() << endl;
			cout << setprecision(23) << "Points: " << "(" << points[i].get_first_point().getX() << "," << points[i].get_first_point().getY() << ") (" << points[i].get_second_point().getX() << "," << points[i].get_second_point().getY() << ")" << " Distance: " << points[i].get_distance() << " Time (seconds): " << times[i].count() << endl;
		}
		else {
			myfile << "Points: " << "(" << points[i].get_first_point().getX() << "," << points[i].get_first_point().getY() << ") (" << points[i].get_second_point().getX() << "," << points[i].get_second_point().getY() << ")" << " Distance: " << points[i].get_distance() << " Time (seconds): " << times[i].count() << endl;
			cout << "Points: " << "(" << points[i].get_first_point().getX() << "," << points[i].get_first_point().getY() << ") (" << points[i].get_second_point().getX() << "," << points[i].get_second_point().getY() << ")" << " Distance: " << points[i].get_distance() << " Time (seconds): " << times[i].count() << endl;
		}
	}
	myfile.close();
}

int main()
{
	TwoPoints res3 = part3();
	TwoPoints res4 = part4();
	create_results_file({ res3, res4 });
}