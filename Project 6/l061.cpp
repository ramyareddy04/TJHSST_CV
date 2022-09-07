// Ramya Reddy
// Period 4
// 1/31/2022

#include<iostream>
#include <iomanip> 
#include <fstream>
#include<math.h>
#include <string>
#include <vector>
#include <list>
#include <cstring>
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
	int width, height;
	double(*arr)[10000][3];
	list<Point> coordinates;
public:
	Graph();
	Graph(int w, int h)
	{
		width = w;
		height = h;
	}
	// helper method that instantiates graphs with white pixels
	void fill_graph()
	{
		arr = new double[height][10000][3];
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				arr[i][j][0] = 0;
				arr[i][j][1] = 0;
				arr[i][j][2] = 0;
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
			scaled_points.push_back(round(it->getX() * width));
			scaled_points.push_back(round(it->getY() * height));
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
	void illuminate(int x, int y, double colorR, double colorG, double colorB)
	{
		if ((x >= 0 and x <= height) and (y >= 0 and y <= width)) {
			if (colorR == -1) {
				arr[x][y][0] += 1;
				arr[x][y][1] += 1;
				arr[x][y][2] += 1;
			}
			else {
				arr[x][y][0] = colorR;
				arr[x][y][1] = colorG;
				arr[x][y][2] = colorB;
			}
		}
	}
	// helper method that creates line with negative slope and driving x axis
	void negative_x(int x1, int y1, int x2, int y2, int deltaX, int deltaY, int error, int colorR, int colorG, int colorB)
	{
		int j = y1;
		for (int i = x1; i >= x2; i--) {
			illuminate(i, j, colorR, colorG, colorB);
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
	void positive_x(int x1, int y1, int x2, int y2, int deltaX, int deltaY, int error, int colorR, int colorG, int colorB)
	{
		int j = y1;
		for (int i = x1; i <= x2; i++) {
			illuminate(i, j, colorR, colorG, colorB);
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
	void positive_y(int x1, int y1, int x2, int y2, int deltaX, int deltaY, int error, int colorR, int colorG, int colorB)
	{
		int i = x1;
		for (int j = y1; j <= y2; j++) { // if true move down else move up
			illuminate(i, j, colorR, colorG, colorB);
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
	void negative_y(int x1, int y1, int x2, int y2, int deltaX, int deltaY, int error, int colorR, int colorG, int colorB)
	{
		int i = x1;
		for (int j = y1; j >= y2; j--) { // if true move down else move up
			illuminate(i, j, colorR, colorG, colorB);
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
	void create_line(int x1, int y1, int x2, int y2, int colorR, int colorG, int colorB)
	{
		int deltaX = x2 - x1;
		int deltaY = y2 - y1;
		int error = deltaY - deltaX;
		if (abs(deltaX) > abs(deltaY)) { // if true driving axis is x-axis else driving axis is y-axis
			if (deltaX < 0) { // if true move right else move left
				negative_x(x1, y1, x2, y2, deltaX, deltaY, error, colorR, colorG, colorB);
			}
			else {
				positive_x(x1, y1, x2, y2, deltaX, deltaY, error, colorR, colorG, colorB);
			}
		}
		else {
			if (deltaY < 0) {
				negative_y(x1, y1, x2, y2, deltaX, deltaY, error, colorR, colorG, colorB);
			}
			else {
				positive_y(x1, y1, x2, y2, deltaX, deltaY, error, colorR, colorG, colorB);
			}
		}
	}
	// rasterized algorithm that draws circle
	void create_circle(int r, int x0, int y0, int colorR, int colorG, int colorB) {
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
			illuminate(x + x0, y + y0, colorR, colorG, colorB);
			illuminate(x + x0, -y + y0, colorR, colorG, colorB);
			illuminate(-x + x0, y + y0, colorR, colorG, colorB);
			illuminate(-x + x0, -y + y0, colorR, colorG, colorB);
			illuminate(y + x0, x + y0, colorR, colorG, colorB);
			illuminate(-y + x0, x + y0, colorR, colorG, colorB);
			illuminate(y + x0, -x + y0, colorR, colorG, colorB);
			illuminate(-y + x0, -x + y0, colorR, colorG, colorB);
			y2_new -= (2 * x) - 3;
		}
	}
	// creates ppm file with origin on the bottom left
	void create_file(string filename, int max)
	{
		ofstream myfile;
		myfile.open(filename, ios::out);
		myfile << "P3 \n" << width << " " << height << "\n" << max << "\n";
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				myfile << arr[i][j][0] << " " << arr[i][j][1] << " " << arr[i][j][2] << " ";
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
	double get_coord(int i, int j)
	{
		return arr[i][j][0];
	}
	double get_coord(int i, int j, int RGB)
	{
		return arr[i][j][RGB];
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
			arr = NULL;
			delete[] arr;
		}
	}
};

// creates vector that contains RGB values for every single pixel
vector<vector<int>> read_ppm_file(string fileName)
{
	string line;
	int numLine = 1;
	int w, h, i, j, max;
	vector<int> temp;
	vector<vector<int>> values;
	ifstream myfile(fileName);
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			i = 0;
			j = 0;
			if (numLine == 1) {

			}
			else if (numLine == 2) {
				while (i <int(line.size())) {
					if (line[i] == ' ') {
						w = stod(line.substr(0, i));
						h = stod(line.substr(i + 1, line.length()));
						i = int(line.size());
					}
					i += 1;
				}
			}
			else if (numLine == 3) {
				max = stod(line);
				values.push_back({ w, h, max });
			}
			else {
				while (i<int(line.size())) {
					if (line[i] == ' ') {
						temp.push_back(stod(line.substr(j, i)));
						j = i + 1;
					}
					i += 1;
				}
			}
			numLine += 1;
		}
	}

	for (int k = 0; k < (int)temp.size();) {
		values.push_back({ temp[k], temp[k + 1], temp[k + 2] });
		k += 3;
	}
	myfile.close();
	return values;
}

// returns grayscale array of image
void create_greyscale(vector<vector<int>> v)
{
	int w = v[0][0], h = v[0][1], max = v[0][2];
	int row, col;
	v.erase(v.begin());
	Graph g(w, h);
	g.fill_graph();
	int temp;
	for (int i = 0; i < (int)v.size(); i++) {
		row = i / w;
		col = i % w;
		temp = (v[i][0] + v[i][1] + v[i][2]) / 3;
		g.illuminate(row, col, temp, temp, temp);
	}
	g.create_file("imageg.ppm", max);
}

int width, height;
// recusive area fill method
void area_fill(int row, int col, Graph& img)
{
	if (row >= 0 and row < height and col >= 0 and col < width) {
		if (img.get_coord(row, col) == 0 or img.get_coord(row, col) == 1) {
			return;
		}
		img.illuminate(row, col, 1, 1, 1);

		area_fill(row - 1, col, img);
		area_fill(row + 1, col, img);
		area_fill(row, col - 1, img);
		area_fill(row, col + 1, img);
		area_fill(row - 1, col - 1, img);
		area_fill(row - 1, col + 1, img);
		area_fill(row + 1, col - 1, img);
		area_fill(row + 1, col + 1, img);
	}
}

int threshold = 70, threshold2 = 210, threshold3 = 60;
string name = "image.ppm";

int round_to_multiple(double inp)
{
	int temp = inp / 45;
	if (inp - temp * 45 <= 22.5) {
		return temp * 45;
	}
	else {
		return (temp + 1) * 45;
	}
}

// performs sobel operator on image, and returns image using two thresholds and direction of gradient
void sobel_conv(vector<vector<int>> v, vector<vector<int>> org)
{
	vector<vector<int>> xKernel = { {-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1} }, yKernel = { {-1, -2, -1}, {0, 0, 0}, {1, 2, 1} };
	int w = v[0][0], h = v[0][1];
	v.erase(v.begin());
	org.erase(org.begin());

	Graph g(w, h);
	g.fill_graph(); // CONTAINS ULT. EDGES
	Graph mag(w, h);
	mag.fill_graph(); // CONTAINS MAG. FOR EACH PIXEL
	Graph angle(w, h);
	angle.fill_graph(); // CONTAINS EXACT ANGLE FOR EACH PIXEL
	Graph votes(w, h);
	votes.fill_graph(); // CONTAINS VOTES FOR IT TO BE CENTER FOR EACH PIXEL
	Graph centers(w, h);
	centers.fill_graph(); // CONTAINS ULT. EDGES

	int tempX, tempY, row, col;
	double temp;
	threshold = pow(threshold, 2);
	threshold2 = pow(threshold2, 2);
	for (int i = 0; i < (int)v.size(); i++) {
		tempX = 0, tempY = 0;
		row = i / w;
		col = i % w;
		centers.illuminate(row, col, org[i][0], org[i][1], org[i][2]);
		if ((row < h - 1) and (row >= 1) and (col < w - 1) and (col >= 1)) {
			for (int j = -1; j <= 1; j++) {
				for (int k = -1; k <= 1; k++) {
					tempX += xKernel[j + 1][k + 1] * v[i + (w * k) + (j)][0];
					tempY += yKernel[j + 1][k + 1] * v[i + (w * k) + (j)][0];
				}
			}
			if (pow(tempX, 2) + pow(tempY, 2) > threshold2) {
				if (pow(tempX, 2) + pow(tempY, 2) > threshold) {
					g.illuminate(row, col, 3, 3, 3);
				}
				else {
					g.illuminate(row, col, 2, 2, 2);
				}
			}
			temp = atan2(tempY, tempX);
			angle.illuminate(row, col, temp, temp, temp);
			temp = pow(tempX, 2) + pow(tempY, 2);
			mag.illuminate(row, col, temp, temp, temp);
		}
	}

	width = w;
	height = h;
	for (int i = 1; i < h - 1; i++) {
		for (int j = 1; j < w - 1; j++) {
			if (g.get_coord(i, j) == 3) {
				area_fill(i, j, g);
			}
		}
	}

	for (int i = 1; i < h - 1; i++) {
		for (int j = 1; j < w - 1; j++) {
			if (g.get_coord(i, j) == 0 or g.get_coord(i, j) == 2) {
				g.illuminate(i, j, 0, 0, 0);
			}
			else {

				temp = angle.get_coord(i, j) * 180 / M_PI;
				if (temp < 0) {
					temp += 180;
				}
				if (temp != 0) {
					temp = round_to_multiple(temp);
				}
				if (temp == 45) {
					if (mag.get_coord(i, j) < mag.get_coord(i + 1, j - 1) or mag.get_coord(i, j) < mag.get_coord(i - 1, j + 1)) {
						g.illuminate(i, j, 0, 0, 0);
					}
				}
				else if (temp == 135) {
					if (mag.get_coord(i, j) < mag.get_coord(i - 1, j - 1) or mag.get_coord(i, j) < mag.get_coord(i + 1, j + 1)) {
						g.illuminate(i, j, 0, 0, 0);
					}
				}
				else if (temp == 90) {
					if (mag.get_coord(i, j) < mag.get_coord(i - 1, j) or mag.get_coord(i, j) < mag.get_coord(i + 1, j)) {
						g.illuminate(i, j, 0, 0, 0);
					}
				}
				else {
					if (mag.get_coord(i, j) < mag.get_coord(i, j - 1) or mag.get_coord(i, j) < mag.get_coord(i, j + 1)) {
						g.illuminate(i, j, 0, 0, 0);
					}
				}
			}
		}
	}

	double dX, dY;
	for (int i = 1; i < h - 1; i++) {
		for (int j = 1; j < w - 1; j++) {
			if (g.get_coord(i, j) == 1) {
				dX = cos(angle.get_coord(i, j));
				dY = sin(angle.get_coord(i, j));
				votes.create_line(i - w * dX, j - w * dY, i + w * dX, j + w * dY, -1, -1, -1);
			}
		}
	}

	int max = 0;
	for (int i = 1; i < h - 1; i++) {
		for (int j = 1; j < w - 1; j++) {
			if (votes.get_coord(i, j) > max) {
				max = votes.get_coord(i, j);
			}
		}
	}
	votes.create_file("imagev.ppm", max);

	for (int i = 0; i < 30; i++) {
		int row = i * int(h / 30);
		int lastRow = min((i + 1) * int(h / 30), h);
		for (int j = 0; j < 30; j++) {
			int col = j * int(w / 30);
			int lastCol = min((j + 1) * int(w / 30), w);
			double density = 0.0;
			for (int r = row; r < lastRow; r++) {
				for (int c = col; c < lastCol; c++) {
					density += votes.get_coord(r, c);
				}
			}
			density /= ((lastRow - row + 1) * (lastCol - col + 1));
			if (density < threshold3 / 2) {
				density = 0.75 * threshold3;
			}
			else {
				density = threshold3;
			}
			for (int r = row; r < lastRow; r++) {
				for (int c = col; c < lastCol; c++) {
					if (votes.get_coord(r, c) > density) {
						centers.create_circle(4, r, c, 255, 0, 0);
						centers.create_circle(3, r, c, 255, 0, 0);
						centers.create_circle(2, r, c, 255, 0, 0);
						centers.create_circle(1, r, c, 255, 0, 0);
					}
				}
			}
		}
	}
	centers.create_file("imageCC.ppm", 255);
}

void canny_edge()
{
	vector<vector<int>> inp = read_ppm_file(name);
	create_greyscale(inp);
	vector<vector<int>> gscale = read_ppm_file("imageg.ppm");
	sobel_conv(gscale, inp);
}

int main(int argc, char* argv[])
{
	if (argc > 1) {
		for (int i = 1; i < argc - 1;) {
			if (strcmp(argv[i], "-H") == 0) {
				threshold = stoi(argv[i + 1]);
			}
			else if (strcmp(argv[i], "-L") == 0) {
				threshold2 = stoi(argv[i + 1]);
			}
			else if (strcmp(argv[i], "-F") == 0) {
				name = argv[i + 1];
			}
			else {
				threshold3 = stoi(argv[i + 1]);
			}
			i += 2;
		}
	}
	canny_edge();
}