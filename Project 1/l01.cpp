// Ramya Reddy
// Period 4
// 08/26/2021

#include<iostream>
#include <fstream>
#include<math.h>
#include <vector>
using namespace std;

const int size = 800;
int(*arr)[size] = new int[size][size];

// helper method that instantiates graphs with white pixels
void fill_graph()
{
    for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			arr[i][j] = 1;
		}
	}
}

// helper method to generate 6 different values for the three pairs of coordinates for the triangle
vector<double> generate_random_points()
{
    vector<double> points;
    for (int i = 0; i < 6; i++) {
		points.push_back(((double)rand() / (RAND_MAX)));
	}
    return points;
}

// helper method that scales and rounds given points
vector<int> scale(vector<double> points)
{
	vector<int> scaled_points;
	for (int i = 0; i < points.size(); i++)	{
		scaled_points.push_back(int(points[i] * 800));
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

// modified Bresenham's algorithm that draws line for all twelve cases
void create_line(int x1, int y1, int x2, int y2)
{
	int deltaX = x2 - x1;
	int deltaY = y2 - y1;
	int error = deltaY - deltaX;
	if (abs(deltaX) > abs(deltaY)) { // if true driving axis is x-axis else driving axis is y-axis
		int j = y1;
		if (deltaX < 0) { // if true move right else move left
			for (int i = x1; i >= x2; i--) {
				illuminate(i, j);
				if (error >= 0) {
					if (deltaY < 0) { // if true move down else move up
						j -= 1;
					}
					else if (deltaY > 0) {
						j += 1;
					}
					error -= abs(deltaX);
				}
				error += abs(deltaY);
			}
		}
		else {
			for (int i = x1; i <= x2; i++) {
				illuminate(i, j);
				if (error >= 0) {
					if (deltaY < 0) { // if true move down else move up
						j -= 1;
					}
					else if (deltaY > 0) {
						j += 1;
					}
					error -= abs(deltaX);
				}
				error += abs(deltaY);
			}
		}
	}
	else {
		int i = x1;
		if (deltaY < 0) {
			for (int j = y1; j >= y2; j--) { // if true move down else move up
				illuminate(i, j);
				if (error >= 0) {
					if (deltaX < 0) { // if true move left else move right
						i -= 1;
					}
					else if (deltaX > 0) {
						i += 1;
					}
					error -= abs(deltaY);
				}
				error += abs(deltaX);
			}
		}
		else {
			for (int j = y1; j <= y2; j++) {
				illuminate(i, j);
				if (error >= 0) {
					if (deltaX < 0) { // if true move left else move right
						i -= 1;
					}
					else if (deltaX > 0) {
						i += 1;
					}
					error -= abs(deltaY);
				}
				error += abs(deltaX);
			}
		}
	}
}

// rasterized algorithm that draws circle
void create_circle(int r, int x0, int y0) {
	int x, y, xmax, y2, y2_new, ty;
	xmax = (int)(r * 0.70710678);
	y = r;
	y2 = y * y;
	ty = (2 * y) - 1; y2_new = y2;
	for (x = 0; x <= xmax + 1; x++) {
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

// draws triangle based on given vertices
void create_triangle(vector<int> points)
{
	create_line(points[0], points[1], points[2], points[3]);
	create_line(points[2], points[3], points[4], points[5]);
	create_line(points[0], points[1], points[4], points[5]);
}

// finds the coefficients of a line in form Ax+By=C that passes through two given points
vector<double> find_equation_of_line(vector<double> points)
{
	vector<double> coefficients;
	coefficients.push_back(points[3] - points[1]);
	coefficients.push_back(points[0] - points[2]);
	coefficients.push_back(coefficients[0] * points[0] + coefficients[1] * points[1]);
	return coefficients;
}

// finds circumcenter of given triangle
pair<double, double> find_circumcenter(vector<double> points)
{
	vector<double> eq1 = { points[0], points[1], points[2], points[3] };
	vector<double> eq2 = { points[2], points[3], points[4], points[5] };
	vector<double> coeffs1 = find_equation_of_line(eq1);  // currently in the form Ax+By=C
	vector<double> coeffs2 = find_equation_of_line(eq2);

	double c0 = (-coeffs1[1] * ((eq1[0] + eq1[2]) / 2) + coeffs1[0] * ((eq1[1] + eq1[3]) / 2)); // calculates C for the perpendicular bisector of each lines
	double c1 = (-coeffs2[1] * ((eq2[0] + eq2[2]) / 2) + coeffs2[0] * ((eq2[1] + eq2[3]) / 2));

	double det = coeffs1[0] * coeffs2[1] - coeffs1[1] * coeffs2[0]; // finds intersection of perpendicular bisectors of the lines
	double x = (coeffs2[0] * c0 - coeffs1[0] * c1) / det;
	double y = (coeffs2[1] * c0 - coeffs1[1] * c1) / det;
	return make_pair(x,y);
}

// draws the circumcircle, euler line, 9 point circle and incircle
void draw_figures(vector<double> points)
{
	double x0, y0;
	double a = sqrt((pow(points[0] - points[2], 2) + pow(points[1] - points[3], 2)));
	double b = sqrt((pow(points[2] - points[4], 2) + pow(points[3] - points[5], 2)));
	double c = sqrt((pow(points[0] - points[4], 2) + pow(points[1] - points[5], 2)));
	double s = (a + b + c) / 2;

	// creates incircle
	double R = sqrt((s - a) * (s - b) * (s - c) / s);
	x0 = ((points[4] * a) + (points[0] * b) + (points[2] * c)) / (s * 2); // calculates incenter by taking weighted average of coordinates
	y0 = ((points[5] * a) + (points[1] * b) + (points[3] * c)) / (s * 2);
    auto temp = scale({R, x0, y0});
	create_circle(temp[0], temp[1], temp[2]);

	// creates circumcircle
	R = (a * b * c) / (4 * R * s);
	auto pair = find_circumcenter(points);
	x0 = pair.first;
	y0 = pair.second;
    temp = scale({R, x0, y0});
	create_circle(temp[0], temp[1], temp[2]);

	// creates euler_line based on line made by centroid and circumcenter
	vector<double> pts = { (points[0] + points[2] + points[4]) / 3, (points[1] + points[3] + points[5]) / 3, x0, y0 };
	auto eq = find_equation_of_line(pts);
    temp = scale({eq});
	create_line(int( size * temp[2] / temp[0] ), 0, int( (size * eq[2] - (size * eq[1])) / eq[0] ), size);  
    
	// creates 9 point circle as the circumcircle made connecting the midpoints of the original triangle
	R /= 2;
	vector<double> ninePtCircle{ (points[0] + points[2]) / 2, (points[1] + points[3]) / 2, (points[2] + points[4]) / 2, (points[3] + points[5]) / 2, (points[0] + points[4]) / 2, (points[1] + points[5]) / 2 };
	pair = find_circumcenter(ninePtCircle);
	x0 = pair.first;
	y0 = pair.second;
    temp = scale({R, x0, y0});
	create_circle(temp[0], temp[1], temp[2]);
}

// creates ppm file with origin on the bottom left
void create_file()
{
	ofstream myfile;
	myfile.open("triangle.ppm", ios::out);
	myfile << "P3 " << size << " " << size << " 1\n";
	for (int i = size - 1; i >= 0; i--) {
		for (int j = 0; j < size; j++) 	{
			myfile << arr[j][i] << " " << arr[j][i] << " " << arr[j][i] << " ";
		}
		myfile << endl;
	}
	myfile.close();
}

int main()
{
	srand(time(0));
	fill_graph();
	vector<double> points = generate_random_points();
	create_triangle(scale(points));
	draw_figures(points);
	create_file();
	return 0;
}