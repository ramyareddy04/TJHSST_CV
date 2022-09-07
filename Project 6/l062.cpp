// Ramya Reddy
// Period 4
// 2/14/2022

#include<iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>
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
        return (x==p.x and y==p.y);
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
            } else {
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
  double fraction_of_circle(int r, int x0, int y0, Graph other) {
		int x, y, xmax, y2, y2_new, ty;
		xmax = (int)(r / sqrt(2));
		y = r;
		y2 = y * y;
		ty = (2 * y) - 1;
		y2_new = y2;
    double temp = 0, total = 0;
		for (x = 0; x < xmax + 2; x++) {
			if ((y2 - y2_new) >= ty) {
				y2 -= ty;
				y -= 1;
				ty -= 2;
			}
      if (other.get_coord(x + x0, y+y0)==1)
        temp += 1;
      if (other.get_coord(x + x0, -y+y0)==1)
        temp += 1;
      if (other.get_coord(-x + x0, y+y0)==1)
        temp += 1;
      if (other.get_coord(-x + x0, -y+y0)==1)
        temp += 1;
      if (other.get_coord(y + y0, x + x0)==1)
        temp += 1;
      if (other.get_coord(-y + y0, x + x0)==1)
        temp += 1;
      if (other.get_coord(y + y0, -x + x0)==1)
        temp += 1;
      if (other.get_coord(-y + y0, -x + x0)==1)
        temp += 1;
			y2_new -= (2 * x) - 3;
      total += 8;
		}
    return temp/total;
	}
	// creates ppm file with origin on the bottom left
	void create_file(string filename, int max)
	{
		ofstream myfile;
		myfile.open(filename, ios::out);
		myfile << "P3 \n" << width << " " << height << "\n" << max<<"\n";
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
                
            } else if (numLine == 2) {
                while (i <int(line.size())) {
                    if (line[i] == ' ') {
                        w = stod(line.substr(0, i));
                        h = stod(line.substr(i+1, line.length()));
                        i = int(line.size());
                    }
                    i += 1;
                }
            } else if (numLine == 3) {
                max = stod(line);
                values.push_back({w, h, max});
            } else {
                while (i<int(line.size())) {
                    if (line[i] == ' ') {
                        temp.push_back(stod(line.substr(j, i)));
                        j = i+1;
                    }
                    i += 1;
                }
            }
            numLine += 1;
		}
	}
    
    for (int k=0; k<(int)temp.size();) {
        values.push_back({temp[k], temp[k+1], temp[k+2]});
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
    for (int i=0; i<(int)v.size(); i++) {
        row = i/w;
        col = i%w;
        temp = (v[i][0]+v[i][1]+v[i][2])/3;
        g.illuminate(row, col, temp, temp, temp);
    }
    g.create_file("imageg.ppm", max);
}

int width, height;
// recusive area fill method
void area_fill(int row, int col, Graph& img)
{
    if (row >= 0 and row <height and col >=0 and col < width) {
        if (img.get_coord(row, col)==0 or img.get_coord(row, col)==1) {
            return;
        }
        img.illuminate(row, col, 1, 1, 1);
        
        area_fill(row-1, col, img);
        area_fill(row+1, col, img);
        area_fill(row, col-1, img);
        area_fill(row, col+1, img);
        area_fill(row-1, col-1, img);
        area_fill(row-1, col+1, img);
        area_fill(row+1, col-1, img);
        area_fill(row+1, col+1, img);
    }
}

int threshold = 80, threshold2 = 150, threshold3 = 60, s=30, minRadius = 70, maxRadius = 175, minSeparation = 65;
double radiusThreshold = 0;

string name = "image.ppm";

int round_to_multiple(double inp)
{
    int temp = inp/45;
    if (inp-temp*45 <= 22.5) {
        return temp*45;
    } else {
        return (temp+1)*45;
    }
}

void create_results(string filename, int count[5])
{
    ofstream myfile;
    myfile.open(filename, ios::out);
    double sum = 0.0;
    string type[5] = {"Pennies", "Nickels", "Dimes", "Quarters", "Silver Dollars"};
    double val[5] = {0.01, 0.05, 0.10, 0.25, 0.50};
    for (int i=0; i<5; i++) {
        myfile << count[i] << " " << type[i] << endl;
        cout << count[i] << " " << type[i] << endl;
        sum += count[i]*val[i];
    }
    myfile << "Total: $" << sum << endl;
    cout << "Total: $" << sum << endl;
    myfile.close();
}

// performs sobel operator on image, and returns image using two thresholds and direction of gradient, and performs hough transform
void sobel_conv(vector<vector<int>> v, vector<vector<int>> org)
{
    vector<vector<int>> xKernel = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}}, yKernel = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
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
    Graph radius(w, h);
	radius.fill_graph(); // CONTAINS OPTIMAL RADIUS FOR EACH CENTER
    Graph finRadius(w, h);
	finRadius.fill_graph(); // CONTAINS FINAL RADIUS OF EACH CENTER
    Graph percentage(w, h);
	percentage.fill_graph(); // CONTAINS CORRESPONDING PERCENTAGE OF OPTIMAL RADIUS
    Graph centers(w, h);
	centers.fill_graph(); // CONTAINS FINAL IMG W/ CIRCLES

    // DOUBLE THRESHOLDING
    int tempX, tempY, row, col;
    double temp;
    threshold = pow(threshold, 2);
    threshold2 = pow(threshold2, 2);
    for (int i=0; i<(int)v.size(); i++) {
        tempX = 0, tempY = 0;
        row = i/w;
        col = i%w;
        centers.illuminate(row, col, org[i][0], org[i][1], org[i][2]);
        if ((row<h-1) and (row>=1) and (col<w-1) and (col>=1)) {
            for (int j=-1; j<=1; j++) {
                for (int k=-1; k<=1; k++) {
                    tempX += xKernel[j+1][k+1]*v[i+(w*k)+(j)][0];
                    tempY += yKernel[j+1][k+1]*v[i+(w*k)+(j)][0];
                }
            }
            if (pow(tempX, 2) + pow(tempY, 2) > threshold2) {
                if (pow(tempX, 2) + pow(tempY, 2) > threshold) {
                    g.illuminate(row, col, 3, 3, 3);
                } else {
                    g.illuminate(row, col, 2, 2, 2);
                }
            }
            temp = atan2(tempX, tempY);
            angle.illuminate(row, col, temp, temp, temp);
            temp = pow(tempX, 2)+pow(tempY, 2);
            mag.illuminate(row, col, temp, temp, temp);
        } 
    }
    
    // RECURSIVE AREA FILL
    width  = w;
    height = h;
    for (int i = 1; i < h-1; i++) {
        for (int j=1; j < w-1; j++) {
            if (g.get_coord(i, j)==3) {
                area_fill(i, j, g);
            }
        }
    }
    
    // APPLYING NMS
    for (int i=0; i<h; i++) {
        for (int j=0; j<w; j++) {
            if (g.get_coord(i,j)==0 or g.get_coord(i, j)==2) {
                if (g.get_coord(i,j)==2)
                    g.illuminate(i, j, 0, 0, 0);
            }  else {
                temp = angle.get_coord(i, j)*180/M_PI;
                if (temp < 0) {
                    temp += 180;
                }            
                if (temp!=0) {
                    temp = round_to_multiple(temp);
                }
                if (temp==45) {
                    if (mag.get_coord(i, j) < mag.get_coord(i-1, j-1) or mag.get_coord(i, j) < mag.get_coord(i+1, j+1)) {
                        g.illuminate(i, j, 0, 0, 0);
                    } 
                } else if (temp==135) {
                    if (mag.get_coord(i, j) < mag.get_coord(i+1, j-1) or mag.get_coord(i, j) < mag.get_coord(i-1, j+1)) {
                        g.illuminate(i, j, 0, 0, 0);
                    } 
                } else if (temp==90) {
                    if (mag.get_coord(i, j) < mag.get_coord(i-1, j) or mag.get_coord(i, j) < mag.get_coord(i+1, j)) {
                        g.illuminate(i, j, 0, 0, 0);
                    } 
                } else {
                    if (mag.get_coord(i, j) < mag.get_coord(i, j-1) or mag.get_coord(i, j) < mag.get_coord(i, j+1)) {
                        g.illuminate(i, j, 0, 0, 0);
                    } 
                }
            }
        }
    }

    // CREATING VOTING SYSTEM
    double dX, dY;
    int end, begin;
    for (int i=1; i<h-1; i++) {
        for (int j=1; j<w-1; j++) {
            if (g.get_coord(i, j)==1) {
                dX = 1/cos(angle.get_coord(i,j));
                dY = 1/sin(angle.get_coord(i,j));
                end = min(i/dX, j/dY);
                begin = min((h-i)/dX, (w-j)/dY);
                votes.create_line(i - end*dX, j - end*dY, i + begin*dX, j +begin*dY, -1, -1, -1);
            }
        }
    }
    
    int maxV=0;
    for (int i=1; i<h-1; i++) {
        for (int j=1; j<w-1; j++) {
            if (votes.get_coord(i, j) > maxV) {
                maxV = votes.get_coord(i, j);
            }
        }
    }
    votes.create_file("imagev.ppm", maxV);
    
    // DETERMINING POTENTIAL CENTERS AND THEIR OPTIMAL RADIUS
    double tempR, bestR;
    int bad = 0;
    for (int i=0; i<s; i++) {
        int row = i*int(h/s);
        int lastRow = min((i+1)*int(h/s), h);
        for (int j=0; j<s; j++) {
            int col = j*int(w/s);
            int lastCol = min((j+1)*int(w/s), w);
            double density = 0.0;
            for (int r=row; r<lastRow; r++) {
                for (int c=col; c<lastCol; c++) {
                    density += votes.get_coord(r, c);
                }
            }
            density /= ((lastRow-row+1)*(lastCol-col+1));
            density *= 1.2;
            for (int r=row; r<lastRow; r++) {
                for (int c=col; c<lastCol; c++) {
                    int avg = (centers.get_coord(r,c,0)+centers.get_coord(r,c,1)+centers.get_coord(r,c,2))/3;
                    bool isEdge = false;
                    for (int dx=-1; dx<=1;) {
                        for (int dy=-1; dy<=1;) {
                            if (r+dx >= 0 and r+dx <h and c+dy >=0 and c+dy < w) {
                                if (g.get_coord(r+dx, c+dy)==1) {
                                    isEdge = true;
                                    break;
                                }
                            }
                            dy += 2;
                        }
                        dx += 2;
                    }
                    if (isEdge == false and votes.get_coord(r, c) > density and abs(avg-200)>38 and abs(avg-170)>30) { 
                        double posR[maxRadius-minRadius+1];
                        for (int idx=0; idx<=maxRadius-minRadius; idx++) {
                            posR[idx] = 0.0;
                        } 
                        for (int radius=minRadius; radius<=maxRadius; radius++) {
                            if (radius <= 120 or radius >= 165) {
                                int halfRad = radius*0.6;
                                if ((r-halfRad >= 0 and abs(centers.get_coord(r-halfRad, c)-200) <= 40) or (r+halfRad < h and abs(centers.get_coord(r+halfRad, c)-200) <= 40)) {
                                    bad += 1;
                                } else {
                                    posR[radius-minRadius] = centers.fraction_of_circle(radius, r, c, g);
                                }
                            }
                        }
                        bool isCenter = false;
                        tempR = radiusThreshold/100;
                        for (int idx=0; idx<=maxRadius-minRadius;) {
                            if (posR[idx] > tempR) {
                                if (idx+minRadius >= 165) {
                                    if (posR[idx] > 0.2 and centers.get_coord(r,c,0) < 1.25*centers.get_coord(r,c,2)) {
                                        isCenter = true;
                                        tempR = posR[idx];
                                        bestR = idx;
                                    }
                                } else if (idx+minRadius >= 103) {
                                    if (centers.get_coord(r,c,0) < 1.25*centers.get_coord(r,c,2) and ((maxRadius > 120 and posR[idx] > 0.2) or maxRadius <= 120)) {
                                        isCenter = true;
                                        tempR = posR[idx];
                                        bestR = idx;
                                    }
                                } else {
                                    isCenter = true;
                                    tempR = posR[idx];
                                    bestR = idx;
                                }
                            }
                            idx += 1;
                        }
                        if (isCenter) {
                            radius.illuminate(r,c,bestR+minRadius,0,0);
                            finRadius.illuminate(r,c,bestR+minRadius,0,0);
                            percentage.illuminate(r,c,tempR,0,0);
                        }
                    }
                }
            }
        }
    }
    
    // ELIMINATING REDUNDANT CIRCLES
    int l = minSeparation;
    for (int i=0; i<h; i++) {
        for (int j=0; j<w; j++) {
            if (radius.get_coord(i, j)!=0) {
                tempR = percentage.get_coord(i,j)*radius.get_coord(i,j);
                int r0=i, c0=j;
                for (int r=max(0, i-l); r<min(h, i+l); r++) {
                    for (int c=max(0, j-l); c<min(w, j+l); c++) {
                        if (percentage.get_coord(r,c)*radius.get_coord(r,c) > tempR) {
                            tempR = percentage.get_coord(r,c)*radius.get_coord(r,c);
                            r0 = r;
                            c0 = c;
                        }
                    }
                }
                for (int r=max(0, i-l); r<min(h, i+l); r++) {
                    for (int c=max(0, j-l); c<min(w, j+l); c++) {
                        if (percentage.get_coord(r,c)*radius.get_coord(r,c) < tempR or (percentage.get_coord(r,c)*radius.get_coord(r,c) == tempR and (r!=r0 or c!=c0))) {
                            percentage.illuminate(r,c,0,0,0);
                            finRadius.illuminate(r,c,0,0,0);
                            //radius.illuminate(r,c,0,0,0);
                        }
                    }
                }
            }                        
        }
    }

    // TRANSPOSING CIRCLES AND DETERMINING COIN TYPE
    int count[5] = {0, 0, 0, 0, 0};
    int type;
    for (int r=0; r<h; r++) {
        for (int c=0; c<w; c++) {
            if (finRadius.get_coord(r,c)!=0) {
                bestR = finRadius.get_coord(r,c);
                int avgR = centers.get_coord(r, c, 0), avgB = centers.get_coord(r, c, 2), total = 1;
                if ((r - (bestR / 2)) >= 0) {
                    avgR += centers.get_coord(r - (bestR / 2), c, 0);
                    avgB += centers.get_coord(r - (bestR / 2), c, 2);
                    total += 1;
                }
                if ((r + (bestR / 2)) >= 0) {
                    avgR += centers.get_coord(r + (bestR / 2), c, 0);
                    avgB += centers.get_coord(r + (bestR / 2), c, 2);
                    total += 1;
                }
                avgR /= total;
                avgB /= total;
                if (bestR < 85) {
                    if (avgR>1.2*avgB) {
                        type = 0;
                    } else {
                        type = 2;
                    }
                } else if (bestR < 103) {
                    if (avgR>1.3*avgB) {
                        type = 0;
                    } else {
                        type = 1;
                    }
                } else if (bestR < 150) {
                    type = 3;
                } else {
                    type = 4;
                }
                count[type] += 1;
                for (int i=0; i<=3; i++) {
                    if (type==0) {
                        centers.create_circle(i, r, c, 255, 0, 0);
                        centers.create_circle(bestR+i, r, c, 255, 0, 0);
                    } else if (type == 1) {
                        centers.create_circle(i, r, c, 128, 0, 128);
                        centers.create_circle(bestR+i, r, c, 128, 0, 128);
                    } else if (type == 2) {
                        centers.create_circle(i, r, c, 0, 0, 255);
                        centers.create_circle(bestR+i, r, c, 0, 0, 255);
                    } else if (type == 3) {
                        centers.create_circle(i, r, c, 0, 255, 0);
                        centers.create_circle(bestR+i, r, c, 0, 255, 0);
                    } else {
                        centers.create_circle(i, r, c, 255, 255, 0);
                        centers.create_circle(bestR+i, r, c, 255, 255, 0);
                    }
                }
            }                        
        }
    }
    centers.create_file("coins.ppm", 255);
    create_results("results.txt", count);
}

void coin_detection()
{
    vector<vector<int>> inp = read_ppm_file(name);
    create_greyscale(inp);
    vector<vector<int>> gscale = read_ppm_file("imageg.ppm");
    sobel_conv(gscale, inp);
}

int main(int argc,  char* argv[])
{
    if (argc > 1) {
        for (int i=1; i<argc-1;) {
            if (strcmp(argv[i],"-H") == 0) {
                threshold = stoi(argv[i+1]);
            } else if (strcmp(argv[i],"-L") == 0) {
                threshold2 = stoi(argv[i+1]);
            } else if (strcmp(argv[i],"-F") == 0) {
                name = argv[i+1];
            } else if (strcmp(argv[i],"-S") == 0) {
                s = stoi(argv[i+1]);
            } else if (strcmp(argv[i],"-TC") == 0) {
                radiusThreshold = stod(argv[i+1]);
            } else if (strcmp(argv[i],"-MinR") == 0) {
                minRadius = stoi(argv[i+1]);
            } else if (strcmp(argv[i],"-MaxR") == 0) {
                maxRadius = stoi(argv[i+1]);
            } else if (strcmp(argv[i],"-MS") == 0) {
                minSeparation = stoi(argv[i+1]);
            } else {
                threshold3 = stoi(argv[i+1]);
            }
            i+=2;
        }
    }
    coin_detection();
}