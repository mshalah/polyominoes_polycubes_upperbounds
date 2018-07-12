//
//  utils.h
//  mira_twigs
//
//  Created by Mira Shalah on 9/26/17.
//  Copyright (c) 2017 Mira Shalah. All rights reserved.
//

#ifndef mira_twigs_utils_h
#define mira_twigs_utils_h
#include <math.h>
#include "matrix.h"

// Basic operations on 3D points

#define PI 3.14159265

using namespace std;

struct Point3 {
    int x,y,z;
public:
    int X() const { return x;}
    int Y() const { return y;}
    int Z() const { return z;}
    Point3() : x(0), y(0), z(0) {};
    Point3(int _x, int _y, int _z) : x(_x), y(_y), z(_z) {};
    Point3(const Point3& p) : x(p.X()), y(p.Y()), z(p.Z()) {};
    bool operator==(const Point3& p) const {
        return x == p.X() && y == p.Y() && z == p.Z();
    }
};

Point3 operator+( const Point3& p1, const Point3& p2 )
{
    Point3 res(p1.X() + p2.X(), p1.Y() + p2.Y(), p1.Z() + p2.Z());
    return res;
};

std::ostream& operator<<(std::ostream& os, Point3& p)
{
    os << "(" << p.X() << ", " << p.Y() << ", " << p.Z() << ")";
    return os;
}

Point3 operator*(const Matrix<int>& m, const Point3& p)
{
    return Point3(m.m11_*p.X() + m.m12_*p.Y() + m.m13_*p.Z(),
                  m.m21_*p.X() + m.m22_*p.Y() + m.m23_*p.Z(),
                  m.m31_*p.X() + m.m32_*p.Y() + m.m33_*p.Z());
};


Point3 rotateX(Point3 p, int theta) {
    Point3 newp(p.X(),
                  round(p.Y() * cos (theta * PI / 180.0 ) + p.Z() * sin (theta * PI / 180.0) ),
                  round(p.Z() * cos (theta * PI / 180.0 ) - p.Y() * sin (theta * PI / 180.0) ));
    return newp;
}

Point3 rotateY(Point3 p, int theta) {
    return Point3(round(p.X() * cos (theta * PI / 180.0 ) - p.Z() * sin (theta * PI / 180.0) ),
                  p.Y(),
                  round(p.X() * sin (theta * PI / 180.0 ) + p.Z() * cos (theta * PI / 180.0) ));
}

Point3 rotateZ(Point3 p, int theta) {
    return Point3(round(p.X() * cos (theta * PI / 180.0 ) + p.Y() * sin (theta * PI / 180.0) ),
                  round(p.Y() * cos (theta * PI / 180.0 ) - p.X() * sin (theta * PI / 180.0) ),
                  p.Z());
}

Point3 reflectXY(Point3 p) {
    return Point3(p.X(), p.Y(), -p.Z());
}

Point3 reflectXZ(Point3 p) {
    return Point3(p.X(), -p.Y(), p.Z());
}

Point3 reflectYZ(Point3 p) {
    return Point3(-p.X(), p.Y(), p.Z());
}


#endif
