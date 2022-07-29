// vector.h

// Vector structure and its operations for 3-d space vectors.

// Written by: Yanting Wang           12/25/2002

#ifndef __VECTOR__H
#define __VECTOR__H

#include <iostream>
#include <math.h>

struct Vector {
    double x, y, z;

    Vector() : x(0.0), y(0.0), z(0.0) {
    }
    Vector(double x, double y, double z) : x(x), y(y), z(z) {
    }

    void set(double xx, double yy, double zz) {
        x = xx;
        y = yy;
        z = zz;
    }

    // square length of a vector
    double r2() const {
        return x * x + y * y + z * z;
    }

    // normalize a vector
    Vector normalize() const {
        double r = sqrt(r2());
        return Vector(x / r, y / r, z / r);
    }

    const Vector &operator=(const Vector &v) {
        x = v.x;
        y = v.y;
        z = v.z;
        return v;
    }
};

inline Vector operator+(const Vector &v1, const Vector &v2) {
    return Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline Vector operator-(const Vector &v1, const Vector &v2) {
    return Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline Vector operator-(const Vector &v1) {
    return Vector(-v1.x, -v1.y, -v1.z);
}

// scalar product
inline double operator*(const Vector &v1, const Vector &v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// vector product
inline Vector operator^(const Vector &v1, const Vector &v2) {
    return Vector(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
                  v1.x * v2.y - v1.y * v2.x);
}

inline Vector operator*(const Vector &v, const double d) {
    return Vector(v.x * d, v.y * d, v.z * d);
}

inline Vector operator/(const Vector &v, const double d) {
    return Vector(v.x / d, v.y / d, v.z / d);
}

inline std::ostream &operator<<(std::ostream &output, const Vector &p) {
    output << p.x << " " << p.y << " " << p.z;
    return output;
}

inline std::istream &operator>>(std::istream &input, Vector &p) {
    input >> p.x >> p.y >> p.z;
    return input;
}

#endif //__VECTOR__H