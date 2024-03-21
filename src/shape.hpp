#pragma once
#include <cmath>
#include <eigen3/Eigen/Core>

namespace Geometry{

    class Shape{
        //----------------------
        // Pure virtual methods
        //----------------------
    public:
        // parametric representation of a shape
        virtual Eigen::Vector2d x(double t) const = 0;
        // Derivative of x
        virtual Eigen::Vector2d dx(double t) const = 0;
        // Derivative of dx
        virtual Eigen::Vector2d ddx(double t) const = 0;
        // Derivative of ddx
        virtual Eigen::Vector2d dddx(double t) const = 0;

        //---------
        // Methods
        //---------
    public:
        // Distance func.
        double r(double t, double tau) const{
            return (x(t) - x(tau)).norm();
        }

        // surface element at t=t
        double ds(double t) const{
            return dx(t).norm();
        }

        // Unit outward normal
        Eigen::Vector2d normal(double t){
            Eigen::Vector2d a = dx(t);
            Eigen::Vector2d out = {
                a[1],-a[0]
            };

            return out / a.norm();
        };

    };

    class Circle : public Shape{
        //-----------
        // Variables
        //-----------
    public:
        // Radius
        double R;

        //--------------
        // Constructors
        //--------------
    public:
        Circle(double R) : R(R){};

        //---------
        // Methods
        //---------
    public:
        Eigen::Vector2d x(double t) const;
        Eigen::Vector2d dx(double t) const;
        Eigen::Vector2d ddx(double t) const;
        Eigen::Vector2d dddx(double t) const;

    };

    class Ellipse : public Shape{
        //-----------
        // Variables
        //-----------
    public:
        // center
        Eigen::Vector2d center;
        // length of the 1st semi-axis
        double R1;
        // length of the 2nd semi-axis
        double R2;
        // angle of the 1st semi-axis (+x if angle=0)
        double angle;

        Ellipse(const Eigen::Vector2d& center, const double& R1, const double& R2, const double& angle) : center(center), R1(R1), R2(R2), angle(angle){};

        Eigen::Vector2d x(double t) const;
        Eigen::Vector2d dx(double t) const;
        Eigen::Vector2d ddx(double t) const;
        Eigen::Vector2d dddx(double t) const;
    };

}