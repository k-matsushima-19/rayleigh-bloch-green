#include "shape.hpp"

namespace Geometry{
    Eigen::Vector2d Circle::x(double t) const{
        // 0<= t < 2*pi

        Eigen::Vector2d out(R*std::cos(t), R*std::sin(t));

        return out;
    };

    Eigen::Vector2d Circle::dx(double t) const{
        // 0<= t < 2*pi
        Eigen::Vector2d out(-R*std::sin(t), R*std::cos(t));

        return out;
    };

    Eigen::Vector2d Circle::ddx(double t) const{
        // 0<= t < 2*pi

        Eigen::Vector2d out(-R*std::cos(t), -R*std::sin(t));

        return out;
    };

    Eigen::Vector2d Circle::dddx(double t) const{
        // 0<= t < 2*pi
        Eigen::Vector2d out(+R*std::sin(t), -R*std::cos(t));

        return out;
    };

    Eigen::Vector2d Ellipse::x(double t) const{
        // 0<= t < 2*pi

        Eigen::Vector2d out(
            R1*std::cos(angle)*std::cos(t) - R2*std::sin(angle)*std::sin(t), 
            R1*std::sin(angle)*std::cos(t) + R2*std::cos(angle)*std::sin(t)
        );

        return out;
    };

    Eigen::Vector2d Ellipse::dx(double t) const{
        // 0<= t < 2*pi

        Eigen::Vector2d out(
            -R1*std::cos(angle)*std::sin(t) - R2*std::sin(angle)*std::cos(t), 
            -R1*std::sin(angle)*std::sin(t) + R2*std::cos(angle)*std::cos(t)
        );

        return out;
    };

    Eigen::Vector2d Ellipse::ddx(double t) const{
        // 0<= t < 2*pi

        Eigen::Vector2d out(
            -R1*std::cos(angle)*std::cos(t) + R2*std::sin(angle)*std::sin(t), 
            -R1*std::sin(angle)*std::cos(t) - R2*std::cos(angle)*std::sin(t)
        );

        return out;
    };

    Eigen::Vector2d Ellipse::dddx(double t) const{
        // 0<= t < 2*pi

        Eigen::Vector2d out(
            +R1*std::cos(angle)*std::sin(t) + R2*std::sin(angle)*std::cos(t), 
            +R1*std::sin(angle)*std::sin(t) - R2*std::cos(angle)*std::cos(t)
        );

        return out;
    };
}