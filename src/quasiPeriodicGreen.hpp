#pragma once
#include "mycomplex.hpp"
#include <set>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <iostream>



//---------------------------------
// Quasi-periodic Green's function
//---------------------------------
class QuasiPeriodicGreen{
    //-----------
    // Variables
    //-----------
public:
    // freq.
    double k;
    // Spacing
    double L;
    // Floquet wavenumber
    std::complex<double> beta;
    // Index set
    std::set<int> I;

    //-----------------------
    // Kummer transformation
    //-----------------------
    double* qs;
    double* cs;
    double* sqrt_qs;
    // Number of Kummer terms (approx. less than 8 for numerical stability)
    const int M = 6;

    //--------
    // Consts
    //--------
    // (qmin + |Im[beta]|/(k*L))^2 <= q <= (qmax + |Im[beta]|/(k*L))^2
    
    // constexpr static double qmax = M*M*0.5; //16.0;
    constexpr static int NMAX = 10000; //100000;
    

    //--------------
    // Constructors
    //--------------
    QuasiPeriodicGreen(double k, double L, std::complex<double> beta, std::set<int> I, int M);
    QuasiPeriodicGreen(double k, double L, std::complex<double> beta, std::set<int> I, int M, double beta_imag_max);
    
    //------------
    // Destructor
    //------------
    ~QuasiPeriodicGreen(){
        delete [] qs;
        delete [] cs;
        delete [] sqrt_qs;
    }

    //---------
    // Methods
    //---------
    std::complex<double> regular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y) const;
    std::complex<double> singular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y) const;
    std::complex<double> calc(const Eigen::Vector2d& x, const Eigen::Vector2d& y) const{return regular_part(x,y)+singular_part(x,y);};
    // vec dot grad(x) of the regular part
    std::complex<double> grad_x_regular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y, const Eigen::Vector2d& vec) const;
    std::complex<double> grad_x_singular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y, const Eigen::Vector2d& vec) const;
    std::complex<double> laplace_regular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y) const;
    // grad_x dot grad_y of the regular part
    std::complex<double> grad_xy_regular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y, const Eigen::Vector2d& a, const Eigen::Vector2d& b) const;

    // New
    std::complex<double> grad_x(const Eigen::Vector2d& x, const Eigen::Vector2d& y, const Eigen::Vector2d& vec) const;
    std::complex<double> grad_xy(const Eigen::Vector2d& x, const Eigen::Vector2d& y, const Eigen::Vector2d& a, const Eigen::Vector2d& b) const;

    // F
    std::complex<double> F_conv(const int& m, const std::complex<double> beta, const double& alpha) const;
    std::complex<double> F(const int& m, const std::complex<double> beta, const double& alpha) const;
    std::complex<double> F_term(const int& m, const std::complex<double> beta, const double& alpha, const int& s) const;
    std::complex<double> F_tilde(const int& m, const std::complex<double> beta, const double& alpha) const;
    // std::complex<double> F_hat(const int& m, const std::complex<double> beta, const double& alpha) const;
    // std::complex<double> F_term_tmp(const int& m, const std::complex<double> beta, const double& alpha, const int& s) const;

    // H
    std::complex<double> H(const int& m, const std::complex<double> beta, const double& alpha) const;
    std::complex<double> H_term(const int& m, const std::complex<double> beta, const double& alpha, const int& s) const;
    std::complex<double> H_tilde(const int& m, const std::complex<double> beta, const double& alpha) const;

    // P
    std::complex<double> P(const int& m, const std::complex<double> beta, const double& alpha) const;
    std::complex<double> P_term(const int& m, const std::complex<double> beta, const double& alpha, const int& s) const;
    std::complex<double> P_tilde(const int& m, const std::complex<double> beta, const double& alpha) const;

    //------------------
    // Static functions
    //------------------
    static double bessel_jn(int n, double x){
        return std::cyl_bessel_j(n, x);
    }

    static double bessel_yn(int n, double x){
        return std::cyl_neumann(n, x);
    }

    static std::complex<double> bessel_hn1(int n, double x){
        return std::complex<double>(bessel_jn(n, x), bessel_yn(n, x));
    }

    static double bessel_kn(int n, double x){
        // std::cout << x << std::endl;
        return std::cyl_bessel_k(n, x);
    }

    static double bessel_in(int n, double x){
        // std::cout << x << std::endl;
        return std::cyl_bessel_i(n, x);
    }

};