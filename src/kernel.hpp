#pragma once
#include "mycomplex.hpp"
#include <set>
#include <cmath>
#include "shape.hpp"
//#include "complex_bessel.h"
#include "quasiPeriodicGreen.hpp"






template <typename T> // codomain of the kernel
class Kernel{
    //-----------
    // Variables
    //-----------
public:
    Geometry::Shape* shape;

    //----------------------
    // Pure virtual methods
    //----------------------
public:
    virtual T X(double t, double tau) const = 0;
    // Singular part of X
    virtual T X1(double t, double tau) const = 0;
    // Diagonal of X1
    virtual T X1_diag(double t) const = 0;
    // Diagonal of X2
    virtual T X2_diag(double t) const = 0;
    // Diagonal function
    virtual T a(double t) const = 0;

    //---------
    // Methods
    //---------
public:
    // Regular part of X
    T X2(double t, double tau) const{
        return X(t,tau) - X1(t,tau)*std::log(4*std::pow(std::sin((t-tau)/2),2));
    }
};

// Kernel N_s
class N : public Kernel<double>{
    //-----------
    // Variables
    //-----------
public:
    // frequency (wavenumber)
    double k;
    // Kummer parameter q > (Im[beta]/(k*L))^2
    double q;
    double sqrt_q; // = sqrt(q)

    //--------------
    // Constructors
    //--------------
    N(Geometry::Shape* shape, double k, double sqrt_q) : k(k), sqrt_q(sqrt_q){
        this->shape = shape; 
        this->q = sqrt_q*sqrt_q;
    };

    //---------------
    // Deconstructor
    //---------------
public:
    ~N(){
        shape = nullptr;
    }

    //---------
    // Methods
    //---------
public:
    /**
    @return std::complex<double> function value   \f$(x_1,y_1)\f$
    **/
    double X(double t, double tau) const;
    double X1(double t, double tau) const;
    double X1_diag(double t) const;
    double X2_diag(double t) const;

    double a(double t) const{return 0;};
    
};


// Kernel K_beta
class K : public Kernel<std::complex<double>>{
    //-----------
    // Variables
    //-----------
public:
    // frequency (wavenumber)
    double k;
    // Spacing
    double L;
    // Floquet wavenumber beta
    std::complex<double> beta;
    // Sheet
    std::set<int> I;
    // coupling parameter
    std::complex<double> eta;

    // Quasi-periodic Green
    QuasiPeriodicGreen* green;

    //--------------
    // Constructors
    //--------------
public:
    K(Geometry::Shape* shape, double k, double L, std::complex<double> beta, std::set<int> I, std::complex<double> eta, int M) : k(k), L(L), beta(beta), I(I), eta(eta){
        this->shape = shape;
        green = new QuasiPeriodicGreen(k, L, beta, I, M);
    }

    //------------
    // Destructor
    //------------
public:
    ~K(){
        delete green;
    }

    //---------
    // Methods
    //---------
    std::complex<double> X(double t, double tau) const;
    std::complex<double> X1(double t, double tau) const;
    std::complex<double> X1_diag(double t) const;
    std::complex<double> X2_diag(double t) const;

    std::complex<double> a(double t) const{return 0;};

};

// -N + 2*dd(G^R)
class Brakhage : public Kernel<std::complex<double>>{
    //-----------
    // Variables
    //-----------
public:
    // frequency (wavenumber)
    double k;
    // Spacing
    double L;
    // Floquet wavenumber beta
    std::complex<double> beta;
    // Sheet
    std::set<int> I;
    // coupling parameter
    std::complex<double> eta;

    QuasiPeriodicGreen* green;

    //--------------
    // Constructors
    //--------------
public:
    Brakhage(Geometry::Shape* shape, double k, double L, std::complex<double> beta, std::set<int> I, std::complex<double> eta, int M);
    Brakhage(Geometry::Shape* shape, double k, double L, std::complex<double> beta, std::set<int> I, std::complex<double> eta, int M, double beta_imag_max);

    //------------
    // Destructor
    //------------
public:
    ~Brakhage(){
        delete green;
    }

    //---------
    // Methods
    //---------
    std::complex<double> X(double t, double tau) const;
    std::complex<double> X1(double t, double tau) const;
    std::complex<double> X1_diag(double t) const;
    std::complex<double> X2_diag(double t) const;

    std::complex<double> a(double t) const{
        const static std::complex<double> ione(0, 1);
        return ione*eta*shape->ds(t);
    };

};

template <typename T> // codomain of the kernel
class KernelPotential{
    //-----------
    // Variables
    //-----------
public:
    Geometry::Shape* shape;

    //----------------------
    // Pure virtual methods
    //----------------------
public:
    virtual T X(Eigen::Vector2d x, double tau) const = 0;
};

class BrakhagePotential : public KernelPotential<std::complex<double>>{
public:
    //-----------
    // Variables
    //-----------
    // frequency (wavenumber)
    double k;
    // Spacing
    double L;
    // Floquet wavenumber beta
    std::complex<double> beta;
    // Sheet
    std::set<int> I;
    // coupling parameter
    std::complex<double> eta;

    QuasiPeriodicGreen* green;

    //--------------
    // Constructors
    //--------------
public:
    BrakhagePotential(Geometry::Shape* shape, double k, double L, std::complex<double> beta, std::set<int> I, std::complex<double> eta, int M);

    std::complex<double> X(Eigen::Vector2d x, double tau) const;
};