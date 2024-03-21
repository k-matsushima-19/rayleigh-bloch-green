#include "kernel.hpp"

//-----------
// Constants
//-----------
const static double EULER = 0.577215664901532861;
const static std::complex<double> ione(0, 1);
const static double PI = M_PI;
static const Eigen::Vector2d e1 = Eigen::Vector2d(1.0,0.0);

//------------------
// Bessel functions
//------------------
static double bessel_jn(int n, double x){
    return std::cyl_bessel_j(n, x);
}

// static std::complex<double> bessel_jn(int n, std::complex<double> x){
//     return sp_bessel::besselJ(n, x);
// }

static double bessel_yn(int n, double x){
    return std::cyl_neumann(n, x);
}

// static std::complex<double> bessel_yn(int n, std::complex<double> x){
//     return sp_bessel::besselY(n, x);
// }

static std::complex<double> bessel_hn1(int n, double x){
    return std::complex<double>(bessel_jn(n, x), bessel_yn(n, x));
}

// static std::complex<double> bessel_hn1(int n, std::complex<double> x){
//     return sp_bessel::hankelH1(n, x);
// }

static double bessel_kn(int n, double x){
    return std::cyl_bessel_k(n, x);
}

static double bessel_in(int n, double x){
    return std::cyl_bessel_i(n, x);
}

//---
// N
//---
double N::X(double t, double tau) const{

    double r = shape->r(t,tau);
    double zeta = ( shape->dx(t).dot(shape->x(t)-shape->x(tau)) ) * ( shape->dx(tau).dot(shape->x(t)-shape->x(tau)) ) / std::pow(r,2);

    return 
        -sqrt_q*k*0.5 * (
            2*sqrt_q*k/PI * bessel_kn(0,sqrt_q*k*r) * zeta 
            -2/PI * bessel_kn(1,sqrt_q*k*r)/r * ( shape->dx(t).dot(shape->dx(tau)) -2*zeta )
        )
        + 1 / ( 4*PI*std::pow(std::sin((t-tau)/2),2) )
    ;
}

double N::X1(double t, double tau) const{
    double r = shape->r(t,tau);
    double zeta = ( shape->dx(t).dot(shape->x(t)-shape->x(tau)) ) * ( shape->dx(tau).dot(shape->x(t)-shape->x(tau)) ) / std::pow(r,2);

    return
        sqrt_q*k/(2*PI) * (
            sqrt_q*k*bessel_in(0,sqrt_q*k*r) * zeta
            + bessel_in(1,sqrt_q*k*r)/r * (shape->dx(t).dot(shape->dx(tau)) -2*zeta)
        )
    ;
}

double N::X1_diag(double t) const{
    return q*std::pow(k*shape->ds(t),2)/(4*PI);
}

double N::X2_diag(double t) const{
    double ds = shape->ds(t);

    return 
        (1+2*EULER+2*std::log(sqrt_q*k*ds*0.5)) * q*std::pow(k*ds,2)/(4*PI)
        + 1/(12*PI) 
        + std::pow(shape->dx(t).dot(shape->ddx(t)),2) / (2*PI*std::pow(ds,4))
        - shape->ddx(t).dot(shape->ddx(t)) / (4*PI*std::pow(ds,2))
        - shape->dx(t).dot(shape->dddx(t)) / (6*PI*std::pow(ds,2))
    ;
}

//---
// K
//---
std::complex<double> K::X(double t, double tau) const{
    // double r = (shape->x(t)-shape->x(tau)).norm();

    // std::complex<double> term1 = 0;
    // {
    //     for(int s=0; s<green->M; s++){
    //         term1 += green->cs[s] * green->qs[s] * bessel_kn(0,green->sqrt_qs[s]*k*r);
    //     }

    //     term1 *= k*k/M_PI * shape->normal(t).dot(shape->normal(tau)) * shape->ds(t)*shape->ds(tau);
    // }

    // std::complex<double> term2 = 0;
    // {
    //     for(int s=0; s<green->M; s++){
    //         term2 += green->cs[s] * green->sqrt_qs[s] * bessel_kn(1,green->sqrt_qs[s]*k*r);
    //     }

    //     term2 *= -ione*k*eta/M_PI * shape->normal(t).dot(shape->x(t)-shape->x(tau))/r * shape->ds(t)*shape->ds(tau);
    // }

    // std::complex<double> term3 = 2 * green->laplace_regular_part(shape->x(t),shape->x(tau)) * shape->normal(t).dot(shape->normal(tau)) * shape->ds(t)*shape->ds(tau);

    // std::complex<double> term4 = -2*ione*eta * green->grad_x_regular_part(shape->x(t),shape->x(tau),shape->normal(t)) * shape->ds(t)*shape->ds(tau);

    // return term1 + term2 + term3 + term4;
    return 0;
     
}

std::complex<double> K::X1(double t, double tau) const{
    double r = (shape->x(t)-shape->x(tau)).norm();

    std::complex<double> term1 = 0;
    {
        for(int s=0; s<green->M; s++){
            term1 += green->cs[s] * green->qs[s] * bessel_in(0,green->sqrt_qs[s]*k*r);
        }

        term1 *= -k*k/(2*M_PI) * shape->normal(t).dot(shape->normal(tau)) * shape->ds(t)*shape->ds(tau);
    }

    std::complex<double> term2 = 0;
    {
        for(int s=0; s<green->M; s++){
            term2 += green->cs[s] * green->sqrt_qs[s] * bessel_in(1,green->sqrt_qs[s]*k*r);
        }

        term2 *= -ione*k*eta/(2*M_PI) * shape->normal(t).dot(shape->x(t)-shape->x(tau))/r * shape->ds(t)*shape->ds(tau);
    }

    return term1 + term2;

}

std::complex<double> K::X1_diag(double t) const{

    double sum = 0;
    for(int s=0; s<green->M; s++){
        sum += green->cs[s] * green->qs[s];
    }

    return -std::pow(k*shape->ds(t),2)/(2*M_PI) * sum;
}

std::complex<double> K::X2_diag(double t) const{

    double term1 = 0;
    {
        for(int s=0; s<green->M; s++){
            term1 += green->cs[s] * green->qs[s] * (EULER + std::log(green->sqrt_qs[s]*k*shape->ds(t)/2) );
        }

        term1 *= -std::pow(k*shape->ds(t),2)/M_PI;
    }

    std::complex<double> term2 = 0;
    {
        for(int s=0; s<green->M; s++){
            term2 += green->cs[s];
        }

        term2 *= ione*eta/(2*M_PI)*shape->normal(t).dot(shape->ddx(t));
    }

    std::complex<double> term3 = 0; // 2*green->laplace_regular_part(shape->x(t),shape->x(t)) * std::pow(shape->ds(t),2);

    std::complex<double> term4 = -2*ione*eta*green->grad_x_regular_part(shape->x(t),shape->x(t),shape->normal(t)) * std::pow(shape->ds(t),2);

    return term1 + term2 + term3 + term4;
    
}

//----------
// Brakhage
//----------
Brakhage::Brakhage(Geometry::Shape* shape, double k, double L, std::complex<double> beta, std::set<int> I, std::complex<double> eta, int M) : k(k), L(L), beta(beta), I(I), eta(eta){
    this->shape = shape;
    
    green = new QuasiPeriodicGreen(k, L, beta, I, M);
    
};

Brakhage::Brakhage(Geometry::Shape* shape, double k, double L, std::complex<double> beta, std::set<int> I, std::complex<double> eta, int M, double beta_imag_max) : k(k), L(L), beta(beta), I(I), eta(eta){
    this->shape = shape;
    
    green = new QuasiPeriodicGreen(k, L, beta, I, M, beta_imag_max);
    
};

std::complex<double> Brakhage::X(double t, double tau) const{

    double r = shape->r(t,tau);
    double zeta = ( shape->dx(t).dot(shape->x(t)-shape->x(tau)) ) * ( shape->dx(tau).dot(shape->x(t)-shape->x(tau)) ) / std::pow(r,2);

    //-----------
    // Compute N
    //-----------
    double N = 0;
    for(int s=0; s<green->M; s++){
        N += green->cs[s] *
        (
                -green->sqrt_qs[s]*k*0.5 * (
                2*green->sqrt_qs[s]*k/PI * bessel_kn(0,green->sqrt_qs[s]*k*r) * zeta 
                -2/PI * bessel_kn(1,green->sqrt_qs[s]*k*r)/r * ( shape->dx(t).dot(shape->dx(tau)) -2*zeta )
            )
            + 1 / ( 4*PI*std::pow(std::sin((t-tau)/2),2) )
        );
    }

    //-----------------
    // Compute dd(G^R)
    //-----------------
    std::complex<double> dd_GR = green->grad_xy_regular_part(shape->x(t), shape->x(tau), shape->dx(t), shape->dx(tau));

    //-----------
    // Compute Q
    //-----------
    std::complex<double> Q = 2*std::pow(k,2)*(
        green->regular_part(shape->x(t), shape->x(tau)) + green->singular_part(shape->x(t), shape->x(tau))
    ) * shape->normal(t).dot(shape->normal(tau)) * shape->ds(t)*shape->ds(tau);

    //-----------
    // Compute Y
    //-----------
    std::complex<double> Y = -2*ione*eta * (
          green->grad_x_regular_part(shape->x(t), shape->x(tau), shape->normal(t))
        + green->grad_x_singular_part(shape->x(t), shape->x(tau), shape->normal(t)) 
    ) * shape->ds(t)*shape->ds(tau);

    return N - 2*dd_GR + Q + Y;

}

std::complex<double> Brakhage::X1(double t, double tau) const{
    double r = shape->r(t,tau);
    double zeta = ( shape->dx(t).dot(shape->x(t)-shape->x(tau)) ) * ( shape->dx(tau).dot(shape->x(t)-shape->x(tau)) ) / std::pow(r,2);

    //---
    // N
    //---
    double N = 0;
    for(int s=0; s<green->M; s++){
        double c = green->cs[s];
        double sqrt_q = green->sqrt_qs[s];

        N += c*sqrt_q * (
            sqrt_q*k*bessel_in(0,sqrt_q*k*r) * zeta
            + bessel_in(1,sqrt_q*k*r)/r * (shape->dx(t).dot(shape->dx(tau)) -2*zeta)
        );
    }
    N *= k/(2*M_PI);

    //---
    // Q
    //---
    std::complex<double> Q = 0;
    {
        for(int s=0; s<green->M; s++){
            double c = green->cs[s];
            double sqrt_q = green->sqrt_qs[s];

            Q += c*bessel_in(0,sqrt_q*k*r);
        }

        Q *= std::pow(k,2)/(2*M_PI) * shape->normal(t).dot(shape->normal(tau)) * shape->ds(t)*shape->ds(tau);
    }

    //---
    // Y
    //---
    std::complex<double> Y = 0;
    {
        for(int s=0; s<green->M; s++){
            double c = green->cs[s];
            double sqrt_q = green->sqrt_qs[s];

            Y += c*sqrt_q*bessel_in(1,sqrt_q*k*r);
        }

        Y *= -ione*k*eta/(2*M_PI) * shape->normal(t).dot(shape->x(t)-shape->x(tau))/r  * shape->ds(t)*shape->ds(tau);

    }

    return N + Q + Y;
}

std::complex<double> Brakhage::X1_diag(double t) const{
    return -std::pow(k*shape->ds(t),2)/(4*PI);
}

std::complex<double> Brakhage::X2_diag(double t) const{
    double ds = shape->ds(t);

    //---
    // N
    //---
    double N = 
        + (1+2*EULER)*std::pow(k*ds,2)/(4*PI) 
        - 1/(12*PI)
        - std::pow(shape->dx(t).dot(shape->ddx(t)),2) / (2*PI*std::pow(ds,4))
        + shape->ddx(t).dot(shape->ddx(t)) / (4*PI*std::pow(ds,2))
        + shape->dx(t).dot(shape->dddx(t)) / (6*PI*std::pow(ds,2))
    ;
    double sum = 0;
    for(int s=0; s<green->M; s++){
        double c = green->cs[s];
        double sqrt_q = green->sqrt_qs[s];
        double q = green->qs[s];

        sum += c*q*std::log(sqrt_q*k*ds/2);
            
    }
    N += +sum * std::pow(k*ds,2)/(2*PI);

    //-----------------
    // Compute dd(G^R)
    //-----------------
    std::complex<double> dd_GR = green->grad_xy_regular_part(shape->x(t), shape->x(t), shape->dx(t), shape->dx(t));

    //---
    // Q
    //---
    double sum_log = 0;
    for(int s=0; s<green->M; s++){
        double c = green->cs[s];
        double sqrt_q = green->sqrt_qs[s];
        double q = green->qs[s];

        sum_log += c*std::log(sqrt_q*k*ds/2);
            
    }

    std::complex<double> Q = +std::pow(k*ds,2) * (
        2*green->regular_part(shape->x(t), shape->x(t)) - EULER/M_PI + sum_log/M_PI
    );

    //---
    // Y
    //---
    std::complex<double> Y = 
        -ione*eta/(2*M_PI) * shape->normal(t).dot(shape->ddx(t))
        -2*ione*eta * std::pow(ds,2) * green->grad_x_regular_part(shape->x(t), shape->x(t), shape->normal(t))
    ;

    return N - 2*dd_GR + Q + Y;
}

//----------
// Brakhage
//----------
BrakhagePotential::BrakhagePotential(Geometry::Shape* shape, double k, double L, std::complex<double> beta, std::set<int> I, std::complex<double> eta, int M) : k(k), L(L), beta(beta), I(I), eta(eta){
    this->shape = shape;
    
    green = new QuasiPeriodicGreen(k, L, beta, I, M);
    
};

// (-nu(tau)*grad_x_G(x,x(tau)) - i*eta*G(x,x(tau))) * ds(tau)
std::complex<double> BrakhagePotential::X(Eigen::Vector2d x, double tau) const{
    return (
        - (
               green->grad_x_regular_part(x, shape->x(tau), shape->normal(tau))
            + green->grad_x_singular_part(x, shape->x(tau), shape->normal(tau))
          )
        - ione*eta*
          (
            green->calc(x, shape->x(tau))
          )
    ) * shape->ds(tau);
}