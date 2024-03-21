#include <vector>
#include <fstream>
#include "quasiPeriodicGreen.hpp"

static const Eigen::Vector2d e1 = Eigen::Vector2d(1.0,0.0);
static const std::complex<double> ione = std::complex<double>(0, 1);

static const bool PRINT_ITER = false;

static const double THRES = 1e-14;

// Perhaps need some special treats for large |z| to avoid overflow ...
static std::complex<double> sigma(const std::complex<double>& z){
    if(std::abs(z) < 5.0) return std::exp(-z) * std::sinh(z);
    else                  return 0.5*(1 - std::exp(-2*z));
}

QuasiPeriodicGreen::QuasiPeriodicGreen(double k, double L, std::complex<double> beta, std::set<int> I, int M) : k(k), L(L), beta(beta), I(I), M(M){
    qs = new double[M];
    sqrt_qs = new double[M];
    cs = new double[M];

    // Set q
    // double _qmin = std::pow((qmin + abs(beta.imag()))/(k*L),2);
    // double _qmax = std::pow((qmax + abs(beta.imag()))/(k*L),2);
    double qmin = 0.1;
    double qmax = 1.0;//M*M*0.5; //16.0;

    double _qmin = std::pow((qmin + abs(beta.imag())/(k*L)),2);
    double _qmax = std::pow((qmax + abs(beta.imag())/(k*L)),2);
    for(int i=0; i<M; i++){
        // !!! qsをbetaの関数にするのはだめ
        qs[i] = _qmin + i*(_qmax-_qmin)/(M-1);
        // qs[i] = std::pow((std::pow(1.1,i) + abs(beta.imag())/(k*L)),2);

        sqrt_qs[i] = std::sqrt(qs[i]);
    }

    // Set c
    for(int s=0; s<M; s++){
        cs[s] = -1.0;
        for(int i=0; i<M; i++){
            if(s != i) cs[s] *= (qs[i]+1)/(qs[i]-qs[s]);
        }

        // std::cout << cs[s] << std::endl;
    }
};

// For fixed parameters q and c
QuasiPeriodicGreen::QuasiPeriodicGreen(double k, double L, std::complex<double> beta, std::set<int> I, int M, double beta_imag_max) : k(k), L(L), beta(beta), I(I), M(M){
    qs = new double[M];
    sqrt_qs = new double[M];
    cs = new double[M];

    // Set q
    // double _qmin = std::pow((qmin + abs(beta.imag()))/(k*L),2);
    // double _qmax = std::pow((qmax + abs(beta.imag()))/(k*L),2);
    double qmin = 0.1;
    double qmax = 1.0;//M*M*0.5; //16.0;

    double _qmin = std::pow((qmin + abs(beta_imag_max)/(k*L)),2);
    double _qmax = std::pow((qmax + abs(beta_imag_max)/(k*L)),2);
    for(int i=0; i<M; i++){
        // !!! qsをbetaの関数にするのはだめ
        qs[i] = _qmin + i*(_qmax-_qmin)/(M-1);
        // qs[i] = std::pow((std::pow(1.1,i) + abs(beta.imag())/(k*L)),2);

        sqrt_qs[i] = std::sqrt(qs[i]);
    }

    // Set c
    for(int s=0; s<M; s++){
        cs[s] = -1.0;
        for(int i=0; i<M; i++){
            if(s != i) cs[s] *= (qs[i]+1)/(qs[i]-qs[s]);
        }

        // std::cout << cs[s] << std::endl;
    }
};



// std::complex<double> QuasiPeriodicGreen::regular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y) const{

//     //-----------------------
//     // Sum in physical space
//     //-----------------------
//     std::complex<double> sum_phys = 0.0;
//     for(int s=0; s<M; s++){
//         // from -infty to -1, and 1 to +infty
//         for(int sign=0; sign<2; sign++){
//             std::complex<double> summand = 0.0;
//             for(int _n=1; _n<NMAX; _n++){
//                 int n = _n * (2*sign-1); // +1,+2,... if sign=1 , -1,-2,... if sign=0

//                 double r = (x-y-n*L*e1).norm();

//                 std::complex<double> tmp = bessel_kn(0,sqrt_qs[s]*k*r) * std::exp(ione*n*beta);

//                 summand += tmp;

//                 // std::cout << s << " " << n << " " << abs(tmp/summand) << std::endl;

//                 bool flag = abs(tmp) < THRES * abs(summand) || (abs(tmp) == 0);

//                 if(flag){
//                     // std::cout << "# converged at " << _n << std::endl;
//                     break;
//                 }

//             }
//             sum_phys += -cs[s]/(2*M_PI) * summand;
//         }
//     }

//     std::complex<double>* exps = new std::complex<double>[M+1];

//     //------------------
//     // Reciprocal space
//     //------------------
//     std::complex<double> sum_reciprocal = 0;
//     for(int _m=0; _m<NMAX; _m++){
//         // std::complex<double> summand = 0.0;
//         std::complex<double> prev = sum_reciprocal;

//         for(int sign=0; sign<2; sign++){
//             // +0,+2,... if sign=1 , -1,-2,... if sign=0
//             int m = _m * (2*sign-1); 
//             if(m == 0 && sign == 1) continue;

//             std::complex<double> km = (beta+2*m*M_PI)/L;
//             std::complex<double> exp1 = std::exp(ione*km*(x[0]-y[0]));

//             std::complex<double> tmp = 0;

//             // compute exp (x2-y2) for each s
//             for(int s=0; s<M; s++){
//                 std::complex<double> k2 = std::sqrt(km*km + qs[s]*k*k);
//                 exps[s] = cs[s] * std::exp(-k2 * abs(x[1]-y[1])) / (2*L*k2);

//                 tmp += exps[s];
//             }

//             // if m is in I
//             std::complex<double> k2;
//             if (I.count(m)) {
//                 // std::cout << m << std::endl;
//                 k2 = -std::sqrt(km*km - k*k);
//             }
//             else{
//                 k2 = std::sqrt(km*km - k*k);
//             }
//             exps[M] = std::exp(-k2 * abs(x[1]-y[1])) / (2*L*k2);
//             tmp += exps[M];

//             // Finally add exps term by term (to avoid loss of significant digits)
//             for(int s=0; s<M+1; s++){
//                 sum_reciprocal += exp1 * exps[s];
//             }

//         }

//         // std::cout << _m << " " << std::abs(sum_reciprocal - prev) / std::abs(sum_reciprocal) << std::endl;

//         if(std::abs(sum_reciprocal - prev) < 1e-14 * std::abs(sum_reciprocal)){
//             if(PRINT_ITER) std::cout << "# converged at " << _m << std::endl;
//             break;
//         }

//     }

//     delete [] exps;

//     return sum_phys + sum_reciprocal;

// }


std::complex<double> QuasiPeriodicGreen::regular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y) const{

    //-----------------------
    // Sum in physical space
    //-----------------------
    std::complex<double> sum_phys = 0.0;
    for(int s=0; s<M; s++){
        // from -infty to -1, and 1 to +infty
        for(int sign=0; sign<2; sign++){
            std::complex<double> summand = 0.0;
            for(int _n=1; _n<NMAX; _n++){
                int n = _n * (2*sign-1); // +1,+2,... if sign=1 , -1,-2,... if sign=0

                double r = (x-y-n*L*e1).norm();

                if(-n*beta.imag() > 500) break;

                std::complex<double> tmp = bessel_kn(0,sqrt_qs[s]*k*r) * std::exp(ione*n*beta);

                summand += tmp;

                // std::cout << s << " " << n << " " << abs(tmp/summand) << std::endl;

                bool flag = abs(tmp) < THRES * abs(summand) || (abs(tmp) == 0);

                if(flag){
                    // std::cout << "# converged at " << _n << std::endl;
                    break;
                }

            }
            sum_phys += -cs[s]/(2*M_PI) * summand;
        }
    }

    std::complex<double>* exps = new std::complex<double>[M+1];

    //------------------
    // Reciprocal space
    //------------------
    std::complex<double> sum_reciprocal = 0;
    for(int _m=0; _m<NMAX; _m++){
        // std::complex<double> summand = 0.0;
        std::complex<double> prev = sum_reciprocal;

        for(int sign=0; sign<2; sign++){
            // +0,+2,... if sign=1 , -1,-2,... if sign=0
            int m = _m * (2*sign-1); 
            if(m == 0 && sign == 1) continue;

            std::complex<double> km = (beta+2*m*M_PI)/L;

            std::complex<double> f;
            // if m is in I
            if (I.count(m)) f = F_tilde(m, beta, std::abs(x[1]-y[1]));
            else            f =       F(m, beta, std::abs(x[1]-y[1]));

            sum_reciprocal += std::exp(ione*km*(x[0]-y[0])) * f / (2*L);

        }

        if(std::abs(sum_reciprocal - prev) < 1e-14 * std::abs(sum_reciprocal)){
            if(PRINT_ITER) std::cout << "# converged at " << _m << std::endl;
            break;
        }

    }

    delete [] exps;

    return sum_phys + sum_reciprocal;

}

std::complex<double> QuasiPeriodicGreen::singular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y) const{
    //-----------------------
    // Sum in physical space
    //-----------------------
    std::complex<double> sum_phys = 0.0;
    for(int s=0; s<M; s++){
        double r = (x-y).norm();

        std::complex<double> tmp = bessel_kn(0,sqrt_qs[s]*k*r);
        
        sum_phys += -cs[s]/(2*M_PI) * tmp;
        
    }

    return sum_phys;
}


// std::complex<double> QuasiPeriodicGreen::grad_x_regular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y, const Eigen::Vector2d& vec) const{

//     int sgn = 1;
//     if(x[1]-y[1] < 0) sgn = -1;

//     //-----------------------
//     // Sum in physical space
//     //-----------------------
//     std::complex<double> sum_phys = 0.0;
//     for(int s=0; s<M; s++){
//         // from -infty to -1, and 1 to +infty
//         for(int sign=0; sign<2; sign++){
//             std::complex<double> summand = 0.0;
//             for(int _n=1; _n<NMAX; _n++){
//                 int n = _n * (2*sign-1); // +1,+2,... if sign=1 , -1,-2,... if sign=0
                

//                 double r = (x-y-n*L*e1).norm();

//                 // Be careful: exp may overflow for large n!
//                 if(std::exp(-n*beta.imag()) == HUGE_VAL) break;

//                 std::complex<double> tmp = bessel_kn(1,sqrt_qs[s]*k*r) * std::exp(ione*n*beta) * vec.dot(x-y-n*L*e1)/r;

//                 summand += tmp;

//                 // std::cout << s << " " << n << " " << abs(tmp/summand) << std::endl;

//                 // std::cout << _n << " " << abs(tmp)/abs(summand) << " " << sqrt_qs[s]*k*r << " " << ione*n*beta << std::endl;
//                 bool flag = abs(tmp) < THRES * abs(summand) || (abs(tmp) == 0);

//                 // if(_n > 20){
//                 //     std::cout << _n  << " " << flag << " " << abs(k/(2*M_PI)*cs[s]*sqrt_qs[s]*summand + sum_phys) << " " << abs(summand) << " " << flag << std::endl;
//                 //     exit(2);
//                 // }
                

//                 // if(abs(tmp) < 1e-18 * abs(summand)) break;
//                 if(flag){
//                     // std::cout << "# converged at " << _n << std::endl;
//                     break;
//                 }

//             }
//             sum_phys += k/(2*M_PI)*cs[s]*sqrt_qs[s] * summand;
//         }
//     }

//     std::complex<double>* exps = new std::complex<double>[M+1];

//     //------------------
//     // Reciprocal space
//     //------------------
//     std::complex<double> sum_reciprocal = 0;
//     for(int _m=0; _m<NMAX; _m++){
//         // std::complex<double> summand = 0.0;
//         std::complex<double> prev = sum_reciprocal;

//         for(int sign=0; sign<2; sign++){
//             // +0,+2,... if sign=1 , -1,-2,... if sign=0
//             int m = _m * (2*sign-1); 
//             if(m == 0 && sign == 1) continue;

//             std::complex<double> km = (beta+2*m*M_PI)/L;
//             std::complex<double> exp1 = std::exp(ione*km*(x[0]-y[0]));

//             std::complex<double> tmp = 0;

            
//             // compute exp (x2-y2) for each s
//             for(int s=0; s<M; s++){
//                 std::complex<double> k2 = std::sqrt(km*km + qs[s]*k*k);
//                 exps[s] = cs[s] * (ione*km*vec[0] - sgn*vec[1]*k2) * std::exp(-k2 * abs(x[1]-y[1])) / (2*L*k2);

//                 tmp += exps[s];
//             }

//             // if m is in I
//             std::complex<double> k2 = std::sqrt(km*km - k*k);
//             if (I.count(m)) {
//                 exps[M] = -(ione*km*vec[0] + sgn*vec[1]*k2) * std::exp(+k2 * abs(x[1]-y[1])) / (2*L*k2);
//             }
//             else{
//                 exps[M] = +(ione*km*vec[0] - sgn*vec[1]*k2) * std::exp(-k2 * abs(x[1]-y[1])) / (2*L*k2);
//             }
            
//             tmp += exps[M];


//             // Finally add exps term by term (to avoid loss of significant digits)
//             for(int s=0; s<M+1; s++){
//                 sum_reciprocal += exp1 * exps[s];
//             }

//             // // Finally add exps (to avoid loss of significant digits)
//             // std::complex<double> prev = summand;
//             // std::complex<double> tmp2 = 0;
//             // for(int s=0; s<M+1; s++){
//             //     tmp2 += exp1 * exps[s];
//             //     summand += exp1 * exps[s];
//             // }

//         }

//         // sum_reciprocal += summand;

//         if(std::abs(sum_reciprocal - prev) < 1e-14 * std::abs(sum_reciprocal)){
//             if(PRINT_ITER) std::cout << "# converged at " << _m << std::endl;
//             break;
//         }

//         // if(abs(summand) < 1e-14*abs(sum_reciprocal)){
//         //     if(PRINT_ITER) std::cout << "# converged at " << _m << std::endl;

//         //     break;
//         // }

//     }

//     delete [] exps;

//     return sum_phys + sum_reciprocal;
// }

std::complex<double> QuasiPeriodicGreen::grad_x_regular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y, const Eigen::Vector2d& vec) const{

    int sgn = 1;
    if(x[1]-y[1] < 0) sgn = -1;

    //-----------------------
    // Sum in physical space
    //-----------------------
    std::complex<double> sum_phys = 0.0;
    for(int s=0; s<M; s++){
        // from -infty to -1, and 1 to +infty
        for(int sign=0; sign<2; sign++){
            std::complex<double> summand = 0.0;
            for(int _n=1; _n<NMAX; _n++){
                int n = _n * (2*sign-1); // +1,+2,... if sign=1 , -1,-2,... if sign=0
                

                double r = (x-y-n*L*e1).norm();

                // Be careful: exp may overflow for large n!
                // if(std::exp(-n*beta.imag()) == HUGE_VAL) break;
                if(-n*beta.imag() > 500) break;

                std::complex<double> tmp = bessel_kn(1,sqrt_qs[s]*k*r) * std::exp(ione*n*beta) * vec.dot(x-y-n*L*e1)/r;

                summand += tmp;

                // std::cout << s << " " << n << " " << abs(tmp/summand) << std::endl;

                // std::cout << _n << " " << abs(tmp)/abs(summand) << " " << sqrt_qs[s]*k*r << " " << ione*n*beta << std::endl;
                bool flag = abs(tmp) < THRES * abs(summand) || (abs(tmp) == 0);

                // if(_n > 20){
                //     std::cout << _n  << " " << flag << " " << abs(k/(2*M_PI)*cs[s]*sqrt_qs[s]*summand + sum_phys) << " " << abs(summand) << " " << flag << std::endl;
                //     exit(2);
                // }
                

                // if(abs(tmp) < 1e-18 * abs(summand)) break;
                if(flag){
                    // std::cout << "# converged at " << _n << std::endl;
                    break;
                }

            }
            sum_phys += k/(2*M_PI)*cs[s]*sqrt_qs[s] * summand;
        }
    }

    std::complex<double>* exps = new std::complex<double>[M+1];

    //------------------
    // Reciprocal space
    //------------------
    std::complex<double> sum_reciprocal = 0;
    for(int _m=0; _m<NMAX; _m++){
        // std::complex<double> summand = 0.0;
        std::complex<double> prev = sum_reciprocal;

        for(int sign=0; sign<2; sign++){
            // +0,+2,... if sign=1 , -1,-2,... if sign=0
            int m = _m * (2*sign-1); 
            if(m == 0 && sign == 1) continue;



            std::complex<double> km = (beta+2*m*M_PI)/L;
            std::complex<double> f       = F(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> f_tilde = F_tilde(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> h       = H(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> h_tilde = H_tilde(m, beta, std::abs(x[1]-y[1]));

            // if m is in I
            if (I.count(m)) sum_reciprocal += std::exp(ione*km*(x[0]-y[0])) / (2*L) * (
                ione*km*vec[0]*f_tilde - sgn*vec[1]*h_tilde
            );
            else            sum_reciprocal += std::exp(ione*km*(x[0]-y[0])) / (2*L) * (
                ione*km*vec[0]*f       - sgn*vec[1]*h
            );

        }

        // sum_reciprocal += summand;

        if(std::abs(sum_reciprocal - prev) < 1e-14 * std::abs(sum_reciprocal)){
            if(PRINT_ITER) std::cout << "# converged at " << _m << std::endl;
            break;
        }

    }

    delete [] exps;

    return sum_phys + sum_reciprocal;
}

std::complex<double> QuasiPeriodicGreen::grad_x_singular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y, const Eigen::Vector2d& vec) const{


    //-----------------------
    // Sum in physical space
    //-----------------------
    std::complex<double> sum_phys = 0.0;
    for(int s=0; s<M; s++){
        
        double r = (x-y).norm();

        std::complex<double> summand = bessel_kn(1,sqrt_qs[s]*k*r) * vec.dot(x-y)/r;
        
        sum_phys += k/(2*M_PI)*cs[s]*sqrt_qs[s] * summand;
        
    }

    return sum_phys;
}

std::complex<double> QuasiPeriodicGreen::grad_x(const Eigen::Vector2d& x, const Eigen::Vector2d& y, const Eigen::Vector2d& vec) const{

    int sgn = 1;
    if(x[1]-y[1] < 0) sgn = -1;

    //-----------------------
    // Sum in physical space
    //-----------------------
    std::complex<double> sum_phys = 0.0;
    for(int s=0; s<M; s++){
        // from -infty to -1, and 1 to +infty
        for(int sign=0; sign<2; sign++){
            std::complex<double> summand = 0.0;
            for(int _n=0; _n<NMAX; _n++){
                int n = _n * (2*sign-1); // +1,+2,... if sign=1 , -1,-2,... if sign=0
                
                if(n == 0 && sign == 1) continue;

                double r = (x-y-n*L*e1).norm();

                // Be careful: exp may overflow for large n!
                // if(std::exp(-n*beta.imag()) == HUGE_VAL) break;
                if(-n*beta.imag() > 500) break;

                std::complex<double> tmp = bessel_kn(1,sqrt_qs[s]*k*r) * std::exp(ione*n*beta) * vec.dot(x-y-n*L*e1)/r;

                summand += tmp;

                // std::cout << s << " " << n << " " << abs(tmp/summand) << std::endl;

                // std::cout << _n << " " << abs(tmp)/abs(summand) << " " << sqrt_qs[s]*k*r << " " << ione*n*beta << std::endl;
                bool flag = abs(tmp) < THRES * abs(summand) || (abs(tmp) == 0);

                // if(_n > 20){
                //     std::cout << _n  << " " << flag << " " << abs(k/(2*M_PI)*cs[s]*sqrt_qs[s]*summand + sum_phys) << " " << abs(summand) << " " << flag << std::endl;
                //     exit(2);
                // }
                

                // if(abs(tmp) < 1e-18 * abs(summand)) break;
                if(flag){
                    // std::cout << "# converged at " << _n << std::endl;
                    break;
                }

            }
            sum_phys += k/(2*M_PI)*cs[s]*sqrt_qs[s] * summand;
        }
    }

    std::complex<double>* exps = new std::complex<double>[M+1];

    //------------------
    // Reciprocal space
    //------------------
    std::complex<double> sum_reciprocal = 0;
    for(int _m=0; _m<NMAX; _m++){
        // std::complex<double> summand = 0.0;
        std::complex<double> prev = sum_reciprocal;

        for(int sign=0; sign<2; sign++){
            // +0,+2,... if sign=1 , -1,-2,... if sign=0
            int m = _m * (2*sign-1); 
            if(m == 0 && sign == 1) continue;



            std::complex<double> km = (beta+2*m*M_PI)/L;
            std::complex<double> f       = F(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> f_tilde = F_tilde(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> h       = H(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> h_tilde = H_tilde(m, beta, std::abs(x[1]-y[1]));

            // if m is in I
            if (I.count(m)) sum_reciprocal += std::exp(ione*km*(x[0]-y[0])) / (2*L) * (
                ione*km*vec[0]*f_tilde - sgn*vec[1]*h_tilde
            );
            else            sum_reciprocal += std::exp(ione*km*(x[0]-y[0])) / (2*L) * (
                ione*km*vec[0]*f       - sgn*vec[1]*h
            );

        }

        // sum_reciprocal += summand;

        if(std::abs(sum_reciprocal - prev) < 1e-14 * std::abs(sum_reciprocal)){
            if(PRINT_ITER) std::cout << "# converged at " << _m << std::endl;
            break;
        }

    }

    delete [] exps;

    return sum_phys + sum_reciprocal;
}

// std::complex<double> QuasiPeriodicGreen::laplace_regular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y) const{

//     //-----------------------
//     // Sum in physical space
//     //-----------------------
//     std::complex<double> sum_phys = 0.0;
//     // bool converged[2] = {false,false};
//     for(int s=0; s<M; s++){
//         // from -infty to -1, and 1 to +infty
//         for(int sign=0; sign<2; sign++){
//             std::complex<double> summand = 0.0;
//             for(int _n=1; _n<NMAX; _n++){
//                 int n = _n * (2*sign-1); // +1,+2,... if sign=1 , -1,-2,... if sign=0
                
                

//                 double r = (x-y-n*L*e1).norm();

//                 // Be careful: exp may overflow for large n!
//                 if(std::exp(-n*beta.imag()) == HUGE_VAL) break;

//                 std::complex<double> tmp = bessel_kn(0,sqrt_qs[s]*k*r) * std::exp(ione*n*beta);

//                 summand += tmp;

//                 // std::cout << s << " " << n << " " << abs(tmp/summand) << std::endl;
//                 // bool flag = abs(tmp) < 1e-18 * abs(summand) || (abs(tmp) == 0);
//                 bool flag = abs(tmp) < THRES * abs(summand) || (abs(tmp) == 0);

//                 // if(_n < 10) std::cout << _n  << " " << abs(tmp) << " " << abs(summand) << " " << flag << std::endl;

                
//                 if(flag){
//                     if(PRINT_ITER) std::cout << "# converged at " << _n << std::endl;
                    
//                     break;
//                 }

//             }
//             sum_phys += cs[s]*qs[s]*k*k/(2*M_PI) * summand;

//             // if(converged[0] && converged[1]) break;
//         }
//     }

//     //------------------
//     // Reciprocal space
//     //------------------
//     std::complex<double>* exps = new std::complex<double>[M+1];
//     std::complex<double> sum_reciprocal = 0;
//     for(int _m=0; _m<NMAX; _m++){
//         std::complex<double> summand = 0.0;
//         for(int sign=0; sign<2; sign++){
//             // +0,+2,... if sign=1 , -1,-2,... if sign=0
//             int m = _m * (2*sign-1); 
//             if(m == 0 && sign == 1) continue;

//             std::complex<double> km = (beta+2*m*M_PI)/L;
//             std::complex<double> exp1 = std::exp(ione*km*(x[0]-y[0]));

//             std::complex<double> tmp = 0;

//             // compute exp (x2-y2) for each s
//             for(int s=0; s<M; s++){
//                 std::complex<double> k2 = std::sqrt(km*km + qs[s]*k*k);
//                 exps[s] = -cs[s]*qs[s]*k*k * std::exp(-k2 * abs(x[1]-y[1])) / (2*L*k2);

//                 tmp += exps[s];
//             }

//             // if m is in I
//             std::complex<double> k2;
//             if (I.count(m)) {
//                 // std::cout << m << std::endl;
//                 k2 = -std::sqrt(km*km - k*k);
//             }
//             else{
//                 k2 = std::sqrt(km*km - k*k);
//             }
//             exps[M] = k*k*std::exp(-k2 * abs(x[1]-y[1])) / (2*L*k2);
//             tmp += exps[M];


//             // if(m > 1) std::cout << m << " " << tmp*exp1 / sum_reciprocal << std::endl;

//             // Finally add exps (to avoid loss of significant digits)
//             std::complex<double> prev = summand;
//             std::complex<double> tmp2 = 0;
//             for(int s=0; s<M+1; s++){
//                 tmp2 += exp1 * exps[s];
//                 summand += exp1 * exps[s];
//             }

//         }

//         sum_reciprocal += summand;

//         if(abs(summand) < 1e-14*abs(sum_reciprocal)){
//             if(PRINT_ITER) std::cout << "# converged at " << _m << std::endl;

//             break;
//         }

//     }

//     delete [] exps;

//     return sum_phys + sum_reciprocal;

// }

// std::complex<double> QuasiPeriodicGreen::grad_xy_regular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y, const Eigen::Vector2d& a, const Eigen::Vector2d& b) const{

//     //-----------------------
//     // Sum in physical space
//     //-----------------------
//     std::complex<double> sum_phys = 0.0;
//     // bool converged[2] = {false,false};
//     for(int s=0; s<M; s++){
//         // from -infty to -1, and 1 to +infty
//         for(int sign=0; sign<2; sign++){
//             std::complex<double> summand = 0.0;
//             for(int _n=1; _n<NMAX; _n++){
//                 int n = _n * (2*sign-1); // +1,+2,... if sign=1 , -1,-2,... if sign=0

//                 double r = (x-y-n*L*e1).norm();
//                 double zeta = ( a.dot(x-y-n*L*e1) ) * ( b.dot(x-y-n*L*e1)  ) / std::pow(r,2);

//                 // Be careful: exp may overflow for large n!
//                 if(std::exp(-n*beta.imag()) == HUGE_VAL) break;

//                 std::complex<double> tmp = (
//                       bessel_kn(0,sqrt_qs[s]*k*r)   * sqrt_qs[s]*k * zeta
//                     + bessel_kn(1,sqrt_qs[s]*k*r)/r *               (2*zeta - a.dot(b))
//                  ) * std::exp(ione*n*beta);

//                 summand += tmp;

//                 bool flag = abs(tmp) < THRES * abs(summand) || (abs(tmp) == 0);
                
//                 if(flag){
//                     if(PRINT_ITER) std::cout << "# converged at " << _n << std::endl;
                    
//                     break;
//                 }

//             }
//             sum_phys += cs[s]*sqrt_qs[s]*k/(2*M_PI) * summand;

//             // if(converged[0] && converged[1]) break;
//         }
//     }

//     //------------------
//     // Reciprocal space
//     //------------------
//     int mu;
//     if(x[1]-y[1]>0) mu = +1;
//     else mu = -1;

//     std::ofstream ff("./tmp.dat");

//     std::complex<double>* exps = new std::complex<double>[M+1];
//     std::complex<double> sum_reciprocal = 0;
//     for(int _m=0; _m<NMAX; _m++){
//         // std::complex<double> summand = 0.0;
//         std::complex<double> prev = sum_reciprocal;

//         // Use long double for exp1 * exps[s] to avoid loss of digits
//         for(int sign=0; sign<2; sign++){
//             // +1,+2,... if sign=1 , 0,-1,-2,... if sign=0
//             int m = _m * (2*sign-1); 
//             if(m == 0 && sign == 1) continue;

//             std::complex<double> km = (beta+2*m*M_PI)/L;
//             std::complex<double> exp1 = std::exp(ione*km*(x[0]-y[0]));

//             std::complex<double> tmp = 0;

//             // if m is in I
//             int xi;
//             if (I.count(m)) {
//                 // std::cout << m << std::endl;
//                 xi = +1;
//             }
//             else{
//                 xi = -1;
//             }

//             // compute exp (x2-y2) for each s
//             for(int s=0; s<M; s++){
//                 std::complex<double> k2;
//                 double c;
//                 int xi;
//                 k2 = std::sqrt(km*km + qs[s]*k*k);

//                 exps[s] = cs[s] * (-ione*km*b[0] + mu*k2*b[1]) * (+ione*km*a[0] - mu*a[1]*k2) / (2*L*k2) * std::exp(-k2*std::abs(x[1]-y[1]));
//                 //exps[s] = -cs[s]*qs[s]*k*k * std::exp(-k2 * abs(x[1]-y[1])) / (2*L*k2);

//                 tmp += exps[s];
//             }

//             std::complex<double> k2 = std::sqrt(km*km - k*k);        
//             exps[M] = (-ione*km*b[0] - xi*mu*k2*b[1]) * (-ione*km*xi*a[0] - mu*a[1]*k2) / (2*L*k2) * std::exp(xi*k2*std::abs(x[1]-y[1]));
            
//             // exps[M] = k*k*std::exp(-k2 * abs(x[1]-y[1])) / (2*L*k2);
//             tmp += exps[M];


//             // Finally add exps term by term (to avoid loss of significant digits)
//             for(int s=0; s<M+1; s++){
//                 sum_reciprocal += exp1 * exps[s];
//             }

//         }

//         ff << _m << " " << std::abs(sum_reciprocal - prev) / std::abs(sum_reciprocal) << "\n";

//         if(std::abs(sum_reciprocal - prev) < 1e-14 * std::abs(sum_reciprocal)){
//             if(PRINT_ITER) std::cout << "# converged at " << _m << std::endl;
//             break;
//         }

//     }

//     delete [] exps;

//     // std::cout << sum_phys << " " << sum_reciprocal << std::endl;

//     return sum_phys + sum_reciprocal;
    
// }

std::complex<double> QuasiPeriodicGreen::grad_xy_regular_part(const Eigen::Vector2d& x, const Eigen::Vector2d& y, const Eigen::Vector2d& a, const Eigen::Vector2d& b) const{

    //-----------------------
    // Sum in physical space
    //-----------------------
    std::complex<double> sum_phys = 0.0;
    // bool converged[2] = {false,false};
    for(int s=0; s<M; s++){
        // from -infty to -1, and 1 to +infty
        for(int sign=0; sign<2; sign++){
            std::complex<double> summand = 0.0;
            for(int _n=1; _n<NMAX; _n++){
                int n = _n * (2*sign-1); // +1,+2,... if sign=1 , -1,-2,... if sign=0

                double r = (x-y-n*L*e1).norm();
                double zeta = ( a.dot(x-y-n*L*e1) ) * ( b.dot(x-y-n*L*e1)  ) / std::pow(r,2);

                // Be careful: exp may overflow for large n!
                // if(std::exp(-n*beta.imag()) == HUGE_VAL) break; // THIS DOES NOT WORK!!!
                if(-n*beta.imag() > 500) break;

                std::complex<double> tmp = (
                      bessel_kn(0,sqrt_qs[s]*k*r)   * sqrt_qs[s]*k * zeta
                    + bessel_kn(1,sqrt_qs[s]*k*r)/r *               (2*zeta - a.dot(b))
                 ) * std::exp(ione*n*beta);

                summand += tmp;

                bool flag = abs(tmp) < THRES * abs(summand) || (abs(tmp) == 0);
                
                if(flag){
                    if(PRINT_ITER) std::cout << "# converged at " << _n << std::endl;
                    
                    break;
                }

            }
            sum_phys += cs[s]*sqrt_qs[s]*k/(2*M_PI) * summand;

            // if(converged[0] && converged[1]) break;
        }
    }

    //------------------
    // Reciprocal space
    //------------------
    int mu;
    if(x[1]-y[1]>0) mu = +1;
    else mu = -1;

    // std::ofstream ff("./tmp.dat");

    std::complex<double>* exps = new std::complex<double>[M+1];
    std::complex<double> sum_reciprocal = 0;
    for(int _m=0; _m<NMAX; _m++){
        // std::complex<double> summand = 0.0;
        std::complex<double> prev = sum_reciprocal;

        // Use long double for exp1 * exps[s] to avoid loss of digits
        for(int sign=0; sign<2; sign++){
            // +1,+2,... if sign=1 , 0,-1,-2,... if sign=0
            int m = _m * (2*sign-1); 
            if(m == 0 && sign == 1) continue;

            std::complex<double> km = (beta+2*m*M_PI)/L;
            std::complex<double> f       = F(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> f_tilde = F_tilde(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> h       = H(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> h_tilde = H_tilde(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> p       = P(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> p_tilde = P_tilde(m, beta, std::abs(x[1]-y[1]));

            // if m is in I
            if (I.count(m)) sum_reciprocal += std::exp(ione*km*(x[0]-y[0])) / (2*L) * (
                //ione*km*vec[0]*f_tilde - sgn*vec[1]*h_tilde
                + a[0]*b[0] * std::pow(km,2)*f_tilde
                + ione*mu*(a[0]*b[1]+a[1]*b[0]) * km*h_tilde
                - a[1]*b[1] * p_tilde
            );
            else            sum_reciprocal += std::exp(ione*km*(x[0]-y[0])) / (2*L) * (
                //ione*km*vec[0]*f_tilde - sgn*vec[1]*h_tilde
                + a[0]*b[0] * std::pow(km,2)*f
                + ione*mu*(a[0]*b[1]+a[1]*b[0]) * km*h
                - a[1]*b[1] * p
            );

            // std::complex<double> exp1 = std::exp(ione*km*(x[0]-y[0]));

            // std::complex<double> tmp = 0;

            // // if m is in I
            // int xi;
            // if (I.count(m)) {
            //     // std::cout << m << std::endl;
            //     xi = +1;
            // }
            // else{
            //     xi = -1;
            // }

            // // compute exp (x2-y2) for each s
            // for(int s=0; s<M; s++){
            //     std::complex<double> k2;
            //     double c;
            //     int xi;
            //     k2 = std::sqrt(km*km + qs[s]*k*k);

            //     exps[s] = cs[s] * (-ione*km*b[0] + mu*k2*b[1]) * (+ione*km*a[0] - mu*a[1]*k2) / (2*L*k2) * std::exp(-k2*std::abs(x[1]-y[1]));
            //     //exps[s] = -cs[s]*qs[s]*k*k * std::exp(-k2 * abs(x[1]-y[1])) / (2*L*k2);

            //     tmp += exps[s];
            // }

            // std::complex<double> k2 = std::sqrt(km*km - k*k);        
            // exps[M] = (-ione*km*b[0] - xi*mu*k2*b[1]) * (-ione*km*xi*a[0] - mu*a[1]*k2) / (2*L*k2) * std::exp(xi*k2*std::abs(x[1]-y[1]));
            
            // // exps[M] = k*k*std::exp(-k2 * abs(x[1]-y[1])) / (2*L*k2);
            // tmp += exps[M];


            // // Finally add exps term by term (to avoid loss of significant digits)
            // for(int s=0; s<M+1; s++){
            //     sum_reciprocal += exp1 * exps[s];
            // }

        }

        // ff << _m << " " << std::abs(sum_reciprocal - prev) / std::abs(sum_reciprocal) << "\n";

        if(std::abs(sum_reciprocal - prev) < 1e-14 * std::abs(sum_reciprocal)){
            if(PRINT_ITER) std::cout << "# converged at " << _m << std::endl;
            break;
        }

    }

    delete [] exps;

    // std::cout << sum_phys << " " << sum_reciprocal << std::endl;

    return sum_phys + sum_reciprocal;
    
}


std::complex<double> QuasiPeriodicGreen::grad_xy(const Eigen::Vector2d& x, const Eigen::Vector2d& y, const Eigen::Vector2d& a, const Eigen::Vector2d& b) const{

    //-----------------------
    // Sum in physical space
    //-----------------------
    std::complex<double> sum_phys = 0.0;
    // bool converged[2] = {false,false};
    for(int s=0; s<M; s++){
        // from -infty to -1, and 1 to +infty
        for(int sign=0; sign<2; sign++){
            std::complex<double> summand = 0.0;
            for(int _n=0; _n<NMAX; _n++){
                int n = _n * (2*sign-1); // +1,+2,... if sign=1 , -1,-2,... if sign=0
                if(n == 0 && sign == 1) continue;


                double r = (x-y-n*L*e1).norm();
                double zeta = ( a.dot(x-y-n*L*e1) ) * ( b.dot(x-y-n*L*e1)  ) / std::pow(r,2);

                // Be careful: exp may overflow for large n!
                // if(std::exp(-n*beta.imag()) == HUGE_VAL) break;
                if(-n*beta.imag() > 500) break;

                std::complex<double> tmp = (
                      bessel_kn(0,sqrt_qs[s]*k*r)   * sqrt_qs[s]*k * zeta
                    + bessel_kn(1,sqrt_qs[s]*k*r)/r *               (2*zeta - a.dot(b))
                 ) * std::exp(ione*n*beta);

                summand += tmp;

                bool flag = abs(tmp) < THRES * abs(summand) || (abs(tmp) == 0);
                
                if(flag){
                    if(PRINT_ITER) std::cout << "# converged at " << _n << std::endl;
                    
                    break;
                }

            }
            sum_phys += cs[s]*sqrt_qs[s]*k/(2*M_PI) * summand;

            // if(converged[0] && converged[1]) break;
        }
    }

    //------------------
    // Reciprocal space
    //------------------
    int mu;
    if(x[1]-y[1]>0) mu = +1;
    else mu = -1;

    // std::ofstream ff("./tmp.dat");

    std::complex<double>* exps = new std::complex<double>[M+1];
    std::complex<double> sum_reciprocal = 0;
    for(int _m=0; _m<NMAX; _m++){
        // std::complex<double> summand = 0.0;
        std::complex<double> prev = sum_reciprocal;

        // Use long double for exp1 * exps[s] to avoid loss of digits
        for(int sign=0; sign<2; sign++){
            // +1,+2,... if sign=1 , 0,-1,-2,... if sign=0
            int m = _m * (2*sign-1); 
            if(m == 0 && sign == 1) continue;

            std::complex<double> km = (beta+2*m*M_PI)/L;
            std::complex<double> f       = F(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> f_tilde = F_tilde(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> h       = H(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> h_tilde = H_tilde(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> p       = P(m, beta, std::abs(x[1]-y[1]));
            std::complex<double> p_tilde = P_tilde(m, beta, std::abs(x[1]-y[1]));

            // if m is in I
            if (I.count(m)) sum_reciprocal += std::exp(ione*km*(x[0]-y[0])) / (2*L) * (
                //ione*km*vec[0]*f_tilde - sgn*vec[1]*h_tilde
                + a[0]*b[0] * std::pow(km,2)*f_tilde
                + ione*mu*(a[0]*b[1]+a[1]*b[0]) * km*h_tilde
                - a[1]*b[1] * p_tilde
            );
            else            sum_reciprocal += std::exp(ione*km*(x[0]-y[0])) / (2*L) * (
                //ione*km*vec[0]*f_tilde - sgn*vec[1]*h_tilde
                + a[0]*b[0] * std::pow(km,2)*f
                + ione*mu*(a[0]*b[1]+a[1]*b[0]) * km*h
                - a[1]*b[1] * p
            );


        }

        // ff << _m << " " << std::abs(sum_reciprocal - prev) / std::abs(sum_reciprocal) << "\n";

        if(std::abs(sum_reciprocal - prev) < 1e-14 * std::abs(sum_reciprocal)){
            if(PRINT_ITER) std::cout << "# converged at " << _m << std::endl;
            break;
        }

    }

    delete [] exps;

    // std::cout << sum_phys << " " << sum_reciprocal << std::endl;

    return sum_phys + sum_reciprocal;
    
}

std::complex<double> QuasiPeriodicGreen::F_conv(const int& m, const std::complex<double> beta, const double& alpha) const{

    

    std::complex<double> km = (beta+2*m*M_PI)/L;

    std::complex<double> exp;
    std::complex<double> sqrt;

    std::complex<double> result = 0;

    // 0
    sqrt = std::sqrt(km*km - k*k);
    exp  = std::exp(-alpha*sqrt);
    result += exp / sqrt;

    for(int s=0; s<M; s++){
        sqrt = std::sqrt(km*km + qs[s]*k*k);
        exp  = std::exp(-alpha*sqrt);

        result += cs[s] * exp / sqrt;
    }

    return result;

}

// // Compute F^{s-1,s}
// std::complex<double> QuasiPeriodicGreen::F_term_tmp(const int& m, const std::complex<double> beta, const double& alpha, const int& s) const{
//     std::complex<double> km = (beta+2*m*M_PI)/L;

//     std::complex<double> exp;
//     std::complex<double> sqrt;

//     if(s>=0) sqrt = std::sqrt(km*km + qs[s]*k*k);
//     else     sqrt = std::sqrt(km*km -       k*k);
//     exp  = std::exp(-alpha*sqrt);

//     return exp / sqrt;
// }

// Compute F^{s-1,s}
std::complex<double> QuasiPeriodicGreen::F_term(const int& m, const std::complex<double> beta, const double& alpha, const int& s) const{

    std::complex<double> km = (beta+2*m*M_PI)/L;

    std::complex<double> exp;

    // q and sqrt for s-1
    double q_1;
    std::complex<double> sqrt_1;
    if(s == 0){
        sqrt_1 = std::sqrt(km*km - k*k);
        q_1    = -1;
    }
    else{
        sqrt_1 = std::sqrt(km*km + qs[s-1]*k*k);
        q_1    = qs[s-1];
    }

    // q and sqrt for s
    double q_2;
    std::complex<double> sqrt_2;
    sqrt_2 = std::sqrt(km*km + qs[s]*k*k);
    q_2    = qs[s];

    std::complex<double> tmp = alpha*std::pow(k,2)*(q_2-q_1)/(sqrt_1+sqrt_2);

    if(std::abs(
        sqrt_1 * sqrt_2 * (sqrt_2*std::exp(-alpha*sqrt_1) + sqrt_1*std::exp(-alpha*sqrt_2))
    ) < 1e-50){
        return 0.0;
    }
    else{
        return
            (
                2*std::pow(km,2) * std::exp(-2*alpha*sqrt_1) * sigma(tmp) + std::pow(k,2)*(q_2*std::exp(-2*alpha*sqrt_1) - q_1*std::exp(-2*alpha*sqrt_2))
            ) / 
            (
                sqrt_1 * sqrt_2 * (sqrt_2*std::exp(-alpha*sqrt_1) + sqrt_1*std::exp(-alpha*sqrt_2))
            );
    }

    // return
    //     (
    //         2*std::pow(km,2) * std::exp(-2*alpha*sqrt_1) * sigma(tmp) + std::pow(k,2)*(q_2*std::exp(-2*alpha*sqrt_1) - q_1*std::exp(-2*alpha*sqrt_2))
    //     ) / 
    //     (
    //         sqrt_1 * sqrt_2 * (sqrt_2*std::exp(-alpha*sqrt_1) + sqrt_1*std::exp(-alpha*sqrt_2))
    //     );
}

std::complex<double> QuasiPeriodicGreen::F(const int& m, const std::complex<double> beta, const double& alpha) const{
    std::complex<double> result = F_term(m, beta, alpha, 0);
    //F_term_tmp(m, beta, alpha, -1) - F_term_tmp(m, beta, alpha, 0);

    double sum_cs = 1;
    for(int s=0; s<M-1; s++){
        sum_cs += cs[s];
        result += sum_cs * F_term(m, beta, alpha, s+1);

        // result += sum_cs * (F_term_tmp(m, beta, alpha, s) - F_term_tmp(m, beta, alpha, s+1));
    }

    return result;
}

std::complex<double> QuasiPeriodicGreen::F_tilde(const int& m, const std::complex<double> beta, const double& alpha) const{

    std::complex<double> km = (beta+2*m*M_PI)/L;

    std::complex<double> exp;
    std::complex<double> sqrt;

    std::complex<double> result = 0;

    // 0
    sqrt = -std::sqrt(km*km - k*k);
    exp  = std::exp(-alpha*sqrt);
    result += exp / sqrt;

    for(int s=0; s<M; s++){
        sqrt = std::sqrt(km*km + qs[s]*k*k);
        exp  = std::exp(-alpha*sqrt);

        result += cs[s] * exp / sqrt;
    }

    return result;

}

// std::complex<double> QuasiPeriodicGreen::F_hat(const int& m, const std::complex<double> beta, const double& alpha) const{

//     std::complex<double> km = (beta+2*m*M_PI)/L;

//     std::complex<double> exp;
//     std::complex<double> sqrt;

//     std::complex<double> result = 0;

//     // 0
//     sqrt = -std::sqrt(km*km - k*k);
//     exp  = std::exp(-alpha*sqrt);
//     result += -exp / sqrt;

//     for(int s=0; s<M; s++){
//         sqrt = std::sqrt(km*km + qs[s]*k*k);
//         exp  = std::exp(-alpha*sqrt);

//         result += cs[s] * exp / sqrt;
//     }

//     return result;

// }

//---
// H
//---
std::complex<double> QuasiPeriodicGreen::H(const int& m, const std::complex<double> beta, const double& alpha) const{
    std::complex<double> result = H_term(m, beta, alpha, 0);

    double sum_cs = 1;
    for(int s=0; s<M-1; s++){
        sum_cs += cs[s];
        result += sum_cs * H_term(m, beta, alpha, s+1);
    }

    return result;
}

std::complex<double> QuasiPeriodicGreen::H_term(const int& m, const std::complex<double> beta, const double& alpha, const int& s) const{

    std::complex<double> km = (beta+2*m*M_PI)/L;

    std::complex<double> exp;

    // q and sqrt for s-1
    double q_1;
    std::complex<double> sqrt_1;
    if(s == 0){
        sqrt_1 = std::sqrt(km*km - k*k);
        q_1    = -1;
    }
    else{
        sqrt_1 = std::sqrt(km*km + qs[s-1]*k*k);
        q_1    = qs[s-1];
    }

    // q and sqrt for s
    double q_2;
    std::complex<double> sqrt_2;
    sqrt_2 = std::sqrt(km*km + qs[s]*k*k);
    q_2    = qs[s];

    std::complex<double> tmp = (alpha/2)*std::pow(k,2)*(q_2-q_1)/(sqrt_1+sqrt_2);

    return
        2 * std::exp(-alpha*sqrt_1) * sigma(tmp)
    ;
}

std::complex<double> QuasiPeriodicGreen::H_tilde(const int& m, const std::complex<double> beta, const double& alpha) const{

    std::complex<double> km = (beta+2*m*M_PI)/L;

    std::complex<double> exp;
    std::complex<double> sqrt;

    std::complex<double> result = 0;

    // 0
    sqrt = -std::sqrt(km*km - k*k);
    exp  = std::exp(-alpha*sqrt);
    result += exp;

    for(int s=0; s<M; s++){
        sqrt = std::sqrt(km*km + qs[s]*k*k);
        exp  = std::exp(-alpha*sqrt);

        result += cs[s] * exp;
    }

    return result;

}

//---
// P
//---
std::complex<double> QuasiPeriodicGreen::P(const int& m, const std::complex<double> beta, const double& alpha) const{
    std::complex<double> result = P_term(m, beta, alpha, 0);

    double sum_cs = 1;
    for(int s=0; s<M-1; s++){
        sum_cs += cs[s];
        result += sum_cs * P_term(m, beta, alpha, s+1);
    }

    return result;
}

std::complex<double> QuasiPeriodicGreen::P_term(const int& m, const std::complex<double> beta, const double& alpha, const int& s) const{

    std::complex<double> km = (beta+2*m*M_PI)/L;

    std::complex<double> exp;

    // q and sqrt for s-1
    double q_1;
    std::complex<double> sqrt_1;
    if(s == 0){
        sqrt_1 = std::sqrt(km*km - k*k);
        q_1    = -1;
    }
    else{
        sqrt_1 = std::sqrt(km*km + qs[s-1]*k*k);
        q_1    = qs[s-1];
    }

    // q and sqrt for s
    double q_2;
    std::complex<double> sqrt_2;
    sqrt_2 = std::sqrt(km*km + qs[s]*k*k);
    q_2    = qs[s];

    std::complex<double> tmp = alpha*std::pow(k,2)*(q_2-q_1)/(sqrt_1+sqrt_2);

    if(std::abs(
        sqrt_1*std::exp(-alpha*sqrt_1) + sqrt_2*std::exp(-alpha*sqrt_2)
    ) < 1e-50){
        return 0;
    }
    else{
        return
            (
                2*std::pow(km,2) * std::exp(-2*alpha*sqrt_1) * sigma(tmp) + std::pow(k,2)*(q_1*std::exp(-2*alpha*sqrt_1) - q_2*std::exp(-2*alpha*sqrt_2))
            ) / 
            (
                sqrt_1*std::exp(-alpha*sqrt_1) + sqrt_2*std::exp(-alpha*sqrt_2)
            );
    }

    // return
    //     (
    //         2*std::pow(km,2) * std::exp(-2*alpha*sqrt_1) * sigma(tmp) + std::pow(k,2)*(q_1*std::exp(-2*alpha*sqrt_1) - q_2*std::exp(-2*alpha*sqrt_2))
    //     ) / 
    //     (
    //         sqrt_1*std::exp(-alpha*sqrt_1) + sqrt_2*std::exp(-alpha*sqrt_2)
    //     );
}

std::complex<double> QuasiPeriodicGreen::P_tilde(const int& m, const std::complex<double> beta, const double& alpha) const{

    std::complex<double> km = (beta+2*m*M_PI)/L;

    std::complex<double> exp;
    std::complex<double> sqrt;

    std::complex<double> result = 0;

    // 0
    sqrt = -std::sqrt(km*km - k*k);
    exp  = std::exp(-alpha*sqrt);
    result += exp * sqrt;

    for(int s=0; s<M; s++){
        sqrt = std::sqrt(km*km + qs[s]*k*k);
        exp  = std::exp(-alpha*sqrt);

        result += cs[s] * exp * sqrt;
    }

    return result;

}