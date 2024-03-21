#pragma once
#include "kernel.hpp"

template <typename T>
class Collocation{
    //-----------
    // Variables
    //-----------
public:
    Kernel<T>* kernel;
    // 2n: number of knots
    int n;
    // knots t_i
    double* knots;
    // Coefficient for ths linear system
    double* Rs;
    // Coefficient for ths linear system (hyper-singular)
    double* Ts;
    

    //--------------
    // Constructors
    //--------------
public:
    Collocation(Kernel<T>* kernel, int n) : n(n){
        this->kernel = kernel;



        //---------
        // Knots
        //---------
        knots = new double[2*n];
        for(int k=0; k<2*n; k++){
            knots[k] = M_PI * k / n;
        }

        //---------
        // Coeff R
        //---------
        Rs = new double[2*n];

        for(int k=0; k<2*n; k++){
            Rs[k] = -std::pow(-1,k)*M_PI / std::pow(n,2);

            for(int m=1; m<=n-1; m++){
                Rs[k] += -2*M_PI/n * std::cos(m*k*M_PI/n) / m;
            }
        }

        //---------
        // Coeff T
        //---------
        Ts = new double[2*n];
        for(int k=0; k<2*n; k++){
            if(k==0){
                Ts[k] = -n*0.5;
            }
            else if(k%2 != 0){
                Ts[k] = 1 / (2*n*std::pow(std::sin(knots[k]/2),2));
            }
            else{
                Ts[k] = 0;
            }
        }

    };

    //------------
    // Destructor
    //------------
public:
    ~Collocation(){
        kernel = nullptr;
        delete [] knots;
        delete [] Rs;
        delete [] Ts;
    }

    //---------
    // Methods
    //---------
public:
    // 係数行列の(i,k)成分を計算
    T calc_entry(const int& i, const int&k) const{

        assert(0 <= i && i < 2*n);
        assert(0 <= k && k < 2*n);

        T out;

        if(i != k){
            out = Ts[abs(i-k)] + Rs[abs(i-k)] * kernel->X1(knots[i],knots[k]) + M_PI/n * kernel->X2(knots[i],knots[k]);
        }
        else{
            out = Ts[abs(i-k)] + Rs[abs(i-k)] * kernel->X1_diag(knots[i]) + M_PI/n * kernel->X2_diag(knots[i]);

            out += kernel->a(knots[i]);
        }

        return out;
    };

    T calc_potential(const KernelPotential<T>& kernelPotential, const Eigen::Vector2d& x, const Eigen::Matrix< T , Eigen::Dynamic , 1>& density) const{
        
        T* tmp = new T[2*n];

        int j;
        #pragma omp parallel for private(j)
        for(j=0; j<2*n; j++){
            // out += M_PI/n * kernelPotential.X(x,knots[j]) * density[j];
            tmp[j] = M_PI/n * kernelPotential.X(x,knots[j]) * density[j];
        }

        T out = 0;
        for(j=0; j<2*n; j++){
            out += tmp[j];
        }
        
        delete [] tmp;

        return out;
    }

    std::complex<double> calc_amplitude_farfield(const BrakhagePotential& kernelPotential, const Eigen::Matrix< std::complex<double> , Eigen::Dynamic , 1>& density, const int& m, const std::string& side) const{

        assert(side == "top" || side == "bottom");
        
        const std::complex<double> ione(0.0,1.0);
        const double k = kernelPotential.k;
        const double L = kernelPotential.L;
        const std::complex<double> beta = kernelPotential.beta;
        const std::set<int> I = kernelPotential.I;
        const std::complex<double> eta = kernelPotential.eta;

        bool unphysical = I.find(m) != I.end(); // if m is in I

        int sign = 1;
        if(side == "bottom") sign *= -1;
        if(unphysical) sign *= -1; 

        std::complex<double> km = (beta + 2*m*M_PI) / L;
        // sqrt(km^2-k^2)
        std::complex<double> k2 = std::sqrt(km*km - k*k);

        // integrate
        std::complex<double> out = 0;
        for(int j=0; j<2*n; j++){
            // quadrature point
            Eigen::Vector2d y = kernel->shape->x(knots[j]);
            Eigen::Vector2d normal = kernel->shape->normal(knots[j]);

            std::complex<double> ex = std::exp(-ione*km*y[0] + sign*k2*y[1]);

            out += M_PI/n * (
                -ione*km* ex * normal[0]
                +sign*k2* ex * normal[1]
                -ione*eta*ex
            ) * density[j];

        }

        if(unphysical) out /= -2*L*k2;
        else           out /= +2*L*k2;

        return out;

    }

};