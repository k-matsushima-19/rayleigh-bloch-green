#include <iostream>
#include <set>
#include <fstream>
#include <iomanip>

#include "shape.hpp"
#include "kernel.hpp"
#include "collocation.hpp"
#include "blockSSM.hpp"
#include "path.hpp"

static const std::complex<double> ione = std::complex<double>(0, 1);

// Compute an eigenvalue beta = (pi,0.400) at k=pi with I=\emptyset, reported by Bennetts and Peter (doi:10.1017/jfm.2022.247)
void example1(void){
    const double k = M_PI; //2.78662;
    const double R = 0.25;
    const std::complex<double> eta = +k;
    std::set<int> I = {};
    const double L = 1.0;
    const int N = 40; //20;
    const int M = 5;
    const std::complex<double> exact(M_PI,0.400204084483313);

    const double beta_imag_max = 0.4;
    

    Geometry::Circle circle(R);

    //-----------------------------
    // Call-back functions for SSM
    //-----------------------------
    std::function<Eigen::MatrixXcd(std::complex<double> beta, const Eigen::MatrixXcd& vec)> f = [&circle,&k,&L,&I,&eta,&N,&beta_imag_max,&M](std::complex<double> beta, const Eigen::MatrixXcd& vec){

        Brakhage kernel(&circle, k, L, beta, I, eta, M, beta_imag_max);
        Collocation c(&kernel, N);

        Eigen::MatrixXcd A(2*N,2*N);
        int i,j;

        #pragma omp parallel for private(i,j)
        for(j=0; j<2*N; j++){
            for(i=0; i<2*N; i++){
                A(i,j) = c.calc_entry(i,j);
            }
        }

        Eigen::MatrixXcd tmp = A.partialPivLu().solve(vec);
        return tmp;
    };

    // Params for SSM
    const int ssm_nskbn = 50;
    const double ssm_rskbn = 0.02;
    const std::complex<double> ssm_center = exact + ssm_rskbn*0.2;

    BlockSSM ssm(2, 2, 1e-4, 1e-4, 2*N);
    std::vector<EigenPair> pairs = ssm.run(f, ssm_nskbn, ssm_rskbn, ssm_center);

    std::cout << "# Eigenvalues" << std::endl;
    for(auto& pair : pairs){
        std::cout << pair.eigenvalue << std::endl;
    }
}

// Compute some eigenvalues and corresponding eigenmodes listed in the paper
void example2(void){
    const double R = 0.25;
    const int N = 40;
    const double L = 1;
    const double rskbn = 1e-3; // 0.01;
    const int nskbn = 60; // 40;

    const double XMIN = -0.5;
    const double XMAX = +0.5;
    const double YMIN = -2.50;
    const double YMAX = +2.50;

    const int NX = 50;
    const int NY = 250;//250;

    Geometry::Circle circle(R);

    
    auto SSM_run = [&circle,&N,&L,&rskbn](const double& k, const std::complex<double>& ssm_center, const std::set<int>& sheet_center) -> std::vector<EigenPair>{
        Path path(k*L, ssm_center, rskbn, sheet_center);

        //-----------------------------
        // Call-back functions for SSM
        //-----------------------------
        std::function<Eigen::MatrixXcd(const double& theta, const Eigen::MatrixXcd& vec)> f = [&circle,&N,&L,&k,&path](const double& theta, const Eigen::MatrixXcd& vec){

            

            // Quadrature point
            std::complex<double> beta = path.center + path.radius*std::exp(ione*theta);
            std::set<int> sheet = path.calc_sheet(theta);

            Brakhage kernel(&circle, k, L, beta, sheet, k, 5, std::abs(path.center)+path.radius);
            Collocation c(&kernel, N);

            Eigen::MatrixXcd A(2*N,2*N);
            int i,j;

            #pragma omp parallel for private(i,j)
            for(j=0; j<2*N; j++){
                for(i=0; i<2*N; i++){
                    A(i,j) = c.calc_entry(i,j);
                }
            }

            Eigen::MatrixXcd tmp = A.partialPivLu().solve(vec);
            return tmp;
        };

        BlockSSM ssm(2, 2, 1e-4, 1e-4, 2*N);
        std::vector<EigenPair> pairs = ssm.run(f, nskbn, path.radius, path.center);

        // return std::abs(pairs[0].eigenvalue.imag());
        for(const auto& pair : pairs){
            std::cout << "# eigenvalue: " << pair.eigenvalue << std::endl;
        }
        return pairs;
        

    };

    // name,k,beta,I
    typedef std::tuple<std::string,double,std::complex<double>,std::set<int>> conf;
    std::vector<conf> variables = 
    {
        // char, k, beta, sheet
        {"B_1", 2.782626896330247, {M_PI, +1e-4}, {}}, 
        {"B_2", 2.782626896330247, {M_PI, +1e-4}, {}}, 
        {"E_1", 6.143443758927157, {+1e-4, -1e-4}, {}},
        {"E_2", 6.143443758927157, {+1e-4, +1e-4}, {}},

        {"I_1", 3.506894388029359, {M_PI, +1e-4}, {}}, 
        {"I_2", 3.506894388029359, {M_PI, -1e-4}, {}}, 
        {"N_1", 12.06329680000105, {5.502543944050176,+1e-4}, {}},
        {"N_2", 12.06329680000105, {5.502543944050176,-1e-4}, {}}, 
    };

    for(auto& var : variables){
        std::cout << "#" << std::endl;
        std::cout << "# " << std::get<0>(var) << std::endl;
        std::cout << "#" << std::endl;

        std::string filename = std::get<0>(var);
        double k = std::get<1>(var);
        std::complex<double> center = std::get<2>(var);
        std::set<int> I = std::get<3>(var);

        std::vector<EigenPair> pairs = SSM_run(k, center, I);

        int ind = 0;
        for(const auto& pair : pairs){
            // std::set<int> sheet(pair.sheet.begin(), pair.sheet.end());
            std::set<int> sheet = I;
            assert(sheet == pair.sheet);

            //------------------------
            // Calc Farfield
            //------------------------
            {
                Brakhage kernel(&circle, k, L, pair.eigenvalue, sheet, k, 5, std::abs(center.imag())+rskbn);
                BrakhagePotential kp(&circle, k, L, pair.eigenvalue, sheet, k, 5);
                Collocation c(&kernel, N);

                std::vector<std::complex<double>> Bps;
                std::vector<std::complex<double>> Bms;
                double scale = 0;
                double phase = 0;
                for(int m=-20; m<20; m++){
                    // calc B^+ (top side)
                    std::complex<double> Bp = c.calc_amplitude_farfield(kp, pairs[0].eigenvector, m, "top");
                    // calc B^- (bottom side)
                    std::complex<double> Bm = c.calc_amplitude_farfield(kp, pairs[0].eigenvector, m, "bottom");

                    scale += std::pow(std::abs(Bp),2);
                    if(m == 0) phase = std::arg(Bp);

                    Bps.push_back(Bp);
                    Bms.push_back(Bm);

                }

                std::ofstream ff("../output/"+filename+"_"+std::to_string(ind)+"_farfield.dat");
                for(int m=-20; m<20; m++){
                    std::complex<double> Bp = Bps[m+20]/std::sqrt(scale) * std::exp(-ione*phase);
                    std::complex<double> Bm = Bms[m+20]/std::sqrt(scale) * std::exp(-ione*phase);

                    ff << m << " " << std::setprecision(16) << 
                        Bp.real() << " " << Bp.imag() << " " <<
                        Bm.real() << " " << Bm.imag() << "\n";

                }

            }

            //------------------------
            // Calc eigenmode
            //------------------------
            {
                
                Brakhage kernel(&circle, k, L, pair.eigenvalue, sheet, k, 5, std::abs(center.imag())+rskbn);
                BrakhagePotential kp(&circle, k, L, pair.eigenvalue, sheet, k, 5);
                Collocation c(&kernel, N);

                // Top
                std::ofstream ff("../output/"+filename+"_"+std::to_string(ind)+"_eigenmode.dat");

                std::complex<double> scale = 1.0;
                for(int iy=0;iy<NY;iy++){
                    for(int ix=0;ix<NX;ix++){
                        const Eigen::Vector2d x(
                            XMIN + (XMAX-XMIN)*ix/NX,
                            YMIN + (YMAX-YMIN)*iy/NY
                        );

                        std::complex<double> u = c.calc_potential(kp, x, pair.eigenvector);

                        ff << std::setprecision(16) << x[0] << " " << x[1] << " " << u.real() << " " << u.imag() << "\n";

                    }
                    ff.flush();
                }
            }
            
            ind++;
        }

    }    
}

int main(void){
    example2();
    example1();
    return 0;
}