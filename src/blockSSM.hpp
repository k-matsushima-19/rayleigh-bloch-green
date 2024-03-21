#pragma once
#include <cassert>
#include <complex>
#include <eigen3/Eigen/Dense>
// #include <eigen3/Eigen/Eigenvalues>
#include <math.h>
#include "progressbar.hpp"

const double PI = M_PI;
const std::complex<double> I(0.0,1.0);

class EigenPair{
public:
    // eigenvalue
    std::complex<double> eigenvalue;
    // eigenvector
    Eigen::VectorXcd eigenvector;
    // sheet on which the eigenvalue lies, if available
    // std::vector<int> sheet;
    std::set<int> sheet;
};

class BlockSSM{
    //-----------
    // Variables
    //-----------
public:
    // block SSMの列数 (計算可能な重複固有値の数の上限) L
    int L;
    // 計算可能な固有値の数の上限 = K*L
    int K;
    // 特異値の切り捨て
    double delta;
    // 固有値の有無の判定基準
    double eps;
    // 行列のサイズ
    int n;

    // column-fullrank random matrices
    Eigen::MatrixXcd V, U;

    //--------------
    // Constructors
    //--------------
public:
    BlockSSM(const int& L, const int& K, const double& delta, const double& eps, const int& n) : L(L), K(K), delta(delta), eps(eps), n(n){
        assert(L <= n);

        //-----------------------------------
        // V = column-fullrank random matrix
        //-----------------------------------
        // repeat until V is full-rank
        while(true){
            V = Eigen::MatrixXcd::Random(n,L); // dynamic, complex, double, random values
            
            Eigen::JacobiSVD<Eigen::MatrixXcd> svd(V);
            Eigen::VectorXd s = svd.singularValues();
            if(s(s.size()-1)/s(0) > 1e-4) break;
        }
        
        //-----------------------------------
        // U = column-fullrank random matrix
        //-----------------------------------
        // repeat until U is full-rank
        while(true){
            U = Eigen::MatrixXcd::Random(n,L); // dynamic, complex, double, random values
            
            Eigen::JacobiSVD<Eigen::MatrixXcd> svd(U);
            Eigen::VectorXd s = svd.singularValues();
            if(s(s.size()-1)/s(0) > 1e-4) break;
        }

        

    };

    //--------------
    // Methods
    //--------------
public:
    // std::vector<EigenPair> run(std::function<std::complex<double>(std::complex<double> z, int i, int j)> g, int N, double rho, std::complex<double> z0) const{
        
    //     std::function<Eigen::MatrixXcd(std::complex<double> z, Eigen::MatrixXcd vec)> f = [&g,this](std::complex<double> z, const Eigen::MatrixXcd& vec){
    //         Eigen::MatrixXcd A(n,n);

            

    //         for(int l=0; l<n; l++){
    //             for(int k=0; k<n; k++){
    //                 A(k,l) = g(z,k,l);
    //             }
    //         }

    //         Eigen::MatrixXcd ans = A.partialPivLu().solve(vec);
    //         return ans;
    //     };

    //     return run(f,N,rho,z0);

    // };

    // f = A^{-1}(z)*ベクトルvecを計算する関数 を渡して，det A(z)=0となるzを計算するメソッド
    // f: complex<double>,Eigen::MatrixXcd -> Eigen::MatrixXcd を返す関数
    // N: 台形積分の点数
    // rho: 積分円の半径
    // z0: 積分円の中心 
    // int n -> 値渡し (コピーして渡す)
    // const int& n -> 参照渡し (コピーせず，参照（ポインタ)を渡す)
    std::vector<EigenPair> run(std::function<Eigen::MatrixXcd(std::complex<double> z, const Eigen::MatrixXcd& vec)> f, int N, double rho, std::complex<double> z0) const{

        

        // 複素平面の積分円上の積分点のリスト
        std::vector<std::complex<double>> zs;
        for(int i=0; i<N; i++){
            double theta = i * 2*PI / N;
            zs.push_back(z0 + rho*std::exp(I*theta));
        }

        // For each z, 行列Yを作る
        // Y = A^{-1}(z) * V
        progressbar bar(N);
        std::vector<Eigen::MatrixXcd> Ys;
        for(const std::complex<double>& z : zs){
            // Eigen::MatrixXcd Y;// = f(z, V);
            Ys.push_back(f(z,V));
            bar.update();
        }

        //-----------
        // Compute S
        //-----------
        // (n,L)のゼロ行列を2*K個作る
        std::vector<Eigen::MatrixXcd> S; // = Eigen::MatrixXcd::Zero(n,L);
        for(int k=0; k<2*K; k++){
            S.push_back(Eigen::MatrixXcd::Zero(n,L));
        }

        // Ymaxを計算する
        // Ymax = 行列Yの成分の絶対値の最大値
        double Ymax = 0.0;
        for(int i=0; i<N; i++){
            std::complex<double> z = zs[i]; // <-
            Eigen::MatrixXcd* Y = &Ys[i];

            // 行列(*Y)の成分の絶対値の最大値を計算
            for (int i = 0; i < Y->rows(); ++i) {
                for (int j = 0; j < Y->cols(); ++j) {
                    double abs_value = std::abs((*Y)(i,j));
                    Ymax = std::max(Ymax, abs_value);
                }
            }

            // 行列Sを計算する
            std::complex<double> zc = (z - z0)/rho;
            for(int k=0; k<2*K; k++){
                // S[k]に zc*Yを足す
                // std::cout << S[k].rows() << " " << Y->rows() << std::endl;
                // std::cout << S[k].cols() << " " << Y->cols() << std::endl;

                S[k] += zc * (*Y);
                zc = zc * (z-z0)/rho;
            }

        }

        //----------------------------------------------
        // Check if eigenvalues exist inside the circle
        //----------------------------------------------
        {
            double max_S0 = 0.0;
            for (int i = 0; i < S[0].rows(); ++i) {
                for (int j = 0; j < S[0].cols(); ++j) {
                    double abs_value = std::abs((S[0])(i,j));
                    max_S0 = std::max(max_S0, abs_value);
                }
            }

            std::cout << "# S[0]: " << max_S0 << std::endl;
            std::cout << "# max(|Y|): " << Ymax/N << std::endl;

            if(max_S0 < Ymax/N * eps){
                std::cout << "# No eigenvalue found" << std::endl;

                // 何も格納していないstd::vectorを返す
                std::vector<EigenPair> empty;
                return empty;
            }

        }

        //-----------
        // Compute M
        //-----------
        std::vector<Eigen::MatrixXcd> M;
        for(const Eigen::MatrixXcd& s : S){
            M.push_back(U.adjoint() * s); // U.H * s
        }
        
        //-----------------
        // Hankel matrices
        //-----------------
        Eigen::MatrixXcd H  = Eigen::MatrixXcd::Zero(L*K,L*K);
        Eigen::MatrixXcd Hs = Eigen::MatrixXcd::Zero(L*K,L*K);

        

        int jj=0;
        for(int j=0; j<K; j++){
            int ii=0;
            for(int i=0; i<K; i++){
                H.block(ii, jj, L, L) = M[i+j];
                Hs.block(ii, jj, L, L) = M[i+j+1];

                ii += L;
            }
            jj += L;
        }

        //----------
        // SVD of H
        //----------
        Eigen::BDCSVD<Eigen::MatrixXcd> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::VectorXd s = svd.singularValues();

        std::cout << "# Singular values: "<< std::endl;
        for(int i=0; i<s.size(); i++){
            std::cout << s[i]/s[0] << " ";
        }
        std::cout << std::endl;

        // 積分円の中の固有値の数 = 特異値がdelta以上の数
        int nev = 0;
        for(int i=0; i<s.size(); i++){
            if(s[i]/s[0] > delta) nev++;
        }

        std::cout << "# Number of eigenvalues: " << nev << std::endl;

        

        // Truncate the matrices H and Hs
        Eigen::MatrixXcd H_2  =  H.block(0,0,nev,nev);
        Eigen::MatrixXcd Hs_2 = Hs.block(0,0,nev,nev);

        

        // 一般化固有値問題を解く
        // Eigen::GeneralizedEigenSolver<Eigen::MatrixXcd> ges(Hs,H);
        // GeneralizedEigenSolverはMatrixXcdに対応していない！
        // Ax=cBxをB^{-1}Ax=cxにして解く
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eig(H_2.fullPivLu().solve(Hs_2));
        Eigen::VectorXcd ev = eig.eigenvalues();
        Eigen::MatrixXcd vr = eig.eigenvectors();

        

        

        //------------------------------------
        // check if all the eigenvalues |e|<1
        //------------------------------------
        bool flag = false;
        for(int i=0; i<ev.size(); i++){
            flag = abs(ev[i]) > 1;
            if(flag) break;
        }
        if(flag){
            std::cout << "# Warning: Eigenvalue found outside the path" << std::endl;
        }

        if(2*K < nev){
            std::cout << "# Error: too small K, nev = " << nev << std::endl;

            // 何も格納していないstd::vectorを返す
            std::vector<EigenPair> empty;
            return empty;
        }

        // S[k]を並べた行列
        Eigen::MatrixXcd Sh(n,L*S.size());
        // 各行列を新しい行列にコピー
        for (size_t i = 0; i < S.size(); ++i) {
            Sh.block(0, i * L, n, L) = S[i];
        }

        Eigen::MatrixXcd Sh_2 = Sh.block(0,0,n,nev);
        Eigen::MatrixXcd vectors = Sh_2 * vr;

        // 格納してreturn
        std::vector<EigenPair> ans;
        for(int i=0; i<nev; i++){
            ans.push_back(EigenPair());
            ans.back().eigenvalue = rho*ev[i]+z0;
            ans.back().eigenvector = vectors.col(i);


            // normalize
            int argmax = 0;
            for(int j=0; j<n; j++){
                if(std::abs(ans.back().eigenvector[j]) > std::abs(ans.back().eigenvector[0])) argmax = j;
            }
            std::complex<double> scale = ans.back().eigenvector[argmax];

            for(int j=0; j<n; j++){
                ans.back().eigenvector[j] /= scale;
            }

        }

        return ans;

    };

    // f = A^{-1}(z)*ベクトルvecを計算する関数 を渡して，det A(z)=0となるzを計算するメソッド
    // f: complex<double>,Eigen::MatrixXcd -> Eigen::MatrixXcd を返す関数
    // N: 台形積分の点数
    // rho: 積分円の半径
    // z0: 積分円の中心 
    // int n -> 値渡し (コピーして渡す)
    // const int& n -> 参照渡し (コピーせず，参照（ポインタ)を渡す)
    std::vector<EigenPair> run(std::function<Eigen::MatrixXcd(const double& theta, const Eigen::MatrixXcd& vec)> f, int N, double rho, std::complex<double> z0) const{

        

        // 複素平面の積分円上の積分点のリスト
        std::vector<std::complex<double>> zs;
        for(int i=0; i<N; i++){
            double theta = i * 2*PI / N;
            zs.push_back(z0 + rho*std::exp(I*theta));
        }

        // For each z, 行列Yを作る
        // Y = A^{-1}(z) * V
        progressbar bar(N);
        std::vector<Eigen::MatrixXcd> Ys;
        for(const std::complex<double>& z : zs){
            // Eigen::MatrixXcd Y;// = f(z, V);
            double theta = std::arg(z - z0);
            Ys.push_back(f(theta,V));
            bar.update();
        }

        //-----------
        // Compute S
        //-----------
        // (n,L)のゼロ行列を2*K個作る
        std::vector<Eigen::MatrixXcd> S; // = Eigen::MatrixXcd::Zero(n,L);
        for(int k=0; k<2*K; k++){
            S.push_back(Eigen::MatrixXcd::Zero(n,L));
        }

        // Ymaxを計算する
        // Ymax = 行列Yの成分の絶対値の最大値
        double Ymax = 0.0;
        for(int i=0; i<N; i++){
            std::complex<double> z = zs[i]; // <-
            Eigen::MatrixXcd* Y = &Ys[i];

            // 行列(*Y)の成分の絶対値の最大値を計算
            for (int i = 0; i < Y->rows(); ++i) {
                for (int j = 0; j < Y->cols(); ++j) {
                    double abs_value = std::abs((*Y)(i,j));
                    Ymax = std::max(Ymax, abs_value);
                }
            }

            // 行列Sを計算する
            std::complex<double> zc = (z - z0)/rho;
            for(int k=0; k<2*K; k++){
                // S[k]に zc*Yを足す
                // std::cout << S[k].rows() << " " << Y->rows() << std::endl;
                // std::cout << S[k].cols() << " " << Y->cols() << std::endl;

                S[k] += zc * (*Y);
                zc = zc * (z-z0)/rho;
            }

        }

        //----------------------------------------------
        // Check if eigenvalues exist inside the circle
        //----------------------------------------------
        {
            double max_S0 = 0.0;
            for (int i = 0; i < S[0].rows(); ++i) {
                for (int j = 0; j < S[0].cols(); ++j) {
                    double abs_value = std::abs((S[0])(i,j));
                    max_S0 = std::max(max_S0, abs_value);
                }
            }

            std::cout << "# S[0]: " << max_S0 << std::endl;
            std::cout << "# max(|Y|): " << Ymax/N << std::endl;

            if(max_S0 < Ymax/N * eps){
                std::cout << "# No eigenvalue found" << std::endl;

                // 何も格納していないstd::vectorを返す
                std::vector<EigenPair> empty;
                return empty;
            }

        }

        //-----------
        // Compute M
        //-----------
        std::vector<Eigen::MatrixXcd> M;
        for(const Eigen::MatrixXcd& s : S){
            M.push_back(U.adjoint() * s); // U.H * s
        }
        
        //-----------------
        // Hankel matrices
        //-----------------
        Eigen::MatrixXcd H  = Eigen::MatrixXcd::Zero(L*K,L*K);
        Eigen::MatrixXcd Hs = Eigen::MatrixXcd::Zero(L*K,L*K);

        

        int jj=0;
        for(int j=0; j<K; j++){
            int ii=0;
            for(int i=0; i<K; i++){
                H.block(ii, jj, L, L) = M[i+j];
                Hs.block(ii, jj, L, L) = M[i+j+1];

                ii += L;
            }
            jj += L;
        }

        //----------
        // SVD of H
        //----------
        Eigen::BDCSVD<Eigen::MatrixXcd> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::VectorXd s = svd.singularValues();

        std::cout << "# Singular values: "<< std::endl;
        for(int i=0; i<s.size(); i++){
            std::cout << s[i]/s[0] << " ";
        }
        std::cout << std::endl;

        // 積分円の中の固有値の数 = 特異値がdelta以上の数
        int nev = 0;
        for(int i=0; i<s.size(); i++){
            if(s[i]/s[0] > delta) nev++;
        }

        std::cout << "# Number of eigenvalues: " << nev << std::endl;

        

        // Truncate the matrices H and Hs
        Eigen::MatrixXcd H_2  =  H.block(0,0,nev,nev);
        Eigen::MatrixXcd Hs_2 = Hs.block(0,0,nev,nev);

        

        // 一般化固有値問題を解く
        // Eigen::GeneralizedEigenSolver<Eigen::MatrixXcd> ges(Hs,H);
        // GeneralizedEigenSolverはMatrixXcdに対応していない！
        // Ax=cBxをB^{-1}Ax=cxにして解く
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eig(H_2.fullPivLu().solve(Hs_2));
        Eigen::VectorXcd ev = eig.eigenvalues();
        Eigen::MatrixXcd vr = eig.eigenvectors();

        

        

        //------------------------------------
        // check if all the eigenvalues |e|<1
        //------------------------------------
        bool flag = false;
        for(int i=0; i<ev.size(); i++){
            flag = abs(ev[i]) > 1;
            if(flag) break;
        }
        if(flag){
            std::cout << "# Warning: Eigenvalue found outside the path" << std::endl;
        }

        if(2*K < nev){
            std::cout << "# Error: too small K, nev = " << nev << std::endl;

            // 何も格納していないstd::vectorを返す
            std::vector<EigenPair> empty;
            return empty;
        }

        // S[k]を並べた行列
        Eigen::MatrixXcd Sh(n,L*S.size());
        // 各行列を新しい行列にコピー
        for (size_t i = 0; i < S.size(); ++i) {
            Sh.block(0, i * L, n, L) = S[i];
        }

        Eigen::MatrixXcd Sh_2 = Sh.block(0,0,n,nev);
        Eigen::MatrixXcd vectors = Sh_2 * vr;

        // 格納してreturn
        std::vector<EigenPair> ans;
        for(int i=0; i<nev; i++){
            ans.push_back(EigenPair());
            ans.back().eigenvalue = rho*ev[i]+z0;
            ans.back().eigenvector = vectors.col(i);


            // normalize
            int argmax = 0;
            for(int j=0; j<n; j++){
                if(std::abs(ans.back().eigenvector[j]) > std::abs(ans.back().eigenvector[0])) argmax = j;
            }
            std::complex<double> scale = ans.back().eigenvector[argmax];

            for(int j=0; j<n; j++){
                ans.back().eigenvector[j] /= scale;
            }

        }

        return ans;

    };


};