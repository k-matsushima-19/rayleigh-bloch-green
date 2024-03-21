#include "path.hpp"
#include <cassert>
#include <iostream>
#include <algorithm>

const static std::complex<double> ione(0, 1);

// Function to calculate the intersection angles
static std::vector<double> calculateIntersectionAngles(double R, std::complex<double> c) {
    std::vector<double> angles;

    // Line Im[z] = 0
    double delta = R * R - c.imag() * c.imag();
    if (delta >= 0) {
        double sqrt_delta = std::sqrt(delta);
        std::complex<double> z1(sqrt_delta + c.real(), 0);
        std::complex<double> z2(-sqrt_delta + c.real(), 0);
        angles.push_back(std::arg(z1 - c));
        angles.push_back(std::arg(z2 - c));
    }

    // Line Re[z] = 0
    delta = R * R - c.real() * c.real();
    if (delta >= 0) {
        double sqrt_delta = std::sqrt(delta);
        std::complex<double> z1(0, sqrt_delta + c.imag());
        std::complex<double> z2(0, -sqrt_delta + c.imag());
        angles.push_back(std::arg(z1 - c));
        angles.push_back(std::arg(z2 - c));
    }

    // Line Re[z] = 2*pi
    double a = 2 * M_PI;
    delta = R * R - (a - c.real()) * (a - c.real());
    if (delta >= 0) {
        double sqrt_delta = std::sqrt(delta);
        std::complex<double> z1(a, sqrt_delta + c.imag());
        std::complex<double> z2(a, -sqrt_delta + c.imag());
        angles.push_back(std::arg(z1 - c));
        angles.push_back(std::arg(z2 - c));
    }

    std::sort(angles.begin(), angles.end());

    return angles;
}

Path::Path(const double& kL, const std::complex<double>& center, const double& radius, const std::set<int>& sheet) : kL(kL), center(center), radius(radius), sheet(sheet){

    assert(0 < radius && radius < M_PI);
    assert(0 < center.real() &&  center.real() < 2*M_PI);

    std::vector<double> angles = calculateIntersectionAngles(radius, center);
    // for(auto& e : angles) std::cout << e << std::endl;

    assert(angles.size() != 1);

    // The circle does not cross any lines
    if(angles.size() == 0){
        segments.push_back(Segment(-M_PI,+M_PI,sheet));
    }
    else{
        int ii = -1;
        
        if(center.imag() > 0){
            if(center.real() < 0){
                // Find the arc with Re<0 and Im>0
                for(int i=0; i<angles.size(); i++){
                    int j;
                    if(i == angles.size()-1) j = 0;
                    else j = i+1;

                    // center of the arc
                    // std::complex<double> c = center + radius*std::exp(ione*(angles[i] + angles[j])/2);
                    std::complex<double> c;
                    {
                        double angle1 = angles[i];
                        double angle2 = angles[j];
                        if(angle1 > angle2) angle2 += 2*M_PI;

                        c = center + radius*std::exp(ione*(angle1 + angle2)/2);
                    }

                    if(c.real() <= 0 && c.imag() >= 0){
                        ii = i;
                        break;
                    }
                }
            }
            else if(center.real() < 2*M_PI){
                // Find the arc with 0<Re<2pi and Im>0
                for(int i=0; i<angles.size(); i++){
                    int j;
                    if(i == angles.size()-1) j = 0;
                    else j = i+1;

                    // center of the arc
                    std::complex<double> c;
                    {
                        double angle1 = angles[i];
                        double angle2 = angles[j];
                        if(angle1 > angle2) angle2 += 2*M_PI;

                        c = center + radius*std::exp(ione*(angle1 + angle2)/2);
                    }

                    if(0<=c.real() && c.real()<=2*M_PI && c.imag() >= 0){
                        ii = i;
                        break;
                    }
                }
            }
            else{
                // Find the arc with 2pi<Re and Im>0
                for(int i=0; i<angles.size(); i++){
                    int j;
                    if(i == angles.size()-1) j = 0;
                    else j = i+1;

                    // center of the arc
                    std::complex<double> c;
                    {
                        double angle1 = angles[i];
                        double angle2 = angles[j];
                        if(angle1 > angle2) angle2 += 2*M_PI;

                        c = center + radius*std::exp(ione*(angle1 + angle2)/2);
                    }

                    if(c.real() >= 2*M_PI && c.imag() >= 0){
                        ii = i;
                        break;
                    }
                }
            }
        }
        else{
            if(center.real() < 0){
                // Find the arc with Re<0 and Im<0
                for(int i=0; i<angles.size(); i++){
                    int j;
                    if(i == angles.size()-1) j = 0;
                    else j = i+1;

                    // center of the arc
                    std::complex<double> c;
                    {
                        double angle1 = angles[i];
                        double angle2 = angles[j];
                        if(angle1 > angle2) angle2 += 2*M_PI;

                        c = center + radius*std::exp(ione*(angle1 + angle2)/2);
                    }

                    if(c.real() <= 0 && c.imag() <= 0){
                        ii = i;
                        break;
                    }
                }
            }
            else if(center.real() < 2*M_PI){
                // Find the arc with 0<Re<2pi and Im<0
                for(int i=0; i<angles.size(); i++){
                    int j;
                    if(i == angles.size()-1) j = 0;
                    else j = i+1;

                    // center of the arc
                    std::complex<double> c;
                    {
                        double angle1 = angles[i];
                        double angle2 = angles[j];
                        if(angle1 > angle2) angle2 += 2*M_PI;

                        c = center + radius*std::exp(ione*(angle1 + angle2)/2);
                    }

                    if(0<=c.real() && c.real()<=2*M_PI && c.imag() <= 0){
                        ii = i;
                        break;
                    }
                }
            }
            else{
                // Find the arc with 2pi<Re and Im<0
                for(int i=0; i<angles.size(); i++){
                    int j;
                    if(i == angles.size()-1) j = 0;
                    else j = i+1;

                    // center of the arc
                    std::complex<double> c;
                    {
                        double angle1 = angles[i];
                        double angle2 = angles[j];
                        if(angle1 > angle2) angle2 += 2*M_PI;

                        c = center + radius*std::exp(ione*(angle1 + angle2)/2);
                    }
                    
                    if(c.real() >= 2*M_PI && c.imag() <= 0){
                        ii = i;
                        break;
                    }
                }
            }
        }

        // // Find the starting angle ! BUG
        // int ii = -1;
        // if(center.real() < M_PI){
        //     for(int i=0; i<angles.size(); i++){
        //         int j;
        //         if(i == angles.size()-1) j = 0;
        //         else j = i+1;

        //         // center of the arc
        //         std::complex<double> c = center + radius*std::exp(ione*(angles[i] + angles[j])/2);
        //         if(c.real()*center.real()>0 && c.imag()*center.imag()>0){
        //             ii = i;
        //             break;
        //         }
        //     }
        // }
        // else{
        //     for(int i=0; i<angles.size(); i++){
        //         int j;
        //         if(i == angles.size()-1) j = 0;
        //         else j = i+1;

        //         // center of the arc
        //         std::complex<double> c = center + radius*std::exp(ione*(angles[i] + angles[j])/2);
        //         if(c.real()*center.real()<2*M_PI && c.imag()*center.imag()>0){
        //             ii = i;
        //             break;
        //         }
        //     }
        // }

        assert(ii >= 0);

        std::rotate(angles.begin(), angles.begin() + ii, angles.end());

        segments.push_back(Segment(angles[0], angles[1], sheet));

        

        std::set<int> sheet_prev = sheet;
        for(int i=1; i<angles.size(); i++){
            // crossing point
            std::complex<double> cp = center + radius*std::exp(ione*angles[i]);

            if(std::abs(cp.imag()) > 1e-10*radius){
                int ind = -(cp.real()/(2*M_PI));
                add_remove_index(sheet_prev, ind);
            }
            else{
                std::set<int> jumps = jump(cp.real());
                for(auto& m : jumps) add_remove_index(sheet_prev, m);
            }

            if(i+1 == angles.size()){
                segments.push_back(Segment(angles[i], angles[0], sheet_prev));
            }
            else{
                segments.push_back(Segment(angles[i], angles[i+1], sheet_prev));
            }
            

            
        }

        // Check
        // crossing point
        std::complex<double> cp = center + radius*std::exp(ione*angles[0]);

        if(std::abs(cp.imag()) > 1e-10*radius){
            int ind = -(cp.real()/(2*M_PI));
            add_remove_index(sheet_prev, ind);
        }
        else{
            std::set<int> jumps = jump(cp.real());
            for(auto& m : jumps) add_remove_index(sheet_prev, m);
        }

        valid = sheet_prev == segments[0].sheet;

        
    }

    for(auto& seg : segments){
        // std::cout << seg.s1 << " " << seg.s2 << std::endl;
        // for(auto& e : seg.sheet) std::cout << e << " "; std::cout << std::endl;
    }


    // // c := 2*pi*argmin |2*pi*m-center|
    // int ind = -std::round(center.real() / (2*M_PI));
    // double c = 2*M_PI * (-ind);

    // //---------------
    // // Four segments
    // //---------------
    // if(std::abs(center-c) < radius){
    //     nseg = 4;
    //     // Sheets on each quadrant
    //     segments = { {}, {}, {}, {} };
        
    //     segments[0].second = 0;
    //     segments[1].second = M_PI/2;
    //     segments[2].second = M_PI;
    //     segments[3].second = 3*M_PI/2;

    //     if(center.imag() > 0){
    //         if(center.real() > c){
    //             // First quadrant
    //             segments[0].first = sheet;

    //             // Second quadrant
    //             {
    //                 segments[1].first = segments[0].first;
    //                 add_remove_index(segments[1].first, ind);
    //             }
                

    //             // Third quadrant
    //             {
    //                 segments[2].first = segments[1].first;
    //                 std::set<int> jumps = jump(center.real()-radius);
    //                 for(auto& m : jumps) add_remove_index(segments[2].first, m);
    //             }
                

    //             // Fourth quadrant
    //             {
    //                 segments[3].first = segments[2].first;
    //                 add_remove_index(segments[3].first, ind);
    //             }
                

    //             // Check
    //             {
    //                 std::set<int> tmp = segments[3].first;
    //                 std::set<int> jumps = jump(center.real()+radius);
    //                 for(auto& m : jumps) add_remove_index(tmp, m);
    //                 valid = tmp == segments[0].first;
    //             }

    //             // for(auto& e : segments[0]) std::cout << e << " "; std::cout << std::endl;
    //             // for(auto& e : segments[1]) std::cout << e << " "; std::cout << std::endl;
    //             // for(auto& e : segments[2]) std::cout << e << " "; std::cout << std::endl;
    //             // for(auto& e : segments[3]) std::cout << e << " "; std::cout << std::endl;


    //         }
    //         else{

    //             // Second quadrant
    //             segments[1].first = sheet;
                

    //             // Third quadrant
    //             {
    //                 segments[2].first = segments[1].first;
    //                 std::set<int> jumps = jump(center.real()-radius);
    //                 for(auto& m : jumps) add_remove_index(segments[2].first, m);
    //             }
                

    //             // Fourth quadrant
    //             {
    //                 segments[3].first = segments[2].first;
    //                 add_remove_index(segments[3].first, ind);
    //             }

    //             // First quadrant
    //             segments[0].first = segments[3].first;
    //             std::set<int> jumps = jump(center.real()+radius);
    //             for(auto& m : jumps) add_remove_index(segments[0].first, m);
                

    //             // Check
    //             {
    //                 std::set<int> tmp = segments[0].first;
    //                 add_remove_index(tmp, ind);
    //                 valid = tmp == segments[1].first;
    //             }

    //             // for(auto& e : segments[0]) std::cout << e << " "; std::cout << std::endl;
    //             // for(auto& e : segments[1]) std::cout << e << " "; std::cout << std::endl;
    //             // for(auto& e : segments[2]) std::cout << e << " "; std::cout << std::endl;
    //             // for(auto& e : segments[3]) std::cout << e << " "; std::cout << std::endl;
    //         }
            
    //     }
    //     else{
    //         if(center.real() < c){
    //             // Third quadrant
    //             {
    //                 segments[2].first = sheet;
    //             }
                

    //             // Fourth quadrant
    //             {
    //                 segments[3].first = segments[2].first;
    //                 add_remove_index(segments[3].first, ind);
    //             }

    //             // First quadrant
    //             segments[0].first = segments[3].first;
    //             std::set<int> jumps = jump(center.real()+radius);
    //             for(auto& m : jumps) add_remove_index(segments[0].first, m);
                
    //             // Second quadrant
    //             {
    //                 segments[1].first = segments[0].first;
    //                 add_remove_index(segments[1].first, ind);
    //             }

    //             // Check
    //             {
    //                 std::set<int> tmp = segments[1].first;
    //                 std::set<int> jumps = jump(center.real()-radius);
    //                 for(auto& m : jumps) add_remove_index(tmp, m);
    //                 valid = tmp == segments[2].first;
    //             }

    //             // for(auto& e : segments[0]) std::cout << e << " "; std::cout << std::endl;
    //             // for(auto& e : segments[1]) std::cout << e << " "; std::cout << std::endl;
    //             // for(auto& e : segments[2]) std::cout << e << " "; std::cout << std::endl;
    //             // for(auto& e : segments[3]) std::cout << e << " "; std::cout << std::endl;
    //         }
    //         else{                    

    //             // Fourth quadrant
    //             {
    //                 segments[3].first = sheet;
    //             }

    //             // First quadrant
    //             segments[0].first = segments[3].first;
    //             std::set<int> jumps = jump(center.real()+radius);
    //             for(auto& m : jumps) add_remove_index(segments[0].first, m);
                
    //             // Second quadrant
    //             {
    //                 segments[1].first = segments[0].first;
    //                 add_remove_index(segments[1].first, ind);
    //             }

    //             // Third quadrant
    //             {
    //                 segments[2].first = segments[1].first;
    //                 std::set<int> jumps = jump(center.real()-radius);
    //                 for(auto& m : jumps) add_remove_index(segments[2].first, m);
    //             }

    //             // Check
    //             {
    //                 std::set<int> tmp = segments[2].first;
    //                 add_remove_index(tmp, ind);
    //                 valid = tmp == segments[3].first;
    //             }

    //             // for(auto& e : segments[0]) std::cout << e << " "; std::cout << std::endl;
    //             // for(auto& e : segments[1]) std::cout << e << " "; std::cout << std::endl;
    //             // for(auto& e : segments[2]) std::cout << e << " "; std::cout << std::endl;
    //             // for(auto& e : segments[3]) std::cout << e << " "; std::cout << std::endl;
    //         }
            
    //     }

    // }
        
    // //---------------------------------------
    // // Two segments, Crossing imaginary axis
    // //---------------------------------------
    // else if(std::abs(center.real() - c) < radius){
    //     // Find two crossing points on re=c
    //     std::complex<double> z1, z2;
    //     {
    //         double delta = radius * radius - (c - center.real()) * (c - center.real());
    //         if (delta < 0) {
    //             std::cout << "No intersection points." << std::endl;
    //             assert(false);
    //         }

    //         double sqrt_delta = std::sqrt(delta);

    //         z1 = std::complex<double>(c, sqrt_delta + center.imag());
    //         z2 = std::complex<double>(c, -sqrt_delta + center.imag());
    //     }

    //     // Angles corresponding to z1, z2 [-pi,pi]
    //     double theta1 = std::arg(z1-center);
    //     double theta2 = std::arg(z2-center);

    //     // center is left
    //     if(theta1 < M_PI/2){
    //         nseg = 2;
    //         segments = { {}, {}};

    //         segments[0].second = theta1;
    //         segments[1].second = theta2;

    //         segments[0].first = sheet;

    //         segments[1].first = segments[0].first;
    //         add_remove_index(segments[1].first, ind);
    //     }
    //     else{
    //         nseg = 2;
    //         segments = { {}, {}};

    //         segments[0].second = theta2;
    //         segments[1].second = theta1;

    //         segments[0].first = sheet;

    //         segments[1].first = segments[0].first;
    //         add_remove_index(segments[1].first, ind);
    //     }
    // }

    // //---------------------------------------
    // // Two segments, Crossing real axis
    // //---------------------------------------
    // else if(std::abs(center.imag()) < radius){
    //     // Find two crossing points on Im=0
    //     std::complex<double> z1, z2;
    //     {
    //         double delta = radius * radius - (center.imag()) * (center.imag());
    //         if (delta < 0) {
    //             std::cout << "No intersection points." << std::endl;
    //             assert(false);
    //         }

    //         double sqrt_delta = std::sqrt(delta);

    //         z1 = std::complex<double>(sqrt_delta + center.real(), 0);
    //         z2 = std::complex<double>(-sqrt_delta + center.real(), 0);
    //     }

    //     // Angles corresponding to z1, z2 [-pi,pi]
    //     double theta1 = std::arg(z1-center);
    //     double theta2 = std::arg(z2-center);

    // }



}

std::set<int> Path::jump(const double& beta) const{
    // Find MMAX such that -2*pi*MMAX + kL < beta
    // MMAX > (kL-beta)/(2*pi)
    int MMAX = std::ceil((kL-beta)/(2*M_PI));

    // Find MMIN such that -2*pi*MMIN - kL > beta
    // MMIN < (-kL-beta)/(2*pi)
    int MMIN = std::floor((-kL-beta)/(2*M_PI));

    std::set<int> out;

    // For each branch cut
    for(int m=MMIN; m<=MMAX; m++){
        double start = -2*M_PI*m - kL;
        double end = -2*M_PI*m + kL;

        // check if beta is in [start,end]
        if(start <= beta && beta <= end) out.insert(m);
    }

    return out;

}

// std::set<int> Path::jump_old(const double& beta) const{
//     int ind = -std::round(beta / (2*M_PI));
//     double beta_0 = 2*M_PI * (-ind);
//     double beta_1 = beta_0 + 2*M_PI;

//     std::set<int> out = {};

//     // Check if the beta is on the branch cut ind, ind+1, ind+2,...
//     int m = 0;
//     while(true){
//         if(beta - beta_0+2*M_PI*m < kL) out.insert(ind+m);
//         else break;

//         m++;
//     }
    
//     // Check if the beta is on the branch cut ind-1, ind-2, ind-3,...
//     m = +1;
//     while(true){
//         if(beta_1 - beta +2*M_PI*m < kL) out.insert(ind-m);
//         else break;

//         m++;
//     }

//     return out;

// }

static double normalize_angle(double angle) {
    const double TWO_PI = 2.0 * M_PI;
    angle = fmod(angle, TWO_PI);
    if (angle <= -M_PI) {
        angle += TWO_PI;
    } else if (angle > M_PI) {
        angle -= TWO_PI;
    }
    return angle;
}

// Return sheet at the given point
std::set<int> Path::calc_sheet(const double& _theta) const{

    if(segments.size() == 1) return sheet;

    // [-pi,pi]
    double theta = normalize_angle(_theta);
    

    // Check which sheet the point belongs to
    bool flag = false;
    std::set<int> tmp;
    for(const auto& seg : segments){
        if(seg.s1 > seg.s2){
            tmp = seg.sheet;
            flag = true;
        }

        if(seg.s1 <= theta && theta < seg.s2) return seg.sheet;
    }

    assert(flag);
    return tmp;

    // std::cout << theta << std::endl;

    // std::cerr << "# Error: failed at Path::calc_sheet" << std::endl;
    // assert(false);
}

// Return sheet at given point within the path
std::set<int> Path::calc_sheet_inside(const std::complex<double>& beta) const{

    
    if(segments.size() == 1){
        return sheet;
    }
    else if(segments.size() == 2){
        if(std::abs(center.imag()) < radius){
            // if the path is cut along the real axis

            if(center.imag()*beta.imag() > 0){
                return sheet;
            }
            else{
                if(center.imag() < 0){
                    return calc_sheet(M_PI/2);
                }
                else{
                    return calc_sheet(-M_PI/2);
                }
            }
        }
        else if(center.real() < M_PI){
            // if the path is cut along the Re=0 axis
            
            if(beta.real() > 0){
                return sheet;
            }
            else{
                return calc_sheet(M_PI);
            }
        }
        else{
            // if the path is cut along the Re=2pi axis
            
            if(beta.real() < 2*M_PI){
                return sheet;
            }
            else{
                return calc_sheet(0);
            }
        }
    }
    else if(segments.size() == 3){
        if(center.imag() < 0 && center.real() < M_PI){
            if(beta.imag() < 0 && beta.real() > 0){
                return sheet;
            }
            else if(beta.imag()>=0){
                return calc_sheet(M_PI/2);
            }
            else{
                return calc_sheet(M_PI);
            }
        }
        else if(center.imag() >= 0 && center.real() < M_PI){
            if(beta.imag() >= 0 && beta.real() > 0){
                return sheet;
            }
            else if(beta.imag() < 0){
                return calc_sheet(-M_PI/2);
            }
            else{
                return calc_sheet(M_PI);
            }
        }
        else if(center.imag() < 0 && center.real() >= M_PI){
            if(beta.imag() < 0 && beta.real() < 2*M_PI){
                return sheet;
            }
            else if(beta.imag()>=0){
                return calc_sheet(M_PI/2);
            }
            else{
                return calc_sheet(0);
            }
        }
        else{
            if(beta.imag() >= 0 && beta.real() < 2*M_PI){
                return sheet;
            }
            else if(beta.imag() < 0){
                return calc_sheet(-M_PI/2);
            }
            else{
                return calc_sheet(0);
            }
        }

    }
    else{
        if(center.imag() < 0 && center.real() < M_PI){
            if(beta.imag() < 0 && beta.real() > 0){
                return sheet;
            }
            else if(beta.imag()>=0 && beta.real() > 0){
                return calc_sheet(M_PI/2);
            }
            else if(beta.imag() < 0 && beta.real() <= 0){
                return calc_sheet(M_PI);
            }
            else{
                return calc_sheet(std::arg(0-center));
            }
        }
        else if(center.imag() >= 0 && center.real() < M_PI){
            if(beta.imag() >= 0 && beta.real() > 0){
                return sheet;
            }
            else if(beta.imag() < 0 && beta.real() > 0){
                return calc_sheet(-M_PI/2);
            }
            else if(beta.imag() >= 0 && beta.real() <= 0){
                return calc_sheet(M_PI);
            }
            else{
                return calc_sheet(std::arg(0-center));
            }
        }
        else if(center.imag() < 0 && center.real() >= M_PI){
            if(beta.imag() < 0 && beta.real() < 2*M_PI){
                return sheet;
            }
            else if(beta.imag()>=0 && beta.real() < 2*M_PI){
                return calc_sheet(M_PI/2);
            }
            else if(beta.imag() < 0 && beta.real() >= 2*M_PI){
                return calc_sheet(0);
            }
            else{
                return calc_sheet(std::arg(2*M_PI-center));
            }
        }

        else{
            if(beta.imag() >= 0 && beta.real() < 2*M_PI){
                return sheet;
            }
            else if(beta.imag() < 0 && beta.real() < 2*M_PI){
                return calc_sheet(-M_PI/2);
            }
            else if(beta.imag() >= 0 && beta.real() >= 2*M_PI){
                return calc_sheet(0);
            }
            else{
                return calc_sheet(std::arg(2*M_PI-center));
            }
        }
    }
}