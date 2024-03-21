#pragma once
#include <set>
#include <vector>
#include "mycomplex.hpp"

// Add index m to sheet if m is not in sheet
// Remove    m          otherwise
static void add_remove_index(std::set<int>& sheet, const int& m){
    if(sheet.count(m) > 0) {
        // if found
        sheet.erase(m);
    } 
    else{
        sheet.insert(m);
    }
}

class Segment{
public:
    double s1, s2;
    std::set<int> sheet;

    Segment(const double& s1, const double& s2, const std::set<int>& sheet) : s1(s1), s2(s2), sheet(sheet){};
};

class Path{
public:
    //-----------
    // Variables
    //-----------
    // Freq. times unit length
    const double kL;
    // Center of the circle
    const std::complex<double> center;
    // Radius
    const double radius;
    // Sheet at the center
    const std::set<int> sheet;

    // 
    // int nseg;
    std::vector<Segment> segments;

    // valid = false if the path contains a branch point
    bool valid = true;

    //--------------
    // Constructors
    //--------------
    Path(const double& kL, const std::complex<double>& center, const double& radius, const std::set<int>& sheet);

    //---------
    // Methods
    //---------
private:
    // jumps of the sheet at beta on the real axis
    std::set<int> jump(const double& beta) const;
public:
    // Return sheet at the given point
    std::set<int> calc_sheet(const double& _theta) const;

    std::set<int> calc_sheet_inside(const std::complex<double>& beta) const;

};