//
//  main.cpp
//  EigenTutorial
//
//  Created by Drinor Cacaj on 10.10.22.
//

#include "Eigen/Dense"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include "Configuration.hpp"

using namespace std;

double Configuration::q = 6;
std::vector<double> Configuration::rhos = linspace(0, M_PI, 100);
std::vector<double> Configuration::lambdas = {0.41025641, 0.4312952 , 0.4534129 , 0.47666485, 0.5011092 ,
    0.52680711, 0.55382285, 0.58222403, 0.61208167, 0.64347047,
    0.67646896, 0.71115967, 0.7476294 , 0.78596937, 0.82627549,
    0.86864859, 0.91319468, 0.96002517, 1.00925723, 1.06101401,
    1.11542499, 1.17262627, 1.23276095, 1.29597946, 1.36243995,
    1.43230866, 1.50576039, 1.58297887, 1.66415727, 1.74949867,
    1.82161807};
std::vector<double> Configuration::hom_spectrum = std::vector<double>(Configuration::lambdas.size(),1.);
std::vector<double> Configuration::zero_spectrum = std::vector<double>(Configuration::lambdas.size(),0.);

int main(int argc, const char * argv[]) {
    std::cout << std::setprecision(15);
    std::cout << "loading input ..";
    
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_csv_file>" << std::endl;
        return 1;
    }

    const std::string csvFilePath = argv[1];
    std::ifstream file(csvFilePath);

    if (!file) {
        std::cerr << "Error opening file: " << csvFilePath << std::endl;
        return 1;
    }
    
    std::cout << " done" << std::endl;
    
    
    std::vector<Configuration> Ps;
    std::vector<double> data;
    std::string line, value;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        data = vector<double>();
        while (std::getline(iss, value, ',')) {
            try {
                double number = std::stod(value);
                data.push_back(number);
            } catch (const std::exception& e) {
                std::cerr << "Error converting value to double: " << value << std::endl;
            }
        }
        unsigned long point_sources_n = data.size()/(2 + Configuration::lambdas.size());
        
        
        vector<vector<double>> point_sources = vector<vector<double>>( point_sources_n , vector<double>( 2 , 0) );
        vector<vector<double>> spectra = vector<vector<double>>( point_sources_n , vector<double>( Configuration::lambdas.size() , 0) );

        for (size_t j = 0; j < point_sources_n; ++j){
            point_sources[j][0] = data[(2+Configuration::lambdas.size())*j];
            point_sources[j][1] = data[(2+Configuration::lambdas.size())*j+1];
            for (size_t i = 0; i < Configuration::lambdas.size(); ++i){
                spectra[j][i] = data[(2+Configuration::lambdas.size())*j+ i+2];
            }
        }
        Ps.push_back(Configuration( point_sources , spectra ));
    }
    file.close();
    
    std::cout << "computing ..";
    std::vector< std::pair<Configuration, double> > Qs_Dlosses(Ps.size());
    std::transform(Ps.begin(), Ps.end(), Qs_Dlosses.begin(), Configuration::get_best_Qstar);
    std::cout << " done" << std::endl;
    
    std::cout << "loading output ..";
    std::ofstream myfile;
    myfile.open ("output.csv");
    
    for (size_t i = 0; i < Qs_Dlosses.size(); ++i){
        for( size_t k = 0; k < Qs_Dlosses[i].first.point_sources.size(); ++k ){
            myfile << Qs_Dlosses[i].first.point_sources[k][0] << "," << Qs_Dlosses[i].first.point_sources[k][1] << ",";
            for (size_t j = 0; j < Configuration::lambdas.size(); ++j){
                myfile << Qs_Dlosses[i].first.spectra[k][j] << ",";
            }
        }
        myfile << Qs_Dlosses[i].second << std::endl;
        /*cout << endl;
        cout << "P = " << Ps[i].to_str() << std::endl;
        cout << "Q = " << Qs_Dlosses[i].first.to_str() << std::endl;
        std::cout << "minsublossQ - lossPQ : " << Qs_Dlosses[i].second << std::endl;*/
    }
    myfile.close();
    std::cout << " done" << std::endl;
    
    return 0;
}
