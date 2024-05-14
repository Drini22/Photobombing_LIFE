//
//  Airy.hpp
//  Photobombing
//
//  Created by Drinor Cacaj on 30.07.23.
//

#ifndef Airy_hpp
#define Airy_hpp

#include "Eigen/Dense"
#include <stdio.h>
#include <vector>
#include <cmath>
#include <iostream>



class Airy {
private:
    void set_ts(){
        std::vector<std::vector<double>> result(delta_xs.size(), std::vector<double>(delta_ys.size(), 0.));
        double dist;
        for (size_t k = 0; k < this->point_sources.size(); ++k){
            for (int i = 0; i < delta_xs.size(); i++){
                for (int j = 0; j < delta_xs.size(); j++){
                    dist = std::sqrt(pow(delta_xs[i] - this->point_sources[k][0],2) - pow(delta_ys[j] - this->point_sources[k][1],2));
                    result[i][j] += pow(2*::j1( dist )/dist,2);
                }
            }
        }
    }
    
    
    
public:
    std::vector<std::vector<double>> point_sources;
    std::vector<double> spectra;
    std::vector<std::vector<double>> ts;
    static std::vector<double> delta_xs;
    static std::vector<double> delta_ys;
    
    Airy(  )
    {
        
    }
    
    Airy( std::vector<std::vector<double>> point_sources , std::vector<double> spectra)
    : point_sources(point_sources), spectra(spectra)
    {
        this->set_ts();
    }
    
    
    static double loss_func(const Airy& P,const Airy& Q){
        double result = 0;
        for (size_t j = 0; j < delta_xs.size()-1; ++j){
            for (size_t i = 0; i < delta_ys.size()-1; ++i){
                result += (pow(P.ts[i][j]-Q.ts[i][j],2) + pow(P.ts[i+1][j]-Q.ts[i+1][j],2) + pow(P.ts[i][j+1]-Q.ts[i][j+1],2) + pow(P.ts[i+1][j+1]-Q.ts[i+1][j+1],2))/4;
            }
        }
        return sqrt(result/delta_xs.size()/delta_ys.size()*( delta_xs.back() - delta_xs[0] )*( delta_ys.back()-delta_ys[0] ));
    }
    
    
    static double loss_func_empty(const Airy& Q){
        double result = 0;
        for (size_t j = 0; j < delta_xs.size()-1; ++j){
            for (size_t i = 0; i < delta_ys.size()-1; ++i){
                result += (pow(Q.ts[i][j],2) + pow(Q.ts[i+1][j],2) + pow(Q.ts[i][j+1],2) + pow(Q.ts[i+1][j+1],2))/4;
            }
        }
        return sqrt(result/delta_xs.size()/delta_ys.size()*( delta_xs.back() - delta_xs[0] )*( delta_ys.back()-delta_ys[0] ));
    }
    
    
    
    std::string to_str(){
        std::string text;
        for(int i = 0; i < this->point_sources.size(); ++i){
            text += "[" + std::to_string(this->point_sources[i][0]) + ","
            + std::to_string(this->point_sources[i][1]);
        }
        return text;
    }
    
    
    
    static Airy get_Q_star_N2(Airy& P){
        std::vector<Airy> P_tests;
        std::vector<double> losses_pl;
        
        
        
        
        int minindex = 0;
        for( int i = 1; i < losses_pl.size(); ++i ){
            if (losses_pl[minindex] > losses_pl[i] ){
                minindex = i;
            }
        }
        
        //std::cout << "minindex : " << minindex << std::endl;
        return P_tests[minindex];
    }
    
    
    static std::pair<double, double> hess_grad_N2(Airy& P, Airy& Q, double q = 6){
        
        std::vector<double> grad_vals( delta_xs.size() , 0 );
        std::vector<double> hess_vals( delta_xs.size() , 0 );
        /*
        for (size_t j = 0; j < delta_xs.size() ; ++j){
            grad_vals[j] = (-2*Q.spectra[0]*(::j0( Q.point_sources[0][0] - delta_xs[j]) - ::j2( Q.point_sources[0][0] - delta_xs[j]))*::j1( Q.point_sources[0][0] - delta_xs[j])/pow(Q.point_sources[0][0] - delta_xs[j], 2) + 4*Q.spectra[0]*pow(::j1( Q.point_sources[0][0] - delta_xs[j]), 2)/pow(Q.point_sources[0][0] - delta_xs[j], 3))*(P.spectra[0]*pow(::j1( P.point_sources[0][0] - delta_xs[j]), 2)/pow(P.point_sources[0][0] - delta_xs[j], 2) + P.spectra[1]*pow(::j1( P.point_sources[1][0] - delta_xs[j]), 2)/pow(P.point_sources[1][0] - delta_xs[j], 2) - Q.spectra[0]*pow(::j1( Q.point_sources[0][0] - delta_xs[j]), 2)/pow(Q.point_sources[0][0] - delta_xs[j], 2));
            
            hess_vals[j] = (-2*Q.spectra[0]*(besselj(0, Q.point_sources[0][0] - delta_xs[j]) - besselj(2, Q.point_sources[0][0] - delta_xs[j]))*besselj(1, Q.point_sources[0][0] - delta_xs[j])/pow(Q.point_sources[0][0] - delta_xs[j], 2) + 4*Q.spectra[0]*pow(besselj(1, Q.point_sources[0][0] - delta_xs[j]), 2)/pow(Q.point_sources[0][0] - delta_xs[j], 3))*(-Q.spectra[0]*(besselj(0, Q.point_sources[0][0] - delta_xs[j]) - besselj(2, Q.point_sources[0][0] - delta_xs[j]))*besselj(1, Q.point_sources[0][0] - delta_xs[j])/pow(Q.point_sources[0][0] - delta_xs[j], 2) + 2*Q.spectra[0]*pow(besselj(1, Q.point_sources[0][0] - delta_xs[j]), 2)/pow(Q.point_sources[0][0] - delta_xs[j], 3)) + (P.spectra[0]*pow(besselj(1, P.point_sources[0][0] - delta_xs[j]), 2)/pow(P.point_sources[0][0] - delta_xs[j], 2) + P.spectra[1]*pow(besselj(1, P.point_sources[1][0] - delta_xs[j]), 2)/pow(P.point_sources[1][0] - delta_xs[j], 2) - Q.spectra[0]*pow(besselj(1, Q.point_sources[0][0] - delta_xs[j]), 2)/pow(Q.point_sources[0][0] - delta_xs[j], 2))*(-2*Q.spectra[0]*((1.0/2.0)*besselj(0, Q.point_sources[0][0] - delta_xs[j]) - 1.0/2.0*besselj(2, Q.point_sources[0][0] - delta_xs[j]))*(besselj(0, Q.point_sources[0][0] - delta_xs[j]) - besselj(2, Q.point_sources[0][0] - delta_xs[j]))/pow(Q.point_sources[0][0] - delta_xs[j], 2) - 2*Q.spectra[0]*(-3.0/2.0*besselj(1, Q.point_sources[0][0] - delta_xs[j]) + (1.0/2.0)*besselj(3, Q.point_sources[0][0] - delta_xs[j]))*besselj(1, Q.point_sources[0][0] - delta_xs[j])/pow(Q.point_sources[0][0] - delta_xs[j], 2) + 8*Q.spectra[0]*(besselj(0, Q.point_sources[0][0] - delta_xs[j]) - besselj(2, Q.point_sources[0][0] - delta_xs[j]))*besselj(1, Q.point_sources[0][0] - delta_xs[j])/pow(Q.point_sources[0][0] - delta_xs[j], 3) - 12*Q.spectra[0]*pow(besselj(1, Q.point_sources[0][0] - delta_xs[j]), 2)/pow(Q.point_sources[0][0] - delta_xs[j], 4));
        }
        */
        double grad_val = 0;
        double hess_val = 0;
        
        //trapeze
        for (size_t j = 0; j < delta_xs.size()-1 ; ++j){
            grad_val += ( grad_vals[j] + grad_vals[j+1] )/2;
            hess_val += ( hess_vals[j] + hess_vals[j+1] )/2;
        }
        //std::cout << hess_val << std::endl;
        //std::cout << grad_val << std::endl;
        
        return std::make_pair(grad_val, hess_val );
    }
};
    

#endif /* Airy_hpp */
