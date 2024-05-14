//
//  Configuration.hpp
//  advphotobombing
//
//  Created by Drinor Cacaj on 28.06.23.

//
//  Configuration.hpp
//  advphotobombing
//
//  Created by Drinor Cacaj on 28.06.23.
//

#ifndef Configuration_hpp
#define Configuration_hpp

#include "Eigen/Dense"
#include <stdio.h>
#include <vector>
#include <cmath>
#include <iostream>


constexpr std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; i++) {
        result[i] = start + i * step;
    }
    return result;
}

class Configuration {
private:
    void set_ts(){
        std::vector<std::vector<double>> result(lambdas.size(), std::vector<double>(rhos.size(), 0.));
        
        for (size_t j = 0; j < rhos.size(); ++j){
            for (size_t i = 0; i < lambdas.size(); ++i){
                for (size_t k = 0; k < this->point_sources.size(); ++k){
                    result[i][j] += this->spectra[k][i]*pow(sin(M_PI*this->point_sources[k][0]*cos(this->point_sources[k][1]+rhos[j])*lambdas[i]),2)*sin(2*M_PI*q*this->point_sources[k][0]*sin(this->point_sources[k][1]+rhos[j])*lambdas[i]);
                }
            }
        }
        this->ts = result;
    }
    
    
    
public:
    std::vector<std::vector<double>> point_sources;
    std::vector<std::vector<double>> spectra;
    std::vector<std::vector<double>> ts;
    static double q;
    static std::vector<double> rhos;
    static std::vector<double> lambdas;
    static std::vector<double> hom_spectrum;
    static std::vector<double> zero_spectrum;
    
    Configuration(  )
    {
        
    }
    
    Configuration( std::vector<std::vector<double>> point_sources , std::vector<std::vector<double>> spectra)
    : point_sources(point_sources), spectra(spectra)
    {
        this->set_ts();
    }
    
    double get_signal_intensity(){
        double summ = 0;
        for (size_t i = 0; i < lambdas.size() - 1; ++i){
            for (size_t j = 0; j < rhos.size() - 1; ++j){
                summ += (abs(this->ts[i][j]) + abs(this->ts[i][j+1]) + abs(this->ts[i+1][j]) + abs(this->ts[i+1][j+1]))/4;
            }
        }
        return summ/rhos.size()/lambdas.size()*( rhos.back() - rhos[0] )*( lambdas.back()-lambdas[0] );
        
    }
    
    static std::vector<double> get_rt(double semi, double theta, double phi) {
        double x = semi*cos(theta);
        double y = semi*sin(theta)*cos(phi);
        double r = sqrt(pow(x,2)+pow(y,2));
        double t = acos(x/r);
        if (y < 0) t = -t;
        std::vector<double> rt = { r*0.716/1.25 ,t };
        return rt;
    }
    
    Configuration step_opt( Configuration& P, double q = 6){
        std::pair<Eigen::MatrixXd, Eigen::MatrixXd> diffs;
        std::vector<std::vector<double>> new_point_sources = this->point_sources;
        
        
        diffs = hess_grad_N2( P , *this , q);
        
        
        Eigen::MatrixXd step = diffs.second.inverse()* diffs.first;
        
        for (int i = 0; i < new_point_sources.size(); ++i){
            new_point_sources[i][0] += -step(i); //update delta_r
            new_point_sources[i][1] += -step(new_point_sources.size()+i); //update delta_theta
        }
        
        Configuration P_return = Configuration(new_point_sources, this->spectra);
        
        /*if (loss_func(P_return, P) > 1.0001*loss_func(*this, P)){ // if grad/hess does not work use 1step gradient descent
         std::vector<std::vector<double>> newnew_point_sources = this->point_sources;
         Eigen::MatrixXd grad = diffs.first; // first step
         // Find biggest step size
         double pres = 1e-11;
         double step_size = 1e-2;
         double lastloss = loss_func(*this, P);
         double newloss;
         P_return = Configuration(this->point_sources, this->spectra);
         Configuration P_temp = Configuration(this->point_sources, this->spectra);
         int counter = 0;
         do { // do gradient descent until better position
         if (step_size < pres) break;
         newnew_point_sources = P_return.point_sources;
         //std::cout << step_size << std::endl;
         for (int i = 0; i < newnew_point_sources.size(); ++i){
         newnew_point_sources[i][0] += -step_size*grad(i); //update delta_r
         newnew_point_sources[i][1] += -step_size*grad(new_point_sources.size()+i); //update delta_theta
         }
         P_temp = Configuration(newnew_point_sources, this->spectra);
         newloss = loss_func(P_temp, P);
         if(newloss < lastloss){
         lastloss = newloss;
         counter++;
         P_return = P_temp;
         grad = grad_func(P, P_return);
         //std::cout << loss_func(P_return, P) << "\t" << step_size << std::endl;
         step_size = step_size*2;
         } else {
         step_size = step_size/2;
         }
         } while ( counter < 300 );
         //std::cout << loss_func(P_return, P) - loss_func(*this, P) << std::endl;
         //std::cout << "finish gd" << std::endl;
         //std::cout << P_return.to_str() << std::endl;
         //std::cout << "loss : " << loss_func(P_return, P) << " : " << loss_func(*this, P)   << std::endl;
         //std::cout << "non-hess" << std::endl <<  diffs.second << std::endl;
         }*/
        //std::cout << "hess" << std::endl;
        return P_return;
    }
    
    
    
    static double loss_func(const Configuration& P,const Configuration& Q){
        double result = 0;
        for (size_t j = 0; j < rhos.size()-1; ++j){
            for (size_t i = 0; i < lambdas.size()-1; ++i){
                result += (pow(P.ts[i][j]-Q.ts[i][j],2) + pow(P.ts[i+1][j]-Q.ts[i+1][j],2) + pow(P.ts[i][j+1]-Q.ts[i][j+1],2) + pow(P.ts[i+1][j+1]-Q.ts[i+1][j+1],2))/4;
            }
        }
        return sqrt(result/rhos.size()/lambdas.size()*( rhos.back() - rhos[0] )*( lambdas.back()-lambdas[0] ));
    }
    
    static std::vector<std::vector<double>> eraseRow(std::vector<std::vector<double>> vec2D, size_t rowIndex) {
        vec2D.erase(vec2D.begin() + rowIndex);
        return vec2D;
    }
    
    static std::vector<std::vector<double>> keepRow(std::vector<std::vector<double>> vec2D, size_t rowIndex) {
        return {vec2D[rowIndex]};
    }
    
    static double loss_func_empty(const Configuration& Q){
        double result = 0;
        for (size_t j = 0; j < rhos.size()-1; ++j){
            for (size_t i = 0; i < lambdas.size()-1; ++i){
                result += (pow(Q.ts[i][j],2) + pow(Q.ts[i+1][j],2) + pow(Q.ts[i][j+1],2) + pow(Q.ts[i+1][j+1],2))/4;
            }
        }
        return sqrt(result/rhos.size()/lambdas.size()*( rhos.back() - rhos[0] )*( lambdas.back()-lambdas[0] ));
    }
    
    Configuration lin_regr( const Configuration& P ) {
        //std::cout << "pssize begin : " << this->point_sources.size();
        size_t point_sources_n = this->point_sources.size();
        std::vector<std::vector<std::vector<double>>> F_space;
        for (size_t k = 0; k < point_sources_n; ++k){
            F_space.push_back( Configuration({this->point_sources[k]}, {hom_spectrum}).ts );
        }
        
        //update spectra
        std::vector<std::vector<double>> coeffs( point_sources_n, std::vector<double>(lambdas.size(),0.) );
        for (int k = 0; k < lambdas.size() ; ++k){
            Eigen::MatrixXd eigenMatrix(point_sources_n, rhos.size());
            
            for (int i = 0; i < point_sources_n; ++i) {
                for (int j = 0; j < rhos.size(); ++j) {
                    eigenMatrix(i, j) = F_space[i][k][j];
                }
            }
            Eigen::VectorXd eigenVector(rhos.size());
            for (int i = 0; i < rhos.size(); ++i) {
                eigenVector(i) = P.ts[k][i];
            }
            Eigen::MatrixXd transposedMatrix = eigenMatrix.transpose();
            Eigen::MatrixXd PhitPhi = eigenMatrix* transposedMatrix;
            //std::cout << PhitPhi << std::endl;
            //std::cout << "det : " << abs(PhitPhi.determinant() - 0.) << std::endl;
            if ( PhitPhi.determinant() == 0 ){ // here because this occurs very rarely
                std::cout << "det = 0" << std::endl;
                std::vector<int> subpointsources = getIndep_point_sources(PhitPhi);
                std::vector<std::vector<double>> newpoint_sources(subpointsources.size(),std::vector<double>(2,0));
                std::vector<std::vector<double>> newspectra(subpointsources.size(),std::vector<double>(lambdas.size(),0));
                for (int i = 0; i < subpointsources.size(); ++i) {
                    newpoint_sources[i] = P.point_sources[subpointsources[i]];
                    newspectra[i] = P.spectra[subpointsources[i]];
                }
                return Configuration( newpoint_sources, newspectra ).lin_regr( P );
            }
            Eigen::VectorXd result = ( PhitPhi ).inverse() * (eigenMatrix * eigenVector);
            
            //std::cout << eigenMatrix* transposedMatrix << std::endl;
            //std::cout << "  " << std::endl;
            
            for (int i = 0; i < point_sources_n; ++i) {
                coeffs[i][k] = result(i);
            }
        }
        
        // update point sources
        double summ;
        std::vector<std::vector<double>> point_sources = this->point_sources;
        for (int i = 0; i < point_sources_n; ++i) {
            summ = 0;
            for (int k = 0; k < lambdas.size() ; ++k){
                summ += std::signbit(coeffs[i][k]);
            }
            //std::cout << "summ : " << summ << std::endl;
            if (summ == 0){
                //std::cout << "phys" << std::endl;
            }
            else if ( summ > 0 && summ < lambdas.size() ){
                //std::cout << "non-phys" << std::endl;
                return *this; // non physical solution
            }
            else if ( summ == lambdas.size() ){
                //std::cout << "antiphys" << std::endl;
                point_sources[i][1] = fmod(point_sources[i][1] + M_PI, 2*M_PI); // opposite side solution
                for (int k = 0; k < lambdas.size() ; ++k){
                    coeffs[i][k] = -coeffs[i][k]; //positive flux
                    //std::cout << coeffs[i][k] << std::endl;
                }
            }
        }
        
        return Configuration( point_sources, coeffs );
    }
    
    
    static std::vector<int> getIndep_point_sources(const Eigen::MatrixXd& matrix) {
        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(matrix);
        std::vector<int> keptRows = {0};
        
        for (int i = 1; i < matrix.rows(); ++i){
            Eigen::MatrixXd submatrix(keptRows.size(),keptRows.size());
            for (int n = 1; n < keptRows.size(); ++n) {
                for (int m = 1; m < keptRows.size(); ++m)  {
                    submatrix(n,m) = matrix(keptRows[n],keptRows[m]);
                }
                Eigen::FullPivLU<Eigen::MatrixXd> sub_lu_decomp(matrix);
                if (sub_lu_decomp.rank() == sub_lu_decomp.rows()){
                    keptRows.push_back(i);
                }
            }
        }
        return keptRows;
    }
    
    std::string to_str(){
        std::string text;
        for(int i = 0; i < this->point_sources.size(); ++i){
            text += "[" + std::to_string(this->point_sources[i][0]) + ","
            + std::to_string(this->point_sources[i][1]) + ","
            + std::to_string(this->spectra[i][24]) + "..."
            + std::to_string(this->spectra[i].back()) + "]";
        }
        return text;
    }
    
    static Configuration concatenate(const Configuration& A, const Configuration& B){
        std::vector<std::vector<double>> point_sources= A.point_sources;
        std::vector<std::vector<double>> spectra = A.spectra;
        
        for(int i =0; i < B.point_sources.size(); i++){
            point_sources.push_back(B.point_sources[i]);
            spectra.push_back(B.spectra[i]);
        }
        return Configuration( point_sources, spectra );
    }
    
    static double min_loss_Q_subP(Configuration& P, Configuration& Q_star ){
        std::vector<std::vector<int>> subsets;
        int count = pow(2, P.point_sources.size())-1;
        for (int i = 1; i < count; i++) {
            std::vector<int> subset;
            for (int j = 0; j < P.point_sources.size(); j++) {
                if ((i & (1 << j)) != 0) subset.push_back(j);
            }
            subsets.push_back(subset);
        }
        
        std::vector<double> losses;
        
        losses.push_back(loss_func_empty(Q_star));
        //std::cout << "\n--start" << std::endl;
        for(auto subset : subsets){
            std::vector<std::vector<double>> sub_point_sources;
            std::vector<std::vector<double>> sub_spectra;
            for(int i :subset){
                sub_point_sources.push_back(P.point_sources[i]);
                sub_spectra.push_back(P.spectra[i]);
            }
            Configuration subP = Configuration(sub_point_sources, sub_spectra);
            
            //std::cout << "subP = " << subP.to_str() << std::endl;
            //std::cout << "Q_star = " << Q_star.to_str() << std::endl;
            losses.push_back(loss_func( subP , Q_star));
            //std::cout << "loss : " << losses.back() << std::endl;
        }
        //std::cout << "--end" << std::endl;
        
        int minindex = 0;
        for( int i = 1; i < losses.size(); ++i ){
            if (losses[minindex] > losses[i] ){
                minindex = i;
            }
        }
        //std::cout << "minindex : " << minindex << std::endl;
        //std::cout << "minloss : " << losses[minindex] << std::endl;
        return losses[minindex];
    }
    
    static double min_loss_Q_subP_Q(Configuration& P, Configuration& Q_star){
        return pow(min_loss_Q_subP( P, Q_star ),2) - pow(loss_func(P, Q_star),2);
    }
    
    static std::pair<Configuration, double> get_best_Qstar(Configuration& P){
        
        if ( P.point_sources.size() == 1 ){
            return std::make_pair(P, -1);
        }
        
        std::vector<std::vector<int>> subsets;
        int count = pow(2, P.point_sources.size());
        for (int i = 1; i < count; i++) {
            std::vector<int> subset;
            for (int j = 0; j < P.point_sources.size(); j++) {
                if ((i & (1 << j)) != 0) subset.push_back(j);
            }
            //if (subset.size() > 1) subsets.push_back(subset); // all check
            if (subset.size() == 2) subsets.push_back(subset); // 2by2 check
        }
        
        std::vector<double> Dlosses;
        std::vector<Configuration> Q_stars;
        
        //std::cout << "\n--start" << std::endl;
        for(auto subset : subsets){
            std::vector<std::vector<double>> sub_point_sources;
            std::vector<std::vector<double>> sub_spectra;
            for(int i :subset){
                sub_point_sources.push_back(P.point_sources[i]);
                sub_spectra.push_back(P.spectra[i]);
            }
            Configuration subP = Configuration(sub_point_sources, sub_spectra);
            //std::cout << "subP : " << subP.to_str() << std::endl;
            Q_stars.push_back(get_Q_star(subP));
            //std::cout << "Q_star : " << Q_stars.back().to_str() << std::endl;
            //std::cout << "subP = " << subP.to_str() << std::endl;
            Dlosses.push_back(min_loss_Q_subP_Q( subP , Q_stars.back()));
            //std::cout << "Dloss : " << Dlosses.back() << std::endl;
        }
        
        int maxindex = 0;
        for( int i = 1; i < Dlosses.size(); ++i ){
            if (Dlosses[maxindex] < Dlosses[i] ){
                maxindex = i;
            }
        }
        std::cout << Dlosses[maxindex] << std::endl;
        return std::make_pair(Q_stars[maxindex], Dlosses[maxindex]);
    }
    
    static Configuration get_Q_star(Configuration& P){
        std::vector<Configuration> P_tests;
        std::vector<double> losses_pl;
        
        //std::cout << std::endl << P.to_str() << std::endl;
        if( P.point_sources.size() == 2){
            for (size_t i = 0; i < P.point_sources.size(); ++i){//N=2
                Configuration P_test = Configuration({ eraseRow(P.point_sources, i) }, eraseRow(P.spectra, i));
                
                for (size_t j = 0; j < 7; ++j){
                    //std::cout << "step" << j << std::endl;
                    P_test = P_test.lin_regr(P);
                    //std::cout << "lin done" << std::endl;
                    P_test = P_test.step_opt(P);
                    //std::cout << "step done" << std::endl;
                }
                //std::cout << "P_test 2 : " << P_test.to_str() << std::endl;
                P_tests.push_back(P_test);
                //std::cout << P_test.to_str() << std::endl;
                losses_pl.push_back(loss_func(P, P_test));
                //std::cout << "loss_pl : " << losses_pl.back() << std::endl;
            }
        } else { //recursive choice for Q_0 computing subPs
            for (size_t i = 0; i < P.point_sources.size(); ++i){
                Configuration P_test = Configuration({ eraseRow(P.point_sources, i) }, eraseRow(P.spectra, i));
                Configuration P_rest = Configuration({ keepRow(P.point_sources, i) }, keepRow(P.spectra, i));
                //std::cout << "P_test : " << P_test.to_str() << std::endl;
                //std::cout << "P_rest : " << P_rest.to_str() << std::endl;
                P_test = concatenate(get_Q_star(P_test), P_rest );
                //std::cout << "P_test 3 : " << P_test.to_str() << std::endl;
                
                for (size_t j = 0; j < 5; ++j){
                    //std::cout << "step" << j << std::endl;
                    P_test = P_test.lin_regr(P);
                    //std::cout << "lin done" << std::endl;
                    P_test = P_test.step_opt(P);
                    //std::cout << "step done" << std::endl;
                }
                P_tests.push_back(P_test);
                //std::cout << P_test.to_str() << std::endl;
                losses_pl.push_back(loss_func(P, P_test));
                //std::cout << "loss_pl : " << losses_pl.back() << std::endl;
            }
        }
        
        
        int minindex = 0;
        for( int i = 1; i < losses_pl.size(); ++i ){
            if (losses_pl[minindex] > losses_pl[i] ){
                minindex = i;
            }
        }
        //std::cout << losses_pl[minindex] << std::endl;
        //std::cout << "minindex : " << minindex << std::endl;
        return P_tests[minindex];
    }
    
    static Eigen::MatrixXd grad_func(Configuration& P, Configuration& Q, double q = 6){
        Eigen::MatrixXd grad = grad_N2(P, Q, q);
        return grad;
    }
    
    static Eigen::MatrixXd grad_N2(Configuration& P, Configuration& Q, double q = 6){
        std::vector<std::vector<std::vector<double>>> grad_vals(2, std::vector<std::vector<double>>(rhos.size(),std::vector<double>( lambdas.size(), 0)));
        
        for (size_t j = 0; j < rhos.size() ; ++j){
            double A_1=cos(Q.point_sources[0][1]+rhos[j]);
            double A_2=sin(Q.point_sources[0][1]+rhos[j]);
            for (size_t i = 0; i < lambdas.size(); ++i){
                
                double A_0=1.0/lambdas[i];
                double A_3=M_PI*A_0*A_1*Q.point_sources[0][0];
                double A_4=2*M_PI*A_0*A_2*Q.point_sources[0][0]*q;
                double A_5=sin(A_3);
                double A_6=sin(A_4);
                double A_7=pow(A_5, 2);
                double A_8=cos(A_3);
                double A_9=cos(A_4);
                double A_10=pow(sin(M_PI*A_0*P.point_sources[0][0]*cos(P.point_sources[0][1]+rhos[j])), 2);
                double A_11=sin(2*M_PI*A_0*P.point_sources[0][0]*q*sin(P.point_sources[0][1]+rhos[j]));
                double A_12=pow(sin(M_PI*A_0*P.point_sources[1][0]*cos(P.point_sources[1][1]+rhos[j])), 2);
                double A_13=sin(2*M_PI*A_0*P.point_sources[1][0]*q*sin(P.point_sources[1][1]+rhos[j]));
                double A_19=A_10*A_11*P.spectra[0][i] + A_12*A_13*P.spectra[1][i] - Q.spectra[0][i]*A_6*A_7;
                double A_30=-4*M_PI*A_0*A_1*Q.spectra[0][i]*A_5*A_6*A_8 - 4*M_PI*A_0*Q.spectra[0][i]*A_2*q*A_7*A_9;
                
                
                grad_vals[0][j][i] = A_19*A_30;
                grad_vals[1][j][i] = 4*Q.spectra[0][i]*A_19*(M_PI*A_0*A_2*Q.point_sources[0][0]*A_5*A_6*A_8 - q*A_3*A_7*A_9);
                
            }
        }
        
        Eigen::MatrixXd grad_val( 2, 1 );
        
        //trapeze
        for (size_t j = 0; j < rhos.size()-1 ; ++j){
            for (size_t i = 0; i < lambdas.size()-1; ++i){
                grad_val(0, 0) += ( grad_vals[0][j][i] + grad_vals[0][j][i+1] + grad_vals[0][j+1][i] + grad_vals[0][j+1][i+1] )/4;
                grad_val(1, 0) += ( grad_vals[1][j][i] + grad_vals[1][j][i+1] + grad_vals[1][j+1][i] + grad_vals[1][j+1][i+1] )/4;
            }
        }
        
        return grad_val ;
    }
    
    static std::pair<Eigen::MatrixXd, Eigen::MatrixXd > hess_grad_N2(Configuration& P, Configuration& Q, double q = 6){
        std::vector<std::vector<std::vector<double>>> grad_vals(2, std::vector<std::vector<double>>(rhos.size(),std::vector<double>( lambdas.size(), 0)));
        std::vector<std::vector<std::vector<std::vector<double>>>> hess_vals(2, std::vector<std::vector<std::vector<double>>>(2, std::vector<std::vector<double>>(rhos.size(),std::vector<double>(lambdas.size(),0.))));
        
        for (size_t j = 0; j < rhos.size() ; ++j){
            double A_1=cos(Q.point_sources[0][1]+rhos[j]);
            double A_2=sin(Q.point_sources[0][1]+rhos[j]);
            for (size_t i = 0; i < lambdas.size(); ++i){
                
                double A_0=1.0/lambdas[i];
                double A_3=M_PI*A_0*A_1*Q.point_sources[0][0];
                double A_4=2*M_PI*A_0*A_2*Q.point_sources[0][0]*q;
                double A_5=sin(A_3);
                double A_6=sin(A_4);
                double A_7=pow(A_5, 2);
                double A_8=cos(A_3);
                double A_9=cos(A_4);
                double A_10=pow(sin(M_PI*A_0*P.point_sources[0][0]*cos(P.point_sources[0][1]+rhos[j])), 2);
                double A_11=sin(2*M_PI*A_0*P.point_sources[0][0]*q*sin(P.point_sources[0][1]+rhos[j]));
                double A_12=pow(sin(M_PI*A_0*P.point_sources[1][0]*cos(P.point_sources[1][1]+rhos[j])), 2);
                double A_13=sin(2*M_PI*A_0*P.point_sources[1][0]*q*sin(P.point_sources[1][1]+rhos[j]));
                double A_14=pow(M_PI, 2);
                double A_15=pow(A_0, 2);
                double A_19=A_10*A_11*P.spectra[0][i] + A_12*A_13*P.spectra[1][i] - Q.spectra[0][i]*A_6*A_7;
                double A_20=pow(A_8, 2);
                double A_21=pow(A_2, 2);
                double A_22=pow(q, 2);
                double A_25=pow(A_1, 2);
                double A_26=pow(Q.point_sources[0][0], 2);
                double A_29=2*M_PI*A_0*Q.spectra[0][i]*A_2*Q.point_sources[0][0]*A_5*A_6*A_8 - 2*Q.spectra[0][i]*q*A_3*A_7*A_9;
                double A_30=-4*M_PI*A_0*A_1*Q.spectra[0][i]*A_5*A_6*A_8 - 4*M_PI*A_0*Q.spectra[0][i]*A_2*q*A_7*A_9;
                
                
                grad_vals[0][j][i] = A_19*A_30;
                
                grad_vals[1][j][i] = 4*Q.spectra[0][i]*A_19*(M_PI*A_0*A_2*Q.point_sources[0][0]*A_5*A_6*A_8 - q*A_3*A_7*A_9);
                
                hess_vals[0][0][j][i] = -2*Q.spectra[0][i]*(M_PI*A_0*A_1*A_30*A_5*A_6*A_8 + M_PI*A_0*A_2*q*A_30*A_7*A_9 + 8*A_1*A_14*A_15*A_19*A_2*q*A_5*A_8*A_9 + 2*A_14*A_15*A_19*A_20*A_25*A_6 - 4*A_14*A_15*A_19*A_21*A_22*A_6*A_7 - 2*A_14*A_15*A_19*A_25*A_6*A_7);
                
                hess_vals[1][0][j][i] = -8*M_PI*A_0*A_1*Q.spectra[0][i]*A_19*q*A_3*A_5*A_8*A_9 - 4*M_PI*A_0*A_1*Q.spectra[0][i]*A_19*q*A_7*A_9 + 4*M_PI*A_0*Q.spectra[0][i]*A_19*A_2*A_20*A_3*A_6 + 8*M_PI*A_0*Q.spectra[0][i]*A_19*A_2*A_22*A_3*A_6*A_7 - 4*M_PI*A_0*Q.spectra[0][i]*A_19*A_2*A_3*A_6*A_7 + 4*M_PI*A_0*Q.spectra[0][i]*A_19*A_2*A_4*A_5*A_8*A_9 + 4*M_PI*A_0*Q.spectra[0][i]*A_19*A_2*A_5*A_6*A_8 + A_29*A_30;
                
                hess_vals[1][1][j][i] = 2*Q.spectra[0][i]*(2*M_PI*A_0*A_2*Q.point_sources[0][0]*A_29*A_5*A_6*A_8 - 2*A_14*A_15*A_19*A_20*A_21*A_26*A_6 + 2*A_14*A_15*A_19*A_21*A_26*A_6*A_7 + 4*A_19*A_22*pow(A_3, 2)*A_6*A_7 + 4*A_19*A_3*A_4*A_5*A_8*A_9 + 2*A_19*A_3*A_5*A_6*A_8 + A_19*A_4*A_7*A_9 - 2*q*A_29*A_3*A_7*A_9);
                
            }
        }
        
        Eigen::MatrixXd grad_val( 2, 1 );
        Eigen::MatrixXd hess_val( 2, 2 );
        
        //trapeze
        for (size_t j = 0; j < rhos.size()-1 ; ++j){
            for (size_t i = 0; i < lambdas.size()-1; ++i){
                grad_val(0, 0) += ( grad_vals[0][j][i] + grad_vals[0][j][i+1] + grad_vals[0][j+1][i] + grad_vals[0][j+1][i+1] )/4;
                grad_val(1, 0) += ( grad_vals[1][j][i] + grad_vals[1][j][i+1] + grad_vals[1][j+1][i] + grad_vals[1][j+1][i+1] )/4;
                hess_val(0, 0) += ( hess_vals[0][0][j][i] + hess_vals[0][0][j][i+1] + hess_vals[0][0][j+1][i] + hess_vals[0][0][j+1][i+1] )/4;
                hess_val(1, 0) +=( hess_vals[1][0][j][i] + hess_vals[1][0][j][i+1] + hess_vals[1][0][j+1][i] + hess_vals[1][0][j+1][i+1] )/4;
                hess_val(1, 1) +=( hess_vals[1][1][j][i] + hess_vals[1][1][j][i+1] + hess_vals[1][1][j+1][i] + hess_vals[1][1][j+1][i+1] )/4;
            }
        }
        hess_val(0,1) = hess_val(1,0);
        
        return std::make_pair(grad_val, hess_val );
    }
};
    
#endif
