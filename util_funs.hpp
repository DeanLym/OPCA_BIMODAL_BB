/*
 * util_funs.hpp
 *
 *  Created on: Mar 5, 2016
 *      Author: yiminliu
 */

#ifndef UTIL_FUNS_HPP_
#define UTIL_FUNS_HPP_

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include "opca_bimodal.h"
#include <exception>
#include <stdexcept>
#include <boost/math/special_functions/erf.hpp>

using namespace::boost::math;


void LogOpcaError(const char* file_name, vector<double> xi, int n){
	ofstream ofs;
	ofs.open(file_name,std::ofstream::out | std::ofstream::app);
	ofs << "O-PCA returns false at: " << std::endl;
	for(int i=0;i<n;i++)
		ofs << xi[i] << std::endl;
	ofs<<endl;
	ofs.close();
}

void TranformUniform2Normal(int l , vector<double> &m){
	for(int i = 0; i < l ; i++){
		m[i] = 1.4142135623 * erf_inv( 2*m[i] - 1 );
	}
}

void TranformNormal2Uniform(int l , vector<double> &m){
	for(int i = 0; i < l ; i++){
		m[i] = 0.5*(erf(m[i]/1.4142135623)+1);
	}
}

double TransformUniform2Normal_2(double u,double mu, double sigma){
	double x;
	x = mu + 1.4142135623 * sigma * erf_inv( 2*u - 1 );
	return x;
}

void GeneratePerm(int Nc, double *x, double *perm){
	for(int i=0; i < Nc ; i++){
		perm[i] = exp(x[i]);
	}
}

void CheckBound(double ub, double lb, int l , vector<double> &m){
	for(int i=0; i < l ; i++){
		m[i] = m[i]>ub?ub:m[i];
		m[i] = m[i]<lb?lb:m[i];
	}
}

string num2str(int x){
	stringstream ss;
	string res;
	ss << x;
	ss >> res;
	return res;
}

int str2num(string x){
	stringstream ss;
	int res;
	ss << x;
	ss >> res;
	return res;
}

int char2num(char* x){
	string y(x);
	stringstream ss;
	int res;
	ss << y;
	ss >> res;
	return res;
}

void ReadData(string fn, int n, double* data){
	//cout << "Saving data to " << fn << endl;
	ifstream in;
	in.open(fn.c_str());

	if (in.is_open()){
		for (int i = 0; i < n; i++){
			in >> data[i];
		}
		in.close();
	}else{
		throw runtime_error("Can not write file while saving data.");
	}
}


void SaveData(string fn, int n, const double* data){
	//cout << "Saving data to " << fn << endl;
	ofstream out;
	out.open(fn.c_str());

	if (out.is_open()){
		for (int i = 0; i < n; i++){
			out << data[i] << endl;
		}
		out.close();
	}else{
		throw runtime_error("Can not write file while saving data.");
	}
}

#endif /* UTIL_FUNS_HPP_ */
