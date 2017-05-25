#include <stdio.h>
#include <stdlib.h>
#include <Rcpp.h>
#include <vector>
#include <math.h>
#include <string>
#include <fstream>
#include <omp.h>
#include "tables_chi_correlation.h"
#include <boost/math/special_functions/gamma.hpp>

void run_test (int** cross_table, size_t phen_rows_num, size_t gen_cols_num, 
				std::vector<double> *stat_vec, std::vector<double> *p_value_vec, std::vector <int> *degree_freedom_vec, 
				std::vector<std::string> *type_vec, std::vector<bool> *warning, 
				std::vector <double> *correlation_coefficient, int iteration){
	
	double result_stat;
	int deg_freedom;
	double p_value;
	double correlation;	
	std::vector <int> row_sums;
	std::vector <int> col_sums;
	//get sums of cols and rows, the sum of all elements	
	row_sums = get_row_sums(cross_table, phen_rows_num, gen_cols_num);
	col_sums = get_col_sums(cross_table, phen_rows_num, gen_cols_num);
	int total_sum = get_total_sum(row_sums);	
	
	bool smaller5 = lessThan5(&row_sums, &col_sums, total_sum);
	
												
	/*if table has expected values < 5, number of columns and rows == 2, then we calculate fisher exact test's value (p-value)
		else we start chi square test, get chi2-test-staticstic and then p-value from chi2-statistic and calculated degree of freedom */
	if(smaller5 && phen_rows_num == 2 && gen_cols_num ==2){
		result_stat = cross_table[0][0];
		p_value =  fisher_exact_test(cross_table);
		deg_freedom = -1;
		type_vec->at(iteration) = "Fisher";
	}
	else {
		result_stat =  chi2_test(cross_table, phen_rows_num, gen_cols_num, &deg_freedom, &row_sums, &col_sums, total_sum);
		p_value =  get_p_value(result_stat, deg_freedom);
		type_vec->at(iteration) = "Chi2";
	}
	correlation = get_corr_coefficient(cross_table, phen_rows_num, gen_cols_num, &row_sums, &col_sums, total_sum);
	
	//push results 
	stat_vec->at(iteration) = result_stat;
	p_value_vec->at(iteration) = p_value;
	warning->at(iteration) = smaller5;
	degree_freedom_vec->at(iteration) = deg_freedom;
	correlation_coefficient->at(iteration) = correlation;
}

double chi2_test(int** cross_table, size_t phen_rows_num, size_t gen_cols_num, int* deg_freedom, std::vector <int> *row_sums, std::vector <int> *col_sums, int total_sum){
	
	//counters of zero-columns and zero-rows
	int zero_cols = 0;
	int zero_rows = 0;
	
	for (auto i = col_sums->begin(); i != col_sums->end(); i++){
		if (*i == 0) 
			zero_cols++;
	}
	
	for (auto i = row_sums->begin(); i != row_sums->end(); i++){
		if(*i == 0)
			zero_rows++;
	}
		/*DEGREE OF FREEDOM*/
	if (zero_cols != 0){
		*deg_freedom = (phen_rows_num - 1) * (gen_cols_num - zero_cols-1);
	} 
	else {
		*deg_freedom = (phen_rows_num-1) * (gen_cols_num-1);
	}
	
	//if we have table 1xN or Mx1 =>  degree of freedom == 0 => calculating of chi2-statistic is not correct 
	if (gen_cols_num - zero_cols == 1 || phen_rows_num - zero_rows == 1)
		return 1.0/0.0;
		
		/*CALCULATING CHI2-TEST-STATICTIC*/
		
	double expectation = 0;
	double difference;
	double chisquare = 0;	
	
	
	for (size_t i = 0; i < phen_rows_num; ++i){
		for (size_t j = 0; j < gen_cols_num; ++j){
			expectation = (double)row_sums->at(i) * col_sums->at(j) / (double)total_sum;
			if (expectation != 0){
				difference = (double)cross_table[i][j] - expectation;
				chisquare += difference*difference / expectation;
				
			
			}
		}
	}
	
	
	return chisquare;
}

double fisher_exact_test(int** cross_table){
	/*table |a|b|
			|c|d|	*/	
	int a = cross_table[0][0];
	int b = cross_table[0][1];
	int c = cross_table[1][0];
	int d = cross_table[1][1];
	
	int n = a+b+c+d;
	
	
	double* logFactorials = new double[n+1];
	//get log(i!)
	get_logFactorials(logFactorials, n);
	
	//log(p) of the original table
	double log_p_firstTable = logHypergeometricProb(logFactorials, a,b,c,d);
	
	double p_value = 0;
	//calculate p-values of all table's versions, and then sum all of them
	for (int v = 0; v <= n; v++){
		if(	a+b-v >= 0 && a+c-v >= 0 && d-a+v >=0){
			double log_p_x = logHypergeometricProb(logFactorials, v, a+b-v, a+c-v, d-a+v);
			if (log_p_x <= log_p_firstTable)
				p_value += exp( log_p_x - log_p_firstTable) ;
				
		}
	}
	
	double log_pValue = log_p_firstTable + log(p_value);
	delete[] logFactorials;
	
	return exp(log_pValue);
	
	
}

void get_logFactorials(double* logFactorials, int n){

	logFactorials[0] = 0;
	
	for (size_t i = 1; i < n+1; i++){
		logFactorials[i] = logFactorials[i-1] + log( (double) i );
	}
}

double logHypergeometricProb(double* logFactorials, int a, int b, int c, int d){
	//hypergeometric distribution formula 	 
	double result = logFactorials[a+b] + logFactorials[c+d] + logFactorials[a+c] + logFactorials[b+d] - logFactorials[a] - logFactorials[b] - logFactorials[c] - logFactorials[d] - logFactorials[a+b+c+d];
	return result;
}

int** generate_cd_crosstable(std::vector <int> phenotype_vector, std::vector<int> cur_genotype_vector, size_t vectors_length, size_t phen_levels_number, size_t gen_levels_number){
	
		
	int** cross_table = new int* [phen_levels_number];
	
	for (size_t k = 0; k < gen_levels_number; k++){
		cross_table[k] = new int [gen_levels_number];
		std::fill_n(cross_table[k], gen_levels_number, 0);
	}
	//we take genotypes 0, 1, 2 and phenotypes 1(-1),2(-1), ... n-1 as indexes 
	for (size_t m = 0; m < vectors_length; m++){
			int cur_gen = cur_genotype_vector[m];
			int cur_phen = phenotype_vector[m];
			if (cur_gen != 3)
				cross_table[cur_phen - 1][cur_gen]++;
	}
	
	return cross_table;
	
}

int** generate_d_crosstable(int** cd_crosstable, size_t phen_rows_num, size_t gen_cols_num){
		
	//table with genotypes-columns 0, 1, 2 and phenotypes 0, 1, 2, ..., n  
	int** cross_table = new int* [phen_rows_num];
	
	for (size_t i = 0; i < phen_rows_num; ++i){
		cross_table[i] = new int[gen_cols_num];
		std::fill_n(cross_table[i], gen_cols_num, 0);
	}
	
	
	for (size_t i = 0; i < phen_rows_num; ++i){
		cross_table[i][0] = cd_crosstable[i][0] + cd_crosstable[i][1];
		cross_table[i][1] = cd_crosstable[i][2];
	}
		
	return cross_table;
}

int** generate_r_crosstable(int** cd_crosstable, size_t phen_rows_num, size_t gen_cols_num){
	//table with genotypes-columns 0, 1, 2 and phenotypes 0, 1, 2, ..., n  
	int** cross_table = new int*[phen_rows_num];
	
	for (size_t i = 0; i < phen_rows_num; ++i){
		cross_table[i] = new int[gen_cols_num];
		std::fill_n(cross_table[i], gen_cols_num, 0);
	}
	
	for (size_t i = 0; i < phen_rows_num; ++i){
		cross_table[i][0] = cd_crosstable[i][0];
		cross_table[i][1] = cd_crosstable[i][1] + cd_crosstable[i][2];
	}
	
	
	
	return cross_table;
}

int** generate_allele_crosstable(int** cd_crosstable, size_t phen_rows_num, size_t gen_cols_num){
	//table with genotypes-columns 0, 1, 2 and phenotypes 0, 1, 2, ..., n  
	int** cross_table = new int*[phen_rows_num];
	
	for (size_t i = 0; i < phen_rows_num; ++i){
		cross_table[i] = new int[gen_cols_num];
		std::fill_n(cross_table[i], gen_cols_num, 0);
	}
	
	for (size_t i = 0; i < phen_rows_num; ++i){
		cross_table[i][0] = 2*cd_crosstable[i][0] + cd_crosstable[i][1];
		cross_table[i][1] = 2*cd_crosstable[i][2] + cd_crosstable[i][1];
		}
	return cross_table;
}

std::vector <int> get_row_sums(int** table, size_t phen_rows_num, size_t gen_cols_num) { 
	std::vector <int> row_sums_vec;
	int gensum = 0;	
	
	for (size_t i = 0; i < phen_rows_num; ++i){
		for (size_t j = 0; j < gen_cols_num; ++j){
			gensum += table[i][j];
		}
		row_sums_vec.push_back(gensum);
		gensum = 0;
	}
	
	return row_sums_vec;
}

std::vector <int> get_col_sums (int** table, size_t phen_rows_num, size_t gen_cols_num){
	std::vector <int> col_sums_vec;
	int gensum = 0;	
	
	for (size_t j = 0; j < gen_cols_num; ++j){
		for (size_t i = 0; i < phen_rows_num; ++i){
			gensum += table[i][j];
		}
		col_sums_vec.push_back(gensum);
		gensum = 0;
	}
	
	return col_sums_vec;
}

int get_total_sum(std::vector<int> row_sums){
	int total = 0;
	for (auto i = row_sums.begin(); i != row_sums.end(); i++)
		total += *i;
		
	return total;
}
	
double get_corr_coefficient(int** table, size_t rows, size_t cols, std::vector <int> *row_sums, std::vector <int> *col_sums, int total_sum ){
	
	double sum_Xi_Yi = 0;
	for (size_t i = 0; i < rows; i ++){
		for (size_t j = 0; j < cols; j++){
			sum_Xi_Yi += table[i][j]*i*j;
		}
	}
	
	
	sum_Xi_Yi = sum_Xi_Yi / total_sum;
	
	double rows_sum = 0;
	for (auto i = 0; i < rows; i++){
		rows_sum += row_sums->at(i)*i;
	}

	double cols_sum = 0;
	for (size_t j = 0; j < cols; j++){
		cols_sum += col_sums->at(j)*j;
		}
		
	double x_average = rows_sum / total_sum;
	double y_average = cols_sum / total_sum;
	
	
	double multiplier1 = 0;
	double multiplier2 = 0;
	
	for (size_t i = 0; i < rows; i++){
		multiplier1 += row_sums->at(i) * (i - x_average) * (i - x_average) / total_sum;
	}
	
	for (size_t j = 0; j < cols; j++){
		multiplier2 += col_sums->at(j) * (j - y_average) * (j - y_average) / total_sum;
	}
	
	double corr_coefficient = ( sum_Xi_Yi - x_average*y_average ) / sqrt(multiplier1*multiplier2);
	
	return corr_coefficient;
}


bool lessThan5(std::vector <int> *row_sums, std::vector <int> *col_sums, int total_sum){
	
	
	for (auto k = row_sums->begin(); k  != row_sums->end(); k++)
		for (auto m = col_sums->begin(); m != col_sums->end(); m++ )
			if(  (double)(*k) * (*m) / (double)total_sum < 5 ){
				return true;
				break;
				}
	
	return false;	
}

double get_p_value(double chi2_stat, size_t deg_freedom){
	if (std::isinf(chi2_stat))
		return 1.0/0.0;
	else
		return ( 1.0 -  boost::math::tgamma_lower(deg_freedom*0.5, chi2_stat*0.5) / boost::math::tgamma(deg_freedom*0.5) );
}
