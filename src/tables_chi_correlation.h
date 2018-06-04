#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <Rcpp.h>
#include <vector>
#include <math.h>
#include <string>
#include <thread>
#ifdef _OPENMP
	#include <omp.h>
#endif
#include <boost/math/special_functions/gamma.hpp>
#include <fstream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <vector>
//#include <math.h>
//#include <string>

void run_test (int** cross_table, size_t phen_rows_num, size_t gen_cols_num, 
				std::vector<double> *stat, std::vector<double> *p_value, std::vector <int> *degree_freedom_vec, 
				std::vector<std::string> *type_vec, std::vector<bool> *warning, 
				std::vector <double> *correlation_coefficient,  int iteration);

double chi2_test(int** cross_table, size_t phen_rows_num, size_t gen_cols_num, int* deg_freedom, std::vector <int> *row_sums, std::vector <int> *col_sums, int total_sum);

double fisher_exact_test(int** cross_table);
void get_logFactorials(double* logFactorials, int n);
double logHypergeometricProb(double* logFactorials, int a, int b, int c, int d);

int** generate_cd_crosstable(std::vector <int> phenotype_vector, std::vector<int> cur_genotype_vector, size_t vectors_length, size_t phen_levels_number, size_t gen_levels_number);
int** generate_d_crosstable(int** cd_crosstable, size_t phen_rows_num, size_t gen_cols_num);
int** generate_r_crosstable(int** cd_crosstable, size_t phen_rows_num, size_t gen_cols_num);
int** generate_allele_crosstable(int** cd_crosstable, size_t phen_rows_num, size_t gen_cols_num);

std::vector <int> get_row_sums(int** table, size_t phen_rows_num, size_t gen_cols_num);
std::vector <int> get_col_sums (int** table, size_t phen_rows_num, size_t gen_cols_num);
int get_total_sum(std::vector<int> row_sums);

double get_corr_coefficient(int** table, size_t rows, size_t cols, std::vector <int> *row_sums, std::vector <int> *col_sums, int total_sum );
bool lessThan5(std::vector <int> *row_sums, std::vector <int> *col_sums, int total_sum);

double get_p_value(double chi_square_stat, size_t deg_freedom);