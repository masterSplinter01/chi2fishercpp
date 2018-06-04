#include "tables_chi_correlation.h"


// [[Rcpp::plugins(openmp)]]

// [[Rcpp::plugins(cpp11)]]  
// [[Rcpp::export]]



Rcpp::List run_chisquare_test(int phenotype_length, int phen_levels_number, int genotype_num_cols, Rcpp::NumericVector phenotype_vector_r, Rcpp::NumericMatrix genotype_matrix_r, int threads){
	
	unsigned threadsSupported = threads; //std::thread::hardware_concurrency() - 1;
	const int gen_levels_number = 3; //genotypes: 0, 1, 2. 3 - no data
	
	std::vector <double> stat_cd(genotype_num_cols);
	std::vector <double> stat_d(genotype_num_cols);
	std::vector <double> stat_r(genotype_num_cols);
	std::vector <double> stat_allele(genotype_num_cols);

	std::vector <int> degree_freedom_cd(genotype_num_cols);
	std::vector <int> degree_freedom_d(genotype_num_cols);
	std::vector <int> degree_freedom_r(genotype_num_cols);
	std::vector <int> degree_freedom_allele(genotype_num_cols);
	
	std::vector <double> p_value_cd(genotype_num_cols);
	std::vector <double> p_value_d(genotype_num_cols);
	std::vector <double> p_value_r(genotype_num_cols);
	std::vector <double> p_value_allele(genotype_num_cols);
	
	std::vector <std::string> typeTest_cd_vec(genotype_num_cols);
	std::vector <std::string> typeTest_d_vec(genotype_num_cols);
	std::vector <std::string> typeTest_r_vec(genotype_num_cols);
	std::vector <std::string> typeTest_allele_vec(genotype_num_cols);
	
	std::vector <bool> warning_cd(genotype_num_cols);
	std::vector <bool> warning_d(genotype_num_cols);
	std::vector <bool> warning_r(genotype_num_cols);
	std::vector <bool> warning_allele(genotype_num_cols);
	
	std::vector <double> correlation_coefficient_cd(genotype_num_cols);
	std::vector <double> correlation_coefficient_d(genotype_num_cols);
	std::vector <double> correlation_coefficient_r(genotype_num_cols);
	std::vector <double> correlation_coefficient_allele(genotype_num_cols);
	
	//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	
	std::vector<int> phenotype_vector = Rcpp::as < std::vector<int> > (phenotype_vector_r); //r vector to cpp vector
	
	int** genotype_matrix = new int*[genotype_num_cols];
	for (size_t i = 0; i < genotype_num_cols; ++i)
		genotype_matrix[i] = new int[phenotype_length];
	
	//r matrix -> r vector -> cpp matrix
	Rcpp::NumericVector gen_vec(genotype_matrix_r);
	
	int k = 0;
	for (int i = 0; i < genotype_num_cols; i++){
		for(int j = 0; j < phenotype_length; j++){
			genotype_matrix[i][j] = gen_vec[k];
			k++;
		}
	}
	
	#pragma omp parallel for num_threads(threadsSupported)
	
	for (size_t j = 0; j < genotype_num_cols; j++){
		int gen;
		std::vector<int> cur_genotype;
		for (int i = 0; i < phenotype_length; i++){
			gen = genotype_matrix[j][i];
			cur_genotype.push_back(gen);
		} //cur_genotype vector - one column of genotype matrix
		
		
		//cd analysis
		int** cd_crosstable = generate_cd_crosstable (phenotype_vector, cur_genotype, phenotype_length, phen_levels_number,
														gen_levels_number); 
		run_test(cd_crosstable, phen_levels_number, gen_levels_number,&stat_cd, &p_value_cd, &degree_freedom_cd, &typeTest_cd_vec, &warning_cd, &correlation_coefficient_cd, j);
		
		
		//d analysis
		int** d_crosstable = generate_d_crosstable(cd_crosstable, phen_levels_number, gen_levels_number-1);
		run_test(d_crosstable, phen_levels_number, gen_levels_number-1, &stat_d, &p_value_d, &degree_freedom_d, &typeTest_d_vec, &warning_d, &correlation_coefficient_d, j);
		
		
		//r analysis
		int** r_crosstable = generate_r_crosstable(cd_crosstable, phen_levels_number, gen_levels_number-1);
		run_test(r_crosstable, phen_levels_number, gen_levels_number-1, &stat_r, &p_value_r, &degree_freedom_r, &typeTest_r_vec, &warning_r, &correlation_coefficient_r, j);
		
				
		//allele analysis
		int** allele_crosstable = generate_allele_crosstable(cd_crosstable,phen_levels_number, gen_levels_number-1);
		run_test(allele_crosstable, phen_levels_number, gen_levels_number-1, &stat_allele, &p_value_allele, &degree_freedom_allele, &typeTest_allele_vec, &warning_allele, &correlation_coefficient_allele, j);
		
	}
	
	
	//std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	
	//Rf_PrintValue(Rcpp::wrap(std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count()));
	

	for (size_t i = 0; i < genotype_num_cols; ++i)
		delete[] genotype_matrix[i];
	delete[] genotype_matrix;
	
	
	return Rcpp::List::create(	Rcpp::List::create(	
									Rcpp::Named("stat_cd") = stat_cd, 
									Rcpp::Named("deg_freedom_cd") = degree_freedom_cd,
									Rcpp::Named("p_value_cd") = p_value_cd,
									Rcpp::Named("corr_cd") = correlation_coefficient_cd,
									Rcpp::Named("typetest_cd") = typeTest_cd_vec,
									Rcpp::Named("elements_lower_than5_cd") = warning_cd
									),
								Rcpp::List::create(
									Rcpp::Named("stat_d") = stat_d, 
									Rcpp::Named("deg_freedom_d") = degree_freedom_d, 
									Rcpp::Named("p_value_d") = p_value_d,
									Rcpp::Named("corr_d") = correlation_coefficient_d,
									Rcpp::Named("typeTest_d") = typeTest_d_vec,
									Rcpp::Named("elements_lower_than5_d") = warning_d
								),
								Rcpp::List::create(
									Rcpp::Named("stat_r") = stat_r, 
									Rcpp::Named("deg_freedom_r") = degree_freedom_r,
									Rcpp::Named("p_value_r") = p_value_r,
									Rcpp::Named("corr_r") = correlation_coefficient_r,
									Rcpp::Named("typeTest_r") = typeTest_r_vec,
									Rcpp::Named("elements_lower_than5_r") = warning_r
								),
								Rcpp::List::create(
									Rcpp::Named("stat_allele") = stat_allele, 
									Rcpp::Named("deg_freedom_allele") = degree_freedom_allele,
									Rcpp::Named("p_value_allele") = p_value_allele,	
									Rcpp::Named("corr_allele") = correlation_coefficient_allele,
									Rcpp::Named("typeTest_allele") = typeTest_allele_vec,	
									Rcpp::Named("elements_lower_than5_allele") = warning_allele
								)		
							);
}


