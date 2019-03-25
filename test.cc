#include "matrix.h"

#include <iostream>
#include <chrono>
#include <complex>


int main() {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_seconds;









	//test 1 - construction, destruction;
	std::vector<std::vector<int>> init1 = 	{{1,2,3,4,5},
											 {5,4,3,2,1},
											 {0,0,0,0,0}};

	start = std::chrono::system_clock::now();
	matrix<int> mat1(init1);
	end = std::chrono::system_clock::now(); 
	elapsed_seconds = end - start; 
	std::cout << "Test 1 completed. Time elapsed: " << elapsed_seconds.count() << std::endl;










	//test 2 - matrix addition
	std::vector<std::vector<int>> init2 = 	{{1,2,3,4,5},
											 {5,4,3,2,1},
											 {0,0,0,3,0},
											 {0,2,1,0,0},
											 {1,0,3,0,5}};

	std::vector<std::vector<int>> init3 = 	{{4,2,3,4,5},
											 {5,1,3,2,1},
											 {4,0,1,3,1},
											 {0,2,1,0,3},
											 {1,1,3,0,5}};

	matrix<int> mat2(init2);
	matrix<int> mat3(init3);

	start = std::chrono::system_clock::now();
	matrix<int> matSum = mat2 + mat3;
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	std::cout << "Test 2 completed. Time elapsed: " << elapsed_seconds.count() << std::endl; 











	//test 3 - scalar multiplication
	std::vector<std::vector<double>> init4 = 	{{1.4,2.2,3.3,4.1,5.6},
											 	 {5.0,4.1,3.2,2.5,1.2},
											 	 {0.0,0.1,0.3,0.0,0.1}};

	double scalar1 = 10;
	matrix<double> mat4(init4); 

	start = std::chrono::system_clock::now();
	mat4 *= scalar1;
	end = std::chrono::system_clock::now();

	elapsed_seconds = end - start;
	std::cout << "Test 3 completed. Time elapsed: " << elapsed_seconds.count() << std::endl; 









	//test 4 - matrix multiplication 1
	std::vector<std::vector<int>> init5 = 	{{1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5}};


	std::vector<std::vector<int>> init6 = 	{{1,2,3,4,5},
											 {5,4,3,2,1},
											 {0,0,0,0,0},
											 {1,3,5,2,1},
											 {1,6,5,2,1},
											 {1,2,3,4,5},
											 {5,4,3,2,1},
											 {0,0,0,0,0},
											 {1,3,5,2,1},
											 {1,6,5,2,1}};

	matrix<int> mat5(init5);
	matrix<int> mat6(init6);

	start = std::chrono::system_clock::now();
	matrix<int> resMat1 = mat5 * mat6;
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	std::cout << "Test 4 completed. Time elapsed: " << elapsed_seconds.count() << std::endl; 











	//test 5 - matrix chain multiplication
	std::vector<std::vector<int>> init7 = 	{{1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5}};

	std::vector<std::vector<int>> init8 = 	{{1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5}};

	std::vector<std::vector<int>> init9 = 	{{1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5}};

	matrix<int> mat7(init7);
	matrix<int> mat8(init8);
	matrix<int> mat9(init9);

	start = std::chrono::system_clock::now();
	matrix<int> resMat2 = (mat7 * mat8) * mat9; 
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	std::cout << "Test 5 completed. Time elapsed: " << elapsed_seconds.count() << std::endl; 








	// test 6 - linear combination of multiplications
	std::vector<std::vector<int>> init10 = 	{{1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5}};

	std::vector<std::vector<int>> init11 = 	{{1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5}};

	std::vector<std::vector<int>> init12 = 	{{1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,5,4,3,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,0,0,0,0,0,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,3,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5},
											 {1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,6,5,2,1,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5}};

	matrix<int> mat10(init10);
	matrix<int> mat11(init11);
	matrix<int> mat12(init12);

	int scalar2 = 1;
	int scalar3 = 2;
	int scalar4 = 3; 

	start = std::chrono::system_clock::now();
	mat10 *= scalar2; 
	mat11 *= scalar3;
	mat12 *= scalar4;

	matrix<int> resMat3 = (mat10 * mat11) * mat12;
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	std::cout << "Test 6 completed. Time elapsed: " << elapsed_seconds.count() << std::endl; 
	







	//Test 7 - stress testing (to be implemented)

	return 0;
}