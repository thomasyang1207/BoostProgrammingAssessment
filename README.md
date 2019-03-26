# BoostProgrammingAssessment

This is a simple implementation of a generic matrix container. All data types are supported as long as addition, subtraction, multiplication, and division are defined for them. Thus, integers, floats, and complex numbers are supported. 

####Matrix representation
Matrices are represented as a 2D array within the matrix class, with a row-major configuration. Thus, for a matrix containing m rows and n columns, the element at row i and column j is A[i*n + j]. 

####Implementation of matrix multiplication:
Given the way I represented the matrix, I interpreted matrix multiplication as multiple linear combinations of rows of the SECOND matrix, using the elements of each row in the first matrix as weights.
