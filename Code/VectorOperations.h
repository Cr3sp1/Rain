#ifndef __VOper_h__
#define __VOper_h__

#include <iostream>
#include <vector>
#include <cassert>
#include <numeric>
#include "VectorStat.h"


// ===============================================================================
// sum of two vectors: sum component by component
// ===============================================================================

template <typename T> 
inline std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
  
  assert(a.size() == b.size());  
  std::vector<T> result(a.size());
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = a[i] + b[i];    
  
  return result;
  
}

// ===============================================================================
// difference of two vectors component by component
// ===============================================================================

template <typename T> 
 inline std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b) {
  
  assert(a.size() == b.size());  
  std::vector<T> result(a.size());
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = a[i] - b[i];    
  
  return result;
  
}

// ===============================================================================  
//  scalar product between two vectors 
// ===============================================================================

template <typename T>
inline T operator*(const std::vector<T> &a, const std::vector<T> &b) {
    assert(a.size() == b.size());
    return std::inner_product(a.begin(), a.end(), b.begin(), T(0));
}

// ===============================================================================  
// cross product between two vectors 
// ===============================================================================
template <typename T> 
inline std::vector<T> CrossProduct(const std::vector<T> &a, const std::vector<T> &b) {
  assert(a.size() == b.size() && a.size() == 3);

  std::vector<T> result(3);
  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = a[2] * b[0] - a[0] * b[2];
  result[2] = a[0] * b[1] - a[1] * b[0];

  return result;
}

// ===============================================================================
// matrix-vector multiplication
// ===============================================================================
template <typename T>
inline std::vector<T> operator*(const std::vector<std::vector<T>>& matrix, const std::vector<T>& vec) {
    if(matrix.size() <= 0 or matrix[0].size() != vec.size()){
      cout << "Error: tried applying matrix: " << endl;
      Print(matrix);
      cout << "To vector : " << endl;
      Print(vec);
    }
    assert(matrix.size() > 0 && matrix[0].size() == vec.size());

    std::vector<T> result(matrix.size(), T(0));  // Initialize result vector with zeros

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < vec.size(); ++j) {
            result[i] += matrix[i][j] * vec[j];
        }
    }

    return result;
}

// ===============================================================================  
// norm of a vector 
// ===============================================================================
template <typename T> 
inline T Norm(const std::vector<T> &a) {
  return sqrt(a*a);
}

// ===============================================================================
// product between a scalar and a vector
// ===============================================================================

template <typename T> 
inline std::vector<T> operator*( T c , const std::vector<T> &a) {
  
  std::vector<T> result(a.size());
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = c * a[i];  

  
  return result;
  
}

// ===============================================================================
// product between a scalar and a vector
// ===============================================================================

template <typename T> 
inline std::vector<T> operator*( const std::vector<T> &a , T c) {
  
  std::vector<T> result(a.size());
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = c * a[i];  

  return result;
  
}


// ===============================================================================
// division between a vector and a scalar  
// ===============================================================================

template <typename T> 
inline std::vector<T> operator/( const std::vector<T> &a , T c) {
  std::vector<T> result(a.size());
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = a[i] / c ;  
  return result;
}

// ===============================================================================
// add to a vector a, vector b, and the result is stored in a
// ===============================================================================

template <typename T> 
inline std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& b) {
  
  assert(a.size() == b.size());  
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) a[i] += b[i];    

  return a;
}

// ===============================================================================
// subtract from a vector a, vector b, and the result is stored in a
// ===============================================================================

template <typename T> std::vector<T>& operator-=(std::vector<T>& a, const std::vector<T>& b) {
  
  assert(a.size() == b.size());  
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) a[i] -= b[i];    
  
  return a;
}

// ===============================================================================
// Transpose of a matrix
// ===============================================================================

template <typename T> std::vector<std::vector<T>> Transpose(const std::vector<std::vector<T>>& matrix) {
  if (matrix.empty()) return {};

  size_t rows = matrix.size();
  size_t cols = matrix[0].size();
  std::vector<std::vector<T>> transposed(cols, std::vector<T>(rows));

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      transposed[j][i] = ( j <= matrix[i].size() ) ? matrix[i][j] : 0;
    }
  }

  return transposed;
}




#endif