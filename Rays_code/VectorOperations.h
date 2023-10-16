#ifndef __VOper_h__
#define __VOper_h__

#include <iostream>
#include <vector>
#include <cassert>


// ===============================================================================
// sum of two vectors: sum component by component
// ===============================================================================

template <typename T> std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
  
  assert(a.size() == b.size());  
  std::vector<T> result(a.size());
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = a[i] + b[i];    
  
  return result;
  
}

// ===============================================================================
// difference of two vectors component by component
// ===============================================================================

template <typename T> std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b) {
  
  assert(a.size() == b.size());  
  std::vector<T> result(a.size());
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = a[i] - b[i];    
  
  return result;
  
}

// ===============================================================================  
//  scalar product between two vectors 
// ===============================================================================

template <typename T> T operator*(const std::vector<T> &a, const std::vector<T> &b) {
  
  assert(a.size() == b.size());  
  T sum = 0 ;
  for (int i = 0; i < static_cast<int>(a.size()); i++) sum += a[i] * b[i];  
  return sum;
  
}

// ===============================================================================  
// cross product between two vectors 
// ===============================================================================
template <typename T> std::vector<T> CrossProduct(const std::vector<T> &a, const std::vector<T> &b) {
  assert(a.size() == b.size() && a.size() == 3);

  std::vector<T> result(3);
  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = a[2] * b[0] - a[0] * b[2];
  result[2] = a[0] * b[1] - a[1] * b[0];

  return result;
}

// ===============================================================================  
// norm of a vector 
// ===============================================================================
template <typename T> T Norm(const std::vector<T> &a) {
  return sqrt(a*a);
}

// ===============================================================================
// product between a scalar and a vector
// ===============================================================================

template <typename T> std::vector<T> operator*( T c , const std::vector<T> &a) {
  
  std::vector<T> result(a.size());
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = c * a[i];  

  
  return result;
  
}

// ===============================================================================
// product between a scalar and a vector
// ===============================================================================

template <typename T> std::vector<T> operator*( const std::vector<T> &a , T c) {
  
  std::vector<T> result(a.size());
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = c * a[i];  

  return result;
  
}


// ===============================================================================
// division between a vector and a scalar  
// ===============================================================================

template <typename T> std::vector<T> operator/( const std::vector<T> &a , T c) {
  std::vector<T> result(a.size());
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = a[i] / c ;  
  return result;
}

// ===============================================================================
// add to a vector a, vector b, and the result is stored in a
// ===============================================================================

template <typename T> std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& b) {
  
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


#endif