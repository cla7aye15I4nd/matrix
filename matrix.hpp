// Copyright (c) 2018 dataisland
 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <bits/stdc++.h>

namespace matrix{
#define ASSERT(expr, message) if (!(expr)) { throw message; }
  template <typename T>
  class Matrix{
    T *x;
    size_t n, m;
    using msg = std::pair<size_t, size_t>;
          
  public:
    T* operator [] (size_t u) { return x + u * m; }
    const T* operator [] (size_t u) const{ return x + u * m; }
    T& operator() (size_t i, size_t j) {
      ASSERT(i >= 0 && i < n && j >= 0 && j < m, "acccess error");
      return (*this)[i][j];
    }
    const T& operator() (size_t i, size_t j) const{
      ASSERT(i >= 0 && i < n && j >= 0 && j < m, "acccess error");
      return (*this)[i][j];
    }
    Matrix<T> row(size_t i) const{
      ASSERT(i >= 0 && i < n, "row error");
      Matrix<T> res(1, m);
      for (size_t j = 0; j < m; ++j)
        res[0][j] = (*this)[i][j];
      return res;
    }

    Matrix<T> column(size_t i) const{
      ASSERT(i >= 0 && i < m, "column error");
      Matrix<T> res(n, 1);
      for (size_t j = 0; j < n; ++j)
        res[j][0] = (*this)[j][i];
      return res;
    }

  public:
    Matrix () = default;
    Matrix (size_t n, size_t m, T v = T()):
      x(new T [n * m]), n(n) , m(m) {
      for (T *p = x, *e = x + n * m; p != e; *p++ = v);
    }
    explicit Matrix(msg sz, T v = T()) 
      :x(new T [sz.first * sz.second]), n(sz.first) , m(sz.second) {
      for (T *p = x, *e = x + n * m; p != e; *p++ = v);
    }
    Matrix (const std::initializer_list<std::initializer_list<T>> &matrix) {
      n = matrix.size(); m = matrix.begin() -> size();
    
      T* ptr = x = new T [n * m];
    
      for (auto &list : matrix) {
        ASSERT(m == list.size(), "initiallizer error");
        for (auto &v : list) { *ptr++ = v; }
      }
    }
    Matrix (Matrix &&v) { n = v.n; m = v.m; x = v.x; v.x = nullptr;}
    Matrix (const Matrix &v) {
      n = v.n; m = v.m; x = new T [n * m];
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
          (*this)[i][j] = v[i][j];
    }
    ~Matrix () { delete [] x; }

  public:
    Matrix<T>& operator= (const Matrix<T> &v) {
      n = v.n; m = v.m; x = new T [n * m];
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
          (*this)[i][j] = v[i][j];
      return *this;
    }
    Matrix<T>& operator= (Matrix<T> &&v) {
      n = v.n; m = v.m; x = v.x; v.x = nullptr;
      return *this;
    }
    template <class K>
    Matrix<T>& operator= (const Matrix<K> &v) {
      n = v.n; m = v.m; x = new T [n * m];
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
          (*this)[i][j] = v[i][j];
      return *this;
    }
          
  public:
    void clear() { delete [] x; }
    size_t rowLength() const{ return n; }
    size_t columnLength() const{ return m; }
    msg size() const{ return {rowLength(), columnLength()}; }
    void resize(size_t _n, size_t _m, T v = T()) {
      if (n * m == _n * _m) {
        n = _n; m = _m;
      } else {
        T* temp = x;
        x = new T [_n * _m];
        size_t c = 0;
        for (size_t i = 0; i < _n; ++i)
          for (size_t j = 0; j < _m; ++j, ++c) {
            x[i * _m + j] = c < n * m ? temp[c]: v;
          }
        delete [] temp;
        n = _n; m = _m;
      }
    }
    void resize(msg sz, T v = T()) { resize(sz.first, sz.second, v); }

  public:
    Matrix<T> operator- () const{
      Matrix<T> res(n, m);
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
          res[i][j] = -(*this)[i][j];
      return res;
    }
    template <typename K>
    Matrix<T> operator+= (const Matrix<K> &rhs) {
      ASSERT(n == rhs.n && m == rhs.m, "error: size mismatch(+=)");
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
          (*this)[i][j] += rhs[i][j];
      return *this;
    }
    template <typename K>
    Matrix<T> operator-= (const Matrix<K> &rhs) {
      ASSERT(n == rhs.n && m == rhs.m, "error: size mismatch(-=)");
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
          (*this)[i][j] += rhs[i][j];
      return *this;
    }
    template <typename K>
    Matrix<T> operator*= (const K& v) const{
      Matrix<T> res(n, m);
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
          (*this)[i][j] *= v;
      return *this;
    }
    Matrix<T> tran() const{
      Matrix<T> res(m, n);
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
          res[j][i] = (*this) [i][j];
      return res;
    }

  public:
    bool operator == (const Matrix<T>& rhs) const{
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
          if ((*this)[i][j] != rhs[i][j])
            return false;
      return true;
    }
    bool operator != (const Matrix<T>& rhs) const{
      return !(*this == rhs);
    }

  public:
    template <typename U>
    friend Matrix<decltype(T() + U())> operator+ (const Matrix<T> &lhs, const Matrix<U> &rhs) {
      ASSERT(lhs.n == rhs.n && lhs.m == rhs.m, "error: size mismatch(+)");
      Matrix<decltype(T() + U())> res(lhs.n, lhs.m);
      for (size_t i = 0; i < res.n; ++i)
        for (size_t j = 0; j < res.m; ++j)
          res[i][j] = lhs[i][j] + rhs[i][j];
      return res;
    }

    template <typename U>
    friend Matrix<decltype(T(0) - U(0))> operator- (const Matrix<T> &lhs, const Matrix<U> &rhs) {
      ASSERT(lhs.n == rhs.n && lhs.m == rhs.m, "error: size mismatch(-)");
      Matrix<decltype(T(0) + U(0))> res(lhs.n, lhs.m);
      for (size_t i = 0; i < res.n; ++i)
        for (size_t j = 0; j < res.m; ++j)
          res[i][j] = lhs[i][j] - rhs[i][j];
      return res;
    }

    template <typename U>
    friend Matrix<decltype(T(0) * U(0))> operator* (const Matrix<T> &lhs, const U& v) {
      Matrix<decltype(T(0) * U(0))> res(lhs.n, lhs.m);
      for (size_t i = 0; i < res.n; ++i)
        for (size_t j = 0; j < res.m; ++j)
          res[i][j] = lhs[i][j] * v;
      return res;
    }

    template <typename U>
    friend Matrix<decltype(T(0) * U(0))> operator* (const Matrix<T> &lhs, const Matrix<U> &rhs) {
      ASSERT(lhs.m == rhs.n, "error: size mismatch(*)");
      Matrix<decltype(T(0) * U(0))> res(lhs.n, rhs.m);
      for (size_t i = 0; i < res.n; ++i)
        for (size_t j = 0; j < res.m; ++j)
          for (size_t k = 0; k < rhs.n; ++k)
            res[i][j] += lhs[i][k] * rhs[k][j];
      return res;
    }
          
  public:
    class iterator {
    public:
      using iterator_category = std::random_access_iterator_tag;
      using value_type        = T;
      using pointer           = T *;
      using reference         = T &;
      using size_type         = size_t;
      using difference_type   = std::ptrdiff_t;
			
      iterator() = default;
      iterator(const iterator&) = default;
      iterator &operator=(const iterator &) = default;
			
    private:
      pointer ptr;
               
    public:
      difference_type operator- (const iterator &p) { return ptr - p; }
      iterator &operator+= (difference_type offset) { return ptr += offset; }
      iterator operator+ (difference_type offset) const{ return ptr + offset;; }
      iterator &operator-= (difference_type offset) { return ptr -= offset; }
      iterator operator- (difference_type offset) const{ return ptr - offset; }	
      iterator &operator++ () { return ++ptr; }	
      iterator &operator-- () { return --ptr; }
      reference operator* () const { return *ptr; }
      pointer operator->() const { return ptr; }
      bool operator==(const iterator &o) const{ return ptr == o; }
      bool operator!=(const iterator &o) const{ return ptr != o; }
    };
    
    iterator begin() { return x; }
    iterator end() { return x + n * m; }
		
    std::pair<iterator, iterator> subMatrix(msg l, msg r) {
      return {x + l.first * m + l.second, x + r.first * m + r.second};
    }
  };
}

namespace sjtu{
  using namespace matrix;
}

#endif //MATRIX_HPP
