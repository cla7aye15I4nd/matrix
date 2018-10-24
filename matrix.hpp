#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <stdexcept>
#include <functional>

namespace matrix{
#define ASSERT(expr, message) if (!(expr)) { throw message; }
  template <typename T>
  class Matrix{
    using size_pair = std::pair<size_t, size_t>;
    
    T *x;
    size_t n, m;
          
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
      auto p = this -> find(i, 0);
      for (auto ptr = res.begin(); ptr != res.end(); ++ptr, ++p)
        *ptr = *p;
      return res;
    }
    
    Matrix<T> column(size_t i) const{
      ASSERT(i >= 0 && i < m, "column error");
      Matrix<T> res(n, 1);
      auto p = this -> find(0, i);
      for (auto ptr = res.begin(); ptr != res.end(); ++ptr, p += m)
        *ptr = *p;
      return res;
    }

  public:
    Matrix (size_t _n = 1, size_t _m = 1, T v = T()) {
      n = _n; m = _m;
      x = new T [n * m];
      for (auto ptr = this -> begin(); ptr != this -> end(); ++ptr)
        *ptr = v;
    }
    explicit Matrix(size_pair sz, T v = T()) 
      :x(new T [sz.first * sz.second]), n(sz.first) , m(sz.second) {
      for (auto ptr = this -> begin(); ptr != this -> end(); ++ptr)
        *ptr = v;
    }
    Matrix (const std::initializer_list<std::initializer_list<T>> &matrix) {
      n = matrix.size(); m = matrix.begin() -> size();
      
      T* ptr = x = new T [n * m];
    
      for (auto &list : matrix) {
        ASSERT(m == list.size(), "initiallizer error");
        for (auto &v : list) { *ptr++ = v; }
      }
    }
    Matrix (Matrix &&v) {
      n = v.n;
      m = v.m;
      x = v.x;
      v.x = nullptr;
    }

    template <typename U>
    Matrix (const Matrix<U> &v) {
      n = v.rowLength();
      m = v.columnLength();
      x = new T [n * m];
      auto ptr = this -> begin();
      for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
        *ptr = *p;
    }
    Matrix (const Matrix &v) {
      n = v.rowLength(); m = v.columnLength();
      x = new T [n * m];
      auto ptr = this -> begin();
      for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
        *ptr = *p;
    }
    
    ~Matrix () {
      if (n && m)
        delete [] x;
    }

  public:
    
    Matrix<T>& operator= (const Matrix<T> &v) {
      if (this == &v) return *this;
      n = v.n; m = v.m; x = new T [n * m];
      auto ptr = this -> begin();
      for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
        *ptr = *p;
      return *this;
    }
    Matrix<T>& operator= (Matrix<T> &&v) {
      if (x == v.x) return *this;
      n = v.n; m = v.m; x = v.x; v.x = nullptr;
      return *this;
    }
    template <class K>
    Matrix<T>& operator= (const Matrix<K> &v) {
      n = v.rowLength();
      m = v.columnLength();
      x = new T [n * m];
      auto ptr = this -> begin();
      for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
        *ptr = *p;
      return *this;
    }
          
  public:
    void clear() { n = 0; m = 0; delete [] x; }
    size_t rowLength() const{ return n; }
    size_t columnLength() const{ return m; }
    size_pair size() const{ return {rowLength(), columnLength()}; }
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
    void resize(size_pair sz, T v = T()) { resize(sz.first, sz.second, v); }
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
      ASSERT(n == rhs.n && m == rhs.m, std::invalid_argument("+="));
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
          (*this)[i][j] += rhs[i][j];
      return *this;
    }
    template <typename K>
    Matrix<T> operator-= (const Matrix<K> &rhs) {
      ASSERT(n == rhs.n && m == rhs.m, std::invalid_argument("-="));
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
          (*this)[i][j] -= rhs[i][j];
      return *this;
    }
    template <typename K>
    Matrix<T> operator*= (const K& v) {
      Matrix<T> res(n, m);
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j)
          x[i * m + j] *= v;
      return *this;
    }
    Matrix tran() const{
      Matrix res(m, n);
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
    Matrix<decltype(T() + U())> operator+ (const Matrix<U> &rhs) const{
      ASSERT(n == rhs.rowLength() && m == rhs.columnLength(), std::invalid_argument("+"));
      Matrix<decltype(T() + U())> res(n, m);
      for (size_t i = 0; i < res.rowLength(); ++i)
        for (size_t j = 0; j < res.columnLength(); ++j)
          res[i][j] = (*this)[i][j] + rhs[i][j];
      return res;
    }

    template <typename U>
    Matrix<decltype(T() - U())> operator- (const Matrix<U> &rhs) const{
      ASSERT(n == rhs.n && m == rhs.m, std::invalid_argument("-"));
      Matrix<decltype(T() + U())> res(n, m);
      for (size_t i = 0; i < res.n; ++i)
        for (size_t j = 0; j < res.m; ++j)
          res[i][j] = (*this) [i][j] - rhs[i][j];
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
			
      iterator () = default;
      iterator (const iterator&) = default;
      iterator &operator=(const iterator &) = default;
      iterator (const pointer& ptr):ptr(ptr) {}
    private:
      pointer ptr;
               
    public:
      difference_type operator- (const iterator &p) { return ptr - p.ptr; }      
      iterator operator+ (difference_type offset) const{ return ptr + offset;; }
      iterator operator- (difference_type offset) const{ return ptr - offset; }	
      iterator& operator++ () { ++ptr; return *this; }	
      iterator& operator-- () { --ptr; return *this; }
      iterator& operator+= (difference_type offset) { ptr += offset; return *this; }
      iterator& operator-= (difference_type offset) { ptr -= offset; return *this; }
      
      reference operator* () const { return *ptr; }
      pointer operator->() const { return ptr; }
      bool operator==(const iterator &o) const{ return ptr == o.ptr; }
      bool operator!=(const iterator &o) const{ return ptr != o.ptr; }
    };
    
    iterator begin() const{
      return n == 0 || m == 0 ? nullptr : x;
    }
    iterator end() const{
      return n == 0 || m == 0 ? nullptr : x + n * m;
    }
    iterator find(size_t i, size_t j) const{
      return x + i * m + j;
    }
		
    std::pair<iterator, iterator> subMatrix(size_pair l, size_pair r) const{
      Matrix<T> *mat = new Matrix<T> (r.first - l.first + 1,
                                      r.second - l.second + 1);
      for (size_t i = 0; i < mat->rowLength(); ++i)
        for (size_t j = 0; j < mat->columnLength(); ++j)
          (*mat)(i, j) = (*this) [l.first + i][l.second + j];
      return {mat->begin(), mat->end()};
    }
  };
  
  template <class T, class U>
  Matrix<decltype(T() * U())> operator*(const Matrix<T> &mat, const U &x) {
    Matrix<decltype(T() * U())> res(mat.rowLength(), mat.columnLength());
    for (size_t i = 0; i < res.rowLength(); ++i)
      for (size_t j = 0; j < res.columnLength(); ++j)
        res[i][j] = mat[i][j] * x;
    return res;
  }
  
  template <class T, class U>
  Matrix<decltype(T() * U())> operator*(const U &x, const Matrix<T> &mat) {
    Matrix<decltype(T() * U())> res(mat.rowLength(), mat.columnLength());
    for (size_t i = 0; i < res.rowLength(); ++i)
      for (size_t j = 0; j < res.columnLength(); ++j)
        res[i][j] = mat[i][j] * x;
    return res;
  }

  template <class U, class V>
  Matrix<decltype(U() * V())> operator*(const Matrix<U> &lhs, const Matrix<V> &rhs) {
    ASSERT(lhs.columnLength() == rhs.rowLength(), std::invalid_argument("(* matrix matrix)"));
    Matrix<decltype(U() * V())> res(lhs.rowLength(), rhs.columnLength());
    for (size_t i = 0; i < res.rowLength(); ++i)
      for (size_t j = 0; j < res.columnLength(); ++j)
        for (size_t k = 0; k < rhs.rowLength(); ++k)
          res[i][j] += lhs[i][k] * rhs[k][j];
    return res;
  }
}

namespace sjtu{
  using namespace matrix;
}

#endif //MATRIX_HPP
