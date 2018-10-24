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
      ASSERT(i >= 0 && i < n && j >= 0 && j < m, std::invalid_argument("operator ()"));
      return (*this)[i][j];
    }
    const T& operator() (size_t i, size_t j) const{
      ASSERT(i >= 0 && i < n && j >= 0 && j < m, std::invalid_argument("operator ()"));
      return (*this)[i][j];
    }
    
    Matrix<T> row(size_t i) const{
      ASSERT(i >= 0 && i < n, std::invalid_argument("row"));
      Matrix<T> res(1, m);
      auto p = this -> find(i, 0);
      for (auto ptr = res.begin(); ptr != res.end(); ++ptr, ++p)
        *ptr = *p;
      return res;
    }
    
    Matrix<T> column(size_t i) const{
      ASSERT(i >= 0 && i < m, std::invalid_argument("column"));
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
        ASSERT(m == list.size(), std::invalid_argument("initiallizer_list"));
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
      auto ptr = this -> begin();
      for (auto p = res.begin(); p != res.end(); ++p, ++ptr)
        *p = -*ptr;
      return res;
    }
    template <typename K>
    Matrix<T> operator+= (const Matrix<K> &v) {
      ASSERT(n == v.n && m == v.m, std::invalid_argument("operator+="));
      auto ptr = this -> begin();
      for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
        *ptr += *p;
      return *this;
    }
    template <typename K>
    Matrix<T> operator-= (const Matrix<K> &v) {
      ASSERT(n == v.n && m == v.m, std::invalid_argument("operator-="));
      auto ptr = this -> begin();
      for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
        *ptr -= *p;
      return *this;
    }
    template <typename K>
    Matrix<T> operator*= (const K& v) {
      for (auto ptr = this -> begin(); ptr != this -> end(); ++ptr)
        *ptr *= v;
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
    bool operator == (const Matrix<T>& v) const{
      if (n != v.n || m != v.m) return false;
      auto ptr = this -> begin();
      for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
        if (*ptr != *p) return false;
      return true;
    }
    bool operator != (const Matrix<T>& v) const{
      return !(*this == v);
    }

  public:
    template <typename U>
    Matrix<decltype(T() + U())> operator+ (const Matrix<U> &v) const{
      ASSERT(n == v.rowLength() && m == v.columnLength(), std::invalid_argument("operator+"));
      Matrix<decltype(T() + U())> res(*this);
      auto ptr = v.begin();
      for (auto p = res.begin(); p != res.end(); ++p, ++ptr)
        *p += *ptr;
      return res;
    }

    template <typename U>
    Matrix<decltype(T() - U())> operator- (const Matrix<U> &v) const{
      ASSERT(n == v.n && m == v.m, std::invalid_argument("operator-"));
      Matrix<decltype(T() + U())> res(*this);
      auto ptr = v.begin();
      for (auto p = res.begin(); p != res.end(); ++p, ++ptr)
        *p -= *ptr;
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
      ASSERT(l.first <= r.first && l.second <= r.second, std::invalid_argument("subMatrix"));
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
    auto ptr = mat.begin();
    for (auto p = res.begin(); p != res.end(); ++p, ++ptr)
      *p = *ptr * x;
    return res;
  }
  
  template <class T, class U>
  Matrix<decltype(T() * U())> operator*(const U &x, const Matrix<T> &mat) {
    Matrix<decltype(T() * U())> res(mat.rowLength(), mat.columnLength());
    auto ptr = mat.begin();
    for (auto p = res.begin(); p != res.end(); ++p, ++ptr)
      *p = *ptr * x;
    return res;
  }

  template <class U, class V>
  Matrix<decltype(U() * V())> operator*(const Matrix<U> &u, const Matrix<V> &v) {
    ASSERT(u.columnLength() == v.rowLength(), std::invalid_argument("matrix multiple"));
    Matrix<decltype(U() * V())> res(u.rowLength(), v.columnLength());
    for (size_t i = 0; i < res.rowLength(); ++i)
      for (size_t k = 0; k < v.rowLength(); ++k)
        for (size_t j = 0; j < res.columnLength(); ++j)
          res[i][j] += u[i][k] * v[k][j];
    return res;
  }
}

namespace sjtu{
  using namespace matrix;
}

#endif //MATRIX_HPP
