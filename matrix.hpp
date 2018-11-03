#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cassert>
#include <stdexcept>
#include <functional>

namespace matrix{
#define ASSERT(expr, message) if (!(expr)) { throw std::invalid_argument(message); }
  
    template <typename T>
    class Matrix{
        using size_p = std::pair<size_t, size_t>;
    
        T *x;
        size_t n, m;
        
    public:
    
        T* operator[] (size_t u) { return x + u * m; }
        const T* operator[] (size_t u) const{ return x + u * m; }
    
        T& operator() (size_t i, size_t j) {
            ASSERT(i >= 0 && i < n && j >= 0 && j < m, "operator ()");
            return (*this)[i][j];
        }
        const T& operator() (size_t i, size_t j) const{
            ASSERT(i >= 0 && i < n && j >= 0 && j < m, "operator ()");
            return (*this)[i][j];
        }
        
        Matrix row(size_t i) const{
            ASSERT(i >= 0 && i < n, "row");
            Matrix res(1, m);
            auto p = this -> find(i, 0);
            for (auto ptr = res.begin(); ptr != res.end(); ++ptr, ++p)
                *ptr = *p;
            return res;
        }
        Matrix column(size_t i) const{
            ASSERT(i >= 0 && i < m, "column");
            Matrix res(n, 1);
            auto p = this -> find(0, i);
            for (auto ptr = res.begin(); ptr != res.end(); ++ptr, p += m)
                *ptr = *p;
            return res;
        }

    public:
        Matrix(size_t _n = 1, size_t _m = 1, T v = T())
            :x(new T [_n * _m]), n(_n), m(_m) {
            for (auto ptr = this -> begin(); ptr != this -> end(); ++ptr)
                *ptr = v;
        }
        explicit Matrix(size_p sz, T v = T()) 
            :x(new T [sz.first * sz.second]), n(sz.first) , m(sz.second) {
            for (auto ptr = this -> begin(); ptr != this -> end(); ++ptr)
                *ptr = v;
        }
        Matrix(const std::initializer_list<std::initializer_list<T>> &matrix) {
            n = matrix.size(); m = matrix.begin() -> size();
            
            for (auto &list : matrix) 
                ASSERT(m == list.size(), "initiallizer_list");
            T* ptr = x = new T [n * m];
            for (auto &list : matrix) for (auto &v : list) *ptr++ = v;
        }
        Matrix(Matrix &&v) {
            n = v.n;
            m = v.m;
            x = v.x;
            v.x = nullptr;
        }
        template <typename U>
        Matrix(const Matrix<U> &v) {
            n = v.rowLength();
            m = v.columnLength();
            x = new T [n * m];
            auto ptr = this -> begin();
            for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
                *ptr = *p;
        }
        Matrix(const Matrix &v) {
            n = v.rowLength(); m = v.columnLength();
            x = new T [n * m];
            auto ptr = this -> begin();
            for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
                *ptr = *p;
        }
        
        ~Matrix() {
            if (n || m)
                delete [] x;
        }
        
    public:
        Matrix& operator= (const Matrix &v) {
            if (this == &v) return *this;
            if (x != nullptr) delete [] x;
            n = v.n; m = v.m; x = new T [n * m];
            auto ptr = this -> begin();
            for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
                *ptr = *p;
            return *this;
        }
        Matrix& operator= (Matrix &&v) {
            n = v.n; m = v.m; 
            if (x == v.x) return *this;
            if (x != nullptr) delete [] x;
            x = v.x; v.x = nullptr;
            return *this;
        }
        template <class K>
        Matrix& operator= (const Matrix<K> &v) {
            n = v.rowLength();
            m = v.columnLength();
            if (x != nullptr) delete [] x;
            x = new T [n * m];
            auto ptr = this -> begin();
            for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
                *ptr = *p;
            return *this;
        }
          
    public:
        void clear() { n = m = 0; delete [] x; }
        size_t rowLength() const{ return n; }
        size_t columnLength() const{ return m; }
        size_p size() const{ return {rowLength(), columnLength()}; }
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
        void resize(size_p sz, T v = T()) { resize(sz.first, sz.second, v); }
    public:
        Matrix operator- () const{
            Matrix res(n, m);
            auto ptr = this -> begin();
            for (auto p = res.begin(); p != res.end(); ++p, ++ptr)
                *p = -*ptr;
            return res;
        }
        template <typename K>
        Matrix operator+= (const Matrix<K> &v) {
            ASSERT(n == v.n && m == v.m, "operator+=");
            auto ptr = this -> begin();
            for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
                *ptr += *p;
            return *this;
        }
        template <typename K>
        Matrix operator-= (const Matrix<K> &v) {
            ASSERT(n == v.n && m == v.m, "operator-=");
            auto ptr = this -> begin();
            for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
                *ptr -= *p;
            return *this;
        }
        template <typename K>
        Matrix operator*= (const K& v) {
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
        bool operator == (const Matrix& v) const{
            if (n != v.n || m != v.m) return false;
            auto ptr = this -> begin();
            for (auto p = v.begin(); p != v.end(); ++p, ++ptr)
                if (*ptr != *p) return false;
            return true;
        }
        bool operator != (const Matrix& v) const{
            return !(*this == v);
        }

    public:
        template <typename U>
        Matrix<decltype(T() + U())> operator+ (const Matrix<U> &v) const{
            ASSERT(n == v.rowLength() && m == v.columnLength(), "operator+");
            Matrix<decltype(T() + U())> res(*this);
            auto ptr = v.begin();
            for (auto p = res.begin(); p != res.end(); ++p, ++ptr)
                *p += *ptr;
            return res;
        }

        template <typename U>
        Matrix<decltype(T() - U())> operator- (const Matrix<U> &v) const{
            ASSERT(n == v.n && m == v.m, "operator-");
            Matrix<decltype(T() - U())> res(*this);
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
            iterator (size_type m_origin, size_type m_current,
                      const pointer& beg_origin,
                      const pointer& beg_current,
                      const pointer& ptr)
                :m_origin(m_origin), m_current(m_current),
                 beg_origin(beg_origin), beg_current(beg_current), ptr(ptr) {}
            
        private:
            size_type m_origin, m_current;
            pointer beg_origin, beg_current, ptr;

            difference_type distance() const{
                difference_type diff = beg_current - beg_origin;
                difference_type head = diff % m_origin;
                pointer left = beg_origin + diff - head;
                return (ptr - left) / m_origin * m_current + (ptr - left) % m_origin - head;
            }
        public:
            pointer base() const{ return ptr; }
            difference_type operator- (const iterator &p) const {
                ASSERT(m_origin == p.m_origin && m_current == p.m_current, "iterator error");
                return distance() - p.distance();
            } 
            iterator operator+ (difference_type offset) const{
                difference_type left = m_current - (ptr - beg_current) % m_origin;
                return iterator(m_origin, m_current, beg_origin, beg_current, ptr + offset + 
                                (offset < left ? 0 : (m_origin - m_current) * ((offset - left) / m_current + 1)));
            }
            iterator operator- (difference_type offset) const{
                difference_type left = (ptr - beg_current) % m_origin;
                return iterator(m_origin, m_current, beg_origin, beg_current, ptr - offset +
                                (offset < left ? 0 : (m_current - m_origin) * ((offset - left) / m_current + 1) + left));
            }
            iterator& operator+= (difference_type offset) { return *this = *this + offset; }
            iterator& operator-= (difference_type offset) { return *this = *this - offset; }
            iterator& operator++ () { return *this += 1; }
            iterator operator++(int) { iterator A = *this; ++(*this); return A; }
            iterator& operator-- () { return *this -= 1; }
            iterator operator--(int) { iterator A = *this; --(*this); return A; }

            reference operator* () const { return *ptr; }
            pointer operator->() const { return ptr; }
            bool operator==(const iterator &o) const{
                return
                    m_origin == o.m_origin &&
                    m_current == o.m_current &&
                    beg_origin == o.beg_origin &&
                    beg_current == o.beg_current &&
                    ptr == o.ptr;
            }
            bool operator!= (const iterator &o) const{ return !(*this == o); }
        };

#define matrix_nullptr iterator(0, 0, nullptr, nullptr, nullptr)
        iterator begin() const{
            return n == 0 || m == 0 ? matrix_nullptr : iterator(m, m, x, x, x);
        }
        iterator end() const{
            return n == 0 || m == 0 ? matrix_nullptr : iterator(m, m, x, x, x + n * m);
        }
        iterator find(size_t i, size_t j) const{
            return iterator(m, m, x, x, x + i * m + j);
        }
		
        std::pair<iterator, iterator> subMatrix(size_p l, size_p r) const{
            ASSERT(l.first <= r.first && l.second <= r.second, "subMatrix");
            iterator beg = iterator(m, r.second - l.second + 1,
                                    this->begin().base(), this->begin().base() + l.first * m + l.second,
                                    this->begin().base() + l.first * m + l.second);
            return {beg, beg + (r.second - l.second + 1) * (r.first - l.first + 1)};
        }
    };
  
    template <class U, class V>
    Matrix<decltype(U() * V())> operator*(const Matrix<U> &mat, const V &x) {
        Matrix<decltype(U() * V())> res(mat.rowLength(), mat.columnLength());
        auto ptr = mat.begin();
        for (auto p = res.begin(); p != res.end(); ++p, ++ptr)
            *p = *ptr * x;
        return res;
    }
  
    template <class U, class V>
    Matrix<decltype(U() * V())> operator*(const V &x, const Matrix<U> &mat) {
        Matrix<decltype(U() * V())> res(mat.rowLength(), mat.columnLength());
        auto ptr = mat.begin();
        for (auto p = res.begin(); p != res.end(); ++p, ++ptr)
            *p = *ptr * x;
        return res;
    }

    template <class U, class V>
    Matrix<decltype(U() * V())> operator*(const Matrix<U> &u, const Matrix<V> &v) {
        ASSERT(u.columnLength() == v.rowLength(), "matrix multiple");
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
