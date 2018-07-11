//
//  matrix.h
//
//  Created by Mira Shalah on 9/27/17.
//
// Matrix class for basic matrix operations

#ifndef mira_twigs_matrix_h
#define mira_twigs_matrix_h

#include <iostream>

template<typename T>
class Matrix
{
public:
        T m11_, m12_, m13_, m21_, m22_, m23_, m31_, m32_, m33_;
        
        Matrix( void ){};
        
        Matrix(
               T m11, T m12, T m13,
               T m21, T m22, T m23,
               T m31, T m32, T m33 ):
        m11_( m11 ), m12_( m12 ), m13_( m13 ),
        m21_( m21 ), m22_( m22 ), m23_( m23 ),
        m31_( m31 ), m32_( m32 ), m33_( m33 )
        {
        };
        
        Matrix( const Matrix& m ):
        m11_( m.m11_ ), m12_( m.m12_ ), m13_( m.m13_ ),
        m21_( m.m21_ ), m22_( m.m22_ ), m23_( m.m23_ ),
        m31_( m.m31_ ), m32_( m.m32_ ), m33_( m.m33_ )
        {
        };
        
        ~Matrix( void ){};
        
        Matrix& operator=( const Matrix& other_matrix )
        {
            m11_ = other_matrix.m11_; m12_ = other_matrix.m12_; m13_ = other_matrix.m13_;
            m21_ = other_matrix.m21_; m22_ = other_matrix.m22_; m23_ = other_matrix.m23_;
            m31_ = other_matrix.m31_; m32_ = other_matrix.m32_; m33_ = other_matrix.m33_;
            return *this;
        };
        
        bool operator==( const Matrix& other_matrix ) const
        {
            return m11_ == other_matrix.m11_ && m12_ == other_matrix.m12_ && m13_ == other_matrix.m13_
            && m21_ == other_matrix.m21_ && m22_ == other_matrix.m22_ && m23_ == other_matrix.m23_
            && m31_ == other_matrix.m31_ && m32_ == other_matrix.m32_ && m33_ == other_matrix.m33_;
        };
        
        Matrix operator+( const Matrix& other_matrix ) const
        {
            Matrix out_matrix(
                              m11_ + other_matrix.m11_, m12_ + other_matrix.m12_, m13_ + other_matrix.m13_,
                              m21_ + other_matrix.m21_, m22_ + other_matrix.m22_, m23_ + other_matrix.m23_,
                              m31_ + other_matrix.m31_, m32_ + other_matrix.m32_, m33_ + other_matrix.m33_ );
            return out_matrix;
        };
        
        Matrix& operator+=( const Matrix& other_matrix ) const
        {
            m11_ += other_matrix.m11_; m12_ += other_matrix.m12_; m13_ += other_matrix.m13_;
            m21_ += other_matrix.m21_; m22_ += other_matrix.m22_; m23_ += other_matrix.m23_;
            m31_ += other_matrix.m31_; m32_ += other_matrix.m32_; m33_ += other_matrix.m33_;
            return *this;
        };
        
        Matrix operator-( const Matrix& other_matrix ) const
        {
            Matrix out_matrix(
                              m11_ - other_matrix.m11_, m12_ - other_matrix.m12_, m13_ - other_matrix.m13_,
                              m21_ - other_matrix.m21_, m22_ - other_matrix.m22_, m23_ - other_matrix.m23_,
                              m31_ - other_matrix.m31_, m32_ - other_matrix.m32_, m33_ - other_matrix.m33_ );
            return out_matrix;
        };
        
        Matrix& operator-=( const Matrix& other_matrix ) const
        {
            m11_ -= other_matrix.m11_; m12_ -= other_matrix.m12_; m13_ -= other_matrix.m13_;
            m21_ -= other_matrix.m21_; m22_ -= other_matrix.m22_; m23_ -= other_matrix.m23_;
            m31_ -= other_matrix.m31_; m32_ -= other_matrix.m32_; m33_ -= other_matrix.m33_;
            return *this;
        };
        
        Matrix operator*( const Matrix& other_matrix ) const
        {
            Matrix out_matrix;
            out_matrix.m11_ = m11_ * other_matrix.m11_ + m12_ * other_matrix.m21_ + m13_ * other_matrix.m31_;
            out_matrix.m12_ = m11_ * other_matrix.m12_ + m12_ * other_matrix.m22_ + m13_ * other_matrix.m32_;
            out_matrix.m13_ = m11_ * other_matrix.m13_ + m12_ * other_matrix.m23_ + m13_ * other_matrix.m33_;
            out_matrix.m21_ = m21_ * other_matrix.m11_ + m22_ * other_matrix.m21_ + m23_ * other_matrix.m31_;
            out_matrix.m22_ = m21_ * other_matrix.m12_ + m22_ * other_matrix.m22_ + m23_ * other_matrix.m32_;
            out_matrix.m23_ = m21_ * other_matrix.m13_ + m22_ * other_matrix.m23_ + m23_ * other_matrix.m33_;
            out_matrix.m31_ = m31_ * other_matrix.m11_ + m32_ * other_matrix.m21_ + m33_ * other_matrix.m31_;
            out_matrix.m32_ = m31_ * other_matrix.m12_ + m32_ * other_matrix.m22_ + m33_ * other_matrix.m32_;
            out_matrix.m33_ = m31_ * other_matrix.m13_ + m32_ * other_matrix.m23_ + m33_ * other_matrix.m33_;
            return out_matrix;
        };
        
        Matrix operator*( const T& constant_value ) const
        {
            Matrix out_matrix(
                              m11_ * constant_value, m12_ * constant_value, m13_ * constant_value,
                              m21_ * constant_value, m22_ * constant_value, m23_ * constant_value,
                              m31_ * constant_value, m32_ * constant_value, m33_ * constant_value );
            return out_matrix;
        };
    
        Matrix operator/( const T& constant_value ) const
        {
            Matrix out_matrix(
                              m11_ / constant_value, m12_ / constant_value, m13_ / constant_value,
                              m21_ / constant_value, m22_ / constant_value, m23_ / constant_value,
                              m31_ / constant_value, m32_ / constant_value, m33_ / constant_value );
            return out_matrix;
        };
    };
    
    template<typename T>
    Matrix<T> operator*( const T& constant_value, const Matrix<T>& m )
    {
        Matrix<T> out_matrix( m );
        out_matrix.m11_ *= constant_value;
        out_matrix.m12_ *= constant_value;
        out_matrix.m13_ *= constant_value;
        out_matrix.m21_ *= constant_value;
        out_matrix.m22_ *= constant_value;
        out_matrix.m23_ *= constant_value;
        out_matrix.m31_ *= constant_value;
        out_matrix.m32_ *= constant_value;
        out_matrix.m33_ *= constant_value;
        return out_matrix;
    };

    template<typename T>
    T determinant( const Matrix<T>& m )
    {
        return
        m.m11_ * ( m.m22_ * m.m33_ - m.m23_ * m.m32_ )
        - m.m12_ * ( m.m21_ * m.m33_ - m.m23_ * m.m31_ )
        + m.m13_ * ( m.m21_ * m.m32_ - m.m22_ * m.m31_ );
    };

    template<typename T> Matrix<T> Inverse(const Matrix<T>& m)
    {
        double det = determinant(m);
        double invdet = 1 / det;
        double
        m11 = (m.m22_ * m.m33_ - m.m32_ * m.m23_) * invdet,
        m12 = (m.m13_ * m.m32_ - m.m12_ * m.m33_) * invdet,
        m13 = (m.m12_ * m.m23_ - m.m13_ * m.m22_) * invdet,
        m21 = (m.m23_ * m.m31_ - m.m21_ * m.m33_) * invdet,
        m22 = (m.m11_ * m.m33_ - m.m13_ * m.m31_) * invdet,
        m23 = (m.m21_ * m.m13_ - m.m11_ * m.m23_) * invdet,
        m31 = (m.m21_ * m.m32_ - m.m31_ * m.m22_) * invdet,
        m32 = (m.m31_ * m.m12_ - m.m11_ * m.m32_) * invdet,
        m33 = (m.m11_ * m.m22_ - m.m21_ * m.m12_) * invdet;
        return Matrix<T> (m11, m12, m13, m21, m22, m23, m31, m32, m33);
    }

    template<typename T>
    std::ostream& operator<<( std::ostream& stream, const Matrix<T>& m )
    {
        stream << "["
        << m.m11_ << " " << m.m12_ << " " << m.m13_ << ", "
        << m.m21_ << " " << m.m22_ << " " << m.m23_ << ", "
        << m.m31_ << " " << m.m32_ << " " << m.m33_ << "]";
        return stream;
    };


#endif
