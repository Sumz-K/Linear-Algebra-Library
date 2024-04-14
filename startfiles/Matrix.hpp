#include "Vector.hpp"
#include <stdexcept>
using namespace std;


namespace concepthelper{
    template<typename T>
    concept Arithmetic=std::is_arithmetic_v<T>;

    template<typename T,typename U>
    concept Same=std::is_same_v<T,U>;
}


template<typename T, int Rows, int Cols>
requires concepthelper::Arithmetic<T>

class Matrix {
private:
    Vector<Vector<T>> data;
    int rows;
    int cols;

public:
    Matrix() : data(Rows){
        for (int i = 0; i < Rows; ++i) {
            data[i] = Vector<T>(Cols);
        }
        
    }

    Matrix(const Matrix<T, Rows, Cols>& other) : data(other.data) {}

    Matrix(initializer_list<initializer_list<T>> elems): data(elems.size()){
        int i=0;
        for(auto row:elems){
            if(row.size()!=Cols){
                throw std::invalid_argument("All rows should have same number of elements");

            }
            int j=0;
            for(auto val:row){
                data[i][j]=val;
                j+=1;
            }
            i+=1;
        }
        
    }

    int numRows() const {
        return Rows;
    }

    int numCols() const {
        return Cols;
    }

    Vector<T>& operator[](int index) {
        if (index < 0 || index >= Rows) {
            throw std::out_of_range("Row index out of range");
        }
        return data[index];
    }

    const Vector<T>& operator[](int index) const {
        if (index < 0 || index >= Rows) {
            throw std::out_of_range("Row index out of range");
        }
        return data[index];
    }

    T& operator()(int i, int j) {
        if (i < 0 || i >= Rows || j < 0 || j >= Cols) {
            throw std::out_of_range("Index out of range");
        }
        return data[i][j];
    }

    const T& operator()(int i, int j) const {
        if (i < 0 || i >= Rows || j < 0 || j >= Cols) {
            throw std::out_of_range("Index out of range");
        }
        return data[i][j];
    }

    Matrix<T, Rows, Cols>& operator=(const Matrix<T, Rows, Cols>& other) {
        if (this != &other) {
            data = other.data;
        }
        return *this;
    }

    Matrix<T, Rows, Cols> operator+(const Matrix<T, Rows, Cols>& other) const {
        Matrix<T, Rows, Cols> result;
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result[i][j] = data[i][j] + other[i][j];
            }
        }
        return result;
    }

    Matrix<T, Rows, Cols> operator-(const Matrix<T, Rows, Cols>& other) const {
        Matrix<T, Rows, Cols> result;
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result[i][j] = data[i][j] - other[i][j];
            }
        }
        return result;
    }

    template<typename U>
    Matrix<T,Rows,Cols> operator+(const U&scalar){
        if constexpr (!is_arithmetic_v<U>){
            throw std::invalid_argument("Cant add a scalar of this type\n");
        }
        Matrix<T,Rows,Cols>result;
        for(int i=0;i<Rows;i++){
            result[i]=data[i]+scalar;
        }
        return result;
    }

    template<typename U>
    Matrix<T,Rows,Cols> operator-(const U&scalar){
        if constexpr (!is_arithmetic_v<U>){
            throw std::invalid_argument("Cant subtract a scalar of this type\n");
        }
        Matrix<T,Rows,Cols>result;
        for(int i=0;i<Rows;i++){
            result[i]=data[i]-scalar;
        }
        return result;
    }

    template<typename U>
    Matrix<T, Rows, Cols> operator*(const U& scalar) const {
        if constexpr (!is_arithmetic_v<U>){
            throw std::invalid_argument("Cant multiply a scalar of this type\n");
        }
        Matrix<T,Rows,Cols>result;
        for(int i=0;i<Rows;i++){
            result[i]=data[i]*scalar;
        }
        return result;
    }



    template<typename U>
    Matrix<T, Rows, Cols> operator/(const U& scalar) const {
        if constexpr (!is_arithmetic_v<U>){
            throw std::invalid_argument("Cant divide by a scalar of this type\n");
        }
        Matrix<T,Rows,Cols>result;
        for(int i=0;i<Rows;i++){
            result[i]=data[i]/scalar;
        }
        return result;
    }



    Matrix<T, Cols, Rows> transpose() const {
        Matrix<T, Cols, Rows> result;
        for (int i = 0; i < Cols; ++i) {
            for (int j = 0; j < Rows; ++j) {
                result[i][j] = data[j][i];
            }
        }
        return result;
    }

    void print()  {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                std::cout << (*this)(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }
    T accumulate(){
            T sum=0;
            for(int i=0;i<this->numRows();i++){
                Vector<T> temp=(*this)[i];
                T temp_sum=temp.accumulate();
                sum+=temp_sum;
            }
            return sum;
    }
    
    template <int SlicedRows, int SlicedCols>
    Matrix<T, SlicedRows, SlicedCols> slice(int startrow, int endrow, int startcol, int endcol)
    {
        static_assert(SlicedRows >= 0 && SlicedCols >= 0,
        "Invalid slice size");

        Matrix<T, SlicedRows, SlicedCols> res;
        if(SlicedRows!=endrow-startrow || SlicedCols!=endcol-startcol){
            cout<<"Invalid slicing\n";
            return res;
        }

        

        for (int i = startrow; i < endrow && i - startrow < SlicedRows; ++i)
        {
            for (int j = startcol; j < endcol && j - startcol < SlicedCols; ++j)
            {
                res[i - startrow][j - startcol] = data[i][j];
            }
        }
        return res;
    }

    void sort(){
        
        for(int i=0;i<Rows;i++){
            Vector<int>rowvec=(*this)[i];
            rowvec.sort();
            
            (*this)[i]=rowvec;
        }
        
    }

    

    

};




template<typename T,typename U,int r1,int c1,int r2,int c2>
requires concepthelper::Arithmetic<T> && concepthelper::Arithmetic<U> && concepthelper:: Same<T,U> && (c1==r2)
Matrix<T,r1,c2> multiply(Matrix<T,r1,c1> m1, Matrix<U,r2,c2> m2){
    Matrix<T,r1,c2> res;
    for(int i=0;i<r1;i++){
        for(int j=0;j<c2;j++){
            for(int k=0;k<r2;k++){
                res[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return res;
}



template<typename T,int r,int c>
requires concepthelper::Arithmetic<T>
Matrix<T,r,c> zeros(){
    Matrix<T,r,c> m;
    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            m[i][j]=static_cast<T>(0);
        }
    }
    return m;
}


template<typename T,int r,int c>
requires concepthelper::Arithmetic<T>
Matrix<T,r,c> ones(){
    Matrix<T,r,c> m;
    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            m[i][j]=static_cast<T>(1);
        }
    }
    return m;
}

template<typename T,int r,int c>
requires concepthelper::Arithmetic<T> && (r==c)
Matrix<T,r,c> identitymatrix(){
    Matrix<T,r,c> m;
    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            if(i==j)
                m[i][j]=static_cast<T>(1);
            else 
                m[i][j]=static_cast<T>(0);
        }
    }
    return m;
}


template<typename T>
requires concepthelper::Arithmetic<T>
T determinant(Matrix<T, 2, 2> mat){
    return mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
}

template<typename T, int Rows, int Cols>
requires concepthelper::Arithmetic<T> && (Rows == Cols) && (Rows > 2)
T determinant(const Matrix<T, Rows, Cols>& mat) {
    T det = static_cast<T>(0);

    for (int i = 0; i <Cols; ++i) {
        Matrix<T,Rows-1,Cols-1> minor;
        for(int j=1;j<Rows;j++){
            for(int k=0;k<Cols;k++){
                if(i==k) continue;
                else if(k<i){
                    minor[j-1][k]=mat[j][k];
                }
                else if(k>i){
                    minor[j-1][k-1]=mat[j][k];
                }
            }
        }
        T temp_det=determinant(minor);
        if(i%2==0){
            det+=temp_det*mat[0][i];
        }
        else{
            det-=temp_det*mat[0][i];
        }
    }

    return det;
}
template<typename T, int Rows, int Cols>
requires (concepthelper::Arithmetic<T>)
Matrix<T, Rows - 1, Cols - 1> minorf(const Matrix<T, Rows, Cols>& mat, int row, int col) {
    Matrix<T, Rows - 1, Cols - 1> result;
    for (int i = 0, r = 0; i < Rows; ++i) {
        if (i == row) continue;
        for (int j = 0, c = 0; j < Cols; ++j) {
            if (j == col) continue;
            result[r][c] = mat[i][j];
            ++c;
        }
        ++r;
    }
    return result;
}
template<typename T, int Rows, int Cols>
Matrix<T, Rows, Cols> adjoint(const Matrix<T, Rows, Cols>& mat) {
    Matrix<T, Rows, Cols> adj;
    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            Matrix<T, Rows - 1, Cols - 1> minor;
            int r = 0, c = 0;
            for (int k = 0; k < Rows; ++k) {
                if (k == i) continue;
                for (int l = 0; l < Cols; ++l) {
                    if (l == j) continue;
                    minor[r][c] = mat[k][l];
                    ++c;
                }
                ++r;
                c = 0;
            }
            T cofactor = determinant(minor);
            if ((i + j) % 2 != 0) cofactor = -cofactor;
            adj[j][i] = cofactor;
        }
    }
    return adj;
}
template<typename T, int Rows, int Cols>
Matrix<T, Rows, Cols> inverse(const Matrix<T, Rows, Cols>& mat) {
    Matrix<T, Rows, Cols> adj = adjoint(mat);
    T det =determinant(mat);
    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            adj[i][j] /= det;
        }
    }
    return adj;
}
template<typename T>
requires concepthelper::Arithmetic<T>
Matrix<T, 2, 1> solveQuadraticEquation(T a, T b, T c) {
    Matrix<T, 2, 1> roots;
    T discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        throw std::invalid_argument("Complex roots are not supported");
    }
    roots[0][0] = (-b + std::sqrt(discriminant)) / (2 * a);
    roots[1][0] = (-b - std::sqrt(discriminant)) / (2 * a);
    return roots;
}
template<typename T>
requires concepthelper::Arithmetic<T>
Matrix<T, 3, 1> solveCubicEquation(T a, T b, T c, T d) {
    Matrix<T, 3, 1> roots;

    if (a == 0) {
        throw std::invalid_argument("Not a cubic equation");
    }

    // Normalize coefficients
    T inv_a = 1 / a;
    b *= inv_a;
    c *= inv_a;
    d *= inv_a;

    
    T Q = (3 * b - pow(2 * c, 2)) / 9;
    T R = (9 * b * c - 27 * d - 2 * pow(2 * c, 3)) / 54;
    T D = pow(Q, 3) + pow(R, 2);

    if (D >= 0) {
        // Three real roots
        T sqrt_D = sqrt(D);
        T S = cbrt(R + sqrt_D);
        T t = cbrt(R - sqrt_D);

        roots[0][0] = -b / 3 + (S + t);
        roots[1][0] = -b / 3 - (S + t) / 2;
        roots[2][0] = roots[1][0];
        roots[1][0] += (sqrt(3.0) * (S - t)) / 2;
        roots[2][0] -= (sqrt(3.0) * (S - t)) / 2;
    } else {
        // One real root, two complex roots
        T theta = acos(R / sqrt(-pow(Q, 3)));
        roots[0][0] = -2 * sqrt(Q) * cos(theta / 3) - b / 3;
        roots[1][0] = -2 * sqrt(Q) * cos((theta + 2 * M_PI) / 3) - b / 3;
        roots[2][0] = -2 * sqrt(Q) * cos((theta - 2 * M_PI) / 3) - b / 3;
    }

    return roots;
}
