#include "Vector.hpp"
#include <stdexcept>
using namespace std;



//concepts defined in a separate namespace scope
namespace concepthelper{
    template<typename T>
    concept Arithmetic=std::is_arithmetic_v<T>;

    template<typename T,typename U>
    concept Same=std::is_same_v<T,U>;
}


//all Matrices must be of Arithmetic type

template<typename T, int Rows, int Cols>
requires concepthelper::Arithmetic<T>

class Matrix {
private:
    //data is the container
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

//overload to access a particular row of a matrix which is a vector
    Vector<T>& operator[](int index) {
        if (index < 0 || index >= Rows) {
            throw std::out_of_range("Row index out of range");
        }
        return data[index];
    }
//const qualified
    const Vector<T>& operator[](int index) const {
        if (index < 0 || index >= Rows) {
            throw std::out_of_range("Row index out of range");
        }
        return data[index];
    }

//another way to access a particular element of a matrix is to do mat(i,j)
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


//simple operator overloads
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
    Matrix<T, Rows, Cols> operator*(const U& scalar)  {
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
        Matrix<T,Rows,Cols>result;
        if constexpr (!is_arithmetic_v<U>){
            throw std::invalid_argument("Cant divide by a scalar of this type\n");
        }
        if(scalar==0) {
            throw std::invalid_argument("Cant divide by 0\n");
            
        }
        
        for(int i=0;i<Rows;i++){
            result[i]=data[i]/scalar;
        }
        return result;
    }

    template<typename U>
    Matrix<T, Rows, Cols> operator/(const U& scalar) {
        Matrix<T,Rows,Cols>result;
        if constexpr (!is_arithmetic_v<U>){
            throw std::invalid_argument("Cant divide by a scalar of this type\n");
            
        }
        if(scalar==0) {
            throw std::invalid_argument("Cant divide by 0\n");
            
        }
        
        for(int i=0;i<Rows;i++){
            result[i]=data[i]/scalar;
        }
        return result;
    }

//to find powers of matrix

    Matrix<T,Rows,Cols> operator^(int scalar) const{
        Matrix<T,Rows,Cols> res=(*this);
        while(scalar>1){
            res=multiply(res,res);
            scalar-=1;
        }
        return  res;

    }
    Matrix<T,Rows,Cols> operator^(int& scalar) {
        Matrix<T,Rows,Cols> res=this;
        while(scalar>0){
            res=multiply(res,res);
            scalar-=1;
        }
        return  res;

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
    void print() const {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                std::cout << (*this)(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }


//slicing a matrix to form a submatrix, given number of rows and cols to be sliced
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

//this will sort each individual row of the matrix
    void sort(){
        
        for(int i=0;i<Rows;i++){
            Vector<int>rowvec=(*this)[i];
            rowvec.sort();
            
            (*this)[i]=rowvec;
        }
        
    }

    


    

    

};



//overload to print using cout statement
template<typename T, int Rows, int Cols>
std::ostream& operator<<(std::ostream& os, const Matrix<T, Rows, Cols>& matrix) {
    matrix.print();
    return os;
}


//matrix multiplication
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



//a matrix full of zeroes

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


//a matrix full of ones

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


//identity matrix
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



//partial specialisation for 2x2 will use this as base case for recursion
template<typename T>
requires concepthelper::Arithmetic<T>
T determinant(Matrix<T, 2, 2> mat){
    return mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
}


//the original determinent function, no specialisation, recursive in nature
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

//to find cofactor of a matrix, will use this in multiple other functions

template<typename T,int Rows,int Cols>
requires (Rows==Cols && concepthelper::Arithmetic<T>)
Matrix<T,Rows,Cols> cofactor(Matrix<T,Rows,Cols> mat){
    Matrix<T,Rows,Cols> cf;
    for(int i=0;i<Rows;i++){
        for(int j=0;j<Cols;j++){
            Matrix<T,Rows-1,Cols-1>minor;
            int minorRow=0;
            int minorCol=0;

            for(int k=0;k<Rows;k++){
                if(k==i) continue;
                minorCol=0;
                for(int l=0;l<Cols;l++){
                    if(l==j) continue;
                    minor[minorRow][minorCol]=mat[k][l];
                    minorCol+=1;
                }
                minorRow+=1;
            }
            T minorDet=determinant(minor);

            int sign=pow(-1,i+j+2);
            cf[i][j]=minorDet*sign;

        }
    }
    return cf;
}



//to find adjoint of a matrix

template<typename T, int Rows, int Cols>
requires (Rows==Cols && concepthelper::Arithmetic<T>)
Matrix<T, Rows, Cols> adjoint(Matrix<T, Rows, Cols>& mat) {

    Matrix<T,Rows,Cols> dummy=cofactor(mat);
    Matrix<T,Rows,Cols>adj=dummy;
    for(int i=0;i<Rows;i++){
        for(int j=0;j<Cols;j++){
            adj[i][j]=dummy[j][i];
        }
    }

    return adj;


}



//to find inverse of a matrix
template<typename T,int Rows,int Cols>
requires (Rows==Cols && concepthelper::Arithmetic<T>)
Matrix<double,Rows,Cols> inverse(Matrix<T,Rows,Cols> mat){
    double det=determinant(mat);

    if(static_cast<int>(det)==0){
        throw std::invalid_argument("This matrix is not invertible\n");
    }
    auto adj=adjoint(mat);
    Matrix<double, Rows, Cols> result;
    for(int i=0;i<Rows;i++){
        for(int j=0;j<Cols;j++){
            result[i][j] = static_cast<double>(adj[i][j]) / det;
        }
    }
    return result;
}

//partial specialisation for inverse, 2x2 case the inverse is simple hence specialised it

template<typename T>
Matrix<double,2,2> inverse(Matrix<T,2,2> mat){
    Matrix<double,2,2> res;
    res[0][0]=mat[1][1];
    res[1][1]=mat[0][0];
    res[0][1]=-mat[0][1];
    res[1][0]=-mat[1][0];

    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            res[i][j]=static_cast<double>(res[i][j])/determinant(mat);
        }
    }
    return res;
}



//to find minor of a matrix


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


//to solve a quadratic equation given the coefficients of each of its terms
template<typename T>
requires concepthelper::Arithmetic<T>
Matrix<double, 2, 1> solveQuadraticEquation(T a, T b, T c) {
    Matrix<double, 2, 1> roots;
    T discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        throw std::invalid_argument("Complex roots are not supported");
    }
    roots[0][0] = (-b + sqrt(discriminant)) / (2 * a);
    roots[1][0] = (-b - sqrt(discriminant)) / (2 * a);
    return roots;
}




//to solve a system of linear equations of the form Ax=b
template<typename T,int Rows,int Cols>
requires concepthelper::Arithmetic<T>
Matrix<double,Rows,1> solvelinear(Matrix<T,Rows,Cols> A,Matrix<double,Rows,1> b){
    //Ax=b
    // need to find x
    //double because inverse always returns double
    T det=determinant(A);
    if(det==0){
        throw std::invalid_argument("This system of equations is not solveable\n");
    }
    else{
        auto inv=inverse(A);
        inv.print();
        return multiply(inv,b);
    }

}




