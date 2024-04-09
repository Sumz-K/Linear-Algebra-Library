#include "Vector.hpp"
#include <stdexcept>
using namespace std;
template<typename T, int Rows, int Cols>
class Matrix {
private:
    Vector<Vector<T>> data;

public:
    Matrix() : data(Rows) {
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

    

};


