#include<iostream>
#include"Matrix.h"


template<typename T, int row, int col>
class Compute{
    public:
        Matrix<T,row,col> add(Matrix<T,row,col> m1, Matrix<T,row,col> m2){
            Matrix<T,row,col> res;
            for(int i=0;i<row;i++){
                for(int j=0;j<col;j++){
                    int sum=m1(i,j)+m2(i,j);
                    res.put(sum,i,j);
                }
            }
            return res;
        }
        Matrix<T,row,col> sub(Matrix<T,row,col> m1, Matrix<T,row,col> m2){
            Matrix<T,row,col> res;
            for(int i=0;i<row;i++){
                for(int j=0;j<col;j++){
                    int sum=m1(i,j)-m2(i,j);
                    res.put(sum,i,j);
                }
            }
            return res;
        }


        T accumulate(Matrix<T,row,col> m){
            T sum=0;
            for(int i=0;i<m. numRows();i++){
                Vector_scratch<T> temp=m[i];
                T temp_sum=temp.accumulate();
                sum+=temp_sum;
            }
            return sum;
        }
};


int main(){
    Matrix<int,2,2> m1;
    Matrix<int,2,2> m2;
    int val1=1;
    int val2=4;
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            m1.put(val1++,i,j);
            m2.put(val2--,i,j);
        }
    }
    

    Compute<int,2,2> c;
    Matrix<int,2,2> m3=c.add(m1,m2);
    m3.print();
    Matrix<int,2,2> m4=c.sub(m1,m2);
    //m4.print();


    cout<<c.accumulate(m3);
    return 0;

}