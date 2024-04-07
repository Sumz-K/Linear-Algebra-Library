#include"Vector_scratch.hpp"
#include<type_traits>

#include<math.h>



template<typename T>
class Vector: public Vector_scratch<T>{
    using Vector_scratch<T>::Vector_scratch;

public:
    Vector<T>& operator +(T elem){
        for(int i=0;i<this->current;i++){
            this->container[i]+=elem;
        }
        return (*this);
    }
    Vector<T>& operator -(T elem){
        for(int i=0;i<this->current;i++){
            this->container[i]-=elem;
        }
        return (*this);
    }
    Vector<T>& operator *(T elem){
        for(int i=0;i<this->current;i++){
            this->container[i]*=elem;
        }
        return (*this);
    }
    Vector<T>& operator /(T elem){
        for(int i=0;i<this->current;i++){
            this->container[i]/=elem;
        }
        return (*this);
    }

    Vector<T> reverse();

    double magnitude();







};



template<typename T>
Vector<T> Vector<T>::reverse(){
    for(int i=0;i<this->current/2;i++){
        std::swap(this->container[i], this->container[this->current - i - 1]);
    }
    return (*this);
}



template<typename T>
double Vector<T>::magnitude(){
    double sum=0;
    for(int i=0;i<this->size();i++){
        sum+=this->container[i]*this->container[i];
    }
    return sqrt(sum);
}


template<typename T>
struct IsVector: std::false_type{};

template<typename T>
struct IsVector<Vector<T>>:std::true_type{};

template<typename Vec>
void normalise(Vec &v){
    if constexpr (IsVector<Vec>::value){
        double mag=v.magnitude();
        for(int i=0;i<v.size();i++){
            v[i]=v[i]/mag;
        }

    }
    else{
        throw ("Normalisation not supported for non vector types\n");
    }
}


auto dot=[]<typename T1, typename T2>(Vector<T1> v1,Vector<T2> v2){
    if(v1.size()!=v2.size()){
        throw("The two vectors dont have the same number of elements\n");
    }

    if constexpr (std::is_same_v<T1,T2> && std::is_integral_v<T1>){
        return static_cast<int>(dothelper(v1,v2));
    }
    else{
        return dothelper(v1,v2);
    }


};

template<typename T1,typename T2>
double dothelper(Vector<T1> v1,Vector<T2> v2){
    double res=0;

    for(int i=0;i<v1.size();i++){
        res+=v1[i]*v2[i];
    }
    return res;
}