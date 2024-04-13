#include"Vector_scratch.hpp"
#include<type_traits>
#include"Map.hpp"
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

    const Vector<T>& operator +(T elem) const{
        for(int i=0;i<this->current;i++){
            this->container[i]+=elem;
        }
        return (*this);
    }
    const Vector<T>& operator -(T elem) const{
        for(int i=0;i<this->current;i++){
            this->container[i]-=elem;
        }
        return (*this);
    }
    const Vector<T>& operator *(T elem) const{
        for(int i=0;i<this->current;i++){
            this->container[i]*=elem;
        }
        return (*this);
    }
    const Vector<T>& operator /(T elem) const{
        for(int i=0;i<this->current;i++){
            this->container[i]/=elem;
        }
        return (*this);
    }

    Vector<T> reverse();

    double magnitude();

    template<typename U> friend Vector<T> square();

    template<typename U> friend Vector<T> cube();

    void sort(){
        int left=0;
        int right=this->size()-1;
        sorthelper(left,right);
        
    }

    void sorthelper(int left,int right){
        if(left<right){
            int mid=left+(right-left)/2;
            sorthelper(left,mid);
            sorthelper(mid+1,right);
            merge(left,mid,right);
        }
    }

    void merge(int left, int mid, int right){
        Vector<T> left_subarr(mid - left + 1);
        Vector<T> right_subarr(right - mid);

        for (int i = left; i <= mid; i++){
            left_subarr.push_back(this->container[i]);
        }
        for (int i = mid + 1; i <= right; i++){
            right_subarr.push_back(this->container[i]);
        }

        int i = 0;
        int j = 0;
        int res = left;

        while (i < left_subarr.size() && j < right_subarr.size()){
            if (left_subarr[i] < right_subarr[j]){
                this->container[res++] = left_subarr[i++];
            }
            else{
                this->container[res++] = right_subarr[j++];
            }
        }


        while (i < left_subarr.size()){
            this->container[res++] = left_subarr[i++];
        }


        while (j < right_subarr.size()){
            this->container[res++] = right_subarr[j++];
        }
    }




};



template<typename T>
double mean(Vector<T> vec){
    double sum=0;
    for(int i=0;i<vec.size();i++){
        sum+=vec[i];
    }
    return sum/vec.size();
}

template<typename T>
double median(Vector<T> vec){
    vec.sort();
    int left=0;
    int right=vec.size()-1;
    while(left<right){
        left+=1;
        right-=1;
    }
    if(left==right){
        return static_cast<double>(vec[left]);
    }
    return static_cast<double>((vec[left]+vec[right])/2);
}



template<typename T>

double euclid_distance(Vector<T>v1,Vector<T>v2){
    if(v1.size()!=v2.size()){
        throw std::invalid_argument("Vectors must have the same number of elements to calculate euclidean distance\n");


    }
    else{
        double dist=0.0;
        for(int i=0;i<v1.size();i++){
            dist+=pow((v1[i]-v2[i]),2);
        }
        return sqrt(dist);
    }
}

template<typename T>
T manhattan_distance(Vector<T> v1,Vector<T> v2){
    if(v1.size()!=v2.size()){
        throw std::invalid_argument("Vectors must have the same number of elements to calculate euclidean distance");


    }
    else{
        T dist=0;
        for(int i=0;i<v1.size();i++){
            dist+=abs(v1[i]-v2[i]);
        }
        return static_cast<T>(dist);
    }
}

template<typename T>
Vector<T> square(Vector<T> vec){
    if(!std::is_arithmetic_v<T>){
        throw std::invalid_argument("Type not suitable for exponentiation");
    }
    Vector<T> res;
    for(int i=0;i<vec.size();i++){
        res[i]=vec.container[i]*vec.container[i];
    }
    return res;
}

template<typename T>
Vector<T> cube(Vector<T> vec){
    if(!std::is_arithmetic_v<T>){
        throw std::invalid_argument("Type not suitable for exponentiation");
    }
    Vector<T> res;
    for(int i=0;i<vec.size();i++){
        res[i]=vec.container[i]*vec.container[i]*vec.container[i];
    }
    return res;
}



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
