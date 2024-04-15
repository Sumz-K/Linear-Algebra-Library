#include <iostream>
#include <type_traits>
#define INIT 1
#define MAX 100
using namespace std;


//container for a vector
template <typename T>
class Vector_scratch
{
public:
    T *container;  //the datastore
    int capacity;  //maxsize
    int current;   //how many elements currently in the vector

//getters
    T* getContainer(){
        return this->container;
    }

    int getCapacity(){
        return this->capacity;
    }
    int getCurrent(){
        return this->current;
    }

    Vector_scratch(int size)
    {
        capacity = size;
        current = 0;
        container = new T[capacity];

        for (int i = 0; i < capacity; i++)
        {
            container[i] = 0;
        }
    }

    Vector_scratch()
    {
        capacity = INIT;
        current = 0;
        container = new T[capacity];
        for (int i = 0; i < capacity; i++)
        {
            container[i] = 0;
        }
        
    }

    Vector_scratch(initializer_list<T> list)
    {
        capacity = list.size();
        container = new T[capacity];
        int i = 0;
        for (auto &elem : list)
        {
            container[i++] = elem;
        }
        
        current = capacity;
    }

    Vector_scratch(const Vector_scratch<T>& other) : capacity(other.capacity), current(other.current) {
        container = new T[capacity];
        for (int i = 0; i < current; ++i) {
            container[i] = other.container[i];
        }
    }


    void print()
    {
        if(this->capacity==0){
            cout<<"No memory allocated to the vector\n";
            return ;
        }
        for (int i = 0; i < current; i++)
        {
            cout << container[i] << " ";
        }
        cout << "\n";
    }
    void print() const
    {
        if(this->capacity==0){
            cout<<"No memory allocated to the vector\n";
            return ;
        }
        for (int i = 0; i < current; i++)
        {
            cout << container[i] << " ";
        }
        cout << "\n";
    }


//custom push back implementation
    void push_back(T data)
    {
        if (current < capacity)
        {
            container[current] = data; //if theres extra space in the container just append to it
            current += 1;
        }
        else
        { //otherwise create a new container double the size of the older container(as per stl vector), bring back all the elements of the older container and append the new element
            T *expanded = new T[2 * capacity];
            for (int i = 0; i < capacity; i++)
            {
                expanded[i] = container[i];
            }
            delete[] container;
            capacity = capacity * 2;
            container = expanded;
            container[current] = data;
            current += 1;
        }
    }
//pop back is simply a decrement in number of elements, all accesses of the vector use this variable so popped variables are never accessed
    void pop_back()
    {
        current -= 1;
    }

//similar to pop back,just at the front
    void pop_front()
    {
        for (int i = 0; i < current; i++)
        {
            container[i] = container[i + 1];
        }
        current -= 1;
    }

//operator overload to access the elements of the container easily
    T &operator[](int index)
    {
        if (index >= capacity || index < 0)
        {
            int newSize = (index >= capacity) ? index + 1 : capacity * 2;
            T *expanded = new T[newSize];
            for (int i = 0; i < newSize; i++)
            {
                expanded[i] = (i < capacity) ? container[i] : 0;
            }
            delete[] container;
            capacity = newSize;
            container = expanded;
        }
        if (index >= current)
        {
            current = index + 1;
        }
        return container[index];
    }
    
    const T& operator[](int index) const{
        if (index < 0 || index >= current) {
            throw std::out_of_range("Index out of range");
        }
        return container[index];
    }
    int size()
    {
        return current;
    }

    int size() const{
        return current;
    }
    int maxsize()
    {
        return capacity;
    }


//Using a variadic template to push back multiple elements into the vector
    template <typename... Args>
    void push_back_many(Args &&...args)
    {
        (this->push_back(std::forward<Args>(args)), ...);  //using forward because the elements may be lvalue or rvalue references and we want to handle both
    }

    template<int sz>
    T accumulate() const
    { //to calculate sum of all elements of an array we use a fold expression
        return accumulate_helper(std::make_index_sequence<sz>{}); 

    }

    template <size_t... Is>
    T accumulate_helper(std::index_sequence<Is...>) const //fold expression being used, index sequence to convert the elements of the array to indiviudal elements
    { 
        return (container[Is] + ...);
    }

    

    //For some reason the destructor was buggy hence created a custom function to de-allocate memory from a vector
    void del()
    {
        

        if constexpr (std::is_class_v<T>)
        {
            for (int i = 0; i < current; ++i)
            {
                container[i].~T();
            }
        }
        this->capacity=0;
        delete[] container; 
    }

};
