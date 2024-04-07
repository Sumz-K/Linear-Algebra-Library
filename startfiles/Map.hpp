#include <iostream>
#include <stdexcept> 


using namespace std;

template<typename Key, typename Value>
struct KVPair {
    Key key;
    Value value;
    KVPair* next;
};

template<typename Key, typename Value>
class Unordered_Map {
private:
    int tablesize = 100;
    KVPair<Key, Value>** hashtable;

    int hashfunc(const Key& key) {
        int hashvalue = 0;
        string key_string = std::to_string(key);
        for (char i : key_string) {
            hashvalue = hashvalue * 13 + static_cast<int>(i);
        }
        return hashvalue % tablesize;
    }

public:
    Unordered_Map() {
        hashtable = new KVPair<Key, Value>*[tablesize];
        for (int i = 0; i < tablesize; i++) {
            hashtable[i] = nullptr;
        }
    }

    Unordered_Map(int size) : tablesize(size) {
        hashtable = new KVPair<Key, Value>*[tablesize];
        for (int i = 0; i < tablesize; i++) {
            hashtable[i] = nullptr;
        }
    }

    ~Unordered_Map() {
        for (int i = 0; i < tablesize; i++) {
            KVPair<Key, Value>* current = hashtable[i];
            while (current != nullptr) {
                KVPair<Key, Value>* next = current->next;
                delete current;
                current = next;
            }
        }
        delete[] hashtable;
    }

    void insert(const Key& key, const Value& value) {
        int index = hashfunc(key);
        KVPair<Key, Value>* kv = new KVPair<Key, Value>{key, value, nullptr};

        if (hashtable[index] == nullptr) {
            hashtable[index] = kv;
        } else {
            KVPair<Key, Value>* curr = hashtable[index];
            while (curr->next != nullptr) {
                curr = curr->next;
            }
            curr->next = kv;
        }
    }

    Value get(const Key& key) {
        int index = hashfunc(key);
        KVPair<Key, Value>* kv = hashtable[index];

        while (kv != nullptr) {
            if (kv->key == key) {
                return kv->value;
            }
            kv = kv->next;
        }
        throw std::out_of_range("Key not found in map");

    }

    Value operator[](const Key &key){
        try{
            return get(key);
        }
        catch (std::out_of_range){
            cout<<"Key not in map\n";
            return -1;
        }
    }
};
