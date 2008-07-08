/* 
 * File:   objectfactory.h
 * Author: ruehle
 *
 * Created on July 12, 2007, 11:39 AM
 */

#ifndef _objectfactory_H
#define	_objectfactory_H

#include <map>
#include <list>

using namespace std;

/**
    The REGISTER_OBJECT macro allows to easily register an object in an object factory.
 */
#define REGISTER_OBJECT(factory, object, key) \
    class _register_##object { \
        public: \
        _register_##object() { factory().Register(key, new object(), false); }  \
    }; \
    _register_##object _instance_register_##object;

/**
    \brief template class for object factory

    This class is a template for an object factory. The factory creates an instance of an derived class,
    be giving a key (e.g. a string) for the object which identifies it uniquely. This allows the implementation
    of new features (e.g. new file formats, new mapping algorithms) without touching or recompiling existing bits of code.

*/
template<typename key_t, typename T>
class ObjectFactory
{
public:
    ObjectFactory() {}
    ~ObjectFactory();
    
    /**
        \brief register an object
    
        This function is called to register an object in the factory. After an object is registered,
        an instance of it can be created by calling Create specifying the corresponding key. The object
        must implement the clone command, which creates an instance of the object.
     */
    void Register(const key_t &key, T *obj, bool bSynonym = true);
    /**
       Create an instance of the object identified by key.
    */
    T *Create(const key_t &key);
    
    /**
        get pointer to an object which is stored in the ObjectFactory to e.g. modify creation parameters.
     */
    T *get(const key_t &key);
private:
    map<key_t, T*> _objects;
    list<T*> _delete_these;
};


template<typename key_t, typename T>
ObjectFactory<key_t, T>::~ObjectFactory()
{
    typename list<T*>::iterator iter;
    for(iter=_delete_these.begin();iter!=_delete_these.end(); ++iter)
        delete (*iter);
    _objects.clear();
    _delete_these.clear();
}

template<typename key_t, typename T>
void ObjectFactory<key_t, T>::Register(const key_t &key, T *obj, bool bSynonym)
{
//   cout << "registered: " << key << endl;
       _objects[key] = obj;
    if(!bSynonym) {
        obj->RegisteredAt(*this);
        _delete_these.push_back(obj);
    }
}

template<typename key_t, typename T>
T* ObjectFactory<key_t, T>::Create(const key_t &key)
{
    T *obj;
    obj = get(key);
    if(obj==NULL) return NULL;
    return obj->Clone();
}

template<typename key_t, typename T>
T* ObjectFactory<key_t, T>::get(const key_t &key)
{
    typename map<key_t, T*>::iterator iter;
    iter = _objects.find(key);
    if(iter == _objects.end())
        return NULL;
    return (*iter).second;
}

#endif	/* _objectfactory_H */
