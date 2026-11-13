// Smart Pointer
// par Christophe [Groove] Riccio
// www.g-truc.net
// Modification ULR VALOR accès d'un pointeur tableau

#ifndef __SMART_PTR_H__
#define __SMART_PTR_H__

/**
 *
 *  \brief Classe de gestion de pointeur
 *
 *	Cette classe permet de gérer les pointeurs
 *	L'utilité est de ne pas avoir besoin de supprimer l'objet à la fin de son utilisation.
 *  La fonction est donc de ne plus avoir de fuites mémoire.
 *
 */


template <typename T>
class smart_ptr
{
public:
	smart_ptr();
	smart_ptr(const smart_ptr& SmartPtr);
	smart_ptr(T* Pointer);
	~smart_ptr();

    smart_ptr& operator=(const smart_ptr& SmartPtr);
    smart_ptr& operator=(T* Pointer);

    T& operator*();
    T* operator->();
	T* get() const; //Ajout 27/07/2009 pour accéder à l'élément.
    const T& operator*() const;
    const T* operator->() const;

	bool operator==(const smart_ptr& SmartPtr) const;
	bool operator!=(const smart_ptr& SmartPtr) const;

private:
	void _Release();

    int* _RefCounter;
    T* _Pointer;
};


template <typename T>
class smart_ptr_ar
{
public:
	smart_ptr_ar();
	smart_ptr_ar(const smart_ptr_ar& SmartPtr);
	smart_ptr_ar(T* Pointer);
	~smart_ptr_ar();

    smart_ptr_ar& operator=(const smart_ptr_ar& SmartPtr);
    smart_ptr_ar& operator=(T* Pointer);

	T& operator[](unsigned i); //Ajout 13/07/2007 pour gérer l'utilisation du array à base de smart_ptr

    T& operator*();
    T* operator->();
	T* get() const; //Ajout 27/07/2009 pour accéder à l'élément.
    const T& operator*() const;
    const T* operator->() const;

	bool operator==(const smart_ptr_ar& SmartPtr) const;
	bool operator!=(const smart_ptr_ar& SmartPtr) const;

private:
	void _Release();

    int* _RefCounter;
    T* _Pointer;
};
#include "smart_ptr.inl"

#endif //__SMART_PTR_H__
