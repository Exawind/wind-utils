
/** \file
 *
 * Runtime selection class registration utilities
 *
 * This file defines a few convenience macros for generating the boilerplate
 * required for a class registration mechanisms. The design is inspired by
 * OpenFOAM's runTimeSelection infrastructure but updated to modern C++11
 * standard.
 */

#include <unordered_map>
#include <string>

/** Declare all the static data members necessary in the base class for derived
 * class registration.
 *
 * \param baseCls The name of the Base Class
 * \param argList Argument list for the constructor of derived classes
 * \param parList Parameter list used during invocation of the constructor
 */
#define DECLARE_INHERITANCE_REGISTRY(baseCls, argList, parList)                           \
                                                                                          \
    using baseCls##Reg_create_t = baseCls* (*)argList;                                    \
    using baseCls##Reg_Table_t = std::unordered_map<std::string, baseCls##Reg_create_t>;  \
                                                                                          \
    static baseCls##Reg_Table_t baseCls##Reg_ConstructorTable_;                           \
                                                                                          \
    template<typename T>                                                                  \
    class baseCls##Reg_TableAdder                                                         \
    {                                                                                     \
    public:                                                                               \
        static baseCls* create argList                                                    \
        {                                                                                 \
            return new T parList;                                                         \
        }                                                                                 \
                                                                                          \
        baseCls##Reg_TableAdder(std::string lookup)                                       \
            : lookup_(lookup)                                                             \
        {                                                                                 \
            baseCls##Reg_ConstructorTable_[lookup] = create;                              \
        }                                                                                 \
                                                                                          \
        ~baseCls##Reg_TableAdder()                                                        \
        {                                                                                 \
            baseCls##Reg_ConstructorTable_.erase(lookup_);                                \
        }                                                                                 \
                                                                                          \
    private:                                                                              \
        std::string lookup_;                                                              \
    };

//! Initialize the static members of class registration infrastructure
#define DEFINE_INHERITANCE_REGISTRY(baseCls)                                              \
    baseCls::baseCls##Reg_Table_t baseCls::baseCls##Reg_ConstructorTable_

//! Convenience wrapper to register subclass to base class
#define REGISTER_DERIVED_CLASS(baseCls, subClass, lookup)                                 \
    baseCls::baseCls##Reg_TableAdder<subClass> add##subClass##To##baseCls(lookup)
