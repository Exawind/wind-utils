//  Copyright 2016 National Renewable Energy Laboratory
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//

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
    static baseCls##Reg_Table_t* baseCls##Reg_ConstructorTable_;                          \
                                                                                          \
    static void baseCls##Reg_TableInit() ;                                                \
    static void baseCls##Reg_TableDestory();                                              \
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
            baseCls##Reg_TableInit();                                                     \
            (*baseCls##Reg_ConstructorTable_)[lookup] = create;                           \
        }                                                                                 \
                                                                                          \
        ~baseCls##Reg_TableAdder()                                                        \
        {                                                                                 \
            baseCls##Reg_ConstructorTable_->erase(lookup_);                               \
        }                                                                                 \
                                                                                          \
    private:                                                                              \
        std::string lookup_;                                                              \
    }

//! Initialize the static members of class registration infrastructure
#define DEFINE_INHERITANCE_REGISTRY(baseCls)                                              \
    void baseCls::baseCls##Reg_TableInit()                                                \
    {                                                                                     \
        static bool constructed=false;                                                    \
        if (!constructed) {                                                               \
            baseCls::baseCls##Reg_ConstructorTable_ = new baseCls::baseCls##Reg_Table_t;  \
            constructed = true;                                                           \
        }                                                                                 \
    }                                                                                     \
    void baseCls::baseCls##Reg_TableDestory()                                             \
    {                                                                                     \
        if (baseCls::baseCls##Reg_ConstructorTable_) {                                    \
            delete baseCls::baseCls##Reg_ConstructorTable_;                               \
            baseCls::baseCls##Reg_ConstructorTable_ = nullptr;                            \
        }                                                                                 \
    }                                                                                     \
baseCls::baseCls##Reg_Table_t* baseCls::baseCls##Reg_ConstructorTable_ = nullptr

//! Convenience wrapper to register subclass to base class
#define REGISTER_DERIVED_CLASS(baseCls, subClass, lookup)                                 \
    baseCls::baseCls##Reg_TableAdder<subClass> add##subClass##To##baseCls(lookup)
