#include <unordered_map>
#include <string>

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

#define DEFINE_INHERITANCE_REGISTRY(baseCls)                                              \
    baseCls::baseCls##Reg_Table_t baseCls::baseCls##Reg_ConstructorTable_

#define REGISTER_DERIVED_CLASS(baseCls, subClass, lookup)                                 \
    baseCls::baseCls##Reg_TableAdder<subClass> add##subClass##To##baseCls(lookup)
