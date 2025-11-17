/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

       







namespace dco {

typedef DCO_TAPE_INT_DEFINE DCO_TAPE_INT;
inline DCO_TAPE_INT& stce_9752( DCO_TAPE_INT& stce_9753) { return stce_9753; }
inline const DCO_TAPE_INT& stce_9752(const DCO_TAPE_INT& stce_9753) { return stce_9753; }

typedef DCO_TAPE_INT_DEFINE DCO_INTEGRAL_TAPE_INT;





static const bool stce_9754 = true;





static const bool stce_9755 = true;







static const bool stce_9756 = true;





static const bool stce_9757 = false;



static const bool stce_9758 = true;




template<class stce_9759,
    class stce_9760 = stce_9759,
    class stce_9761 = stce_9759>
struct stce_9762
{
  typedef stce_9759 value_t;
  typedef stce_9760 stce_9763;
  typedef stce_9761 stce_9764;
};

}
       
       
namespace dco
{
typedef unsigned long mem_long_t;
inline size_t get_allocation_size(size_t stce_9765);
struct memory_model {
  enum TYPE
  {
    stce_9766,
    BLOB_TAPE,
    BLOB_TAPE_SPLINT,
    CHUNK_TAPE,
    stce_9767,
    stce_9768,
    stce_9769,
    stce_9770
  };
};
namespace internal
{
template<class stce_9771, class stce_9772>
  struct active_type;
template<class stce_9771, class DATA_TYPE>
  struct stce_9773;
template<class stce_9771, class stce_9774, class stce_9775>
  struct stce_9776;
template<class stce_9771, class stce_9777, class stce_9778, class stce_9775>
  struct stce_9779;
template<class stce_9771, class stce_9777, class stce_9775>
  struct stce_9780;
template<class stce_9771, class stce_9778, class stce_9775>
  struct stce_9781;
template <typename stce_9782, typename enable_if = void>
  struct trait_value;
template <typename stce_9782, typename stce_9783 = void, typename enable_if = void>
  struct stce_9784;
template <typename stce_9782, typename enable_if = void>
  struct stce_9785;
template<typename stce_9782, typename enable_if = void>
  struct stce_9786;
template<typename stce_9782, typename enable_if = void>
  struct stce_9787;
template<enum memory_model::TYPE MEMORY_MODEL,
         class stce_9788,
         enum memory_model::TYPE stce_9789=memory_model::stce_9766>
  class tape;
template<typename stce_9782>
  struct stce_9790 {};
template<memory_model::TYPE MEMORY_MODEL,
         typename stce_9788,
         memory_model::TYPE stce_9789>
  struct stce_9790<tape<MEMORY_MODEL, stce_9788, stce_9789> > {
  typedef typename stce_9788::stce_9763 type;
};
template<class stce_9791, typename stce_9792>
  class stce_9793;
template<class stce_9794>
  struct jacobian_preaccumulator_t;
}
template <class stce_9795> struct mode;
enum stce_9796 {
  stce_9797, stce_9798, stce_9799, stce_9800, stce_9801
};
template <typename stce_9782>
  typename internal::stce_9786<stce_9782>::stce_9802 tape(const stce_9782&);
  template<typename stce_9771>
  struct stce_9803 {
    typedef stce_9771 active_value_t;
    typedef stce_9771 value_t;
    typedef stce_9771 passive_t;
    static const int order = 0;
  };
namespace vector_ns {
template<typename stce_9804, size_t stce_9805> struct alignas( alignof(stce_9804[stce_9805]) ) vector_t;
}
}
extern "C" {
  extern int stce_1000002(int stce_9806=false);
  extern void* stce_1000001(size_t *stce_9807, int);
  extern void* stce_999999(size_t, int);
  extern int stce_1000003(int *stce_9808);
  extern long stce_1000004(long *stce_9808);
  extern void stce_1000000(void* stce_9809);
}
template <typename stce_9782> inline stce_9782 stce_9810(stce_9782 &stce_9808);
template<> inline int stce_9810<int>(int &stce_9808) {
  return stce_1000003(&stce_9808);
}
template<>
  inline long stce_9810<long>(long &stce_9808) {
  return stce_1000004(&stce_9808);
}
namespace dco {
template<memory_model::TYPE MEMORY_MODEL, class stce_9788, memory_model::TYPE stce_9789>
  bool is_null(internal::tape<MEMORY_MODEL, stce_9788, stce_9789> *tape) {
  return (tape == NULL);
}
 static unsigned int stce_9811=0xABCDE;
 static unsigned int stce_9812(unsigned int *stce_9813, unsigned int stce_9814) {
    int stce_9815;
    stce_9815 = *stce_9813 & 1;
    *stce_9813 >>=1;
    if(stce_9815==1)
        *stce_9813=stce_9814;
    return *stce_9813;
 }
 static inline double stce_9816() {
      stce_9812(&stce_9811, 0xB4BCD35C);
      unsigned int stce_9817 = stce_9812(&stce_9811, 0xB4BCD35C);
      return double(stce_9817)/RAND_MAX;
 }
  template <class stce_9782>
  std::string
  type_name()
  {
    typedef typename std::remove_reference<stce_9782>::type stce_9818;
    std::unique_ptr<char, void(*)(void*)> stce_9819
      (
       abi::__cxa_demangle(typeid(stce_9818).name(), nullptr,
               nullptr, nullptr),
       std::free
       );
    std::string stce_9820 = stce_9819 != nullptr ? stce_9819.get() : typeid(stce_9818).name();
    if (std::is_const<stce_9818>::value)
      stce_9820 += " const";
    if (std::is_volatile<stce_9818>::value)
      stce_9820 += " volatile";
    if (std::is_lvalue_reference<stce_9782>::value)
      stce_9820 += "&";
    else if (std::is_rvalue_reference<stce_9782>::value)
      stce_9820 += "&&";
    return stce_9820;
  }
namespace internal {
template <typename stce_9782, bool>
struct stce_9821 { typedef stce_9782 type; };
template <typename stce_9782>
  struct stce_9821<stce_9782, true> { typedef stce_9782& type; };
}
namespace helper {
static bool stce_9822 = true;
template<typename stce_9782>
  struct stce_9823 {
  static const bool stce_9753 = true;
};
template<class stce_9782, class stce_9824>
struct is_same {
  static const bool value = false;
};
template<class stce_9782>
struct is_same<stce_9782, stce_9782> {
  static const bool value = true;
};
template <typename stce_9782>
  struct type_identity {
  typedef stce_9782 type;
};
template<class stce_9825, const int stce_9826=1>
  struct stce_9827 {
  static const int stce_9828 = stce_9826;
  typedef stce_9825 value_t;
  stce_9825 stce_9829[stce_9826];
  void stce_9830(const stce_9825 &stce_9831) {
    for (DCO_TAPE_INT stce_9832 = 0; stce_9832 < stce_9826; ++stce_9832)
      stce_9829[stce_9832] = stce_9831;
  }
  stce_9827(const stce_9825 &stce_9831) {
    stce_9830(stce_9831);
  }
  stce_9827() {
    stce_9830(0);
  }
  stce_9825 &operator[](const DCO_TAPE_INT stce_9808) {
    return stce_9829[stce_9808];
  }
  const stce_9825 &operator[](const DCO_TAPE_INT stce_9808) const {
    return stce_9829[stce_9808];
  }
  stce_9827 &operator= (const stce_9825 &stce_9833) {
    for (int stce_9832 = 0; stce_9832 < stce_9826; ++stce_9832)
      this->stce_9829[stce_9832] = stce_9833;
    return *this;
  }
  stce_9827 &operator += (const stce_9827 &stce_9833) {
    for (int stce_9832 = 0; stce_9832 < stce_9826; ++stce_9832)
      this->stce_9829[stce_9832] += stce_9833.stce_9829[stce_9832];
    return *this;
  }
};
template <typename stce_9825, const int stce_9826>
  bool operator==(stce_9827<stce_9825, stce_9826> const& stce_9834, int const& stce_9820) {
  for (int stce_9832 = 0; stce_9832 < stce_9826; ++stce_9832) {
    if (stce_9834[stce_9832] != stce_9820) {
      return false;
    }
  }
  return true;
}
template <typename stce_9825, typename stce_9835, const int stce_9826>
  bool operator==(stce_9827<stce_9825, stce_9826> const& stce_9834, stce_9827<stce_9835, stce_9826> const& stce_9820) {
  for (int stce_9832 = 0; stce_9832 < stce_9826; ++stce_9832) {
    if (stce_9834[stce_9832] != stce_9820[stce_9832]) {
      return false;
    }
  }
  return true;
}
template<class stce_9782>
struct stce_9836 {
  static const int stce_9837;
  typedef stce_9782 stce_9838;
};
template<class stce_9782>
const int stce_9836<stce_9782>::stce_9837 = 1;
template<class stce_9825, int stce_9826>
struct stce_9836<stce_9827<stce_9825, stce_9826> > {
  static const int stce_9837;
  typedef stce_9825 stce_9838;
};
template<class stce_9825, int stce_9826>
const int stce_9836<stce_9827<stce_9825, stce_9826> > ::stce_9837 = stce_9826;
template<class stce_9825, int stce_9826>
  static inline
  stce_9827<stce_9825, stce_9826>
  operator + (const stce_9827<stce_9825, stce_9826> &stce_9839, const stce_9827<stce_9825, stce_9826> &stce_9833) {
  stce_9827<stce_9825, stce_9826> stce_9840;
  for (int stce_9832 = 0; stce_9832 < stce_9826; ++stce_9832)
    stce_9840.stce_9829[stce_9832] = stce_9839.stce_9829[stce_9832] + stce_9833.stce_9829[stce_9832];
  return stce_9840;
}
template<class stce_9825, int stce_9826>
static inline stce_9827<stce_9825, stce_9826> operator *(const stce_9827<stce_9825, stce_9826> &stce_9839,
        const stce_9827<stce_9825, stce_9826> &stce_9833)
{
    stce_9827<stce_9825, stce_9826> stce_9840;
    for (DCO_TAPE_INT stce_9832 = 0; stce_9832 < stce_9826; ++stce_9832)
        stce_9840.stce_9829[stce_9832] = stce_9839.stce_9829[stce_9832] * stce_9833.stce_9829[stce_9832];
    return stce_9840;
}
template<class stce_9825, class stce_9835, int stce_9826>
static inline stce_9827<stce_9825, stce_9826> operator*(const stce_9835 &stce_9841,
        const stce_9827<stce_9825, stce_9826> &stce_9842)
{
    stce_9827<stce_9825, stce_9826> stce_9840;
    for (DCO_TAPE_INT stce_9832 = 0; stce_9832 < stce_9826; ++stce_9832)
        stce_9840.stce_9829[stce_9832] = stce_9842.stce_9829[stce_9832] * stce_9841;
    return stce_9840;
}
template<class stce_9825, class stce_9835, int stce_9826>
static inline stce_9827<stce_9825, stce_9826> operator *(const stce_9827<stce_9825, stce_9826> &stce_9839,
        const stce_9835 &stce_9833)
{
    stce_9827<stce_9825, stce_9826> stce_9840;
    for (DCO_TAPE_INT stce_9832 = 0; stce_9832 < stce_9826; ++stce_9832)
        stce_9840.stce_9829[stce_9832] = stce_9839.stce_9829[stce_9832] * stce_9833;
    return stce_9840;
}
template<class stce_9825, int stce_9826>
  static inline std::ostream &operator << (std::ostream &out, const stce_9827<stce_9825, stce_9826> &stce_9843) {
  out << "(";
  for (DCO_TAPE_INT stce_9832 = 0; stce_9832 < stce_9826 - 1; ++stce_9832) {
    out << stce_9843.stce_9829[stce_9832] << ";";
  }
  out << stce_9843.stce_9829[stce_9826 - 1];
  out << ")";
  return out;
}
}
template<typename TYPE>
  class stce_9844 {
  const TYPE &stce_9782;
public:
  stce_9844(const TYPE &stce_9845) : stce_9782(stce_9845) {}
  template <typename stce_9846>
    stce_9846 *create() const {
    return new stce_9846(stce_9782);
  }
};
template<>
  class stce_9844<void *> {
public:
  stce_9844(void *stce_9847) {
    (void) stce_9847;
  }
  template <typename stce_9846>
    stce_9846 *create() const {
    return new stce_9846();
  }
};
class stce_9848 {
  std::ofstream stce_9849;
  std::string stce_9850;
  std::string stce_9851;
  std::stringstream stce_9852;
  std::stringstream stce_9853;
public:
  stce_9848(const std::string &filename) : stce_9850("digraph {\n"), stce_9851("\n}\n") {
    stce_9852.clear();
    stce_9853.clear();
    stce_9849.open(filename.c_str());
  }
  void stce_9854(std::vector<DCO_TAPE_INT> &stce_9855, std::vector<DCO_TAPE_INT> &stce_9856, DCO_TAPE_INT stce_9857) {
    stce_9852 << stce_9857 << "[shape=box,label=\"(" << stce_9857 << ", ext. function)\",color=\"red\"]\n";
    for (std::vector<DCO_TAPE_INT>::iterator stce_9832 = stce_9855.begin(); stce_9832 != stce_9855.end(); ++stce_9832) {
      stce_9853 << *stce_9832 << " -> " << stce_9857 << "[label=\"unknown" << "\",color=\"red\"]\n";
    }
    for (std::vector<DCO_TAPE_INT>::iterator stce_9832 = stce_9856.begin(); stce_9832 != stce_9856.end(); ++stce_9832) {
      stce_9853 << stce_9857 << " -> " << *stce_9832 << "[label=\"unknown" << "\",color=\"red\"]\n";
    }
  }
  void stce_9858(DCO_TAPE_INT stce_9857) {
    if(stce_9857 == 0){
      return;
    }
    stce_9852 << stce_9857 << "[shape=box,label=\"" << stce_9857 << "\", color=\"green\"]\n";
  }
  template <typename stce_9859>
  void stce_9858(DCO_TAPE_INT stce_9857, stce_9859 stce_9860) {
    if(stce_9857 == 0){
      return;
    }
    stce_9852 << stce_9857 << "[shape=box,label=\"(" << stce_9857 << ", " << stce_9860 << ")\", color=\"green\"]\n";
  }
  template <typename stce_9859>
  void stce_9861(DCO_TAPE_INT stce_9862, DCO_TAPE_INT stce_9863, stce_9859 stce_9864) {
    if(stce_9862 == 0 || stce_9863 == 0){
      return;
    }
    stce_9853 << stce_9862 << " -> " << stce_9863 << "[label=\"" << stce_9864 << "\"]\n";
  }
  ~stce_9848() {
    stce_9849 << stce_9850;
    stce_9849 << stce_9852.str();
    stce_9849 << stce_9853.str();
    stce_9849 << stce_9851;
    stce_9849.close();
  }
};
template<bool stce_9865, class stce_9782, class stce_9866> struct stce_9867 { typedef stce_9782 type; };
template< class stce_9782, class stce_9866> struct stce_9867<false, stce_9782, stce_9866> { typedef stce_9866 type; };
template<bool stce_9865, typename stce_9782 = void> struct enable_if {};
template< typename stce_9782 > struct enable_if<true, stce_9782> { typedef stce_9782 type; };
template <typename stce_9782> struct type_identity { typedef stce_9782 type; };
template<class stce_9782> struct remove_const { typedef stce_9782 type; };
template<class stce_9782> struct remove_const<const stce_9782> { typedef stce_9782 type; };
template<class stce_9782> struct remove_const<const stce_9782&> { typedef stce_9782& type; };
template<class stce_9782> struct remove_const<const stce_9782&&> { typedef stce_9782&& type; };
template<class stce_9782> struct stce_9868 { typedef stce_9782 type; };
template<class stce_9782> struct stce_9868<stce_9782&> { typedef stce_9782 type; };
template<typename stce_9782> struct remove_all { typedef stce_9782 type; };
template<typename stce_9782> struct remove_all<const stce_9782> { typedef typename remove_all<stce_9782>::type type; };
template<typename stce_9782> struct remove_all<stce_9782 const&> { typedef typename remove_all<stce_9782>::type type; };
template<typename stce_9782> struct remove_all<stce_9782&> { typedef typename remove_all<stce_9782>::type type; };
template<typename stce_9782> struct remove_all<stce_9782 const*> { typedef typename remove_all<stce_9782>::type type; };
template<typename stce_9782> struct remove_all<stce_9782*> { typedef typename remove_all<stce_9782>::type type; };
template<typename stce_9782> struct stce_9869 { enum { value = false }; };
template<typename stce_9782> struct stce_9869<stce_9782*> { enum { value = true }; };
template <typename stce_9870, typename stce_9871>
  struct stce_9872 {
  static const bool value = (stce_9870::order < stce_9871::order);
};
template <>
  struct stce_9872<void, void> {
  static const bool value = false;
};
static inline void *alloc(size_t size, bool stce_9873 = false) {
  return stce_999999(size, static_cast<int>(stce_9873));
}
static inline void dealloc(void* stce_9809) {
  stce_1000000(stce_9809);
}
}
       
namespace dco {
  namespace internal {
    template<class stce_9771, class stce_9874 = stce_9771>
    struct stce_9875 {
      DCO_TAPE_INT stce_9876;
      DCO_TAPE_INT stce_9877;
      DCO_TAPE_INT _progvarcounter;
      std::string stce_9878;
      std::ofstream stce_9879;
      std::ofstream stce_9880;
    private:
      stce_9875(const std::string stce_9881) : stce_9876(0), stce_9877(0), _progvarcounter(0), stce_9878(stce_9881) {
      }
    public:
      static stce_9875 *create(const std::string stce_9881, const std::string filename) {
        const char *stce_9850 = "//=================================================//\n"
                             "// This file is generated by dco/c++.              //\n"
                             "// Further inquiries info@stce.rwth-aachen.de      //\n"
                             "//=================================================//\n\n"
                             "#include <cmath>\n\n"
                             "typedef double DCO_TAPE_REAL;\n\n";
        stce_9875 *stce_9840 = new stce_9875(stce_9881);
        stce_9840->stce_9879.open(filename.c_str());
        stce_9840->stce_9879 << stce_9850 << "void " << stce_9881 << "(const double *x, double *x_a1s, double *y, const double *y_a1s) {\n double va0_a1s=0; \n";
        stce_9840->stce_9880.open("ad_tmp");
        return stce_9840;
      }
      static void remove(stce_9875 *tape) {
        tape->stce_9880.close();
        std::ifstream stce_9882("ad_tmp");
        std::vector<std::string> stce_9883;
        std::string stce_9884;
        while (std::getline(stce_9882, stce_9884)) {
          stce_9883.insert(stce_9883.begin(), stce_9884);
        }
        tape->stce_9879 << "// BEGIN OF ADJOINT SECTION\n";
        for (size_t stce_9832 = 0; stce_9832 < stce_9883.size(); ++stce_9832) {
          tape->stce_9879 << stce_9883[stce_9832] << "\n";
        }
        tape->stce_9879 << "}\n";
        tape->stce_9879.close();
      }
      template<class stce_9885>
      void register_variable(stce_9885 &stce_9886) {
        _progvarcounter++;
        stce_9886._tape_index() = _progvarcounter;
        stce_9879 << "const double &va" << _progvarcounter << "=x[" << stce_9876 << "];\n";
        stce_9879 << "double &va" << _progvarcounter << "_a1s=x_a1s[" << stce_9876 << "];\n";
        stce_9876++;
      }
      template<class stce_9885>
      void register_output_variable(stce_9885 &stce_9886) {
        stce_9879 << "y[" << stce_9877 << "]=va" << stce_9886._tape_index() << ";\n";
        stce_9879 << "va" << stce_9886._tape_index() << "_a1s+=y_a1s[" << stce_9877 << "];\n";
        stce_9877++;
      }
    };
    class stce_9887
    {
      template<typename>
      friend class stce_9888;
      DCO_TAPE_INT stce_9889;
      stce_9887(): stce_9889(0) {}
    public:
      std::string stce_9890() const {
        std::stringstream stce_9809;
        stce_9809 << "va" << stce_9889;
        return stce_9809.str();
      }
      DCO_TAPE_INT& _tape_index() {
        return stce_9889;
      }
      DCO_TAPE_INT _tape_index() const {
        return stce_9889;
      }
      void clear() {
        _tape_index() = 0;
      }
    };
    template<class stce_9891>
    class stce_9888: public stce_9887{
    public:
      typedef stce_9891 mode_t;
      stce_9888(){}
      template<class stce_9892>
      DCO_ALWAYS_INLINE void stce_9893(const stce_9892 &stce_9894) {
        stce_9894.stce_9895(stce_9891::global_tape->stce_9879);
        _tape_index() = stce_9891::global_tape->_progvarcounter;
        std::stringstream stce_9896;
        stce_9896 << "va" << ++stce_9891::global_tape->_progvarcounter;
       stce_9891::global_tape->stce_9880 << stce_9894.stce_9890() << "_a1s+=" <<
            stce_9896.str() << "_a1s;\t";
        stce_9894.stce_9897(stce_9891::global_tape->stce_9880);
        stce_9891::global_tape->stce_9880 << "//EOL\n";
        _tape_index()++;
        stce_9891::global_tape->stce_9879 << "double " << stce_9890()
            << "=" << stce_9894.stce_9890() << ";\n";
        stce_9891::global_tape->stce_9879 << "double "
            << stce_9890() << "_a1s=0.0;\n";
      }
    };
  }
}
namespace dco {
namespace internal {
static inline bool is_character(const std::string &stce_9898, size_t stce_9808) {
  if(stce_9808 >=stce_9898.size()) return false;
  return std::isalpha(stce_9898[stce_9808]) != 0 ;
}
static inline std::string replace_all_in_string(std::string& str, const std::string &stce_9862, const std::string &stce_9863) {
  size_t stce_9899=0;
  while ((stce_9899 = str.find(stce_9862, stce_9899)) != std::string::npos) {
    if((stce_9899 > 0 && is_character(str, stce_9899-1))
      || is_character(str, stce_9899 + stce_9862.length())) {
        stce_9899++;
        continue;
    }
    str.replace(stce_9899, stce_9862.length(), stce_9863);
    stce_9899 += stce_9863.length();
  }
  return str;
}
template<class stce_9771>
class stce_9900 {
  static DCO_TAPE_INT stce_9901;
  DCO_TAPE_INT stce_9902;
public:
  stce_9900() :
      stce_9902(stce_9901++) {
  }
   std::string stce_9890() const {
    std::stringstream stce_9809;
    stce_9809 << "im" << stce_9902;
    return stce_9809.str();
  }
protected:
  void stce_9903(std::ofstream &out, std::string stce_9904) const {
    out << "const double " << stce_9890() << "=" << stce_9904 << ";\n";
    out << "double " << stce_9890() << "_a1s=0.0;\n";
  }
  template<class stce_9782>
  std::string str(stce_9782 arg) const {
    return arg.stce_9890();
  }
  std::string str(stce_9771 stce_9905) const {
    std::stringstream str;
    str << std::scientific << std::setprecision(17) << stce_9905;
    return str.str();
  }
  template <class stce_9906>
  void stce_9907(std::string& stce_9904, const stce_9906& arg, const std::string stce_9908) const {
    replace_all_in_string(stce_9904, stce_9908, str(arg));
  }
  void stce_9907(std::string& stce_9904, const stce_9771& arg, const std::string stce_9908) const {
    replace_all_in_string(stce_9904, stce_9908, str(arg));
  }
  template <class stce_9775, class stce_9909, class stce_9910>
  void stce_9911(std::ofstream &out, const stce_9909& stce_9912, const stce_9910& stce_9913) const {
    std::string stce_9904 = stce_9775::stce_9914();
    this->stce_9907(stce_9904, stce_9912, "arg1");
    this->stce_9907(stce_9904, stce_9913, "arg2");
    out << this->str(stce_9912) << "_a1s +=(" << stce_9904 << ")*" << this->stce_9890() << "_a1s;\t";
    stce_9912.stce_9897(out);
  }
  template <class stce_9775, class stce_9910>
  void stce_9911(std::ofstream&, const stce_9771&, const stce_9910&) const {}
  template <class stce_9775, class stce_9909, class stce_9910>
  void stce_9915(std::ofstream &out, const stce_9909& stce_9912, const stce_9910& stce_9913) const {
    std::string stce_9904 = stce_9775::stce_9916();
    this->stce_9907(stce_9904, stce_9912, "arg1");
    this->stce_9907(stce_9904, stce_9913, "arg2");
    out << this->str(stce_9913) << "_a1s +=(" << stce_9904 << ")*" << this->stce_9890() << "_a1s;\t";
    stce_9913.stce_9897(out);
  }
  template <class stce_9775, class stce_9909>
  void stce_9915(std::ofstream&, const stce_9909&, const stce_9771&) const {}
};
template<class stce_9771>
DCO_TAPE_INT stce_9900<stce_9771>::stce_9901 = 0;
}
}
       
       
       
namespace dco {
  class exception {
  public:
    template <typename stce_9917>
    static stce_9917 create(std::string stce_9918, std::string stce_9919 = "", int stce_9884 = 0) {
      std::stringstream stce_9920;
      stce_9920 << "--- dco/c++ --- " << stce_9918;
      if (stce_9919 != "") stce_9920 << " --- " << stce_9919 << ":" << stce_9884 << ".";
      DCO_LOG(dco::logERROR) << "EXCEPTION thrown: " << stce_9920.str();
      return stce_9917(stce_9920.str());
    }
    template <typename stce_9917>
    static stce_9917 stce_9921(std::string stce_9918, std::string stce_9919 = "", int stce_9884 = 0) {
      std::stringstream stce_9920;
      stce_9920 << "--- dco/c++ --- " << stce_9918;
      if (stce_9919 != "") stce_9920 << " --- " << stce_9919 << ":" << stce_9884 << ".";
      DCO_LOG(dco::logERROR) << "EXCEPTION thrown: " << stce_9920.str();
      return stce_9917();
    }
};
}
       
namespace dco {
namespace folding {
template <typename stce_9782, typename enable_if = void>
  struct is_zero_trait {
  static bool get(const stce_9782& stce_9843) {
    return stce_9843 == 0;
  }
};
template <typename stce_9782>
  bool is_zero(const stce_9782& stce_9843) {
  return is_zero_trait<stce_9782>::get(stce_9843);
}
template <typename stce_9782>
  struct is_zero_trait<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_tangent_type>::type> {
  static bool get(const stce_9782& stce_9843) {
    return is_zero(value(stce_9843)) && is_zero(derivative(stce_9843));
  }
};
template <typename stce_9782>
  struct is_zero_trait<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_adjoint_type>::type> {
  static bool get(const stce_9782& stce_9843) {
    return is_zero(value(stce_9843)) && !tape_index(stce_9843);
  }
};
}
}
namespace dco {
namespace folding {
template<class stce_9922>
  class stce_9923 {
  typedef typename stce_9922::stce_9764 stce_9764;
  typedef typename stce_9922::stce_9763 stce_9763;
  typedef typename stce_9922::stce_9924 stce_9925;
public:
  stce_9923(stce_9922* tape) : stce_9926(1), stce_9927(tape->stce_9928()) {}
  template<class stce_9929>
    DCO_ALWAYS_INLINE stce_9929& stce_9930(stce_9763 const& pval) {
    stce_9926 *= pval;
    return static_cast<stce_9929&>(*this);
  }
  DCO_ALWAYS_INLINE stce_9763& pval() { return stce_9926; }
  DCO_ALWAYS_INLINE DCO_TAPE_INT index() const { return stce_9927.index(); }
  template <class stce_9931>
    DCO_ALWAYS_INLINE void insert(stce_9931 const& stce_9843) {
    const DCO_TAPE_INT stce_9932 = stce_9843._tape_index();
    if (!stce_9755 || stce_9932 != 0) {
      if (!stce_9756 || !is_zero(stce_9926)) {
        stce_9927.insert(stce_9932, stce_9926);
      }
    }
  }
protected:
  stce_9763 stce_9926;
private:
  stce_9923();
  stce_9923(stce_9923& stce_9933);
  stce_9925 stce_9927;
};
template<class stce_9922>
  class stce_9934 : public stce_9923<stce_9922> {
  typedef typename stce_9922::stce_9764 stce_9764;
  typedef typename stce_9922::stce_9924 stce_9925;
public:
  stce_9934(stce_9922* tape) : stce_9923<stce_9922>(tape),
    stce_9935(0) {
  }
  template <class stce_9931>
    DCO_ALWAYS_INLINE void insert(stce_9931 stce_9843) {
    stce_9923<stce_9922>::insert(stce_9843);
    stce_9935 += this->stce_9926 * stce_9843.stce_9936;
  }
  stce_9764 stce_9937() const {
    return stce_9935;
  }
private:
  stce_9934();
  stce_9934(stce_9934& stce_9933);
  stce_9764 stce_9935;
};
template<class stce_9938, class stce_9759, class stce_9939>
  DCO_ALWAYS_INLINE void interpret(const dco::internal::active_type<stce_9759, stce_9939> &stce_9843,
                                   stce_9938& stce_9940) {
  stce_9940.insert(stce_9843);
}
template<class stce_9938, class stce_9759, class stce_9941, class stce_9942, class stce_9943>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9779<stce_9759, stce_9941, stce_9942, stce_9943> &stce_9843,
                                   stce_9938& stce_9940) {
  stce_9759 stce_9944 = stce_9940.pval();
  interpret<stce_9938>(stce_9843.stce_9945, stce_9940.template stce_9930<stce_9938>(stce_9843.stce_9946()));
  stce_9940.pval() = stce_9944;
  interpret<stce_9938>(stce_9843.stce_9947, stce_9940.template stce_9930<stce_9938>(stce_9843.stce_9948()));
  stce_9940.pval() = stce_9944;
}
template<class stce_9938, class stce_9759, class stce_9949, class stce_9943>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9776<stce_9759, stce_9949, stce_9943> &stce_9843,
                                   stce_9938& stce_9940) {
  stce_9759 stce_9944 = stce_9940.pval();
  interpret<stce_9938>(stce_9843.stce_9950, stce_9940.template stce_9930<stce_9938>(stce_9843.pval()));
  stce_9940.pval() = stce_9944;
}
template<class stce_9938, class stce_9759, class stce_9941, class stce_9943>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9780<stce_9759, stce_9941, stce_9943> &stce_9843,
                                   stce_9938& stce_9940) {
  stce_9759 stce_9944 = stce_9940.pval();
  interpret<stce_9938>(stce_9843.stce_9945, stce_9940.template stce_9930<stce_9938>(stce_9843.stce_9946()));
  stce_9940.pval() = stce_9944;
}
template<class stce_9938, class stce_9759, class stce_9942, class stce_9943>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9781<stce_9759, stce_9942, stce_9943> &stce_9843,
                                   stce_9938& stce_9940) {
  stce_9759 stce_9944 = stce_9940.pval();
  interpret<stce_9938>(stce_9843.stce_9947, stce_9940.template stce_9930<stce_9938>(stce_9843.stce_9948()));
  stce_9940.pval() = stce_9944;
}
template<class stce_9922>
  class stce_9951 {
public:
  stce_9951() :
    stce_9952(0) {
  }
  template <class stce_9931>
    DCO_ALWAYS_INLINE void insert(stce_9931 stce_9843) {
    stce_9952 = stce_9843.tape();
  }
  stce_9922* tape() const {
    return stce_9952;
  }
private:
  stce_9951(stce_9951& stce_9933);
  stce_9922* stce_9952;
};
template<class stce_9759, class stce_9941, class stce_9942, class stce_9943, class stce_9922>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9779<stce_9759, stce_9941, stce_9942, stce_9943> &stce_9843,
                                   stce_9951<stce_9922>& stce_9940) {
  interpret(stce_9843.stce_9945, stce_9940);
  if(stce_9940.tape())
    {
      if(IS_DEBUG)
        {
          stce_9922 * tape = stce_9940.tape();
          interpret(stce_9843.stce_9947, stce_9940);
          if (stce_9940.tape() && tape != stce_9940.tape())
            throw dco::exception::create<std::runtime_error>(
                                                             "multiple tape error: different tapes on right-hand side of statement",
                                                             "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_folding.hpp", 159);
        }
      return;
    }
  interpret(stce_9843.stce_9947, stce_9940);
}
template<class stce_9759, class stce_9949, class stce_9943, class stce_9922>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9776<stce_9759, stce_9949, stce_9943> &stce_9843,
                                   stce_9951<stce_9922>& stce_9940) {
  interpret(stce_9843.stce_9950, stce_9940);
}
template<class stce_9759, class stce_9941, class stce_9943, class stce_9922>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9780<stce_9759, stce_9941, stce_9943> &stce_9843,
                                   stce_9951<stce_9922>& stce_9940) {
  interpret(stce_9843.stce_9945, stce_9940);
}
template<class stce_9759, class stce_9942, class stce_9943, class stce_9922>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9781<stce_9759, stce_9942, stce_9943> &stce_9843,
                                   stce_9951<stce_9922>& stce_9940) {
  interpret(stce_9843.stce_9947, stce_9940);
}
}
}
namespace dco {
namespace internal {
struct stce_9953 {
  template<class stce_9954, class stce_9955> static auto stce_9956(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957._value() == stce_9958._value()) { return stce_9957._value() == stce_9958._value(); } template<class stce_9954, class stce_9955> static auto stce_9959(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957._value() == stce_9958) { return stce_9957._value() == stce_9958; } template<class stce_9954, class stce_9955> static auto stce_9960(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957 == stce_9958._value()) { return stce_9957 == stce_9958._value(); }
  template<class stce_9954, class stce_9955> static auto stce_9961(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957._value() != stce_9958._value()) { return stce_9957._value() != stce_9958._value(); } template<class stce_9954, class stce_9955> static auto stce_9962(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957._value() != stce_9958) { return stce_9957._value() != stce_9958; } template<class stce_9954, class stce_9955> static auto stce_9963(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957 != stce_9958._value()) { return stce_9957 != stce_9958._value(); }
  template<class stce_9954, class stce_9955> static auto stce_9964(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957._value() < stce_9958._value()) { return stce_9957._value() < stce_9958._value(); } template<class stce_9954, class stce_9955> static auto stce_9965(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957._value() < stce_9958) { return stce_9957._value() < stce_9958; } template<class stce_9954, class stce_9955> static auto stce_9966(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957 < stce_9958._value()) { return stce_9957 < stce_9958._value(); }
  template<class stce_9954, class stce_9955> static auto stce_9967(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957._value() <= stce_9958._value()) { return stce_9957._value() <= stce_9958._value(); } template<class stce_9954, class stce_9955> static auto stce_9968(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957._value() <= stce_9958) { return stce_9957._value() <= stce_9958; } template<class stce_9954, class stce_9955> static auto stce_9969(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957 <= stce_9958._value()) { return stce_9957 <= stce_9958._value(); }
  template<class stce_9954, class stce_9955> static auto stce_9970(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957._value() > stce_9958._value()) { return stce_9957._value() > stce_9958._value(); } template<class stce_9954, class stce_9955> static auto stce_9971(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957._value() > stce_9958) { return stce_9957._value() > stce_9958; } template<class stce_9954, class stce_9955> static auto stce_9972(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957 > stce_9958._value()) { return stce_9957 > stce_9958._value(); }
  template<class stce_9954, class stce_9955> static auto stce_9973(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957._value() >= stce_9958._value()) { return stce_9957._value() >= stce_9958._value(); } template<class stce_9954, class stce_9955> static auto stce_9974(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957._value() >= stce_9958) { return stce_9957._value() >= stce_9958; } template<class stce_9954, class stce_9955> static auto stce_9975(const stce_9954 &stce_9957, const stce_9955 &stce_9958) -> decltype(stce_9957 >= stce_9958._value()) { return stce_9957 >= stce_9958._value(); }
};
template<class stce_9976, stce_9976 *&stce_9952>
struct stce_9977 {
  template<class stce_9954, class stce_9955> static bool stce_9956(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957._value() == stce_9958._value(); const bool stce_9978 = stce_9957.stce_9979(0) == stce_9958.stce_9979(0); const bool stce_9980 = stce_9957.stce_9979(1) == stce_9958.stce_9979(1); stce_9952->stce_9981("compare_AA", stce_9840, stce_9978,stce_9980); return stce_9840; } template<class stce_9954, class stce_9955> static bool stce_9959(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957._value() == stce_9958; const bool stce_9978 = stce_9957.stce_9979(0) == stce_9958; const bool stce_9980 = stce_9957.stce_9979(1) == stce_9958; stce_9952->stce_9981("compare_AP", stce_9840, stce_9978, stce_9980); return stce_9840; } template<class stce_9954, class stce_9955> static bool stce_9960(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957 == stce_9958._value(); const bool stce_9978 = stce_9957 == stce_9958.stce_9979(0); const bool stce_9980 = stce_9957 == stce_9958.stce_9979(1); stce_9952->stce_9981("compare_PA", stce_9840, stce_9978, stce_9980); return stce_9840; }
  template<class stce_9954, class stce_9955> static bool stce_9961(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957._value() != stce_9958._value(); const bool stce_9978 = stce_9957.stce_9979(0) != stce_9958.stce_9979(0); const bool stce_9980 = stce_9957.stce_9979(1) != stce_9958.stce_9979(1); stce_9952->stce_9981("compare_AA", stce_9840, stce_9978,stce_9980); return stce_9840; } template<class stce_9954, class stce_9955> static bool stce_9962(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957._value() != stce_9958; const bool stce_9978 = stce_9957.stce_9979(0) != stce_9958; const bool stce_9980 = stce_9957.stce_9979(1) != stce_9958; stce_9952->stce_9981("compare_AP", stce_9840, stce_9978, stce_9980); return stce_9840; } template<class stce_9954, class stce_9955> static bool stce_9963(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957 != stce_9958._value(); const bool stce_9978 = stce_9957 != stce_9958.stce_9979(0); const bool stce_9980 = stce_9957 != stce_9958.stce_9979(1); stce_9952->stce_9981("compare_PA", stce_9840, stce_9978, stce_9980); return stce_9840; }
  template<class stce_9954, class stce_9955> static bool stce_9964(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957._value() < stce_9958._value(); const bool stce_9978 = stce_9957.stce_9979(0) < stce_9958.stce_9979(0); const bool stce_9980 = stce_9957.stce_9979(1) < stce_9958.stce_9979(1); stce_9952->stce_9981("compare_AA", stce_9840, stce_9978,stce_9980); return stce_9840; } template<class stce_9954, class stce_9955> static bool stce_9965(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957._value() < stce_9958; const bool stce_9978 = stce_9957.stce_9979(0) < stce_9958; const bool stce_9980 = stce_9957.stce_9979(1) < stce_9958; stce_9952->stce_9981("compare_AP", stce_9840, stce_9978, stce_9980); return stce_9840; } template<class stce_9954, class stce_9955> static bool stce_9966(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957 < stce_9958._value(); const bool stce_9978 = stce_9957 < stce_9958.stce_9979(0); const bool stce_9980 = stce_9957 < stce_9958.stce_9979(1); stce_9952->stce_9981("compare_PA", stce_9840, stce_9978, stce_9980); return stce_9840; }
  template<class stce_9954, class stce_9955> static bool stce_9967(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957._value() <= stce_9958._value(); const bool stce_9978 = stce_9957.stce_9979(0) <= stce_9958.stce_9979(0); const bool stce_9980 = stce_9957.stce_9979(1) <= stce_9958.stce_9979(1); stce_9952->stce_9981("compare_AA", stce_9840, stce_9978,stce_9980); return stce_9840; } template<class stce_9954, class stce_9955> static bool stce_9968(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957._value() <= stce_9958; const bool stce_9978 = stce_9957.stce_9979(0) <= stce_9958; const bool stce_9980 = stce_9957.stce_9979(1) <= stce_9958; stce_9952->stce_9981("compare_AP", stce_9840, stce_9978, stce_9980); return stce_9840; } template<class stce_9954, class stce_9955> static bool stce_9969(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957 <= stce_9958._value(); const bool stce_9978 = stce_9957 <= stce_9958.stce_9979(0); const bool stce_9980 = stce_9957 <= stce_9958.stce_9979(1); stce_9952->stce_9981("compare_PA", stce_9840, stce_9978, stce_9980); return stce_9840; }
  template<class stce_9954, class stce_9955> static bool stce_9970(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957._value() > stce_9958._value(); const bool stce_9978 = stce_9957.stce_9979(0) > stce_9958.stce_9979(0); const bool stce_9980 = stce_9957.stce_9979(1) > stce_9958.stce_9979(1); stce_9952->stce_9981("compare_AA", stce_9840, stce_9978,stce_9980); return stce_9840; } template<class stce_9954, class stce_9955> static bool stce_9971(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957._value() > stce_9958; const bool stce_9978 = stce_9957.stce_9979(0) > stce_9958; const bool stce_9980 = stce_9957.stce_9979(1) > stce_9958; stce_9952->stce_9981("compare_AP", stce_9840, stce_9978, stce_9980); return stce_9840; } template<class stce_9954, class stce_9955> static bool stce_9972(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957 > stce_9958._value(); const bool stce_9978 = stce_9957 > stce_9958.stce_9979(0); const bool stce_9980 = stce_9957 > stce_9958.stce_9979(1); stce_9952->stce_9981("compare_PA", stce_9840, stce_9978, stce_9980); return stce_9840; }
  template<class stce_9954, class stce_9955> static bool stce_9973(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957._value() >= stce_9958._value(); const bool stce_9978 = stce_9957.stce_9979(0) >= stce_9958.stce_9979(0); const bool stce_9980 = stce_9957.stce_9979(1) >= stce_9958.stce_9979(1); stce_9952->stce_9981("compare_AA", stce_9840, stce_9978,stce_9980); return stce_9840; } template<class stce_9954, class stce_9955> static bool stce_9974(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957._value() >= stce_9958; const bool stce_9978 = stce_9957.stce_9979(0) >= stce_9958; const bool stce_9980 = stce_9957.stce_9979(1) >= stce_9958; stce_9952->stce_9981("compare_AP", stce_9840, stce_9978, stce_9980); return stce_9840; } template<class stce_9954, class stce_9955> static bool stce_9975(const stce_9954 &stce_9957, const stce_9955 &stce_9958) { const bool stce_9840 = stce_9957 >= stce_9958._value(); const bool stce_9978 = stce_9957 >= stce_9958.stce_9979(0); const bool stce_9980 = stce_9957 >= stce_9958.stce_9979(1); stce_9952->stce_9981("compare_PA", stce_9840, stce_9978, stce_9980); return stce_9840; }
};
}
}
namespace dco {
namespace internal {
template<class stce_9982>
class local_gradient_t {
private:
  typedef typename dco::mode<stce_9982>::tape_t tape_t;
  typedef typename tape_t::stce_9763 stce_9763;
  typedef typename tape_t::stce_9924 stce_9924;
  stce_9982& stce_9983;
  stce_9924 stce_9927;
  local_gradient_t(const local_gradient_t&);
public:
  local_gradient_t(stce_9982 &stce_9984)
    : stce_9983(stce_9984)
    , stce_9927(stce_9984.tape()->stce_9928())
  {}
  ~local_gradient_t() {
    stce_9983._tape_index() = stce_9927.index();
  }
  template<typename stce_9794>
    void put(const stce_9794 &stce_9843, const stce_9763& stce_9985) {
    if (!stce_9755 || stce_9843._tape_index() != 0) {
      bool stce_9986 = (stce_9927.index() == 0);
      stce_9927.insert(stce_9843._tape_index(), stce_9985);
      if (stce_9986)
        stce_9983.stce_9987(stce_9983._value(), stce_9816());
    }
  }
};
template<class stce_9891>
struct stce_9988 {
  typedef typename stce_9891::tape_t::stce_9764 derivative_t;
protected:
  mutable DCO_TAPE_INT stce_9889;
public:
  stce_9988() : stce_9889(0){}
  DCO_STRONG_INLINE stce_9988 &operator =(const stce_9988 &stce_9843) {
    stce_9889 = stce_9843.stce_9889;
    return *this;
  }
  void clear() {
    stce_9889 = 0;
  }
  derivative_t& stce_9989() const {
    return _adjoint();
  }
  derivative_t& _adjoint() const {
    if (IS_DEBUG && !stce_9891::global_tape)
      throw dco::exception::create<std::runtime_error>("Tape not created yet. Please use tape_t::create() before accessing derivative components of adjoint types.");
    return stce_9891::global_tape->_adjoint(stce_9889);
  }
  void stce_9990(const derivative_t &stce_9847) const {
    if (IS_DEBUG && !stce_9891::global_tape)
      throw dco::exception::create<std::runtime_error>("Tape not created yet. Please use tape_t::create() before accessing derivative components of adjoint types.");
    stce_9891::global_tape->_adjoint(stce_9889) = stce_9847;
  }
  bool _is_registered() const {
    return stce_9889 == 0 ? false : true;
  }
  DCO_TAPE_INT & _tape_index() { return stce_9889; }
  DCO_TAPE_INT const& _tape_index() const { return stce_9889; }
  template<class stce_9991>
  DCO_STRONG_INLINE void stce_9987(const stce_9991 &, const double &) { }
  template <typename stce_9931>
    static DCO_STRONG_INLINE void register_variable(stce_9931& stce_9843, DCO_TAPE_INT index, typename stce_9891::tape_t*) {
    if (IS_DEBUG && !stce_9891::global_tape) {
      throw dco::exception::create<std::runtime_error>("Global tape not yet created (is NULL).");
    }
    stce_9843.stce_9889 = index;
  }
  static typename stce_9891::tape_t *tape() {
    return stce_9891::global_tape;
  }
  void set_tape(typename stce_9891::tape_t *) {}
};
template<typename stce_9891>
  class single_tape_data : public stce_9988<stce_9891> {
  typedef stce_9988<stce_9891> stce_9992;
public:
  typedef stce_9891 mode_t;
  typedef stce_9953 stce_9993;
  typedef typename stce_9891::tape_t tape_t;
  typedef typename tape_t::stce_9764 derivative_t;
  single_tape_data() {}
  DCO_STRONG_INLINE single_tape_data &operator =(const single_tape_data &stce_9843) {
    stce_9992::operator= (static_cast<const stce_9992&>(stce_9843));
    return *this;
  }
  template<class stce_9892>
  DCO_ALWAYS_INLINE void stce_9893(const stce_9892 &stce_9994) {
    if (!stce_9891::global_tape || !stce_9891::global_tape->is_active()) {
      this->clear();
      return;
    }
    using namespace folding;
    stce_9923<tape_t> stce_9995(stce_9891::global_tape);
    interpret(stce_9994, stce_9995);
    this->stce_9889 = stce_9995.index();
  }
};
template<typename stce_9891>
  class stce_9996 {
  typedef typename stce_9891::tape_t tape_t;
  mutable DCO_TAPE_INT stce_9889;
  mutable tape_t *stce_9997;
public:
  typedef stce_9891 mode_t;
  typedef stce_9953 stce_9993;
  typedef typename tape_t::stce_9764 derivative_t;
  stce_9996() : stce_9889(0), stce_9997(0) {}
  DCO_STRONG_INLINE stce_9996 &operator =(const stce_9996 &stce_9843) {
   stce_9889 = stce_9843.stce_9889;
   stce_9997 = stce_9843.stce_9997;
    return *this;
  }
  void clear() {
    stce_9889 = 0;
    stce_9997 = 0;
  }
  derivative_t& stce_9989() const {
    return _adjoint();
  }
  derivative_t &_adjoint() const {
    static derivative_t stce_9998 = 0;
    if (stce_9997 == 0) {
      stce_9998 = 0;
      return stce_9998;
    } else {
      return stce_9997->_adjoint(stce_9889);
    }
  }
  void stce_9990(const derivative_t &stce_9847) const {
    if (stce_9997 != 0)
      stce_9997->_adjoint(stce_9889) = stce_9847;
  }
  bool _is_registered() const {
    return stce_9889 == 0 ? false : true;
  }
  DCO_TAPE_INT &_tape_index() {
    return stce_9889;
  }
  const DCO_TAPE_INT &_tape_index() const {
    return stce_9889;
  }
  template<typename stce_9991>
  DCO_STRONG_INLINE void stce_9987(const stce_9991 &, const double &) { }
  template <typename stce_9931>
    static DCO_STRONG_INLINE void register_variable(stce_9931& stce_9843, DCO_TAPE_INT index, tape_t* tape) {
    stce_9843.stce_9997 = tape;
    stce_9843.stce_9889 = index;
  }
  tape_t *tape() const {
    return stce_9997;
  }
  void set_tape(tape_t *stce_9782) {
    stce_9997 = stce_9782;
  }
  template<class stce_9892>
  DCO_ALWAYS_INLINE void stce_9893(const stce_9892 &stce_9994) {
    using namespace folding;
    stce_9951<tape_t> stce_9951;
    interpret(stce_9994, stce_9951);
    tape_t *stce_9999 = stce_9951.tape();
    if (!stce_9999 || !stce_9999->is_active()) {
       clear();
       return;
     }
    stce_9923<tape_t> stce_9995(stce_9999);
    interpret(stce_9994, stce_9995);
    stce_9889 = stce_9995.index();
    stce_9997 = stce_9999;
  }
};
template<typename stce_9891, size_t stce_10000=2>
class stce_10001: public stce_9988<stce_9891> {
  typedef stce_9988<stce_9891> stce_9992;
public:
  typedef stce_9891 mode_t;
  typedef stce_9977<typename stce_9891::tape_t, stce_9891::global_tape> stce_9993;
  typedef typename stce_9891::tape_t tape_t;
  typedef typename tape_t::stce_9764 derivative_t;
  mutable derivative_t stce_9936;
  mutable derivative_t stce_10002[stce_10000];
  void stce_10003(const derivative_t &value, const double &stce_9937) const {
    if(stce_9891::global_tape->stce_10004 == 0)
      this->stce_10002[0] = this->stce_10002[1] = value;
    stce_9891::global_tape->stce_10005(this, stce_9937);
    stce_9936 = stce_9937;
  }
public:
  stce_10001(): stce_9936(0) {}
  DCO_STRONG_INLINE stce_10001 &operator =(const stce_10001 &stce_9843) {
    stce_9992::operator= (static_cast<const stce_9992&>(stce_9843));
    stce_9936 = stce_9843.stce_9936;
    return *this;
  }
  inline void clear() {
    stce_9992::clear();
    stce_9936 = 0;
  }
  template<class stce_9991>
  DCO_STRONG_INLINE void stce_9987(const stce_9991 &value, const double &stce_9937) {
    stce_10003(value, stce_9937);
  }
  template <typename stce_9931>
    static DCO_STRONG_INLINE void register_variable(stce_9931& stce_9843, DCO_TAPE_INT index, tape_t* tape) {
    stce_9992::register_variable(stce_9843, index, tape);
    stce_9843.stce_10003(trait_value<stce_9931>::value(stce_9843), stce_9816());
  }
  const derivative_t& stce_9979(const derivative_t &value, const int stce_9808) const {
      if(this->stce_9889==0) return value;
      return stce_10002[stce_9808];
  }
  template<class stce_9892>
    DCO_ALWAYS_INLINE void stce_9893(const stce_9892 &stce_9994) {
    if (!stce_9891::global_tape || !stce_9891::global_tape->is_active()) {
      clear();
      return;
    }
    {
      using namespace folding;
      stce_9934<tape_t> stce_9995(stce_9891::global_tape);
      interpret(stce_9994, stce_9995);
      this->stce_9889 = stce_9995.index();
      if (this->stce_9889 == 0)
        return;
      stce_9936 = stce_9995.stce_9937();
    }
    stce_10002[0] = stce_9994.stce_9979(0);
    stce_10002[1] = stce_9994.stce_9979(1);
    derivative_t stce_10006 = stce_10002[1] - stce_10002[0];
    double stce_10007 = stce_10006/(2*stce_9891::global_tape->stce_10004);
    stce_9891::global_tape->stce_10008(stce_9936, stce_10007);
    stce_9891::global_tape->stce_10009();
    if(stce_9891::global_tape->stce_10010) {
      const double stce_10011=3.1415692;
      stce_9891::global_tape->_adjoint(this->stce_9889) = stce_10011;
      stce_9891::global_tape->interpret_adjoint();
      double stce_10012=0;
      for(size_t stce_9832=0;stce_9832<stce_9891::global_tape->stce_10013.size();++stce_9832)
        stce_10012 += stce_9891::global_tape->_adjoint(stce_9891::global_tape->stce_10013[stce_9832]) * stce_9891::global_tape->stce_10014[stce_9832];
      stce_9891::global_tape->zero_adjoints();
      stce_9891::global_tape->stce_10015(stce_10012, stce_9936 * stce_10011, stce_10007 * stce_10011);
    }
    stce_9891::global_tape->stce_10009();
  }
};
template <typename stce_10016>
  struct stce_10017 : public dco::mode<typename stce_10016::binary_t>::data_t {
  typedef stce_10016 mode_t;
  static const int gbcp_order = mode_t::gbcp_order;
  static const int order = mode_t::order;
  template<class stce_9892>
    DCO_ALWAYS_INLINE void stce_9893(stce_9892 const&) {
    }
  void clear() {}
  stce_10017& operator=(stce_10017 const&) {
    return *this;
  }
};
}
}
namespace dco
{
  template <class stce_9782>
  class ga1s_codegen {
    typedef stce_9762<stce_9782> stce_10018;
public:
    static const memory_model::TYPE TAPE_MODEL = memory_model::stce_9767;
    typedef dco::internal::stce_9888<ga1s_codegen> data;
    typedef typename stce_10018::value_t value_t;
    typedef typename stce_10018::value_t active_value_t;
    typedef typename stce_10018::stce_9764 derivative_t;
    typedef typename stce_10018::stce_9763 stce_9763;
    typedef typename mode<value_t>::passive_t passive_t;
    typedef dco::internal::active_type<derivative_t, data> type;
    typedef type active_t;
    static const bool is_dco_type = true;
    static const bool is_adjoint_type = true;
    static const bool is_tangent_type = false;
    static const int order = dco::mode<stce_9782>::order + 1;
    typedef typename internal::stce_9875<derivative_t, stce_9763> tape_t;
    static DCOi_THREAD_LOCAL tape_t *global_tape;
    typedef typename internal::local_gradient_t<type> local_gradient_t;
  };
template <class stce_9782>
  DCOi_THREAD_LOCAL typename ga1s_codegen<stce_9782>::tape_t *ga1s_codegen<stce_9782>::global_tape;
}
       
       
       
namespace dco {
  template <class stce_9759>
  class stce_10019 {
  public:
    typedef stce_9759 active_value_t;
    typedef stce_9759 type;
    typedef type active_t;
    typedef stce_9759 value_t;
    typedef stce_9759 passive_t;
    typedef void derivative_t;
    typedef void stce_10020;
    typedef void tape_t;
    typedef void local_gradient_t;
    typedef void local_gradient_with_activity_t;
    typedef void external_adjoint_object_t;
    typedef internal::jacobian_preaccumulator_t<void> jacobian_preaccumulator_t;
    static const bool is_dco_type = false;
    static const bool is_adjoint_type = false;
    static const bool is_tangent_type = false;
    static const int order = 0;
    static tape_t *global_tape;
  };
template <class stce_9759>
  typename stce_10019<stce_9759>::tape_t *stce_10019<stce_9759>::global_tape = 0;
}
namespace dco {
template <class stce_9795>
  struct mode : public dco::stce_10019<stce_9795> {
  static const bool is_intermediate_type = false;
};
template <typename stce_9771, typename stce_9772>
  struct mode<internal::active_type<stce_9771, stce_9772> > :
    public stce_9772::mode_t {
  static const bool is_intermediate_type = false;
};
template <typename stce_9771, typename stce_9772>
  struct mode<internal::active_type<stce_9771, stce_9772>& > :
    public stce_9772::mode_t {
  static const bool is_intermediate_type = false;
};
template <typename stce_9782, typename stce_10021, typename stce_10022>
  struct mode<internal::stce_9776<stce_9782, stce_10021, stce_10022> > :
    public internal::stce_9776<stce_9782, stce_10021, stce_10022>::data_t::mode_t {
  static const bool is_intermediate_type = true;
};
template<class stce_9782, class stce_10023, class stce_10024, class stce_10022>
  struct mode<internal::stce_9779<stce_9782, stce_10023, stce_10024, stce_10022> > :
    public internal::stce_9779<stce_9782, stce_10023, stce_10024, stce_10022>::data_t::mode_t {
  static const bool is_intermediate_type = true;
};
template<class stce_9782, class stce_10021, class stce_10022>
  struct mode<internal::stce_9780<stce_9782, stce_10021, stce_10022> > :
    public internal::stce_9780<stce_9782, stce_10021, stce_10022>::data_t::mode_t {
  static const bool is_intermediate_type = true;
};
template<class stce_9782, class stce_10021, class stce_10022>
  struct mode<internal::stce_9781<stce_9782, stce_10021, stce_10022> > :
    public internal::stce_9781<stce_9782, stce_10021, stce_10022>::data_t::mode_t {
  static const bool is_intermediate_type = true;
};
template<class stce_9771, class DATA_TYPE>
  struct mode<internal::stce_9773<stce_9771, DATA_TYPE> > :
    public DATA_TYPE::mode_t {
  static const bool is_intermediate_type = true;
};
}
       
namespace dco {
template <typename stce_9782, typename stce_9824>
  stce_9782& create_gbcp(stce_9824& stce_10025, typename enable_if<!stce_9869<stce_9782>::value, stce_9782>::type * = 0) {
      return reinterpret_cast<stce_9782&>(stce_10025);
  }
template <typename stce_9782, typename stce_9824>
  stce_9782 create_gbcp(stce_9824& stce_10025, typename enable_if<stce_9869<stce_9782>::value, stce_9782>::type * = 0) {
      return reinterpret_cast<stce_9782>(&stce_10025);
  }
namespace internal {
template<typename stce_9782>
  struct stce_10026 {
  static const bool value = false;
};
template <typename mode_t>
  struct stce_10026<stce_10017<mode_t> > {
  static const bool value = true;
};
template <typename stce_10027, typename data_t>
  struct stce_10026<active_type<stce_10027, data_t> > {
  static const bool value = stce_10026<data_t>::value;
};
template<typename stce_9782>
  struct is_gbcp_type {
  static const bool value = stce_10026<typename dco::remove_all<stce_9782>::type>::value;
};
template<typename stce_9782, int stce_9843> struct stce_10028;
template<typename stce_9782> struct stce_10028<stce_9782, 1> {
    typedef typename mode<stce_9782>::active_value_t& stce_9802;
    typedef const typename mode<stce_9782>::active_value_t& stce_10029;
    static stce_9802 value(stce_9782& at) { return at._value(); }
    static stce_10029 value(const stce_9782& at) { return at._value(); }
};
template<typename stce_9782, int stce_9843> struct stce_10028 {
    typedef typename stce_10028<typename mode<stce_9782>::active_value_t, stce_9843-1>::stce_9802 stce_9802;
    typedef typename stce_10028<typename mode<stce_9782>::active_value_t, stce_9843-1>::stce_10029 stce_10029;
    static stce_9802 value(stce_9782& at) { return stce_10028<typename mode<stce_9782>::active_value_t, stce_9843-1>::value(at._value()); }
    static stce_10029 value(const stce_9782& at) { return stce_10028<typename mode<stce_9782>::active_value_t, stce_9843-1>::value(at._value()); }
};
}
}
namespace dco {
template <typename stce_9782,typename stce_10030=stce_9782>
  struct dco_type_constructable_from {};
template <typename stce_9782>
  struct dco_explicit_type_cast_to {
  const static bool value = false;
};
}
namespace dco {
namespace internal {
template<class stce_9771>
struct stce_10031 {
  const stce_9771 stce_10032;
  stce_10031(stce_9771 const& value) :
      stce_10032(value) {
  }
  inline const stce_9771 &_value() const { return stce_10032; }
  template <typename TYPE,
            typename stce_10033 = typename enable_if<dco_explicit_type_cast_to<TYPE>::value>::type
            > explicit operator TYPE() { return static_cast<TYPE>(_value()); }
};
template<class stce_9771, class DATA_TYPE>
struct stce_9773 : public stce_10031<stce_9771> {
  stce_9773(const stce_9771& value) :
    stce_10031<stce_9771>(value) {
  }
    template <typename stce_9782 = DATA_TYPE>
    DCO_STRONG_INLINE operator typename dco::enable_if<is_gbcp_type<stce_9782>::value, typename stce_10028<stce_9773, stce_9782::gbcp_order>::stce_10029>::type const&() const {
      return stce_10028<stce_9773, stce_9782::gbcp_order>::value(*this);
    }
};
template<class stce_9771>
  struct stce_9773<stce_9771, typename ga1s_codegen<stce_9771>::data>
: public stce_10031<stce_9771>, public stce_9900<stce_9771> {
  stce_9773(const stce_9771& value) :
    stce_10031<stce_9771>(value) {
  }
};
template<class stce_9771, class stce_9774, class stce_9775>
struct stce_9776: public stce_9773<stce_9771, typename stce_9774::data_t>
  {
      const stce_9774 DCO_TEMPORARY_REFORCOPY stce_9950;
      typedef typename stce_9774::data_t data_t;
      typedef typename data_t::mode_t::value_t VALUE_TYPE;
      typedef typename data_t::mode_t::passive_t stce_10034;
    stce_9776(const stce_9774 &arg)
      :stce_9773<stce_9771, data_t>(stce_9775::stce_10035(arg._value())),
       stce_9950(arg) {
      }
    inline const stce_9771 stce_9979(const int stce_9808) const {
        return stce_9775::stce_10035(stce_9950.stce_9979(stce_9808));
    }
    inline const stce_9771 pval() const {
      return stce_9775::stce_10036(this->_value(), stce_9950._value());
    }
    void stce_9895(std::ofstream &out) const {
      stce_9950.stce_9895(out);
      this->stce_9903(out, stce_9775::stce_10037(this->str(stce_9950)));
    }
    void stce_9897(std::ofstream &out) const {
      std::string stce_9904 = stce_9775::stce_10038();
      this->stce_9907(stce_9904, stce_9950, "x");
      out << this->str(stce_9950) << "_a1s +=(" << stce_9904 << ")*" << this->stce_9890() << "_a1s;\t";
      stce_9950.stce_9897(out);
    }
};
    template<class stce_9771, class stce_9777, class stce_9778, class stce_9775>
    struct stce_9779: public stce_9773<stce_9771, typename stce_9777::data_t> {
      const stce_9777 DCO_TEMPORARY_REFORCOPY stce_9945;
      const stce_9778 DCO_TEMPORARY_REFORCOPY stce_9947;
      typedef typename stce_9777::data_t data_t;
      typedef typename data_t::mode_t::value_t VALUE_TYPE;
      typedef typename data_t::mode_t::passive_t stce_10034;
      stce_9779(const stce_9777 &stce_9912, const stce_9778 &stce_9913)
        :stce_9773<stce_9771, data_t>(stce_9775::stce_10035(stce_9912._value(), stce_9913._value())),
        stce_9945(stce_9912),
        stce_9947(stce_9913) {}
      inline const stce_9771 stce_9946() const {
        return stce_9775::stce_10039(this->_value(), stce_9945._value(), stce_9947._value());
      }
      inline const stce_9771 stce_9948() const {
        return stce_9775::stce_10040(this->_value(), stce_9945._value(), stce_9947._value());
      }
      inline const stce_9771 stce_9979(const int stce_9808) const {
        return stce_9775::stce_10035(stce_9945.stce_9979(stce_9808),stce_9947.stce_9979(stce_9808));
      }
      void stce_9895(std::ofstream &out) const {
        stce_9945.stce_9895(out);
        stce_9947.stce_9895(out);
        this->stce_9903(out, stce_9775::stce_10037(this->str(stce_9945),
                                                                                   this->str(stce_9947)));
      }
      void stce_9897(std::ofstream &out) const {
        this->template stce_9911<stce_9775>(out, stce_9945, stce_9947);
        this->template stce_9915<stce_9775>(out, stce_9945, stce_9947);
      }
    };
    template<class stce_9771, class stce_9777, class stce_9775>
    struct stce_9780 : public stce_9773<stce_9771, typename stce_9777::data_t> {
    const stce_9777 DCO_TEMPORARY_REFORCOPY stce_9945;
    const stce_9771 stce_9947;
    typedef typename stce_9777::data_t data_t;
    typedef typename data_t::mode_t::value_t VALUE_TYPE;
    typedef typename data_t::mode_t::passive_t stce_10034;
      stce_9780(const stce_9777 &stce_9912, const stce_9771 &stce_9913)
      :stce_9773<stce_9771, data_t>(stce_9775::stce_10035(stce_9912._value(), stce_9913)),
        stce_9945(stce_9912), stce_9947(stce_9913) {
      }
    inline const stce_9771 stce_9979(const int stce_9808) const {
      return stce_9775::stce_10035(stce_9945.stce_9979(stce_9808), stce_9947);
    }
    inline const stce_9771 stce_9946() const {
        return stce_9775::stce_10039(this->stce_10032, stce_9945._value(), stce_9947);
    }
    void stce_9895(std::ofstream &out) const {
      stce_9945.stce_9895(out);
      this->stce_9903(out, stce_9775::stce_10037(this->str(stce_9945),
          this->str(stce_9947)));
    }
    void stce_9897(std::ofstream &out) const {
      this->template stce_9911<stce_9775>(out, stce_9945, stce_9947);
      this->template stce_9915<stce_9775>(out, stce_9945, stce_9947);
    }
};
    template<class stce_9771, class stce_9778, class stce_9775>
      struct stce_9781: public stce_9773<stce_9771, typename stce_9778::data_t> {
      const stce_9771 stce_9945;
      const stce_9778 DCO_TEMPORARY_REFORCOPY stce_9947;
      typedef typename stce_9778::data_t data_t;
      typedef typename data_t::mode_t::value_t VALUE_TYPE;
      typedef typename data_t::mode_t::passive_t stce_10034;
      stce_9781(const stce_9771 &stce_9912, const stce_9778 &stce_9913)
        :stce_9773<stce_9771, data_t>(stce_9775::stce_10035(stce_9912, stce_9913._value())),
        stce_9945(stce_9912), stce_9947(stce_9913) {
      }
      inline const stce_9771 stce_9979(const int stce_9808) const {
        return stce_9775::stce_10035(stce_9945, stce_9947.stce_9979(stce_9808));
      }
      inline const stce_9771 stce_9948() const {
        return stce_9775::stce_10040(this->stce_10032, stce_9945, stce_9947._value());
      }
      void stce_9895(std::ofstream &out) const {
        stce_9947.stce_9895(out);
        this->stce_9903(out, stce_9775::stce_10037(this->str(stce_9945),
                                                                                   this->str(stce_9947)));
      }
      void stce_9897(std::ofstream &out) const {
        this->template stce_9911<stce_9775>(out, stce_9945, stce_9947);
        this->template stce_9915<stce_9775>(out, stce_9945, stce_9947);
      }
    };
}
}
namespace dco {
namespace internal {
template <typename stce_9782>
struct passive_value_type_of {
    typedef stce_9782 TYPE;
};
template <typename stce_9782>
  class stce_10041 {
protected:
  stce_9782 stce_10032;
public:
  stce_10041() : stce_10032(0) {}
  template <typename stce_9824>
    stce_10041(const stce_9824& stce_9843) : stce_10032(stce_9843) {}
};
template<class stce_9771, class stce_9772>
struct active_type : private stce_10041<stce_9771>, public stce_9772 {
  typedef stce_10041<stce_9771> stce_10042;
  using stce_10042::stce_10032;
public:
  typedef stce_9772 data_t;
  typedef typename data_t::mode_t::value_t VALUE_TYPE;
  typedef typename data_t::mode_t::passive_t stce_10034;
  DCO_STRONG_INLINE const stce_9771& _value() const {
    return stce_10032;
  }
  stce_9772& data() {
    return *static_cast<stce_9772*>(this);
  }
  DCO_STRONG_INLINE stce_9771 &_value() {
    return stce_10032;
  }
  DCO_STRONG_INLINE const stce_9771 stce_9979(const int stce_9808) const {
        return stce_9772::stce_9979(stce_10032, stce_9808);
    }
  DCO_STRONG_INLINE active_type() :
    stce_10042(static_cast<stce_9771>(0.0)) {
  }
    template <typename stce_9782 = stce_9772>
    DCO_STRONG_INLINE operator typename dco::enable_if<is_gbcp_type<stce_9782>::value, typename stce_10028<active_type, stce_9782::gbcp_order>::stce_10029>::type const&() const {
      return stce_10028<active_type, stce_9782::gbcp_order>::value(*this);
    }
    template <typename stce_9782 = stce_9772>
    DCO_STRONG_INLINE operator typename dco::enable_if<is_gbcp_type<stce_9782>::value, typename stce_10028<active_type, stce_9782::gbcp_order>::stce_9802>::type&() {
      return stce_10028<active_type, stce_9782::gbcp_order>::value(*this);
    }
    template <typename stce_9824, typename stce_9782 = stce_9772>
    DCO_STRONG_INLINE typename dco::enable_if<(dco::mode<stce_9824>::order > 0 && dco::internal::is_gbcp_type<stce_9782>::value &&
                                               dco::helper::is_same<typename stce_9782::mode_t::active_t, typename dco::mode<stce_9824>::active_t>::value), active_type>::type& operator =(const stce_9824 &stce_9994) {
      stce_10032 = stce_9994;
      return *this;
    }
    template<typename stce_9824, typename stce_9782 = stce_9772>
    DCO_STRONG_INLINE active_type(const stce_9824 &stce_9753, typename dco::enable_if<(dco::mode<stce_9824>::order > 0 && dco::internal::is_gbcp_type<stce_9782>::value &&
                                                                         dco::helper::is_same<typename stce_9782::mode_t::active_t, typename dco::mode<stce_9824>::active_t>::value), void*>::type = NULL) :
      stce_10042(stce_9753) { }
  template<typename TYPE>
    DCO_STRONG_INLINE active_type(const TYPE &stce_9753, typename dco_type_constructable_from<TYPE>::type* = 0) :
      stce_10042(stce_9753) { }
  template<typename TYPE>
    DCO_STRONG_INLINE typename dco::enable_if<dco_type_constructable_from<TYPE>::value, active_type>::type& operator=(const TYPE &stce_9753) {
    stce_10032 = stce_9753;
    return *this;
  }
  template <typename TYPE, typename stce_9782 = stce_9772,
            typename stce_10043 = typename enable_if<!is_gbcp_type<stce_9782>::value>::type,
            typename stce_10033 = typename enable_if<dco_explicit_type_cast_to<TYPE>::value>::type
            > explicit operator TYPE() { return static_cast<TYPE>(_value()); }
  template<typename TYPE>
  DCO_STRONG_INLINE active_type(const TYPE& stce_9753, typename dco::enable_if<stce_9872<typename TYPE::data_t::mode_t, typename data_t::mode_t>::value, void*>::type = NULL) :
     stce_10042(stce_9753) { }
  DCO_STRONG_INLINE active_type(const stce_10034 &stce_9753) : stce_10042(stce_9753) { }
  DCO_STRONG_INLINE active_type &operator =(const active_type &stce_9843) {
    stce_9772::operator= (static_cast<const stce_9772&>(stce_9843));
    stce_10032 = stce_9843.stce_10032;
    return *this;
  }
  DCO_STRONG_INLINE active_type &operator =(const stce_10034 &stce_9994) {
    this->stce_10032 = stce_9994;
    this->clear();
    return *this;
  }
  template<class stce_10044, class stce_10045>
  DCO_STRONG_INLINE typename dco::enable_if<stce_9872<typename stce_10045::mode_t, typename data_t::mode_t>::value, active_type>::type& operator =(const active_type<stce_10044, stce_10045> &stce_9994) {
    stce_10032 = stce_9994;
    this->clear();
    return *this;
  }
  template<typename TYPE>
  DCO_STRONG_INLINE typename dco_type_constructable_from<TYPE, active_type>::type& operator=(const TYPE &stce_9753){
    stce_10032 = stce_9753;
    this->clear();
    return *this;
  }
    template<class stce_9941, class stce_9942, class stce_9943 > DCO_ALWAYS_INLINE active_type(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843) { this->stce_9893(stce_9843); this->stce_10032 = stce_9843.stce_10032; } template<class stce_9941, class stce_9942, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator=(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843) { this->stce_9893(stce_9843); this->stce_10032 = stce_9843.stce_10032; return *this; }
    template<class stce_9941, class stce_9943 > DCO_ALWAYS_INLINE active_type(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843) { this->stce_9893(stce_9843); this->stce_10032 = stce_9843.stce_10032; } template<class stce_9941, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator=(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843) { this->stce_9893(stce_9843); this->stce_10032 = stce_9843.stce_10032; return *this; }
    template<class stce_9942, class stce_9943 > DCO_ALWAYS_INLINE active_type(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843) { this->stce_9893(stce_9843); this->stce_10032 = stce_9843.stce_10032; } template<class stce_9942, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator=(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843) { this->stce_9893(stce_9843); this->stce_10032 = stce_9843.stce_10032; return *this; }
    template<class stce_9949, class stce_9943 > DCO_ALWAYS_INLINE active_type(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843) { this->stce_9893(stce_9843); this->stce_10032 = stce_9843.stce_10032; } template<class stce_9949, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator=(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843) { this->stce_9893(stce_9843); this->stce_10032 = stce_9843.stce_10032; return *this; }
    template <typename stce_9824, typename stce_9782 = stce_9772> DCO_ALWAYS_INLINE typename enable_if<(dco::mode<stce_9824>::order > 0 && dco::internal::is_gbcp_type<stce_9782>::value && dco::helper::is_same<typename stce_9782::mode_t::active_t, typename dco::mode<stce_9824>::active_t>::value ), active_type>::type& operator += (const stce_9824 &stce_9843){ *this = *this + stce_9843; return *this; } template<class stce_10045> DCO_ALWAYS_INLINE active_type& operator += (const active_type<stce_9771, stce_10045> &stce_9843){ *this = *this + stce_9843; return *this; } template<class stce_9941, class stce_9942, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator += (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843){ *this = *this + stce_9843; return *this; } template<class stce_9941, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator += (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843){ *this = *this + stce_9843; return *this; } template<class stce_9942, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator += (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843){ *this = *this + stce_9843; return *this; } template<class stce_9949, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator += (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843){ *this = *this + stce_9843; return *this; } DCO_ALWAYS_INLINE active_type& operator += (const VALUE_TYPE &stce_9843) { this->_value() += stce_9843; return *this; }
    template <typename stce_9824, typename stce_9782 = stce_9772> DCO_ALWAYS_INLINE typename enable_if<(dco::mode<stce_9824>::order > 0 && dco::internal::is_gbcp_type<stce_9782>::value && dco::helper::is_same<typename stce_9782::mode_t::active_t, typename dco::mode<stce_9824>::active_t>::value ), active_type>::type& operator -= (const stce_9824 &stce_9843){ *this = *this - stce_9843; return *this; } template<class stce_10045> DCO_ALWAYS_INLINE active_type& operator -= (const active_type<stce_9771, stce_10045> &stce_9843){ *this = *this - stce_9843; return *this; } template<class stce_9941, class stce_9942, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator -= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843){ *this = *this - stce_9843; return *this; } template<class stce_9941, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator -= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843){ *this = *this - stce_9843; return *this; } template<class stce_9942, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator -= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843){ *this = *this - stce_9843; return *this; } template<class stce_9949, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator -= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843){ *this = *this - stce_9843; return *this; } DCO_ALWAYS_INLINE active_type& operator -= (const VALUE_TYPE &stce_9843) { this->_value() -= stce_9843; return *this; }
    template <typename stce_9824, typename stce_9782 = stce_9772> DCO_ALWAYS_INLINE typename enable_if<(dco::mode<stce_9824>::order > 0 && dco::internal::is_gbcp_type<stce_9782>::value && dco::helper::is_same<typename stce_9782::mode_t::active_t, typename dco::mode<stce_9824>::active_t>::value ), active_type>::type& operator *= (const stce_9824 &stce_9843){ *this = *this * stce_9843; return *this; } template<class stce_10045> DCO_ALWAYS_INLINE active_type& operator *= (const active_type<stce_9771, stce_10045> &stce_9843){ *this = *this * stce_9843; return *this; } template<class stce_9941, class stce_9942, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator *= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843){ *this = *this * stce_9843; return *this; } template<class stce_9941, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator *= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843){ *this = *this * stce_9843; return *this; } template<class stce_9942, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator *= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843){ *this = *this * stce_9843; return *this; } template<class stce_9949, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator *= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843){ *this = *this * stce_9843; return *this; } DCO_STRONG_INLINE active_type& operator *= (const VALUE_TYPE &stce_9843) { *this = *this * stce_9843; return *this; }
    template <typename stce_9824, typename stce_9782 = stce_9772> DCO_ALWAYS_INLINE typename enable_if<(dco::mode<stce_9824>::order > 0 && dco::internal::is_gbcp_type<stce_9782>::value && dco::helper::is_same<typename stce_9782::mode_t::active_t, typename dco::mode<stce_9824>::active_t>::value ), active_type>::type& operator /= (const stce_9824 &stce_9843){ *this = *this / stce_9843; return *this; } template<class stce_10045> DCO_ALWAYS_INLINE active_type& operator /= (const active_type<stce_9771, stce_10045> &stce_9843){ *this = *this / stce_9843; return *this; } template<class stce_9941, class stce_9942, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator /= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843){ *this = *this / stce_9843; return *this; } template<class stce_9941, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator /= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843){ *this = *this / stce_9843; return *this; } template<class stce_9942, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator /= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843){ *this = *this / stce_9843; return *this; } template<class stce_9949, class stce_9943 > DCO_ALWAYS_INLINE active_type& operator /= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843){ *this = *this / stce_9843; return *this; } DCO_STRONG_INLINE active_type& operator /= (const VALUE_TYPE &stce_9843) { *this = *this / stce_9843; return *this; }
  DCO_STRONG_INLINE active_type &operator ++() {
    ++this->stce_10032;
    return *this;
  }
  DCO_STRONG_INLINE active_type &operator --() {
    --this->stce_10032;
    return *this;
  }
  DCO_STRONG_INLINE active_type operator ++(int) {
    active_type stce_9840(*this);
    ++this->stce_10032;
    return stce_9840;
  }
  DCO_STRONG_INLINE active_type operator --(int) {
    active_type stce_9840(*this);
    --this->stce_10032;
    return stce_9840;
  }
  std::string stce_9890() const {
    if (stce_9772::_tape_index() == 0) {
      std::stringstream stce_10046;
      return stce_10046.str();
    } else
      return stce_9772::stce_9890();
  }
  void stce_9895(const std::ofstream &) const {
  }
  void stce_9897(const std::ofstream &) const {
  }
};
template<class stce_9771, class stce_9772>
struct passive_value_type_of<active_type<stce_9771, stce_9772> > {
  typedef typename passive_value_type_of<stce_9771>::TYPE TYPE;
};
}
}
namespace dco {
  namespace internal {
  template<typename stce_9891>
    struct ts_data {
    typedef stce_9891 mode_t;
    typedef stce_9953 stce_9993;
    typedef typename mode_t::derivative_t derivative_t;
    mutable derivative_t stce_10047;
    ts_data() : stce_10047(0) {}
    ts_data &operator = (const ts_data &stce_10048) {
      stce_10047 = stce_10048.stce_10047;
      return *this;
    }
    void stce_9990(const derivative_t &stce_9753) {
        stce_10047 = stce_9753;
      }
      void stce_10049(derivative_t &stce_9753) const {
        stce_9753 = stce_10047;
      }
      void clear() {
        stce_10047 = 0;
      }
      derivative_t &stce_9989() const
      {
       	stce_1000002();;
        return stce_10047;
      }
      template<class stce_9892>
      void stce_9893(const stce_9892 &stce_9994) {
        stce_10047 = stce_10050(stce_9994, 1.0);
      }
      template<class stce_9771, class stce_9939 >
      static derivative_t stce_10050(const dco::internal::active_type<stce_9771, stce_9939> &stce_9843, const derivative_t &pval) {
        return stce_9843.stce_10047 * pval;
      }
      template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 >
      static derivative_t stce_10050(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843, const derivative_t &pval) {
        return stce_10050(stce_9843.stce_9945, stce_9843.stce_9946() * pval) + stce_10050(stce_9843.stce_9947, stce_9843.stce_9948() * pval);
      }
      template<class stce_9771, class stce_9949, class stce_9943 >
      static derivative_t stce_10050(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843, const derivative_t &pval) {
        return stce_10050(stce_9843.stce_9950, stce_9843.pval() * pval);
      }
      template<class stce_9771, class stce_9941, class stce_9943 >
      static derivative_t stce_10050(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843, const derivative_t &pval) {
        return stce_10050(stce_9843.stce_9945, stce_9843.stce_9946() * pval);
      }
      template<class stce_9771, class stce_9942, class stce_9943 >
      static derivative_t stce_10050(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843, const derivative_t &pval) {
        return stce_10050(stce_9843.stce_9947, stce_9843.stce_9948() * pval);
      }
      typedef void tape_t;
      typedef DCO_TAPE_INT stce_10051;
      void *tape() const {
        return NULL;
      }
      DCO_TAPE_INT _tape_index() const {
        return 0;
      }
    };
  }
}
namespace dco {
  namespace internal {
  template<typename stce_9891>
    struct stce_10052 {
    typedef stce_9891 mode_t;
    typedef stce_9953 stce_9993;
    typedef typename mode_t::derivative_t derivative_t;
    typedef typename derivative_t::value_t stce_10020;
private:
    typedef typename mode_t::value_t stce_9771;
public:
      static const int stce_9837;
      mutable derivative_t tlms;
      stce_10052() {
        for (int stce_9832 = 0; stce_9832 < stce_9837; ++stce_9832) tlms[stce_9832] = 0;
      }
    derivative_t &stce_9989() const {
        return tlms;
      }
      void stce_9990(const stce_10020 &stce_9753, const int stce_10053) {
        stce_1000002();;
        tlms[stce_10053] = stce_9753;
      }
      void stce_10049(stce_10020 &stce_9753, const int stce_10053) const {
        stce_1000002();;
        stce_9753 = tlms[stce_10053];
      }
      void clear() {
        for (int stce_9832 = 0; stce_9832 < stce_9837; ++stce_9832) {
          tlms[stce_9832] = 0.0;
        }
      }
      template<class stce_9892, class stce_9794>
      static void stce_9893(const stce_9892 &stce_9994, stce_9794 &stce_10054) {
        stce_10052 &data = const_cast<stce_10052 &>(stce_10054);
        for (int stce_9832 = 0; stce_9832 < stce_9837; ++stce_9832) {
          data.tlms[stce_9832] = stce_10050(stce_9994, 1.0, stce_9832);
        }
      }
      template<class stce_9939 >
      static stce_10020 stce_10050(const dco::internal::active_type<stce_9771, stce_9939> &stce_9843, const stce_10020 &pval, const int stce_10053) {
        return stce_9843.tlms[stce_10053] * pval;
      }
      template<class stce_9941, class stce_9942, class stce_9943 >
      static stce_10020 stce_10050(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843, const stce_10020 &pval, const int stce_10053) {
        return stce_10050(stce_9843.stce_9945, stce_9843.stce_9946() * pval, stce_10053) + stce_10050(stce_9843.stce_9947, stce_9843.stce_9948() * pval, stce_10053);
      }
      template<class stce_9949, class stce_9943 >
      static stce_10020 stce_10050(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843, const stce_10020 &pval, const int stce_10053) {
        return stce_10050(stce_9843.stce_9950, stce_9843.pval() * pval, stce_10053);
      }
      template<class stce_9941, class stce_9943 >
      static stce_10020 stce_10050(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843, const stce_10020 &pval, const int stce_10053) {
        return stce_10050(stce_9843.stce_9945, stce_9843.stce_9946() * pval, stce_10053);
      }
      template<class stce_9942, class stce_9943 >
      static stce_10020 stce_10050(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843, const stce_10020 &pval, const int stce_10053) {
        return stce_10050(stce_9843.stce_9947, stce_9843.stce_9948() * pval, stce_10053);
      }
      typedef void stce_10055;
      typedef DCO_TAPE_INT stce_10051;
      void *stce_9952() const {
        return NULL;
      }
      DCO_TAPE_INT _tape_index() const {
        return 0;
      }
    };
  template<typename stce_9891>
    struct stce_10056 {
      typedef stce_9891 mode_t;
      typedef stce_9953 stce_9993;
      typedef typename mode_t::derivative_t derivative_t;
      typedef typename derivative_t::value_t stce_10020;
protected:
      typedef typename mode_t::value_t stce_9771;
public:
      static const int stce_9837 = derivative_t::stce_9828;
      mutable bool stce_10057;
      mutable derivative_t tlms;
      stce_10056() : stce_10057(false) {
        for (int stce_9832 = 0; stce_9832 < stce_9837; ++stce_9832) tlms[stce_9832] = 0;
      }
      void stce_9990(const stce_10020 &stce_9753, const int stce_10053) {
        if (!stce_10057) {
          stce_10057 = true;
        }
        tlms[stce_10053] = stce_9753;
      }
      derivative_t &stce_9989() const {
        stce_1000002();;
        stce_10057 = true;
        return tlms;
      }
      void stce_10049(stce_10020 &stce_9753, const int stce_10053) const {
        stce_1000002();
        if (stce_10057) {
          stce_9753 = tlms[stce_10053];
        } else {
          stce_9753 = 0;
        }
      }
      void clear() {
        for (int stce_9832 = 0; stce_9832 < stce_9837; ++stce_9832) {
          tlms[stce_9832] = 0.0;
        }
      }
      template<class stce_9892>
      DCO_ALWAYS_INLINE void stce_9893(const stce_9892 &stce_9994) {
        bool stce_10058 = this->stce_10059(stce_9994);
        if (stce_10058) {
          for (int stce_9832 = 0; stce_9832 < stce_9837; ++stce_9832) {
            tlms[stce_9832] = stce_10050(stce_9994, 1.0, stce_9832);
          }
          stce_10057 = true;
        } else
          stce_10057 = false;
      }
      template<class stce_9939 >
      static stce_10020 stce_10050(const dco::internal::active_type<stce_9771, stce_9939> &stce_9843, const stce_10020 &pval, const int stce_10053) {
        if (stce_9843.stce_10057) {
          return stce_9843.tlms[stce_10053] * pval;
        } else
          return 0;
      }
      template<class stce_9941, class stce_9942, class stce_9943 >
      static stce_10020 stce_10050(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843, const stce_10020 &pval, const int stce_10053) {
        return stce_10050(stce_9843.stce_9945, stce_9843.stce_9946() * pval, stce_10053) + stce_10050(stce_9843.stce_9947, stce_9843.stce_9948() * pval, stce_10053);
      }
      template<class stce_9949, class stce_9943 >
      static stce_10020 stce_10050(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843, const stce_10020 &pval, const int stce_10053) {
        return stce_10050(stce_9843.stce_9950, stce_9843.pval() * pval, stce_10053);
      }
      template<class stce_9941, class stce_9943 >
      static stce_10020 stce_10050(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843, const stce_10020 &pval, const int stce_10053) {
        return stce_10050(stce_9843.stce_9945, stce_9843.stce_9946() * pval, stce_10053);
      }
      template<class stce_9942, class stce_9943 >
      static stce_10020 stce_10050(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843, const stce_10020 &pval, const int stce_10053) {
        return stce_10050(stce_9843.stce_9947, stce_9843.stce_9948() * pval, stce_10053);
      }
      template<class stce_9939 >
      static bool stce_10059(const dco::internal::active_type<stce_9771, stce_9939> &stce_9843) {
        return stce_9843.stce_10057;
      }
      template<class stce_9941, class stce_9942, class stce_9943 >
      static bool stce_10059(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843) {
        return stce_10059(stce_9843.stce_9945) || stce_10059(stce_9843.stce_9947);
      }
      template<class stce_9949, class stce_9943 >
      static bool stce_10059(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843) {
        return stce_10059(stce_9843.stce_9950);
      }
      template<class stce_9941, class stce_9943 >
      static bool stce_10059(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843) {
        return stce_10059(stce_9843.stce_9945);
      }
      template<class stce_9942, class stce_9943 >
      static bool stce_10059(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843) {
        return stce_10059(stce_9843.stce_9947);
      }
      typedef void stce_10055;
      typedef void tape_t;
      typedef DCO_TAPE_INT stce_10051;
      void *tape() const {
        return NULL;
      }
      DCO_TAPE_INT _tape_index() const {
        return 0;
      }
    };
  }
}
namespace dco {
  namespace internal {
  template<DCO_INTEGRAL_TAPE_INT size, class stce_9771>
    struct stce_10060 {
    typedef stce_9953 stce_9993;
    typedef std::bitset<size> stce_10061;
    typedef stce_9803<stce_9771> mode_t;
      stce_10061 stce_10062;
      stce_10060() {}
      DCO_STRONG_INLINE stce_10060 &operator = (const stce_10060 &stce_10048) {
        stce_10062 = stce_10048.stce_10062;
        return *this;
      }
       void stce_9990(const size_t &stce_9808) {
        stce_10062.set(stce_9808);
      }
       bool stce_10049(const size_t &stce_9808) const {
        return stce_10062[stce_9808];
      }
       void clear() {
        stce_10062.reset();
      }
      template<class stce_9892>
      DCO_ALWAYS_INLINE void stce_9893(const stce_9892 &stce_9994) {
        stce_10060 stce_10063;
        stce_10064(stce_9994, stce_10063.stce_10062);
        stce_10062 = stce_10063.stce_10062;
      }
      template<class stce_9939 >
      static void stce_10064(const dco::internal::active_type<stce_9771, stce_9939> &stce_9843, stce_10061 &stce_10054) {
        stce_10054 |= stce_9843.stce_10062;
      }
      template<class stce_9941, class stce_9942, class stce_9943 >
      static void stce_10064(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843, stce_10061 &stce_10054) {
        stce_10064(stce_9843.stce_9945, stce_10054);
        stce_10064(stce_9843.stce_9947, stce_10054);
      }
      template<class stce_9949, class stce_9943 >
      static void stce_10064(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843, stce_10061 &stce_10054) {
        stce_10064(stce_9843.stce_9950, stce_10054);
      }
      template<class stce_9941, class stce_9943 >
      static void stce_10064(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843, stce_10061 &stce_10054) {
        stce_10064(stce_9843.stce_9945, stce_10054);
      }
      template<class stce_9942, class stce_9943 >
      static void stce_10064(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843, stce_10061 &stce_10054) {
        stce_10064(stce_9843.stce_9947, stce_10054);
      }
    };
  }
}
       
       
namespace dco {
  namespace internal {
struct interpretation_settings {
  bool reset;
  bool stce_10065;
  interpretation_settings(bool stce_10066=false, bool stce_10067=false)
    : reset(stce_10066), stce_10065(stce_10067) {
  }
};
}
}
namespace dco {
    class tape_options {
    private:
      size_t stce_10068, stce_10069;
      bool stce_10070;
      std::string stce_10071;
      bool stce_10072;
      bool stce_10073;
    public:
      tape_options()
 : stce_10068(128*1024*1024),
          stce_10069(0),
   stce_10070(false),
   stce_10071("dco_tape_data.bin"),
   stce_10072(false),
          stce_10073(false) { }
      bool operator==(const tape_options& stce_9933) {
        return
          stce_10068 == stce_9933.stce_10068 &&
          stce_10069 == stce_9933.stce_10069 &&
          stce_10070 == stce_9933.stce_10070 &&
          stce_10072 == stce_9933.stce_10072 &&
          stce_10071 == stce_9933.stce_10071 &&
          stce_10073 == stce_9933.stce_10073;
      }
      bool operator!=(const tape_options& stce_9933) { return !operator==(stce_9933); }
      DCO_TAPE_INT chunk_size() const {
        if (helper::stce_9822) {
          throw dco::exception::create<std::runtime_error>("chunk_size removed; use chunk_size_in_byte() instead");
        }
        return -1;
      }
      size_t chunk_size_in_byte() const { return stce_10068; }
      size_t& chunk_size_in_byte() { return stce_10068; }
      void set_chunk_size_in_byte (mem_long_t size) {
        if (static_cast<double>(size) / static_cast<double>(std::numeric_limits<size_t>::max()) > 1.0) {
          DCO_LOG(logWARNING) << "set_chunk_size: requesting too much memory for your system (size_t too small). Resetting to std::numeric_limits<size_t>::max()";
          stce_10068 = std::numeric_limits<size_t>::max();
        } else {
          stce_10068 = static_cast<size_t>(size);
        }
      }
      void set_chunk_size_in_kbyte(double size) { set_chunk_size_in_byte(1024*static_cast<mem_long_t>(size)); }
      void set_chunk_size_in_mbyte(double size) { set_chunk_size_in_kbyte(1024.*size); }
      void set_chunk_size_in_gbyte(double size) { set_chunk_size_in_mbyte(1024.*size); }
      size_t blob_size_in_byte() const { return stce_10069; }
      size_t& blob_size_in_byte() { return stce_10069; }
      void set_blob_size_in_byte (mem_long_t size) {
        if (static_cast<double>(size) / static_cast<double>(std::numeric_limits<size_t>::max()) > 1.0) {
          DCO_LOG(logWARNING) << "set_blob_size: requesting too much memory for your system (size_t too small). Resetting to std::numeric_limits<size_t>::max()";
          stce_10069 = std::numeric_limits<size_t>::max();
        } else {
          stce_10069 = static_cast<size_t>(size);
        }
      }
      void set_blob_size_in_kbyte(double size) { set_blob_size_in_byte(1024*static_cast<mem_long_t>(size)); }
      void set_blob_size_in_mbyte(double size) { set_blob_size_in_kbyte(1024.*size); }
      void set_blob_size_in_gbyte(double size) { set_blob_size_in_mbyte(1024.*size); }
      bool &deallocation_on_reset() { return stce_10070; }
      const bool &deallocation_on_reset() const { return stce_10070; }
      std::string &filename() { return stce_10071; }
      const std::string &filename() const { return stce_10071; }
      bool& write_to_file() { return stce_10072; }
      bool write_to_file() const { return stce_10072; }
      bool& alloc_use_huge_pages() { return stce_10073; }
      bool alloc_use_huge_pages() const { return stce_10073; }
};
}
       
       
namespace dco {
template <typename stce_9782, typename stce_10074 = void> struct trait_size_of;
template<typename stce_9782> mem_long_t size_of(const stce_9782& stce_9843) {
  return trait_size_of<stce_9782>::get(stce_9843);
}
template<typename stce_9782> mem_long_t size_of(const stce_9782& stce_9843, int stce_10075) {
  return trait_size_of<stce_9782>::get(stce_9843, stce_10075);
}
template <typename stce_9782, typename stce_10074> struct trait_size_of {
  static mem_long_t get(const stce_9782&) {
    if (helper::stce_9822) {
      throw dco::exception::create<std::runtime_error>("You're trying to compute the size of the checkpoint, but a type ("
                                                       + type_name<stce_9782>()
                                                       + ") you checkpointed does not have a correct specialization:\n\n" +
                                                       "template <typename T> struct trait_size_of {\n" +
                                                       "   static mem_long_t get(const T&) { ... }\n" +
                                                       "};\n\n");
    }
    return 0;
  }
};
template<memory_model::TYPE MEMORY_MODEL, class stce_9788, memory_model::TYPE stce_9789>
  struct trait_size_of<dco::internal::tape<MEMORY_MODEL, stce_9788, stce_9789>*> {
  typedef dco::internal::tape<MEMORY_MODEL, stce_9788, stce_9789> tape_t;
  static mem_long_t get(const tape_t *tape,
                        const int stce_10075 = tape_t::stce_10076) {
  return tape->stce_10077(stce_10075);
  }
};
template<typename stce_9782> struct trait_size_of<std::vector<stce_9782> > {
  static mem_long_t get(const std::vector<stce_9782>& stce_9843) { return stce_9843.size() * sizeof(stce_9782) + sizeof(stce_9843); }
};
template<> struct trait_size_of<char> { static mem_long_t get(const char&) { return sizeof(char); } };
template<> struct trait_size_of<wchar_t> { static mem_long_t get(const wchar_t&) { return sizeof(wchar_t); } };
template<> struct trait_size_of<signed char> { static mem_long_t get(const signed char&) { return sizeof(signed char); } };
template<> struct trait_size_of<short> { static mem_long_t get(const short&) { return sizeof(short); } };
template<> struct trait_size_of<int> { static mem_long_t get(const int&) { return sizeof(int); } };
template<> struct trait_size_of<long> { static mem_long_t get(const long&) { return sizeof(long); } };
template<> struct trait_size_of<unsigned char> { static mem_long_t get(const unsigned char&) { return sizeof(unsigned char); } };
template<> struct trait_size_of<unsigned short> { static mem_long_t get(const unsigned short&) { return sizeof(unsigned short); } };
template<> struct trait_size_of<unsigned int> { static mem_long_t get(const unsigned int&) { return sizeof(unsigned int); } };
template<> struct trait_size_of<unsigned long> { static mem_long_t get(const unsigned long&) { return sizeof(unsigned long); } };
template<> struct trait_size_of<long long> { static mem_long_t get(const long long&) { return sizeof(long long); } };
template<> struct trait_size_of<unsigned long long> { static mem_long_t get(const unsigned long long&) { return sizeof(unsigned long long); } };
template<> struct trait_size_of<float> { static mem_long_t get(const float&) { return sizeof(float); } };
template<> struct trait_size_of<double> { static mem_long_t get(const double&) { return sizeof(double); } };
template<> struct trait_size_of<long double> { static mem_long_t get(const long double&) { return sizeof(long double); } };
template<> struct trait_size_of<bool> { static mem_long_t get(const bool&) { return sizeof(bool); } };
}
       
namespace dco {
namespace helper {
template <typename, typename stce_10021>
  struct stce_10078 {
  stce_10078(stce_10021*) {}
  template <typename stce_9782> void register_input(const stce_9782&) { }
  template <typename stce_9782> void register_output(const stce_9782&) { }
};
}
}
       
       
namespace dco {
template<class stce_9761>
  struct stce_10079 {
  virtual ~stce_10079(){};
  virtual stce_9761 stce_10080(const DCO_TAPE_INT stce_9808)=0;
  virtual void stce_10081(const DCO_TAPE_INT stce_9808, const stce_9761 &stce_10082)=0;
};
}
namespace dco {
namespace helper {
template<typename stce_10083>
  class stce_10084 {
  template<memory_model::TYPE, typename, memory_model::TYPE>
    friend class internal::tape;
  typedef typename stce_10083::stce_10085 stce_9764;
protected:
  stce_10083 *stce_10086;
public:
  virtual ~stce_10084() {
  }
  dco::stce_10079<stce_9764> *stce_10087;
public:
  void set_tape(stce_10083 *stce_9905) {
    if (stce_10086 != 0) { return; }
    else { stce_10086 = stce_9905; }
  }
  inline stce_10083 *get_tape() { return stce_10086; }
  stce_10084() : stce_10086(0) {}
  virtual void stce_10088(stce_9848 &stce_9849, DCO_TAPE_INT stce_9857) {
    (void)stce_9849;
    (void)stce_9857;
  }
  virtual mem_long_t size_in_byte() { return sizeof(stce_10083); }
};
}
}
       
namespace dco {
namespace helper {
template<class stce_10089, class stce_10083>
  class stce_10090 : public stce_10084<stce_10083> {
private:
  class stce_10091 {
public:
    virtual ~stce_10091() {};
    virtual mem_long_t size_in_byte() = 0;
    virtual void stce_10092(const std::string&, const std::string&) {};
  };
  template <typename stce_10093>
    class stce_10094 : public stce_10091 {
public:
    const stce_10093 _data;
    stce_10094(const stce_10093 &stce_10095) : _data(stce_10095) {}
    const stce_10093 &stce_10096() const { return _data; }
  };
  template <typename stce_10093>
    class stce_10097 : public stce_10094<stce_10093> {
public:
    stce_10097(stce_10093 data) : stce_10094<stce_10093>(data) {}
    virtual ~stce_10097() { }
    virtual mem_long_t size_in_byte() {
      return size_of(stce_10094<stce_10093>::_data);
    }
    virtual void stce_10092(const std::string& stce_10098, const std::string& stce_10099) {
      if (false && IS_DEBUG && stce_10099 != stce_10098) {
        std::stringstream stce_10063;
        stce_10063 << "read_data error: tried to read " << stce_10098 << " instead of " << stce_10099 << "!";
        throw dco::exception::create<std::runtime_error>(stce_10063.str());
      }
    }
  };
  template <typename stce_10093>
    class stce_10100 : public stce_10094<stce_10093*> {
public:
    const DCO_TAPE_INT stce_10101;
    stce_10100(const stce_10093 *data, DCO_TAPE_INT stce_10102) : stce_10094<stce_10093*>(new stce_10093[size_t(stce_10102)]), stce_10101(stce_10102) {
      for (DCO_TAPE_INT stce_9832 = 0; stce_9832 < stce_10102; ++stce_9832)
        this->_data[stce_9832] = data[stce_9832];
    }
    stce_10100(const stce_10093 *data, const DCO_TAPE_INT stce_10103, const DCO_TAPE_INT stce_10102) :
      stce_10094<stce_10093 *>(new stce_10093[stce_9752(stce_10102)]), stce_10101(stce_10102) {
      for (DCO_TAPE_INT stce_9832 = 0, stce_9808 = 0; stce_9832 < stce_10102; ++stce_9832, stce_9808 += stce_10103)
        this->_data[stce_9832] = data[stce_9808];
    }
    virtual ~stce_10100() {
      delete [] this->_data;
    }
    virtual mem_long_t size_in_byte() {
      return stce_10101 * static_cast<DCO_TAPE_INT>(sizeof(stce_10093));
    }
    virtual void stce_10092(const std::string& stce_10098, const std::string& stce_10099) {
      (void) stce_10098; (void) stce_10099;
      if (false && IS_DEBUG && stce_10099+"*" != stce_10098) {
        std::stringstream stce_10063;
        stce_10063 << "read_data error: tried to read " << stce_10098 << " instead of " << stce_10099 << "!";
        throw dco::exception::create<std::runtime_error>(stce_10063.str());
      }
    }
  };
  DCO_TAPE_INT stce_10104;
  std::vector<stce_10091 *> stce_10105;
  std::vector<std::string> stce_10106;
protected:
  virtual ~stce_10090() {
    for (size_t stce_9832 = 0; stce_9832 < stce_10105.size(); stce_9832++)
      delete stce_10105[stce_9832];
    stce_10105.clear();
  }
public:
  stce_10090(): stce_10084<stce_10083>(), stce_10104(0) {}
  virtual mem_long_t size_in_byte() {
    mem_long_t stce_9920 = stce_10084<stce_10083>::size_in_byte();
    for (size_t stce_9832 = 0; stce_9832 < stce_10105.size(); stce_9832++) {
      stce_9920 += stce_10105[stce_9832]->size_in_byte();
    }
    stce_9920 += stce_10105.size() * sizeof(stce_10091*);
    stce_9920 += stce_10106.size() * sizeof(std::string);
    return stce_9920;
  }
  template<typename stce_10093>
    void write_data(const stce_10093 &stce_10107) {
    stce_10105.push_back(new stce_10097<stce_10093>(stce_10107));
    if (IS_DEBUG)
      stce_10106.push_back(type_name<stce_10093>());
  }
  template<typename stce_10093>
    void write_data(const stce_10093 *const stce_10107, const DCO_TAPE_INT stce_10102) {
    stce_10105.push_back(new stce_10100<stce_10093>(stce_10107, stce_10102));
    if (IS_DEBUG)
      stce_10106.push_back(type_name<stce_10093>());
  }
  template<typename stce_10093>
    void write_data(const stce_10093 *const &stce_10107, const DCO_TAPE_INT stce_10103, const DCO_TAPE_INT stce_10102) {
    stce_10105.push_back(new stce_10100<stce_10093>(stce_10107, stce_10103, stce_10102));
    if (IS_DEBUG) {
      stce_10106.push_back(type_name<stce_10093>());
    }
  }
  template<typename stce_10093>
    const stce_10093 &read_data() {
    const stce_10093 &stce_10107 = static_cast<stce_10094<stce_10093>* >(stce_10105[stce_10104])->stce_10096();
    if (IS_DEBUG)
      stce_10105[stce_10104]->stce_10092(type_name<stce_10093>(), stce_10106.at(stce_10104));
    ++stce_10104;
    if (stce_10104 == static_cast<DCO_TAPE_INT>(stce_10105.size())) stce_10104 = 0;
    return stce_10107;
  }
};
}
}
       
namespace dco {
namespace helper {
template<class stce_10089, class stce_10083>
  class stce_10108 :
           public stce_10090<stce_10089, stce_10083>,
           public stce_10078<stce_10083, stce_10108<stce_10089, stce_10083> > {
  typedef stce_10078<stce_10083, stce_10108> stce_10109;
protected:
  std::vector<DCO_TAPE_INT> stce_9855;
  std::vector<DCO_TAPE_INT> stce_9856;
  DCO_TAPE_INT stce_10110;
  DCO_TAPE_INT stce_10111;
public:
  stce_10083* &tape() { return this->stce_10086; }
  virtual mem_long_t size_in_byte() {
    mem_long_t stce_9920 = stce_10090<stce_10089, stce_10083>::size_in_byte();
    stce_9920 += stce_9855.size() * sizeof(DCO_TAPE_INT);
    stce_9920 += stce_9856.size() * sizeof(DCO_TAPE_INT);
    stce_9920 += 2 * sizeof(DCO_TAPE_INT);
    return stce_9920;
  }
  size_t stce_10112() {
    return stce_9855.size();
  }
  size_t stce_10113() {
    return stce_9856.size();
  }
  void stce_10088(stce_9848 &stce_9849, DCO_TAPE_INT stce_9857) {
    stce_9849.stce_9854(stce_9855, stce_9856, stce_9857 - static_cast<DCO_TAPE_INT>(stce_9856.size()) - 1);
  }
public:
  void stce_10114(const stce_10089 &stce_9843) {
    if ((stce_9843.tape() != 0) && (this->stce_10086 != stce_9843.tape()))
      throw dco::exception::create<std::runtime_error>("impossible binding tape - wrong tape in variable!");
  }
protected:
  ~stce_10108() {
    if (IS_DEBUG && stce_10110 != 0) {
      DCO_LOG(dco::logERROR) << "not all input adjoints incremented in external adjoint object.";
    }
    if (IS_DEBUG && stce_10111 != 0) {
      DCO_LOG(dco::logERROR) << "not all output adjoints read in external adjoint object.";
    }
  }
public:
  stce_10108(const std::pair<DCO_TAPE_INT, DCO_TAPE_INT> &stce_10115):
    stce_10090<stce_10089, stce_10083>(),
    stce_10109(this),
    stce_10110(0),
    stce_10111(0) {
    stce_9855.reserve(stce_10115.first);
    stce_9856.reserve(stce_10115.second);
  }
  stce_10108(): stce_10090<stce_10089, stce_10083>(),
    stce_10109(this),
    stce_10110(0), stce_10111(0) { }
  typename stce_10089::VALUE_TYPE register_input(const stce_10089 &stce_9843) {
    stce_10114(stce_9843);
    const DCO_TAPE_INT stce_9808 = stce_9843._tape_index();
    stce_9855.push_back(stce_9808);
    stce_10109::register_input(stce_9843);
    return stce_9843._value();
  }
  void register_input(const stce_10089 *const stce_9843, typename stce_10089::VALUE_TYPE *stce_10116, const DCO_TAPE_INT stce_10102) {
    stce_9855.reserve(stce_9855.size() + stce_10102);
    for (int stce_9832 = 0; stce_9832 < stce_10102; ++stce_9832) {
      stce_10114(stce_9843[stce_9832]);
      const DCO_TAPE_INT stce_9808 = stce_9843[stce_9832]._tape_index();
      stce_9855.push_back(stce_9808);
      stce_10116[stce_9832] = stce_9843[stce_9832]._value();
      stce_10109::register_input(stce_9843[stce_9832]);
    }
  }
  void register_input(const std::vector<stce_10089> &stce_9843, std::vector<typename stce_10089::VALUE_TYPE> &stce_10116) {
    assert(stce_9843.size() == stce_10116.size());
    register_input(&(stce_9843[0]), &(stce_10116[0]), stce_9843.size());
  }
  std::vector<typename stce_10089::VALUE_TYPE> register_input(const std::vector<stce_10089> &stce_9843) {
    std::vector<typename stce_10089::VALUE_TYPE> stce_10116(stce_9843.size());
    register_input(stce_9843, stce_10116);
    return stce_10116;
  }
  void register_output(stce_10089 *stce_10117, const size_t stce_10102) {
    if (this->stce_10086 == NULL) {
      throw dco::exception::create<std::runtime_error>("impossible binding output - no tape available");
    } else {
      register_output(NULL, stce_10117, stce_10102);
    }
  }
  void register_output(const typename stce_10089::VALUE_TYPE *const stce_10118, stce_10089 *stce_10117, const size_t stce_10102) {
    if (this->stce_10086 == NULL) {
      throw dco::exception::create<std::runtime_error>("impossible binding output - no tape available");
    } else {
      stce_9856.reserve(stce_9856.size() + stce_10102);
      for (size_t stce_9832 = 0; stce_9832 < stce_10102; ++stce_9832) {
        if (stce_10118) {
          stce_10117[stce_9832] = stce_10118[stce_9832];
        }
        this->stce_10086->register_variable(stce_10117[stce_9832]);
        stce_9856.push_back(stce_10117[stce_9832]._tape_index());
        stce_10109::register_output(stce_10117[stce_9832]);
      }
    }
  }
  void register_output(const std::vector<typename stce_10089::VALUE_TYPE> &stce_10118, std::vector<stce_10089> &stce_10117) {
    assert(stce_10118.size() == stce_10117.size());
    register_output(&(stce_10118[0]), &(stce_10117[0]), stce_10118.size());
  }
  std::vector<stce_10089> register_output(const std::vector<typename stce_10089::VALUE_TYPE> &stce_10118) {
    std::vector<stce_10089> stce_10117(stce_10118.size());
    register_output(stce_10118, stce_10117);
    return stce_10117;
  }
  void register_output(std::vector<stce_10089> &stce_10117) {
    register_output(&(stce_10117[0]), stce_10117.size());
  }
  stce_10089 register_output(const typename stce_10089::VALUE_TYPE &stce_10119, stce_10083 *tape = NULL) {
    stce_10089 stce_10120;
    if (tape != NULL) {
      if (this->stce_10086 != NULL && this->stce_10086 != tape) {
        throw dco::exception::create<std::runtime_error>("impossible binding output in external function (register_output) - tape of inputs and outputs differ!");
      }
      this->stce_10086 = tape;
    }
    if (this->stce_10086 != NULL) {
      stce_10120 = stce_10119;
      this->stce_10086->register_variable(stce_10120);
    } else
      throw dco::exception::create<std::runtime_error>("impossible binding output in external function - no tape available");
    stce_9856.push_back(stce_10120._tape_index());
    stce_10109::register_output(stce_10120);
    return stce_10120;
  }
  inline typename dco::mode<stce_10089>::stce_10020 get_output_adjoint() {
    DCO_TAPE_INT stce_9808 = stce_10111;
    stce_10111++;
    if (stce_10111 == static_cast<DCO_TAPE_INT>(stce_9856.size()))
      stce_10111 = 0;
    typename dco::mode<stce_10089>::stce_10020 back = this->stce_10087->stce_10080(stce_9856[static_cast<size_t>(stce_9808)]);
    return back;
  }
  void get_output_adjoint(typename stce_10089::VALUE_TYPE *stce_9809, const size_t stce_10102) {
    DCO_TAPE_INT stce_9808 = stce_10111;
    for (size_t stce_9832 = 0; stce_9832 < stce_10102; ++stce_9832) {
      stce_9809[stce_9832] = get_output_adjoint();
      stce_9808++;
    }
    if (stce_10111 == static_cast<DCO_TAPE_INT>(stce_9856.size())) {
      stce_10111 = 0;
    }
  }
  void get_output_adjoint(std::vector<typename stce_10089::VALUE_TYPE> &stce_9809) {
    assert(stce_9809.size());
    get_output_adjoint(&(stce_9809[0]), stce_9809.size());
  }
  void increment_input_adjoint(const typename stce_10089::VALUE_TYPE *const stce_10082, const DCO_TAPE_INT stce_10102) {
    for (DCO_TAPE_INT stce_9832 = 0; stce_9832 < stce_10102; ++stce_9832) {
      this->stce_10087->stce_10081(stce_9855[stce_10110 + stce_9832], stce_10082[stce_9832]);
    }
    stce_10110 += stce_10102;
    if (stce_10110 == static_cast<DCO_TAPE_INT>(stce_9855.size()))
      stce_10110 = 0;
  }
  void increment_input_adjoint(const std::vector<typename stce_10089::VALUE_TYPE> &stce_10082) {
    assert(stce_10082.size() != 0);
    increment_input_adjoint(&(stce_10082[0]), static_cast<DCO_TAPE_INT>(stce_10082.size()));
  }
  bool stce_10121() {
    if (stce_10110 == 0) return true;
    else return false;
  }
  bool stce_10122() {
    if (stce_10111 == 0) return true;
    else return false;
  }
  void increment_input_adjoint(const typename dco::mode<stce_10089>::stce_10020 &stce_10082) {
    DCO_TAPE_INT stce_9808 = stce_10110;
    stce_10110++;
    if (stce_10110 == static_cast<DCO_TAPE_INT>(stce_9855.size()))
      stce_10110 = 0;
    this->stce_10087->stce_10081(stce_9855[static_cast<size_t>(stce_9808)], stce_10082);
  }
};
}
}
       
namespace dco {
namespace ACM {
template<class stce_10083, class stce_9982>
  struct baseclass : dco::helper::stce_10084<stce_10083> {
private:
  std::vector<stce_9982*> stce_10123;
  std::vector<stce_9982*> stce_10124;
  std::vector<DCO_TAPE_INT> stce_10125;
  std::vector<DCO_TAPE_INT> stce_10126;
  std::vector<void*> _data;
  int stce_10127;
public:
  typedef stce_9982 active_t;
  typedef stce_10083 tape_t;
  baseclass(tape_t *tape) : stce_10127(0) {
    this->set_tape(tape);
  }
  size_t register_input(stce_9982 &stce_9843) {
    size_t stce_9840 = stce_10123.size();
    stce_10123.push_back(&stce_9843);
    stce_10126.push_back(stce_9843._tape_index());
    return stce_9840;
  }
  size_t register_input(std::vector<stce_9982> &stce_10128) {
    size_t stce_9840 = stce_10123.size();
    for(int stce_9832=0;stce_9832<stce_10128.size();++stce_9832) {
      const stce_9982 &stce_10129 = stce_10128[stce_9832];
      this->register_input(stce_10129);
    }
    return stce_9840;
  }
  size_t register_output(stce_9982 &stce_9843) {
    if(stce_9843._tape_index()==0) {
      this->get_tape()->register_variable(stce_9843);
    }
    size_t stce_9840 = stce_10124.size();
    stce_10124.push_back(&stce_9843);
    stce_10125.push_back(stce_9843._tape_index());
    return stce_9840;
  }
  int register_output(std::vector<stce_9982> &stce_10128) {
    int stce_9840 = stce_10124.size();
    for(int stce_9832=0;stce_9832<stce_10128.size();++stce_9832) {
      if(stce_10128[stce_9832]._tape_index()==0) {
        this->get_tape()->register_variable(stce_10128[stce_9832]);
      }
      stce_10125.push_back(stce_10128[stce_9832]._tape_index());
      stce_10124.push_back(&stce_10128[stce_9832]);
    }
    return stce_9840;
  }
  template<typename stce_9782> int write_data(const stce_9782 &data) {
    stce_9782* stce_10130 = new stce_9782(data);
    int stce_9840 = _data.size();
    _data.push_back(stce_10130);
    return stce_9840;
  }
  template<typename stce_9782> const stce_9782& read_data(const int stce_10131=-1) {
    int stce_9808=0;
    if(stce_10131==-1) {
      stce_9808 = stce_10127%_data.size(); stce_10127++;
    }
    void* stce_10132 = _data[stce_9808];
    return *static_cast<stce_9782*>(stce_10132);
  }
  size_t stce_10133() {
    return stce_10123.size();
  }
  size_t stce_10134() {
    return stce_10123.size();
  }
  typename stce_9982::VALUE_TYPE& input_value(const int stce_9808) {
    return stce_10123[stce_9808]->_value();
  }
  typename stce_9982::VALUE_TYPE& output_value(const int stce_9808) {
    return stce_10124[stce_9808]->_value();
  }
  typename dco::mode<stce_9982>::stce_10020& input_adjoint(const int stce_9808) {
    return this->get_tape()->_adjoint( stce_10126[stce_9808] );
  }
  typename dco::mode<stce_9982>::stce_10020& output_adjoint(const int stce_9808) {
    return this->get_tape()->_adjoint( stce_10125[stce_9808] );
  }
  virtual void evaluate_primal() {
    std::cout << "WARNING: you need to implement evaluate_primal!!!!\n";
  }
  virtual void evaluate_augmented_primal() {
    std::cout << "WARNING: you need to implement evaluate_augmented_primal!!!!\n";
  }
  virtual void evaluate_adjoint() {
    std::cout << "WARNING: you need to implement evaluate_adjoint!!!!\n";
  }
};
template<class stce_10083, class stce_9982>
  static void stce_10135( baseclass<stce_10083, stce_9982> *stce_10136 ) {
  stce_10136->evaluate_adjoint();
}
}
}
namespace dco {
namespace stce_10137 {
    template<class stce_9991>
    struct stce_10138 {
        std::vector<stce_9991*> stce_10123;
        std::vector<stce_9991*> stce_10124;
        std::vector<void*> _data;
        int stce_10127;
        stce_10138() : stce_10127(0) {}
        virtual ~stce_10138() {}
        template<typename stce_9782> int write_data(const stce_9782 &data) {
            stce_9782* stce_10130 = new stce_9782(data);
            int stce_9840 = _data.size();
            _data.push_back(stce_10130);
            return stce_9840;
        }
        template<typename stce_9782> const stce_9782& read_data(const int stce_10131=-1) {
            int stce_9808=0;
            if(stce_10131==-1) {
                stce_9808 = stce_10127%_data.size(); stce_10127++;
            }
            void* stce_10132 = _data[stce_9808];
            return *static_cast<stce_9782*>(stce_10132);
        }
        size_t stce_10133() {
            return stce_10123.size();
        }
        size_t stce_10134() {
            return stce_10124.size();
        }
        stce_9991& input_value(const int stce_9808) {
            return *stce_10123[stce_9808];
        }
        stce_9991& output_value(const int stce_9808) {
            return *stce_10124[stce_9808];
        }
        virtual void evaluate_primal() {
          std::cout << "WARNING: you need to implement evaluate_primal!!!!\n";
        }
    };
    template<class stce_9991>
    struct stce_10139 : stce_10138<stce_9991> {
        int register_input(stce_9991 &stce_9843) {
            int stce_9840 = this->stce_10123.size();
            this->stce_10123.push_back(&stce_9843);
            return stce_9840;
        }
        int register_input(std::vector<stce_9991> &stce_10128) {
            int stce_9840 = this->stce_10123.size();
            for(int stce_9832=0;stce_9832<stce_10128.size();++stce_9832) {
                const stce_9991 &stce_10129 = stce_10128[stce_9832];
                this->register_input(stce_10129);
            }
            return stce_9840;
        }
        int register_output(stce_9991 &stce_9843) {
            int stce_9840 = this->stce_10124.size();
            this->stce_10124.push_back(&stce_9843);
            return stce_9840;
        }
        int register_output(std::vector<stce_9991> &stce_10128) {
            int stce_9840 = this->stce_10124.size();
            for(int stce_9832=0;stce_9832<stce_10128.size();++stce_9832) {
                register_output(stce_10128[stce_9832]);
            }
            return stce_9840;
        }
    };
    template<class stce_10083, class stce_9982>
    struct stce_10140 : public dco::helper::stce_10084<stce_10083>, public stce_10138<stce_9982> {
    private:
        std::vector<DCO_TAPE_INT> stce_10125;
        std::vector<DCO_TAPE_INT> stce_10126;
    public:
        typedef stce_9982 active_t;
        typedef stce_10083 tape_t;
        stce_10140(tape_t *tape) {
            this->set_tape(tape);
        }
        stce_9982& stce_10141(const int stce_9808) {
            return *(this->stce_10123[stce_9808]);
        }
        stce_9982& stce_10142(const int stce_9808) {
            return *(this->stce_10124[stce_9808]);
        }
        typename stce_9982::VALUE_TYPE& input_value(const int stce_9808) {
            return this->stce_10123[stce_9808]->_value();
        }
        typename stce_9982::VALUE_TYPE& output_value(const int stce_9808) {
            return this->stce_10124[stce_9808]->_value();
        }
        typename dco::mode<stce_9982>::stce_10020& input_adjoint(const int stce_9808) {
            return this->get_tape()->_adjoint( stce_10126[stce_9808] );
        }
        typename dco::mode<stce_9982>::stce_10020& output_adjoint(const int stce_9808) {
            return this->get_tape()->_adjoint( stce_10125[stce_9808] );
        }
        int register_input(stce_9982 &stce_9843) {
            int stce_9840 = this->stce_10123.size();
            this->stce_10123.push_back(&stce_9843);
            this->stce_10126.push_back(stce_9843._tape_index());
            return stce_9840;
        }
        int register_input(std::vector<stce_9982> &stce_10128) {
            int stce_9840 = this->stce_10123.size();
            for(int stce_9832=0;stce_9832<stce_10128.size();++stce_9832) {
                const stce_9982 &stce_10129 = stce_10128[stce_9832];
                this->register_input(stce_10129);
            }
            return stce_9840;
        }
        int register_output(stce_9982 &stce_9843) {
            if(stce_9843._tape_index()==0) {
                this->get_tape()->register_variable(stce_9843);
            }
            int stce_9840 = this->stce_10124.size();
            this->stce_10124.push_back(&stce_9843);
            this->stce_10125.push_back(stce_9843._tape_index());
            return stce_9840;
        }
        int register_output(std::vector<stce_9982> &stce_10128) {
            int stce_9840 = this->stce_10124.size();
            for(int stce_9832=0;stce_9832<stce_10128.size();++stce_9832) {
                if(stce_10128[stce_9832]._tape_index()==0) {
                    this->get_tape()->register_variable(stce_10128[stce_9832]);
                }
                this->stce_10125.push_back(stce_10128[stce_9832]._tape_index());
                this->stce_10124.push_back(&stce_10128[stce_9832]);
            }
            return stce_9840;
        }
        virtual void evaluate_augmented_primal() {
            std::cout << "WARNING: you need to implement evaluate_augmented_primal!!!!\n";
        }
        virtual void evaluate_adjoint() {
            std::cout << "WARNING: you need to implement evaluate_adjoint!!!!\n";
        }
    };
    template<class stce_10083, class stce_9982>
    static void stce_10135( stce_10140<stce_10083, stce_9982> *stce_10136 ) {
        stce_10136->evaluate_adjoint();
    }
}
}
       
       
       
namespace dco {
namespace internal {
template <typename stce_10143>
  class stce_10144 {
protected:
  stce_10143& stce_10145;
public:
  stce_10144(stce_10143& stce_10146) : stce_10145(stce_10146) { }
  stce_10143& stce_10146() { return stce_10145; }
};
}
}
namespace dco {
namespace internal {
template<memory_model::TYPE memory_model,
         typename stce_10147,
         memory_model::TYPE stce_10148=memory_model::stce_9766> class stce_10149 {};
}
}
namespace dco {
namespace internal {
template<class stce_10147, memory_model::TYPE stce_10150>
  class stce_10149<memory_model::BLOB_TAPE_SPLINT, stce_10147, stce_10150>
  : public stce_10144<stce_10147>, object_logger
{
  typedef stce_10144<stce_10147> stce_10151;
  typedef typename stce_9790<stce_10147>::type stce_10152;
public:
  struct entry {
    stce_10152* stce_10153;
    DCO_TAPE_INT* stce_10154;
    DCO_ALWAYS_INLINE entry(stce_10152* stce_10155,
                            DCO_TAPE_INT* stce_10156)
      : stce_10153(stce_10155),
        stce_10154(stce_10156) {};
    DCO_ALWAYS_INLINE DCO_TAPE_INT arg() {
      return (*stce_10154 < 0) ? std::abs(stce_9752(*stce_10154))-1 : stce_9752(*stce_10154);
    }
    DCO_ALWAYS_INLINE stce_10152& pval() { return *stce_10153; }
    DCO_ALWAYS_INLINE stce_10152 pval() const { return *stce_10153; }
    DCO_ALWAYS_INLINE bool stce_10157() { return *stce_10154 < 0; }
    DCO_ALWAYS_INLINE static bool stce_10157(DCO_TAPE_INT* stce_10158) { return *stce_10158 < 0; }
    DCO_ALWAYS_INLINE static DCO_TAPE_INT stce_10159(DCO_TAPE_INT stce_10054) { return -stce_10054-1; }
    static size_t size_of() { return sizeof(stce_10152) + sizeof(DCO_TAPE_INT); }
  };
  class iterator {
public:
    DCO_TAPE_INT stce_10160;
    entry stce_10161;
    template<class stce_10162>
      void stce_10163(stce_10162 &stce_10164) const {
      assert(sizeof(iterator) <= sizeof(stce_10162)) ;
      stce_10165(&stce_10164, this, sizeof(iterator));
    }
    template<class stce_10162>
      void stce_10166(const stce_10162 &stce_10164) {
      assert(sizeof(iterator) <= sizeof(stce_10162)) ;
      stce_10165(static_cast<void*>(this), &stce_10164, sizeof(iterator));
    }
    DCO_ALWAYS_INLINE const entry* operator->() const { return &stce_10161; }
    DCO_ALWAYS_INLINE entry* operator->() { return &stce_10161; }
    DCO_ALWAYS_INLINE void stce_10167() {
      stce_10161.stce_10153++;
      stce_10161.stce_10154--;
      if (stce_9754 && static_cast<void*>(stce_10161.stce_10153) >= static_cast<void*>(stce_10161.stce_10154)) {
        throw dco::exception::create<std::runtime_error>("Blob tape (splint) container overflow. Allocate bigger blob or use chunk tape container instead.");
      }
    }
    DCO_ALWAYS_INLINE DCO_TAPE_INT stce_10168() { return ++stce_10160; }
    DCO_ALWAYS_INLINE explicit iterator() : stce_10160(-1), stce_10161(NULL, NULL) { }
    DCO_ALWAYS_INLINE explicit iterator(stce_10152* stce_10155, DCO_TAPE_INT* stce_10156, DCO_TAPE_INT index = DCO_TAPE_INT(-1)) :
      stce_10160(index), stce_10161(stce_10155, stce_10156) { }
    DCO_ALWAYS_INLINE DCO_TAPE_INT index() const { return stce_10160; }
    DCO_ALWAYS_INLINE DCO_TAPE_INT& index() { return stce_10160; }
    bool stce_10157() const {
      return entry::stce_10157(stce_10161.stce_10154);
    }
    DCO_ALWAYS_INLINE void operator--() {
      if(entry::stce_10157(stce_10161.stce_10154))
        stce_10160--;
      stce_10161.stce_10153--;
      stce_10161.stce_10154++;
    }
    DCO_ALWAYS_INLINE iterator operator--(int) {
      iterator stce_10063(*this);
      operator--();
      return stce_10063;
    }
    DCO_ALWAYS_INLINE bool operator==(const iterator& stce_9933) const {
      return stce_10161.stce_10153 == stce_9933.stce_10161.stce_10153;
    }
    DCO_ALWAYS_INLINE bool operator!=(const iterator& stce_9933) const {
      return !operator==(stce_9933);
    }
    DCO_ALWAYS_INLINE bool operator<(const iterator& stce_9933) const {
      return stce_10161.stce_10153 < stce_9933.stce_10161.stce_10153;
    }
    DCO_ALWAYS_INLINE bool operator>=(const iterator& stce_9933) const {
      return !operator<(stce_9933);
    }
    DCO_ALWAYS_INLINE bool operator>(const iterator& stce_9933) const {
      return stce_10161.stce_10153 > stce_9933.stce_10161.stce_10153;
    }
    DCO_ALWAYS_INLINE bool operator<=(const iterator& stce_9933) const {
      return !operator>(stce_9933);
    }
  };
  class stce_9924
  {
    iterator& stce_10169;
    DCO_TAPE_INT stce_10160;
public:
    stce_9924(iterator& current) :
      stce_10169(current), stce_10160(0){
    }
    iterator& current() { return stce_10169; }
    stce_9924& operator=(stce_9924 stce_9933) {
      swap(stce_9933);
      return *this;
    }
    DCO_TAPE_INT index() const {
      return stce_10160;
    }
    void swap(stce_9924& stce_9933)
    {
      std::swap(stce_10169, stce_9933.stce_10169);
      std::swap(stce_10160, stce_9933.stce_10160);
    }
    DCO_ALWAYS_INLINE_IF_AGRESSIVE void insert(const DCO_TAPE_INT& stce_10054, const stce_10152& pval) {
      if (IS_DEBUG && stce_10169.index() < stce_10054) {
        throw dco::exception::create<std::runtime_error>("DCO_DEBUG: tape index of read variable bigger than current progvarcounter.");
      }
      if(stce_10160 == 0) {
        stce_10169.stce_10167();
        new (stce_10169.stce_10161.stce_10153) stce_10152(pval);
        new (stce_10169.stce_10161.stce_10154) DCO_TAPE_INT(entry::stce_10159(stce_10054));
        stce_10160 = stce_10169.stce_10168();
        DCO_LOG(logDEBUG2) << "inserting first entry, index now = " << stce_10160 << "; pval = " << pval << ", target = " << stce_10054;
        return;
      }
      if(stce_9757 && stce_10054 == *(stce_10169.stce_10161.stce_10154))
        {
          *(stce_10169.stce_10161.stce_10153) += pval;
          return;
        }
      stce_10169.stce_10167();
      new (stce_10169.stce_10161.stce_10153) stce_10152(pval);
      new (stce_10169.stce_10161.stce_10154) DCO_TAPE_INT(stce_10054);
    }
private:
    stce_9924(const stce_9924& stce_9933);
  };
protected:
  void init(tape_options const& stce_10170, DCO_TAPE_INT stce_10171 = 0) {
    size_t stce_10172;
    if(stce_10170.write_to_file())
      stce_10172 = stce_10170.blob_size_in_byte();
    else
      stce_10172 = get_allocation_size(stce_10170.blob_size_in_byte());
    size_t stce_10173 = stce_10172 / entry::size_of();
    DCO_LOG(logDEBUG1) << "BLOB SPLINT container: user requested allocation size = " << stce_10172
                       << "; size of element = " << entry::size_of()
                       << "; max elements in memory = " << stce_10173;
    if (static_cast<double>(stce_10173) / static_cast<double>(std::numeric_limits<DCO_INTEGRAL_TAPE_INT>::max()) > 1.0) {
      stce_10173 = std::numeric_limits<DCO_TAPE_INT>::max();
      DCO_LOG(logDEBUG1) << "BLOB SPLINT container: too many elements requested for currently used DCO_TAPE_INT. Setting to: " << stce_10173 * 1.e-6 << "M elements";
    }
    DCO_TAPE_INT stce_10174 = static_cast<DCO_TAPE_INT>(stce_10173);
    size_t stce_10175 = stce_10174 * entry::size_of();
    DCO_LOG(logDEBUG) << "BLOB SPLINT container: trying allocation for = " << stce_10175 / 1024 << "K ~=~ " << stce_9752(stce_10174) * 1.e-6 << "M elements";
    stce_10145 = stce_1000001(&stce_10175, stce_10170.alloc_use_huge_pages());
    if (!stce_10145) {
      std::string stce_10176 = "";
      if (stce_10170.alloc_use_huge_pages()) {
        stce_10176 = " You are using huge pages; have you prepared your system correspondingly (see documentation as well)?";
      }
      throw dco::exception::create<std::runtime_error>
        ("Could not allocate memory. Use environment variables (DCO_MEM_RATIO, DCO_MAX_ALLOCATION) and see documentation." + stce_10176);
    }
    stce_10174 = static_cast<DCO_TAPE_INT>(stce_10175 / entry::size_of());
    DCO_LOG(logINFO) << "BLOB SPLINT container: actually allocated size = " << stce_10175 / 1024 << "K ~=~ " << stce_9752(stce_10174) * 1.e-6 << "M elements";
    stce_10177 = static_cast<stce_10152*>(stce_10145);
    stce_10178 = static_cast<DCO_TAPE_INT*>(stce_10145) + stce_10175 / sizeof(DCO_TAPE_INT) - 1;
    stce_10169 = iterator(stce_10177, stce_10178, stce_10171);
    stce_10179 = iterator(stce_10177, stce_10178, stce_10171);
    stce_10180 = iterator(stce_10177 + stce_9752(stce_10174) - 1, stce_10178 - stce_9752(stce_10174) + 1, std::numeric_limits<DCO_TAPE_INT>::max());
    new (stce_10169.stce_10161.stce_10153) stce_10152(0.0);
    new (stce_10169.stce_10161.stce_10154) DCO_TAPE_INT(entry::stce_10159(0));
  }
public:
  explicit stce_10149(tape_options const& stce_10170, stce_10147 &tape) :
    stce_10151(tape),
    object_logger("BLOB SPLINT container"),
    stce_10145(0), stce_10177(0), stce_10178(0),
    stce_10169(), stce_10179(), stce_10180() {
    assert(stce_10150!=memory_model::stce_9770);
    init(stce_10170);
  }
  explicit stce_10149(iterator& stce_9817, stce_10147 &tape) :
    stce_10151(tape),
    object_logger("BLOB SPLINT container"),
    stce_10145(0), stce_10177(0), stce_10178(0),
    stce_10169(stce_9817), stce_10179(), stce_10180() {
    assert(stce_10150==memory_model::stce_9770);
  }
  ~stce_10149() {
    dco::dealloc(stce_10145);
  }
  DCO_ALWAYS_INLINE iterator& stce_9928() { return stce_10169; }
  iterator start() { return stce_10179; }
  iterator start() const { return stce_10179; }
  iterator rend() { return start(); }
  iterator current() const { return stce_10169; }
  iterator current() { return stce_10169; }
  iterator rbegin() const { return current(); }
  iterator capacity() const { return stce_10180; }
  bool empty() const { return stce_10169 == start(); }
  DCO_TAPE_INT size(const iterator& stce_9862, const iterator& stce_9863) const {
    if(stce_9862 < stce_9863)
      throw dco::exception::create<std::runtime_error>("size from < to");
    return static_cast<DCO_TAPE_INT>(stce_9862->stce_10153 - stce_9863->stce_10153);
  }
  DCO_TAPE_INT size(iterator stce_9862) const {
    return size(stce_9862, start());
  }
  DCO_TAPE_INT size() const {
    return size(current(), start());
  }
  DCO_TAPE_INT stce_10181() const {
    return size(stce_10180, start());
  }
  size_t size_in_byte() const { return size(current(), start())*entry::size_of(); }
  size_t stce_10182() const { return size(stce_10180 , start())*entry::size_of(); }
  void erase(const iterator& stce_10158) {
    if(!stce_10177)
      return;
    if(stce_10158 > stce_10180)
      throw dco::exception::create<std::runtime_error>("pos behind stack!");
    if(stce_10158 > stce_10169)
      throw dco::exception::create<std::runtime_error>("reset to position out of tape!");
    stce_10169 = stce_10158;
  }
  void stce_10183(const iterator&) { }
  void stce_10184 ( ) { }
  DCO_TAPE_INT stce_10185() const { return stce_10169.index()+1; }
private:
  stce_10149(const stce_10149& stce_9933);
  void *stce_10145;
  stce_10152 * stce_10177;
  DCO_TAPE_INT * stce_10178;
protected:
  typename stce_9821
    <iterator, stce_10150==memory_model::stce_9770>::type stce_10169;
  iterator stce_10179;
  iterator stce_10180;
};
}
}
       
       
namespace dco {
namespace internal {
template<typename stce_10147, memory_model::TYPE stce_10150>
  class stce_10149<memory_model::BLOB_TAPE, stce_10147, stce_10150> :
    public stce_10144<stce_10147>, public object_logger {
  typedef stce_10144<stce_10147> stce_10151;
  typedef stce_10147 tape_t;
  typedef typename stce_9790<stce_10147>::type stce_10152;
public:
  enum stce_10186 {
    stce_10187,
    stce_10188
  };
      struct entry {
        DCO_TAPE_INT stce_9950;
        stce_10152 stce_9926;
      public:
        static size_t size_of() {
          return sizeof(entry);
        }
        entry(DCO_TAPE_INT arg = 0, stce_10152 const& pval = 0.) :
          stce_9950(arg), stce_9926(pval) {
        }
        DCO_STRONG_INLINE DCO_TAPE_INT arg() const {
          return (stce_9950 < 0)?std::abs(stce_9752(stce_9950)) -1:stce_9752(stce_9950);
        }
        DCO_STRONG_INLINE static DCO_TAPE_INT stce_10159(const DCO_TAPE_INT& arg) {
          return -arg-1;
        }
        DCO_STRONG_INLINE bool stce_10157() const {
          return stce_9950 < 0;
        }
        DCO_STRONG_INLINE void stce_10189(const stce_10152& stce_10103) {
          stce_9926 += stce_10103;
        }
        DCO_STRONG_INLINE stce_10152 pval() const
        {
          return stce_9926;
        }
      }
      ;
      class iterator {
        entry* stce_10190;
        DCO_TAPE_INT stce_10160;
        entry const* stce_10191;
      public:
        template<class stce_10162>
          void stce_10163(stce_10162 &stce_10164) const {
          assert(sizeof(iterator) <= sizeof(stce_10162)) ;
          stce_10165(&stce_10164, this, sizeof(iterator));
        }
        template<class stce_10162>
          void stce_10166(const stce_10162 &stce_10164) {
          assert(sizeof(iterator) <= sizeof(stce_10162)) ;
          stce_10165(static_cast<void*>(this), &stce_10164, sizeof(iterator));
        }
        void stce_10192(DCO_TAPE_INT stce_9808)
        {
          stce_10160 = stce_9808;
        }
        void stce_10167() {
          stce_10190++;
          if (stce_10190 > stce_10191) {
            throw dco::exception::create<std::runtime_error>("Blob tape container overflow. Allocate bigger blob or use chunk tape container instead.");
          }
        }
        DCO_TAPE_INT index() const { return stce_10160; }
        DCO_TAPE_INT& index() { return stce_10160; }
        DCO_TAPE_INT stce_10168() { return ++stce_10160; }
        void stce_10193(entry const* end) { stce_10191 = end; }
        explicit iterator(entry* stce_10158=0, DCO_TAPE_INT index = DCO_TAPE_INT(-1)) : stce_10190(stce_10158), stce_10160(index) { }
        void operator--() {
          if(stce_10190->stce_10157())
            stce_10160--;
          stce_10190--;
        }
        bool stce_10157() const {
          return stce_10190->stce_10157();
        }
        iterator operator--(int) {
          iterator stce_10063(*this);
          operator--();
          return stce_10063;
        }
        entry* operator*() { return stce_10190; }
        entry* operator->() { return stce_10190; }
        entry* operator*() const { return stce_10190; }
        entry* operator->() const { return stce_10190; }
        bool operator==(const iterator& stce_9933) const { return stce_10190 == stce_9933.stce_10190; }
        bool operator!=(const iterator& stce_9933) const { return !operator==(stce_9933); }
        bool operator< (const iterator& stce_9933) const { return stce_10190 < stce_9933.stce_10190; }
        bool operator>=(const iterator& stce_9933) const { return !operator<(stce_9933); }
        bool operator> (const iterator& stce_9933) const { return stce_10190 > stce_9933.stce_10190; }
        bool operator<=(const iterator& stce_9933) const { return !operator>(stce_9933); }
        friend std::ostream& operator<< ( std::ostream &stce_10194, const iterator &stce_10164 ) {
          stce_10194 << "(" << stce_10164.stce_10190 << ":" << stce_10164.stce_10160 << ")";
          return stce_10194;
        }
      };
      class stce_9924 {
        iterator& stce_10169;
        DCO_TAPE_INT stce_10160;
      public:
        stce_9924(iterator& current) : stce_10169(current), stce_10160(0) {}
        iterator& current() { return stce_10169; }
        stce_9924& operator=(stce_9924 stce_9933) {
          swap(stce_9933);
          return *this;
        }
        DCO_TAPE_INT index() const { return stce_10160; }
        void swap(stce_9924& stce_9933) {
          std::swap(stce_10169, stce_9933.stce_10169);
          std::swap(stce_10160, stce_9933.stce_10160);
        }
        DCO_ALWAYS_INLINE_IF_AGRESSIVE void insert(const DCO_TAPE_INT& stce_10054, const stce_10152& pval) {
          if (IS_DEBUG && stce_10169.index() < stce_10054) {
            throw dco::exception::create<std::runtime_error>("DCO_DEBUG: tape index of read variable bigger than current progvarcounter.");
          }
          if(stce_10160 == 0) {
              stce_10169.stce_10167();
              new (*stce_10169) entry(entry::stce_10159(stce_10054),pval);
              stce_10160 = stce_10169.stce_10168();
              return;
            }
          if(stce_9757 && stce_10054 == stce_10169->arg()) {
              stce_10169->stce_10189(pval);
              return;
            }
          stce_10169.stce_10167();
          new (*stce_10169) entry(stce_10054,pval);
        }
      private:
        stce_9924(const stce_9924& stce_9933);
      };
      explicit stce_10149(stce_10149 &stce_9862,
                              stce_10186 stce_10195)
        : stce_10151(stce_9862),
          object_logger("BLOB container"),
          stce_10145(0),
          stce_10169(), stce_10179(), stce_10180(),
          stce_10196(stce_10195),
          stce_10197(stce_9862.stce_10197) {
        if (stce_10195 == stce_10187) {
            stce_10145 = stce_9862.stce_10145;
            stce_10169 = iterator(stce_10145, stce_9862.current().index());
            stce_10179 = stce_10169;
            stce_10180 = stce_9862.stce_10180;
            stce_10169.stce_10193(*stce_10180);
          } else if (stce_10195 == stce_10188) {
          DCO_TAPE_INT stce_10173 = allocate(stce_10197);
          stce_10169 = iterator(stce_10145, (stce_9862.current().index()));
          stce_10179 = iterator(stce_10145 , (stce_9862.current().index()));
          stce_10180 = iterator(stce_10145 + static_cast<DCO_INTEGRAL_TAPE_INT>(stce_10173) - 1, std::numeric_limits<DCO_INTEGRAL_TAPE_INT>::max());
          stce_10169.stce_10193(*stce_10180);
          new (*stce_10169) entry(entry::stce_10159(0));
        }
      }
  DCO_TAPE_INT allocate(const tape_options& stce_10170) {
    size_t stce_10175 = get_allocation_size(stce_10170.blob_size_in_byte());
    stce_10145 = static_cast<entry*>(stce_1000001(&stce_10175, stce_10170.alloc_use_huge_pages()));
    if (!stce_10145) {
      throw dco::exception::stce_9921<std::bad_alloc>
        ("Could not allocate memory. Use environment variables (DCO_MEM_RATIO, DCO_MAX_ALLOCATION) and see documentation.");
    }
    size_t stce_10198 = stce_10175 / entry::size_of();
    if (static_cast<double>(stce_10198) / static_cast<double>(std::numeric_limits<DCO_INTEGRAL_TAPE_INT>::max()) > 1.0) {
      stce_10198 = std::numeric_limits<DCO_TAPE_INT>::max();
    }
    DCO_TAPE_INT stce_10173 = static_cast<DCO_TAPE_INT>(stce_10198);
    DCO_LOG(logINFO) << "BLOB container; size = " << stce_10170.blob_size_in_byte() << "b ~=~ "
                     << stce_10170.blob_size_in_byte() / entry::size_of() << "elements";
    DCO_LOG(logINFO) << "BLOB container; tried allocation size = " << stce_10175
                     << "b ~=~ " << stce_10175 / entry::size_of() << "elements";
    DCO_LOG(logDEBUG) << "sizeof(entry) = " << sizeof(entry) << " == " << entry::size_of();
    return stce_10173;
  }
  void init(tape_options const& stce_10170) {
    DCO_TAPE_INT stce_10173 = allocate(stce_10170);
    stce_10169 = iterator(stce_10145, 0);
    stce_10179 = iterator(stce_10145, 0);
    stce_10180 = iterator(stce_10145 + static_cast<DCO_INTEGRAL_TAPE_INT>(stce_10173) - 1, std::numeric_limits<DCO_INTEGRAL_TAPE_INT>::max());
    stce_10169.stce_10193(*stce_10180);
    new (*stce_10169) entry(entry::stce_10159(0));
  }
  explicit stce_10149(tape_options const& stce_10170,
                          tape_t& tape)
    : stce_10151(tape),
      object_logger("BLOB container"),
      stce_10145(0),
      stce_10169(), stce_10179(), stce_10180(),
      stce_10196(stce_10188),
      stce_10197(stce_10170) {
    assert(stce_10150!=memory_model::stce_9770);
    init(stce_10170);
  }
  explicit stce_10149(iterator& stce_9817,
                          tape_t& tape)
    : stce_10151(tape),
      object_logger("BLOB container"),
      stce_10145(0),
      stce_10169(stce_9817), stce_10179(), stce_10180(),
      stce_10196(stce_10188),
      stce_10197(tape.stce_10170()) {
    assert(stce_10150==memory_model::stce_9770);
  }
      ~stce_10149() {
        if (stce_10196 == stce_10188)
          dco::dealloc(stce_10145);
      }
      DCO_ALWAYS_INLINE iterator& stce_9928() { return stce_10169; }
      entry* &stce_10199() { return stce_10145; }
      iterator current() const { return stce_10169; }
      iterator current() { return stce_10169; }
      iterator& stce_10200() { return stce_10169; }
      iterator start() { return stce_10179; }
      iterator start() const { return stce_10179; }
      iterator rend() { return start(); }
      iterator rbegin() const { return current(); }
      iterator capacity() const { return stce_10180; }
      bool empty() const { return stce_10169 == start(); }
      DCO_TAPE_INT size(const iterator& stce_9862, const iterator& stce_9863) const {
        if(stce_9862 < stce_9863)
          throw dco::exception::create<std::runtime_error>("size from < to");
        return static_cast<DCO_TAPE_INT>(*stce_9862 - *stce_9863);
      }
  tape_options stce_10170() const { return stce_10197; }
  DCO_TAPE_INT size(iterator stce_9862) const { return size(stce_9862, start()); }
  DCO_TAPE_INT size() const { return size(current(), start()); }
  DCO_TAPE_INT stce_10181() const { return size(stce_10180, start()); }
      size_t size_in_byte() const { return size(current(), start())*sizeof(entry); }
      size_t stce_10182() const { return size(stce_10180 , start())*sizeof(entry); }
      void erase(const iterator& stce_10158) {
        if(!stce_10145)
          return;
        if(IS_DEBUG && stce_10158 > stce_10180)
          throw dco::exception::create<std::runtime_error>("pos behind stack!");
        if(stce_10158 > stce_10169)
          throw dco::exception::create<std::runtime_error>("reset to position out of tape!");
        stce_10169 = stce_10158;
        stce_10169.stce_10193(*stce_10180);
      }
        void erase() {
        if(!stce_10145)
          return;
        stce_10169 = stce_10179;
        stce_10169.stce_10193(*stce_10180);
      }
      void stce_10183(const iterator&) { }
      void stce_10184 ( ) { }
      DCO_TAPE_INT stce_10185() { return stce_10169.index()+1; }
private:
  stce_10149(const stce_10149& stce_9933);
  entry * stce_10145;
  typename stce_9821
    <iterator, stce_10150==memory_model::stce_9770>::type stce_10169;
  iterator stce_10179;
  iterator stce_10180;
  stce_10186 stce_10196;
  tape_options stce_10197;
};
}
}
namespace dco {
namespace internal {
template<typename stce_10147, memory_model::TYPE stce_10150>
  class stce_10149<memory_model::CHUNK_TAPE, stce_10147, stce_10150>
  : public stce_10144<stce_10147>, object_logger
{
  typedef stce_10144<stce_10147> stce_10151;
private:
  typedef typename stce_9790<stce_10147>::type stce_10152;
protected:
  typedef stce_10149<memory_model::BLOB_TAPE, stce_10147> stce_10201;
  static void stce_10202(const std::string& stce_10203, char* data, long size)
  {
    DCO_LOG(logDEBUG1) << "save_data_to_file: bytes to write = " << size << ", fname = " << stce_10203;
    std::ofstream out(stce_10203.c_str(), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    if (!out.is_open())
      throw dco::exception::create<std::runtime_error>(
          "could not open file, perhaps index corruption?");
    if (!out.good())
      throw dco::exception::create<std::runtime_error>(
          "file could be opened, but got corrupted afterwards.");
    out.write(data, size);
    if (!out.good())
      throw dco::exception::create<std::runtime_error>(
          "reading the chunk failed, container corrupted.");
    out.close();
  }
  static DCO_INTEGRAL_TAPE_INT stce_10204(const std::string& stce_10203, char* data)
  {
    struct stat stce_10205;
    if (stat(stce_10203.c_str(), &stce_10205) == -1)
      throw dco::exception::create<std::runtime_error>(
          "could not check file size, perhaps index corruption?");
    DCO_INTEGRAL_TAPE_INT size = static_cast<DCO_INTEGRAL_TAPE_INT>(stce_10205.st_size);
    std::ifstream in(stce_10203.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!in.is_open())
      throw dco::exception::create<std::runtime_error>(
          "could not open file, perhaps index corruption?");
    if (!in.good())
      throw dco::exception::create<std::runtime_error>(
          "file could be opened, but got corrupted afterwards.");
    in.read(data, size);
    if (!in.good())
      throw dco::exception::create<std::runtime_error>(
          "reading the chunk failed, container corrupted.");
    return size;
  }
public:
  class stce_10206 : public stce_10201 {
    DCO_TAPE_INT stce_10207;
    tape_options stce_10197;
    bool stce_10208;
    stce_10206 *stce_10209, *stce_10210;
    DCO_TAPE_INT stce_10211;
public:
    typedef typename stce_10201::entry entry;
    typedef typename stce_10201::iterator iterator;
    DCO_TAPE_INT stce_10212() const { return stce_10207; }
    static std::string stce_10213(const stce_10206& stce_10201) {
      std::stringstream stce_10203;
      stce_10203 << stce_10201.stce_10214().filename() << '.' << stce_10201.stce_10212();
      return stce_10203.str();
    }
    stce_10206(stce_10206& stce_9862, typename stce_10201::stce_10186 stce_10215)
      : stce_10201(stce_9862, stce_10215), stce_10207(stce_9862.stce_10212()+1), stce_10197(stce_9862.stce_10214()), stce_10208(true),
        stce_10209(0), stce_10210(0), stce_10211(-1) { }
    stce_10206(tape_options stce_10170, stce_10147 &tape)
      : stce_10201(stce_10170, tape), stce_10207(0), stce_10197(stce_10170), stce_10208(true),
        stce_10209(0), stce_10210(0), stce_10211(-1) { }
    ~stce_10206() {
      if(stce_10214().write_to_file()) std::remove(stce_10213(*this).c_str());
    }
    DCO_TAPE_INT size(const iterator& stce_9862, const iterator& stce_9863) const {
      return stce_10201::size(stce_9862, stce_9863);
    }
    DCO_TAPE_INT size(iterator stce_9862) const {
      if(!*this->start())
        throw dco::exception::create<std::runtime_error>("given iterator does not point "
            "into given tape!");
      return stce_10201::size(stce_9862, this->start());
    }
    DCO_TAPE_INT size() const {
        return stce_10201::size(this->current(), this->start());
    }
    DCO_TAPE_INT stce_10181() const {
      if(!stce_10216())
        return 0;
      return stce_10201::stce_10181();
    }
    size_t size_in_byte() const {
      return stce_10201::size(this->current(), this->start()) * stce_10201::entry::size_of();
    }
    size_t stce_10182() const {
      if(!stce_10216()) { return 0; }
      return stce_10201::stce_10181() * stce_10201::entry::size_of();
    }
    tape_options stce_10214() const {
      return stce_10197;
    }
    bool stce_10216() const { return stce_10208; }
    void stce_10217(bool stce_10218) { stce_10208 = stce_10218; }
    void stce_10219() { stce_10211 = size(); }
    void stce_10220() { stce_10211 = -1; }
    bool stce_10221() { return stce_10211 != size(); }
    stce_10206* &stce_10222() { return stce_10209; }
    stce_10206* &stce_10223() { return stce_10210; }
  };
  typedef typename stce_10201::entry entry;
  typedef stce_10206* stce_10224;
  typedef typename stce_10201::iterator stce_10225;
  static void stce_10226(stce_10206& stce_10201) {
    if(stce_10201.stce_10214().write_to_file()) {
      assert(!stce_10201.stce_10216());
      DCO_LOG(logDEBUG1) << "wake_up: waking up chunk " << stce_10201.stce_10212() << " | " << stce_10201.stce_10199() << " from sleep.";
      DCO_INTEGRAL_TAPE_INT size = stce_10204(stce_10206::stce_10213(stce_10201),
                                                       reinterpret_cast<char*>(stce_10201.stce_10199()));
      (void)size;
      stce_10201.stce_10217(true);
    } else {
      DCO_LOG(logDEBUG1) << "wake_up: waking already awake chunk " << stce_10201.stce_10212() << " | " << stce_10201.stce_10199();
    }
  }
  static void sleep(stce_10206& stce_10201) {
    DCO_LOG(logDEBUG1) << "sleep: chunk " << stce_10201.stce_10212() << " | " << stce_10201.stce_10199() << " dirty = " << stce_10201.stce_10221();
    if(stce_10201.stce_10214().write_to_file() && stce_10201.stce_10221()) {
      assert(stce_10201.stce_10216());
      DCO_LOG(logDEBUG1) << "sleep: writing chunk " << stce_10201.stce_10212() << " | " << stce_10201.stce_10199();
      const long stce_10227 = stce_9752(stce_10201.size()+1) * sizeof(entry);
      DCO_LOG(logDEBUG1) << "sleep: chunk.size() = " << stce_10201.size() << ", bytes to write = " << stce_10227 << ", ptr diff = " << stce_10201.current() << ", " << stce_10201.start() << " = " << *stce_10201.current() - *stce_10201.start() << " " << sizeof(entry);
      entry* data = stce_10201.stce_10199();
      stce_10202(stce_10206::stce_10213(stce_10201),
                        reinterpret_cast<char*>(data), stce_10227);
      stce_10201.stce_10219();
    } else {
      if (!stce_10201.stce_10221()) {
        DCO_LOG(logDEBUG1) << "sleep: not writing data for chunk: " << stce_10201.stce_10212();
      }
    }
  }
  static void stce_10228(stce_10206& stce_10201, stce_10206& stce_10229) {
    if (!stce_10201.stce_10216()) {
      sleep(stce_10229);
      stce_10229.stce_10217(false);
      stce_10226(stce_10201);
    }
  }
  void stce_10228(stce_10206& stce_10201) {
    if (!stce_10201.stce_10216()) {
      DCO_LOG(logDEBUG1) << "activate: chunk " << stce_10201.stce_10212() << " not awake yet";
      stce_10206* stce_10229 = stce_10169.stce_10230();
      while (stce_10229 && !stce_10229->stce_10216()) { stce_10229 = stce_10229->stce_10222(); }
      if (stce_10229) {
        DCO_LOG(logDEBUG1) << "activate: currently awake " << stce_10229->stce_10212();
        stce_10228(stce_10201, *stce_10229);
      } else {
        stce_10226(stce_10201);
      }
    }
  }
  class iterator {
private:
    stce_10224 stce_10231;
    stce_10225 stce_10232;
public:
    template<class stce_10162>
      void stce_10166(const stce_10162 &stce_10164) {
      assert(sizeof(iterator) <= sizeof(stce_10162)) ;
      stce_10165(static_cast<void*>(this), &stce_10164, sizeof(iterator));
    }
    template<class stce_10162>
      void stce_10163(stce_10162 &stce_10164) const {
      assert(sizeof(iterator) <= sizeof(stce_10162)) ;
      stce_10165(&stce_10164, this, sizeof(iterator));
    }
    stce_10224 stce_10230() const {
      return stce_10231;
    }
    stce_10225 stce_10233() const {
      return stce_10232;
    }
    explicit iterator(const stce_10225& stce_10234, stce_10224 stce_10230) :
        stce_10231(stce_10230), stce_10232(stce_10234) {
    }
    explicit iterator(stce_10224 stce_10230) :
        stce_10231(stce_10230), stce_10232(stce_10230->current()) {
    }
    iterator(const iterator& stce_9933) :
        stce_10231(stce_9933.stce_10231), stce_10232(stce_9933.stce_10232) {
    }
    iterator() : stce_10231(), stce_10232() {}
    bool stce_10157() const {
      return stce_10232.stce_10157();
    }
    void operator--() {
      DCO_LOG(logDEBUG4) << stce_10232 << " =?= " << stce_10231->start();
      --stce_10232;
      if (stce_10232 == stce_10231->start()) {
          if(stce_10231->stce_10212()==0) return;
        stce_10231 = stce_10231->stce_10222();
        stce_10228(*stce_10231, *(stce_10231->stce_10223()));
        stce_10232 = stce_10231->current();
      }
    }
    DCO_TAPE_INT index() const { return stce_10232.index(); }
    DCO_TAPE_INT& index() { return stce_10232.index(); }
    DCO_TAPE_INT stce_10168() { stce_10231->stce_10200().stce_10168(); return stce_10232.stce_10168(); }
    stce_10224& stce_10201() { return stce_10231; }
    stce_10225& stce_10235() { return stce_10232; }
    iterator operator--(int) {
      iterator stce_10063(*this);
      operator--();
      return stce_10063;
    }
    entry* operator* () { return *stce_10232; }
    entry* operator* () const { return *stce_10232; }
    entry* operator->() { return *stce_10232; }
    entry* operator->() const { return *stce_10232; }
    bool operator==(const iterator& stce_9933) const {
      return stce_10232 == stce_9933.stce_10232 && stce_10231 == stce_9933.stce_10231;
    }
    bool operator!=(const iterator& stce_9933) const {
      return !operator==(stce_9933);
    }
    bool operator<(const iterator& stce_9933) const {
      return
          (stce_10231 == stce_9933.stce_10231) ?
              stce_10232 < stce_9933.stce_10232 : stce_10232.index() < stce_9933.stce_10232.index();
    }
    bool operator>=(const iterator& stce_9933) const {
      return !operator<(stce_9933);
    }
    bool operator>(const iterator& stce_9933) const {
      return
          (stce_10231 == stce_9933.stce_10231) ?
              stce_10232 > stce_9933.stce_10232 : stce_10232.index() > stce_9933.stce_10232.index();
    }
    bool operator<=(const iterator& stce_9933) const {
      return !operator>(stce_9933);
    }
  };
  class stce_9924
  {
    stce_10224& stce_10236;
    stce_10225& stce_10237;
    stce_10225 stce_10232;
    stce_10225 stce_10238;
    DCO_TAPE_INT stce_10160;
public:
    stce_9924(iterator& stce_9817) :
      stce_10236(stce_9817.stce_10201()),
      stce_10237(stce_9817.stce_10235()),
      stce_10232(stce_10236->current()),
      stce_10238(stce_10236->capacity()),
      stce_10160(0) {
    }
    stce_10225& current() { return stce_10232; }
    ~stce_9924() {
      if(static_cast<DCO_INTEGRAL_TAPE_INT>(stce_10160)) {
        stce_10237 = stce_10236->stce_10200() = stce_10232;
      }
    }
    DCO_TAPE_INT index() const { return stce_10160; }
    DCO_TAPE_INT& index() { return stce_10160; }
    DCO_ALWAYS_INLINE_IF_AGRESSIVE void insert(const DCO_TAPE_INT& stce_10054, const stce_10152& pval) {
      if (IS_DEBUG && stce_10232.index() < stce_10054) {
        throw dco::exception::create<std::runtime_error>("DCO_DEBUG: tape index of read variable bigger than current progvarcounter.");
      }
      if (stce_10160 == 0) {
        stce_10232.stce_10167();
        new (*stce_10232) entry(entry::stce_10159(stce_10054), pval);
        stce_10160 = stce_10232.stce_10168();
        DCO_LOG(logDEBUG2) << "inserting first entry, index now = " << stce_10160 << "; pval = " << pval << ", target = " << stce_10054 << " in chunk " << stce_10236 << " into memory at = " << stce_10232;
        return;
      }
      if(stce_9757 && stce_10232->arg() == stce_10054){
        DCO_LOG(logDEBUG3) << "merging a parallel edge";
        stce_10232->stce_10189(pval);
        return;
      }
      if (*stce_10232 >= *stce_10238) {
        DCO_LOG(logDEBUG1) << "insert: chunk full, updating current and getting next";
        stce_10237 = stce_10236->stce_10200() = stce_10232;
        stce_10239(stce_10236);
        stce_10232 = stce_10236->current();
        stce_10238 = stce_10236->capacity();
      }
      stce_10232.stce_10167();
      new (*stce_10232) entry(stce_10054, pval);
      DCO_LOG(logDEBUG2) << "inserting next entry, index now = " << stce_10160 << "; pval = " << pval << ", target = " << stce_10054 << " in chunk " << stce_10236 << " into memory at = " << stce_10232;
    }
  private:
    stce_9924(const stce_9924& stce_9933);
  };
public:
  void init(tape_options const& stce_10170) {
    tape_options stce_10240 = stce_10170;
    stce_10240.blob_size_in_byte() = stce_10170.chunk_size_in_byte();
    stce_10241 = new stce_10206(stce_10240, stce_10151::stce_10146());
    stce_10169 = iterator(stce_10241);
  }
  explicit stce_10149(tape_options const& stce_10170, stce_10147& tape) :
    stce_10151(tape),
    object_logger("CHUNK container (w/ options)"),
    stce_10241(0),
    stce_10169() {
    assert(stce_10150!=memory_model::stce_9770);
    init(stce_10170);
  }
  explicit stce_10149(iterator &stce_9817, stce_10147& tape) :
    stce_10151(tape),
    object_logger("CHUNK container (w/ options)"),
    stce_10241(0), stce_10169(stce_9817) {
    assert(stce_10150==memory_model::stce_9770);
  }
  ~stce_10149() {
    stce_10206* stce_10164 = stce_10241;
    while (stce_10164) {
      stce_10206* stce_10242 = stce_10164;
      stce_10164 = stce_10164->stce_10223();
      delete stce_10242;
    }
  }
  static void stce_10239(stce_10224& stce_9817)
  {
    DCO_LOG(logDEBUG1) << "current chunk: id = " << stce_9817->stce_10212() << ", index = " << stce_9817->current().index()
                       << ", data pointer = " << stce_9817->stce_10199();
    delete stce_9817->stce_10223();
    if (stce_9817->stce_10214().write_to_file()) {
      sleep(*stce_9817);
      stce_9817->stce_10217(false);
      stce_9817->stce_10223() = new stce_10206(*stce_9817, stce_10206::stce_10187);
    } else {
      stce_9817->stce_10223() = new stce_10206(*stce_9817, stce_10206::stce_10188);
    }
    stce_9817->stce_10223()->stce_10222() = stce_9817;
    stce_9817 = stce_9817->stce_10223();
    DCO_LOG(logDEBUG1) << "next chunk: id = " << stce_9817->stce_10212() << ", index = " << stce_9817->current().index()
                       << ", data pointer = " << stce_9817->stce_10199();
  }
  DCO_ALWAYS_INLINE iterator& stce_9928()
  {
    if (stce_10169.stce_10201()->current() >= stce_10169.stce_10201()->capacity()) {
      DCO_LOG(logDEBUG1) << "insert_ref: chunk full, getting next";
      stce_10239(stce_10169.stce_10201());
    }
    return stce_10169;
  }
  iterator start() { return iterator(stce_10241->start(), stce_10241); }
  iterator start() const { return iterator(stce_10241->start(), stce_10241); }
  iterator rend() { return start(); }
  iterator current() const { return stce_10169; }
  iterator rbegin() const { return current(); }
  bool empty() const { return current() == start(); }
  DCO_TAPE_INT size(const iterator& stce_9862, const iterator& stce_9863) const {
    const stce_10206& stce_10243 = *(stce_9862.stce_10230());
    if(stce_9862.stce_10230() == stce_9863.stce_10230())
      return stce_10243.size(stce_9862.stce_10233(), stce_9863.stce_10233());
    DCO_TAPE_INT size = 0;
    stce_10224 stce_10244 = stce_9862.stce_10230();
    size += stce_10244->size(stce_9862.stce_10233());
    DCO_LOG(logDEBUG2) << "size of last chunk = " << size;
    for (stce_10244 = stce_10244->stce_10222(); stce_10244 != stce_9863.stce_10230(); stce_10244 = stce_10244->stce_10222()) {
      if(stce_10244) {
        size += stce_10244->size();
        DCO_LOG(logDEBUG2) << "adding " << stce_10244->size() << " => size = " << size;
      }
    }
    const stce_10206& stce_10245 = *(stce_9863.stce_10230());
    size += stce_10245.size(stce_10245.current(), stce_9863.stce_10233());
    DCO_LOG(logDEBUG2) << "very first bit: adding " << stce_10245.size(stce_10245.current(), stce_9863.stce_10233()) << " => size = " << size;
    return size;
  }
  DCO_TAPE_INT size(iterator stce_9862) const { return size(stce_9862, start()); }
  DCO_TAPE_INT size() const { return size(current(), start()); }
  mem_long_t size_in_byte() const {
    return size(current(), start())*entry::size_of();
  }
  DCO_TAPE_INT stce_10181() const {
    DCO_TAPE_INT size = 0;
    for(stce_10224 stce_10246 = stce_10241; stce_10246 != NULL; stce_10246 = stce_10246->stce_10223())
      if(stce_10246)
        size += stce_10246->stce_10181();
    return size;
  }
  mem_long_t stce_10182() const {
    mem_long_t size = 0;
    for(stce_10224 stce_10246 = stce_10241; stce_10246 != NULL; stce_10246 = stce_10246->stce_10223())
      if(stce_10246) size += stce_10246->stce_10182();
    return size;
  }
  void erase(const iterator& stce_10158) {
    if (stce_10158 > current())
      throw dco::exception::create<std::invalid_argument>("erase pos out of tape!");
    stce_10224 stce_10247(stce_10241);
    while (stce_10247 != stce_10158.stce_10230())
      stce_10247 = stce_10247->stce_10223();
    stce_10169.stce_10201() = stce_10247;
    stce_10228(*stce_10169.stce_10201());
    stce_10169.stce_10201()->stce_10220();
    stce_10224 stce_10164 = stce_10247->stce_10223();
    stce_10247->stce_10223() = 0;
    while (stce_10164) {
      stce_10224 stce_10242 = stce_10164;
      stce_10164 = stce_10164->stce_10223();
      delete stce_10242;
    }
    stce_10169.stce_10201()->erase(stce_10158.stce_10233());
    stce_10169.stce_10235() = stce_10169.stce_10201()->current();
  }
  void stce_10183(const iterator& stce_9817) {
    DCO_LOG(logDEBUG1) << "chunk tape: prepare_interpretation; activate: " << stce_9817.stce_10230()->stce_10212() << ", is awake = " << stce_9817.stce_10230()->stce_10216();
    stce_10228(*stce_9817.stce_10230());
  }
  void stce_10184() {
    DCO_LOG(logDEBUG1) << "chunk tape: prepare_recording";
    stce_10228(*stce_10169.stce_10201());
  }
  DCO_TAPE_INT stce_10185() { return current().index()+1; }
private:
  stce_10206* stce_10241;
  typename stce_9821
    <iterator, stce_10150==memory_model::stce_9770>::type stce_10169;
};
}
}
       
namespace dco {
namespace internal {
template<class stce_10147>
  class stce_10149<memory_model::stce_9769, stce_10147>
  : public stce_10149<memory_model::CHUNK_TAPE, stce_10147>
{
  typedef stce_10147 tape_t;
  typedef typename stce_9790<stce_10147>::type stce_10152;
public:
  typedef stce_10149<memory_model::CHUNK_TAPE, stce_10147> stce_10248;
  typedef typename stce_10248::entry entry;
  typedef typename stce_10248::stce_9924 stce_9924;
  typedef typename stce_10248::iterator iterator;
  bool stce_10249, stce_10010, stce_10250, stce_10251;
  double stce_10252, stce_10253, stce_10004;
  std::vector<DCO_TAPE_INT> stce_10013;
  std::vector<stce_10152> stce_10014;
  std::vector<const void*> stce_10254;
  DCO_TAPE_INT stce_10255, stce_10256;
  struct stce_10257 {
    static std::string id() { return " gtas "; }
  };
  typedef logging<stce_10257> logger;
  explicit stce_10149(tape_options stce_10170, stce_10147 &tape) :
    stce_10248(stce_10170, tape),
    stce_10249(false),
    stce_10010(false),
    stce_10250(false),
    stce_10251(false),
    stce_10252(1e-6),
    stce_10253(1e-6),
    stce_10004(0.0),
    stce_10255(0),
    stce_10256(-1) {}
  void stce_10258(log_level_enum level) {
    if (stce_10251 && level <= dco::logWARNING) {
      throw dco::exception::create<std::runtime_error>("Warning occured.");
    }
  }
  void stce_10009() {
    if (stce_10256 == stce_10255) {
      stce_10255++;
      throw dco::exception::create<std::runtime_error>("Event counter reached.");
    }
    stce_10255++;
  };
  void stce_10008(stce_10152 stce_9937, double stce_10007)
  {
    if(stce_10250) {
      double stce_10259 = ::fabs(stce_9937 - stce_10007);
      log_level_enum stce_10260 = dco::logWARNING;
      if (stce_10259 < stce_10252) {stce_10260 = dco::logINFO;}
      logger(stce_10260).get_stream() << stce_10255 << ": fd vs. tangent check";
      logger(stce_10260).get_stream() << " delta = " << stce_10259 << " tangent = " << stce_9937
          << " finite differences = " << stce_10007;
      stce_10258(stce_10260);
    }
  }
  void stce_10015(double stce_10012, double stce_10261, double stce_10262)
  {
    const double stce_10263 = ::fabs(stce_10012-stce_10261);
    const double stce_10264 = ::fabs(stce_10012-stce_10262);
    log_level_enum stce_10260 = dco::logWARNING;
    if (stce_10263 < stce_10253) {stce_10260 = dco::logINFO;}
    logger(stce_10260).get_stream() << stce_10255 << ": Tangent/adjoint identity check:";
    logger(stce_10260).get_stream() << " tangent vs. adjoint = " << stce_10263;
    logger(stce_10260).get_stream() << "  <x_(1),x^(1)> = " << stce_10012;
    logger(stce_10260).get_stream() << "  <y_(1),y^(1)> = " << stce_10261;
    stce_10258(stce_10260);
    stce_10260 = dco::logWARNING;
    if (stce_10264 < stce_10252) {stce_10260 = dco::logINFO;}
    logger(stce_10260).get_stream() << " finite differences vs. adjoint = " << stce_10264;
    logger(stce_10260).get_stream() << "  <y_(1),y_{FD}^(1)> = " << stce_10262;
    stce_10258(stce_10260);
  }
  void stce_10265(const iterator &stce_9863) {
    assert(stce_10248::start() == stce_9863);
    stce_10248::stce_10265(stce_9863);
    stce_10013.clear();
    stce_10014.clear();
    stce_10254.clear();
    stce_10004=0;
  }
  void stce_9981(std::string stce_10266, const bool stce_9840, const bool stce_9978, const bool stce_9980)
  {
    if (!stce_10249) {
      stce_10009();
      return;
    }
    log_level_enum stce_10260 = dco::logWARNING;
    if (stce_9840 == stce_9978 && stce_9840 == stce_9980) {
      stce_10260 = dco::logINFO;
    }
    logger(stce_10260).get_stream() << stce_10255 << ": Branch check:" << stce_10266
        << std::boolalpha << " at current value returns " << stce_9840 << ", at minus h " << stce_9978
        << ", and at plus h " << stce_9980;
    stce_10258(stce_10260);
    stce_10009();
  }
  template <class data_t>
  void stce_10005(const data_t* data, const double &stce_9937) {
    if(stce_10004 == 0) {
        stce_10254.push_back(data);
        stce_10013.push_back(data->_tape_index());
        stce_10014.push_back(stce_9937);
    }
  }
  template <typename stce_10267>
  void stce_10268() {
      double stce_10269=0;
      for(size_t stce_9832=0;stce_9832<stce_10254.size();++stce_9832) {
        const stce_10267 &stce_9817 = *static_cast<const stce_10267*>(stce_10254[stce_9832]);
        stce_10269 = std::max(stce_10269, std::fabs(stce_9817.stce_10002[0]));
      }
      const double stce_10270 = 1e-14;
      stce_10004 = std::sqrt(stce_10270) * stce_10269;
      for(size_t stce_9832=0;stce_9832<stce_10254.size();++stce_9832) {
        const stce_10267 &stce_9817 = *static_cast<const stce_10267*>(stce_10254[stce_9832]);
        stce_9817.stce_10002[0] -= stce_10004 * stce_9817.stce_9936;
        stce_9817.stce_10002[1] += stce_10004 * stce_9817.stce_9936;
      }
  }
};
}
}
       
namespace dco {
namespace internal {
template<class stce_10147, memory_model::TYPE stce_10148>
  class stce_10149<memory_model::stce_9770, stce_10147, stce_10148>
  : public stce_10149<stce_10148, stce_10147, memory_model::stce_9770> {
  typedef typename stce_9790<stce_10147>::type stce_10152;
public:
  typedef stce_10149<stce_10148, stce_10147, memory_model::stce_9770> stce_10248;
  explicit stce_10149(tape_options stce_10170, stce_10147 &tape) :
    stce_10248(stce_10169, tape), stce_10169() {
    stce_10248::init(stce_10170);
  }
  struct iterator : public stce_10248::iterator {
    iterator() : stce_10271(1), stce_10272(0) {}
    iterator(const typename stce_10248::iterator& stce_10273) : stce_10248::iterator(stce_10273), stce_10271(-1), stce_10272(-1) {}
    DCO_TAPE_INT stce_10271;
    DCO_TAPE_INT stce_10272;
    const DCO_TAPE_INT& max_dist() const { return stce_10271; }
    void stce_10274(const DCO_TAPE_INT &stce_10275) {
      stce_10271 = stce_10275+1;
      DCO_LOG(logDEBUG3) << " updating max dist: " << stce_10271;
    }
    const DCO_TAPE_INT &stce_10276() const { return stce_10272; }
    DCO_TAPE_INT stce_10168() {
      stce_10272++;
      return stce_10248::iterator::stce_10168();
    }
  };
  class stce_9924 : public stce_10248::stce_9924 {
    iterator &stce_10169;
public:
    stce_9924(iterator& current) : stce_10248::stce_9924(current), stce_10169(current) { }
    iterator& current() { return stce_10169; }
    DCO_ALWAYS_INLINE void insert(const DCO_TAPE_INT& stce_10054, const stce_10152& pval) {
      DCO_LOG(logDEBUG4) << "target = " << stce_10054 << ", parametercounter = " << current().stce_10276() << ", max dist = " << current().max_dist();
      DCO_TAPE_INT stce_10277 = current().index();
      stce_10248::stce_9924::insert(stce_10054, pval);
      if(stce_10054 > current().stce_10276()) {
        DCO_TAPE_INT stce_10278 = stce_10277 - stce_10054+1;
        DCO_LOG(logDEBUG4) << "curdist = " << stce_10278;
        if(stce_10278 - current().max_dist() > 0) {
          current().stce_10274(stce_10278);
        }
      }
    }
private:
    stce_9924(const stce_9924& stce_9933);
  };
  DCO_ALWAYS_INLINE iterator& stce_9928() {
    return static_cast<iterator&>(stce_10248::stce_9928());
  }
  iterator current() const { return stce_10169; }
  iterator current() { return stce_10169; }
  iterator& stce_10200() { return stce_10169; }
  iterator rbegin() const { return stce_10169; }
  DCO_TAPE_INT stce_10185() {
    using ::log;
    using ::pow;
    DCO_TAPE_INT max_dist = stce_10169.max_dist();
    bool stce_10279 = (max_dist & (max_dist - 1)) == 0;
    if (HAS_BITWISE_MODULO && !stce_10279) {
      DCO_TAPE_INT stce_10095 = static_cast<DCO_TAPE_INT>(floor(log(static_cast<double>(max_dist))/log(2.0)))+1;
      stce_10169.stce_10274( (1 << stce_9752(stce_10095)) - 1 );
      DCO_LOG(logDEBUG2) << "rounding adjoint vector size to " << stce_10169.max_dist();
    }
    return stce_10169.stce_10276() + stce_10169.max_dist() + 1;
  }
private:
  iterator stce_10169;
};
}
namespace helper {
template<typename stce_9788, memory_model::TYPE stce_9789, typename external_adjoint_object_t>
  struct stce_10078<internal::tape<memory_model::stce_9770, stce_9788, stce_9789>, external_adjoint_object_t> {
  typedef internal::tape<memory_model::stce_9770, stce_9788, stce_9789> tape_t;
  DCO_TAPE_INT stce_10280;
  DCO_TAPE_INT stce_10281;
  tape_t* &stce_9952;
  stce_10078(external_adjoint_object_t* stce_10282) :
    stce_10280(std::numeric_limits<DCO_TAPE_INT>::max()),
    stce_10281(0),
    stce_9952(stce_10282->tape()) {}
  template <typename stce_9782> void register_input(const stce_9782& stce_9843) {
    DCO_TAPE_INT stce_9808 = stce_9843._tape_index();
    if(stce_9808 > stce_9952->current().stce_10276()) {
      stce_10280 = std::min(stce_10280, stce_9808);
    }
  }
  template <typename stce_9782> void register_output(const stce_9782& stce_10120) {
    stce_10281 = std::max(stce_10281, stce_10120._tape_index());
    if(stce_10280 != std::numeric_limits<DCO_TAPE_INT>::max()) {
        DCO_TAPE_INT max_dist = stce_10281 - stce_10280;
        if(max_dist > stce_9952->stce_10200().max_dist()) {
          stce_9952->stce_10200().stce_10274(max_dist);
        }
    }
  }
};
}
namespace internal {
template<typename stce_9791, typename stce_10147, memory_model::TYPE stce_10148>
  class stce_9793<stce_9791, internal::stce_10149<memory_model::stce_9770, stce_10147, stce_10148> >
  : public stce_9793<stce_9791, internal::stce_10149<stce_10148, stce_10147, memory_model::stce_9770> >
{
  typedef internal::stce_10149<memory_model::stce_9770, stce_10147, stce_10148> stce_10283;
  typedef typename stce_10283::stce_10248 stce_10248;
  typedef stce_9793<stce_9791, stce_10248> stce_10284;
  typedef typename stce_10284::stce_9764 stce_9764;
public:
  stce_10283 &stce_9952;
  stce_9793(DCO_TAPE_INT size, stce_10283 &tape, tape_options stce_10170) :
    stce_10284(size, static_cast<stce_10248&>(tape), stce_10170), stce_9952(tape),stce_10285(true), stce_10286(0) { }
  void init() {
    stce_10287 = stce_9952.current().stce_10276();
    if (HAS_BITWISE_MODULO) {
      stce_10288 = stce_9952.current().max_dist()-1;
    } else {
      stce_10288 = stce_9952.current().max_dist();
    }
    stce_10289 = this->stce_10290 + stce_9752(stce_10287) + 1;
  }
  bool stce_10291(DCO_TAPE_INT stce_10292) const {
    if(stce_10292> stce_9952.current().stce_10276()) {
      stce_10292 -= stce_9952.current().stce_10276();
      stce_10292 %= stce_9952.current().max_dist();
      stce_10292 += stce_9952.current().stce_10276() + 1;
    }
    return static_cast<DCO_TAPE_INT>(this->stce_10293) > stce_10292;
  }
  DCO_TAPE_INT stce_10294(DCO_TAPE_INT stce_9808) const {
    if(stce_9808 > stce_9952.current().stce_10276()) {
      stce_9808 -= stce_9952.current().stce_10276();
      stce_9808 %= stce_9952.current().max_dist();
      stce_9808 += stce_9952.current().stce_10276() + 1;
    }
    return stce_9808;
  }
  typename stce_10284::stce_9764& operator[](DCO_TAPE_INT stce_9808) {
    DCO_LOG(logDEBUG4) << "getting " << stce_9808;
    stce_9808 = stce_10294(stce_9808);
    DCO_LOG(logDEBUG4) << "matched to " << stce_9808;
    return stce_10284::operator[](stce_9808);
  }
  typename stce_10284::stce_9764 operator[](DCO_TAPE_INT stce_9808) const {
    DCO_LOG(logDEBUG4) << "getting " << stce_9808;
    stce_9808 = stce_10294(stce_9808);
    DCO_LOG(logDEBUG4) << "matched to " << stce_9808;
    return stce_10284::operator[](stce_9808);
  }
  stce_9764 at(DCO_TAPE_INT stce_9808) const { return stce_10291(stce_9808)?(*this)[stce_9808]:0; }
  stce_9764 &at(const DCO_TAPE_INT stce_9808) {
    if(!stce_10291(stce_9808))
      this->resize(stce_9952.stce_10185());
    return (*this)[stce_9808];
  }
  void stce_10295(DCO_TAPE_INT stce_9808) {
    if(!stce_10291(stce_9808)) {
      DCO_LOG(logDEBUG2) << "Resizing internal adjoint vector to " << stce_9808;
      this->resize(stce_9952.stce_10185());
    }
  }
  template <typename iterator_t>
  void stce_10296(iterator_t& stce_9817) {
    DCO_TAPE_INT stce_9863 = stce_9817->arg();
    stce_9764 pval = stce_9817->pval();
    DCO_TAPE_INT stce_9862 = stce_9817.index();
    if (stce_10285) {
      DCO_LOG(logDEBUG3) << "setting " << stce_9862 << " to 0 ";
      stce_9764& stce_10297 = HAS_BITWISE_MODULO ?
                           stce_10289[ (stce_9862 - stce_10287) & stce_10288 ] :
                           stce_10289[ (stce_9862 - stce_10287) % stce_10288 ];
      this->stce_10286 = stce_10297;
      stce_10297 = 0;
    }
    if(stce_9863 > stce_10287) {
      stce_9764& stce_10297 = HAS_BITWISE_MODULO ?
                           stce_10289[ (stce_9863 - stce_10287) & stce_10288 ] :
                           stce_10289[ (stce_9863 - stce_10287) % stce_10288 ];
      stce_10297 += pval*this->stce_10286;
    } else {
      (*this)[stce_9863] += pval*this->stce_10286;
    }
    stce_10285 = stce_9817.stce_10157();
  }
private:
  stce_9764* stce_10289;
  DCO_TAPE_INT stce_10287;
  DCO_TAPE_INT stce_10288;
  bool stce_10285;
  stce_9764 stce_10286;
};
}
}
namespace dco {
namespace internal {
template<class stce_9761, typename stce_10143>
  class stce_9793 : public stce_10144<stce_10143>
{
  typedef stce_10144<stce_10143> stce_10151;
public:
  typedef stce_9761 stce_9764;
protected:
  stce_9764 *stce_10290;
  DCO_TAPE_INT stce_10293;
  DCO_TAPE_INT stce_10298;
  DCO_TAPE_INT stce_10299;
  int stce_10300;
  tape_options stce_10197;
public:
  typedef stce_9764 value_t;
  void init() {}
  stce_9793(DCO_TAPE_INT size, stce_10143 &tape, tape_options stce_10170) :
    stce_10151(tape), stce_10290(0), stce_10293(0), stce_10298(0),
    stce_10299(size), stce_10300(0), stce_10197(stce_10170) { }
  ~stce_9793() {
    if (stce_10290) dco::dealloc(stce_10290);
  }
  DCO_TAPE_INT size() const { return stce_10299; }
  bool stce_10291(DCO_TAPE_INT stce_10292) const {
    return stce_10293 > stce_10292;
  }
  mem_long_t stce_10301() const {
    return static_cast<mem_long_t>(stce_10298) * sizeof(stce_9764);
  }
  stce_9764& operator[](DCO_TAPE_INT stce_9808) { return stce_10290[stce_9808]; }
  stce_9764 operator[](DCO_TAPE_INT stce_9808) const { return stce_10290[stce_9808]; }
  void stce_10296(typename stce_10143::iterator& stce_9817) {
    (*this)[stce_9817->arg()] += stce_9817->pval()*(*this)[stce_9817.index()];
  }
  void resize(DCO_TAPE_INT stce_10302) {
    DCO_LOG(logDEBUG2) << "Resizing internal adjoint vector to " << stce_10302;
    if (stce_10302 < stce_10299) {
      stce_10303(stce_10299, stce_10302-1);
    }
    stce_10299 = stce_10302;
  }
  void stce_10304(DCO_TAPE_INT stce_10102) {
    if (++stce_10300 % 100 == 0) {
      DCO_LOG(logWARNING) << "dco_adjoint_vector: resizing adjoint vector quite often (#=" << stce_10300
                          << ") -- might have severe performance impact (reallocation / initialization / copy)";
    }
    if (stce_10298 >= stce_10102) {
      DCO_LOG(logDEBUG1) << "dco_adjoint_vector: resizing adjoint vector to " << stce_10102 << "; no reallocation required.";
      for (DCO_TAPE_INT stce_9832 = stce_10293; stce_9832 < stce_10102; ++stce_9832) {
        stce_10290[stce_9832] = stce_9764();
      }
      stce_10293 = stce_10102;
      return;
    } else {
      DCO_LOG(logDEBUG1) << "dco_adjoint_vector: resizing adjoint vector to " << stce_10102 << "; reallocation required.";
      stce_9764* stce_10305 = static_cast<stce_9764*>(dco::alloc(stce_10102 * sizeof(stce_9764), stce_10197.alloc_use_huge_pages()));
      if (stce_10305 == 0) {
        throw dco::exception::create<std::runtime_error>("Unable to allocate memory for adjoint vector.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_lazy_adjoint_vector.hpp", 101);
      }
      if (stce_10290) {
        for (DCO_TAPE_INT stce_9832 = 0; stce_9832 < stce_10293; ++stce_9832) { stce_10305[stce_9832] = stce_10290[stce_9832]; }
        for (DCO_TAPE_INT stce_9832 = stce_10293; stce_9832 < stce_10102 ; ++stce_9832) { stce_10305[stce_9832] = stce_9764(); }
        dco::dealloc(stce_10290);
      } else {
        for (DCO_TAPE_INT stce_9832 = 0; stce_9832 < stce_10102 ; ++stce_9832) { stce_10305[stce_9832] = stce_9764(); }
      }
      stce_10293 = stce_10102;
      stce_10298 = stce_10102;
      stce_10290 = stce_10305;
    }
  }
  void stce_10295(DCO_TAPE_INT stce_9808) {
    if(!stce_10291(stce_9808)) { stce_10304(stce_9808); }
  }
  stce_9764 at(DCO_TAPE_INT stce_9808) const {
    return stce_10291(stce_9808)?stce_10290[stce_9808]:0;
  }
  stce_9764 &at(const DCO_TAPE_INT stce_9808) {
    if(!stce_10291(stce_9808)) { stce_10304(stce_9808+1); }
    return stce_10290[stce_9808];
  }
  void stce_10303(DCO_TAPE_INT stce_9862, DCO_TAPE_INT stce_9863) {
    if(stce_9862 < stce_9863) { std::swap(stce_9862, stce_9863); }
    if(!stce_10291(stce_9863)) { return; }
    stce_10306();
    for (DCO_TAPE_INT stce_9832 = std::min(stce_9862, stce_10299-1); stce_9832 > stce_9863; --stce_9832) {
      stce_10290[stce_9832] = 0;
    }
  }
  void stce_10306() { stce_10304(stce_10299); }
};
}
}
namespace dco {
template<class stce_10143, class stce_9761, int stce_10307=-1>
  struct adjoint_vector : dco::stce_10079<typename stce_10143::tape_t::stce_9764>, public internal::stce_10144<stce_10143> {
    static const int stce_9837 = stce_10307 < 0 ? -stce_10307 : stce_10307;
  typedef internal::stce_10144<stce_10143> stce_10151;
public:
  typedef dco::helper::stce_9827<stce_9761, stce_9837> value_t;
  typedef stce_9761 adjoint_real_t;
  typedef typename stce_10143::stce_9764 stce_9764;
private:
  value_t* stce_10290;
  DCO_TAPE_INT stce_10308;
public:
  int stce_10309;
  void init() {}
  adjoint_vector(stce_10143 *tape) : stce_10151(*tape), stce_10290(0), stce_10308(0), stce_10309(0) { }
  virtual ~adjoint_vector() {
    if (stce_10290) delete [] stce_10290;
  }
  template<class stce_9794>
    value_t& derivative(const stce_9794 &stce_9886) {
    stce_10310();
    const DCO_TAPE_INT stce_9808 = stce_9886._tape_index();
    return stce_10290[stce_9808];
  }
  void interpret_adjoint() {
    typename stce_10143::iterator stce_9863 = stce_10151::stce_10146().start();
    internal::interpretation_settings stce_10311;
    stce_10312(stce_10151::stce_10146().current(), stce_9863, stce_10311);
  }
  void interpret_adjoint_to(const typename stce_10143::iterator &stce_9863) {
    internal::interpretation_settings stce_10311;
    if (stce_9863 > stce_10151::stce_10146().current())
      throw dco::exception::create<std::runtime_error>("adjoint interpretation: from < to.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_adjoint_vector.hpp", 59);
    else
      stce_10312(stce_10151::stce_10146().current(), stce_9863, stce_10311);
  }
  void interpret_adjoint_from(const typename stce_10143::iterator &stce_9862) {
    typename stce_10143::iterator stce_9863 = stce_10151::stce_10146().start();
    assert(!(stce_9862 < stce_9863));
    typename internal::interpretation_settings stce_10311;
    stce_10312(stce_9862, stce_9863, stce_10311);
  }
  void interpret_adjoint_from_to(const typename stce_10143::iterator &stce_9862, const typename stce_10143::iterator &stce_9863) {
    typename internal::interpretation_settings stce_10311;
    if (stce_9863 > stce_9862)
      throw dco::exception::create<std::runtime_error>("adjoint interpretation: from < to.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_adjoint_vector.hpp", 75);
    else
      stce_10312(stce_9862, stce_9863, stce_10311);
  }
  void interpret_adjoint_and_zero_adjoints_to(const typename stce_10143::iterator &stce_9863) {
    typename stce_10143::iterator stce_9862(stce_10151::stce_10146().current());
    typename internal::interpretation_settings stce_10311;
    stce_10311.reset = false;
    stce_10311.stce_10065 = true;
    stce_10312(stce_9862, stce_9863, stce_10311);
  }
  void interpret_adjoint_and_zero_adjoints_from_to(const typename stce_10143::iterator &stce_9862, const typename stce_10143::iterator &stce_9863) {
    typename internal::interpretation_settings stce_10311;
    stce_10311.reset = false;
    stce_10311.stce_10065 = true;
    stce_10312(stce_9862, stce_9863, stce_10311);
  }
  void zero_adjoints() {
    typename stce_10143::iterator stce_9863 = stce_10151::stce_10146().start();
    stce_10313(stce_10151::stce_10146().current(), stce_9863);
  }
  void zero_adjoints_to(const typename stce_10143::iterator &stce_9863) {
    stce_10313(stce_10151::stce_10146().current(), stce_9863);
  }
  void zero_adjoints_from(const typename stce_10143::iterator &stce_9862) {
    typename stce_10143::iterator stce_9863 = stce_10151::stce_10146().start();
    stce_10313(stce_9862, stce_9863);
  }
  void zero_adjoints_from_to(const typename stce_10143::iterator &stce_9862, const typename stce_10143::iterator &stce_9863) {
    stce_10313(stce_9862, stce_9863);
  }
  value_t& operator[](DCO_TAPE_INT stce_9808) {
    return stce_10290[stce_9808];
  }
  value_t at(DCO_TAPE_INT stce_9808) const { return (stce_9808 < stce_10308)?stce_10290[stce_9808]:0; }
  value_t& at(DCO_TAPE_INT stce_9808) { stce_10295(stce_9808); return stce_10290[stce_9808]; }
  void stce_10295(DCO_TAPE_INT stce_9808) {
    if(stce_10308 < stce_9808) {
      stce_10308 = stce_9808;
    }
    stce_10310();
  }
  void stce_10296(typename stce_10143::iterator& stce_9817) {
    (*this)[stce_9817->arg()] += stce_9817->pval()*(*this)[stce_9817.index()];
  }
  value_t& stce_10049(const DCO_TAPE_INT stce_9808) {
    return stce_10290[stce_9808];
  }
  stce_9761& stce_10049(const DCO_TAPE_INT stce_9808, int stce_10314) {
    return stce_10290[stce_9808][stce_10314];
  }
  template<class stce_10315>
    void stce_10081(const DCO_TAPE_INT &stce_9808, const stce_10315 &stce_9985, const value_t &stce_10082) {
    stce_10290[stce_9808] += (stce_9985 * stce_10082);
  }
  virtual typename stce_10143::stce_9764 stce_10080(const DCO_TAPE_INT stce_9808) {
    return static_cast<typename stce_10143::stce_9764>(stce_10049(stce_9808)[stce_10309]);
  }
  virtual void stce_10081(const DCO_TAPE_INT stce_9808, const typename stce_10143::stce_9764 &stce_10082) {
    stce_10290[stce_9808][stce_10309] += stce_10082;
  }
  void stce_10313(const typename stce_10143::iterator &stce_9862, const typename stce_10143::iterator &stce_9863) {
    stce_10303(stce_9862.index(), stce_9863.index());
  }
  void stce_10303(DCO_TAPE_INT stce_9862, DCO_TAPE_INT stce_9863) {
    if(stce_9862 < stce_9863) { std::swap(stce_9862, stce_9863); }
    if(stce_10290==0) { return; }
    for(DCO_TAPE_INT stce_9832=stce_9863+1;stce_9832<=stce_9862;stce_9832++) {
      stce_10290[stce_9832] = 0;
    }
  }
  void stce_10312(const typename stce_10143::iterator &stce_9862,
                                   const typename stce_10143::iterator &stce_9863,
                                   typename internal::interpretation_settings stce_10311)
  {
    stce_10310();
    stce_10151::stce_10146().stce_10316(*this, stce_9862, stce_9863, stce_10311);
  }
  void stce_10310() {
    if(stce_10290==0) {
      stce_10308 = stce_10151::stce_10146().current().index()+1;
      stce_10290 = new value_t[size_t(stce_10308)]();
    }
    else {
      DCO_TAPE_INT stce_10317 = stce_10308 = stce_10151::stce_10146().current().index()+1;
      if(stce_10317 > stce_10308) {
        value_t *stce_10318 = new value_t[size_t(stce_10317)]();
        for(DCO_TAPE_INT stce_9832=0;stce_9832<stce_10308;++stce_9832)
          stce_10318[stce_9832] = stce_10290[stce_9832];
        delete [] stce_10290;
        stce_10290 = stce_10318;
      }
    }
  }
};
}
       
namespace dco {
template< typename stce_10319, typename adjoint_vector >
  void stce_10320(stce_10319 stce_9817, stce_10319 stce_9863, const adjoint_vector& stce_10321,
                 const std::string &filename = "tape.dot") {
  stce_10320(stce_9817, stce_9863 , stce_10321, stce_9848(filename));
}
template< typename stce_10319, typename adjoint_vector >
  void stce_10320(stce_10319 stce_9817, stce_10319 stce_9863, const adjoint_vector& stce_10321,
                 stce_9848 &stce_9849) {
  if (stce_9817.index() < stce_9863.index()) {
    std::swap(stce_9817, stce_9863);
  }
  DCO_TAPE_INT stce_10322 = stce_9817.index();
  stce_9849.stce_9858(stce_10322,stce_10321.at(stce_10322));
  for (; stce_9817 != stce_9863; --stce_9817) {
    stce_9849.stce_9861(stce_9817->arg(), stce_10322, stce_9817->pval());
    if(stce_9817->stce_10157()) {
      stce_10322--;
      stce_9849.stce_9858(stce_10322,stce_10321.at(stce_10322));
    }
  }
}
template< typename stce_10319, typename adjoint_vector, typename stce_10323 >
  void stce_10320(stce_10319 stce_9817, stce_10319 stce_9863, const adjoint_vector& stce_10321,
                 const stce_10323& stce_10324, stce_9848 &stce_9849) {
  if (stce_9817.index() < stce_9863.index()) {
    std::swap(stce_9817, stce_9863);
  }
  stce_10320(stce_9817, stce_9863, stce_10321, stce_9849);
  typedef typename stce_10323::const_reverse_iterator stce_10325;
  for (stce_10325 stce_10164 = stce_10324.rbegin(); stce_10164 != stce_10324.rend() && stce_9863 <= stce_10164->stce_10326(); ++stce_10164) {
    if (stce_9817 < stce_10164->stce_10326())
      continue;
    stce_10164->stce_10327()->stce_10088(stce_9849, stce_10164->stce_10326().index() + 1);
  }
}
template <typename stce_10319>
  void stce_10328(stce_10319 stce_9817, stce_10319 stce_9863,
                 std::string const& filename) {
  if (stce_9817.index() < stce_9863.index()) { std::swap(stce_9817, stce_9863); }
  std::ofstream stce_9919(filename.c_str());
  stce_9919 << stce_9817.index() << "\n";
  for (; stce_9817 != stce_9863; --stce_9817) {
    stce_9919 << stce_9817->arg() << ", " << stce_9817->pval() << "; " ;
    if(stce_9817->stce_10157())
      stce_9919 << "\n";
  }
}
template <bool USE_SPARSE_INTERPRET, typename stce_10319, typename adjoint_vector>
  DCO_STRONG_INLINE void interpret(stce_10319 stce_9817, const stce_10319& stce_9863,
                                            adjoint_vector& stce_10321) {
  stce_9817.index()--;
  stce_1000002();;
  stce_10321.stce_10295(stce_9810(stce_9752(stce_9817.index())));
  stce_10321.init();
  for (; stce_9817 != stce_9863; --stce_9817) {
    if(USE_SPARSE_INTERPRET && folding::is_zero(stce_10321[stce_9817.index()])) { continue; }
    DCO_LOG(logDEBUG3) << stce_10321[stce_9817->arg()] << " <=> adj[" << stce_9817->arg() << "] += " << stce_9817->pval() * stce_10321[stce_9817.index()] << " from idx " << stce_9817.index();
    stce_10321.stce_10296(stce_9817);
  }
}
template <bool USE_SPARSE_INTERPRET, typename stce_10319, typename adjoint_vector>
  DCO_STRONG_INLINE void stce_10329(stce_10319 stce_9817,
                                           const stce_10319& stce_9863,
                                           adjoint_vector& stce_10321) {
  stce_9817.index()--;
  stce_1000002();;
  stce_10321.stce_10295(stce_9810(stce_9752(stce_9817.index())));
  for (; stce_9817 != stce_9863; --stce_9817) {
    if(USE_SPARSE_INTERPRET && folding::is_zero(stce_10321[stce_9817.index()])) { continue; }
    DCO_LOG(logDEBUG3) << stce_10321[stce_9817->arg()] << " <=> adj[" << stce_9817->arg() << "] += " << stce_9817->pval() * stce_10321[stce_9817.index()] << " from idx " << stce_9817.index();
    stce_10321.stce_10296(stce_9817);
    if(stce_9817->stce_10157()) {
      stce_10321[stce_9817.index()] = 0;
    }
  }
}
template< typename stce_10319, typename adjoint_vector >
  DCO_STRONG_INLINE void stce_10330(stce_10319 stce_9817, const stce_10319& stce_9863, adjoint_vector& stce_10321) {
  stce_10321.stce_10295(stce_9817.index());
  for (; stce_9817 != stce_9863; --stce_9817) {
    std::cout << stce_9817.index() << ": " << stce_9817->arg() << ", " << stce_9817->pval() << ", " << stce_10321[stce_9817.index()] << std::endl;
  }
}
template< typename tape_t >
  DCO_STRONG_INLINE void stce_10330(const tape_t *tape) {
  typename tape_t::iterator stce_9817 = tape->get_position();
  for (; stce_9817 != tape->begin(); --stce_9817) {
    std::cout << stce_9817.index() << ": " << stce_9817->arg() << ", " << stce_9817->pval() << std::endl;
  }
}
}
       
       
       
namespace dco {
namespace internal {
template <typename stce_10331> struct stce_10332 {
  typedef stce_10331 type;
  typedef stce_10331& stce_10333;
  typedef const stce_10331 stce_10334;
};
template <typename stce_9782 > struct stce_10332<std::vector<stce_9782> > {
  typedef std::vector<stce_9782> & type;
  typedef std::vector<stce_9782> & stce_10333;
  typedef std::vector<stce_9782> const& stce_10334;
};
template <stce_9796 stce_10335, typename stce_9782, typename stce_10336 = void>
  struct stce_10337 { };
template <typename stce_9782> struct stce_10337<stce_9797, stce_9782> {
  typedef typename dco::mode<stce_9782>::value_t value_t;
  typedef typename trait_value<stce_9782>::stce_9802 stce_9802;
  typedef typename trait_value<stce_9782>::stce_10029 stce_10029;
  stce_9802 get( stce_9782& stce_9843) { return trait_value<stce_9782>::value(stce_9843); };
  stce_10029 get(const stce_9782& stce_9843) const { return trait_value<stce_9782>::value(stce_9843); };
};
template <typename stce_9782> struct stce_10337<stce_9800, stce_9782> {
  typedef typename dco::mode<stce_9782>::passive_t value_t;
  typedef typename stce_9785<stce_9782>::stce_9802 stce_9802;
  typedef typename stce_9785<stce_9782>::stce_10029 stce_10029;
  stce_9802 get( stce_9782& stce_9843) { return stce_9785<stce_9782>::value(stce_9843); };
  stce_10029 get(const stce_9782& stce_9843) const { return stce_9785<stce_9782>::value(stce_9843); };
};
template <typename stce_9782> struct stce_10337<stce_9798, stce_9782> {
  typedef typename dco::mode<stce_9782>::derivative_t value_t;
  typedef typename stce_9784<stce_9782>::stce_9802 stce_9802;
  typedef typename stce_9784<stce_9782>::stce_9802 stce_10029;
  stce_10029 get(const stce_9782& stce_9843) const { return stce_9784<stce_9782>::value(stce_9843); };
};
template <typename stce_9782> struct stce_10337<stce_9801, stce_9782> {
  typedef typename stce_9787<stce_9782>::stce_10338 stce_9802;
  typedef typename stce_9787<stce_9782>::stce_10029 stce_10029;
  typedef stce_9802 value_t;
  stce_9802 get( stce_9782& stce_9843) { return stce_9787<stce_9782>::value(stce_9843); };
  stce_10029 get(const stce_9782& stce_9843) const { return stce_9787<stce_9782>::value(stce_9843); };
};
template <typename stce_9782, typename stce_10336> struct stce_10337<stce_9799, stce_9782, stce_10336> {
  typedef typename dco::mode<stce_9782>::derivative_t value_t;
  typedef typename stce_9784<stce_9782, stce_10336>::stce_10339 stce_9802;
  typedef typename stce_9784<stce_9782, stce_10336>::stce_10339 stce_10029;
  stce_10336* stce_10290;
  stce_10337(stce_10336* adjoints) : stce_10290(adjoints) {}
  stce_9802 get( stce_9782& stce_9843) { return stce_9784<stce_9782, stce_10336>::value(stce_9843, stce_10290); };
  stce_10029 get(const stce_9782& stce_9843) const { return stce_9784<stce_9782, stce_10336>::value(stce_9843, stce_10290); };
};
template<stce_9796 stce_10335, class stce_9782, typename stce_10331, typename stce_9783 = void>
  struct stce_10340 : stce_10337<stce_10335, stce_9782, stce_9783> {
  typedef typename stce_10332<stce_10331>::type stce_10283;
  typedef typename stce_10332<stce_10331>::stce_10333 stce_10341;
  typedef typename stce_10332<stce_10331>::stce_10334 stce_10342;
  typedef stce_10337<stce_10335, stce_9782, stce_9783> stce_10343;
  typedef typename stce_10343::value_t value_t;
  stce_10283 stce_10344;
  stce_10340(stce_10342 stce_9842) : stce_10344(const_cast<stce_10341>(stce_9842)) {}
  stce_10340(stce_10342 stce_9842, stce_9783* stce_10345) : stce_10343(stce_10345), stce_10344(const_cast<stce_10341>(stce_9842)) {}
  stce_10340 &operator=(const std::vector<value_t> &stce_9894) {
    for (size_t stce_9832 = 0; stce_9832 < stce_10344.size(); ++stce_9832) {
      stce_10343::get(stce_10344[stce_9832]) = stce_9894[stce_9832];
    }
    return *this;
  }
  operator std::vector<value_t>() const {
    std::vector<value_t> stce_9840(stce_10344.size());
    for (size_t stce_9832 = 0; stce_9832 < stce_10344.size(); ++stce_9832) {
      stce_9840[stce_9832] = stce_10343::get(stce_10344[stce_9832]);
    }
    return stce_9840;
  }
  typename stce_10343::stce_9802 operator[](size_t stce_9832) { return stce_10343::get(stce_10344[stce_9832]); }
  typename stce_10343::stce_10029 operator[](size_t stce_9832) const {
    return stce_10343::get(stce_10344[stce_9832]);
  }
  size_t size() const { return stce_10344.size(); }
};
}
}
namespace dco {
template <typename stce_9782> inline typename internal::trait_value<stce_9782>::stce_9802
  value( stce_9782 &stce_9843) { return internal::trait_value<stce_9782>::value(stce_9843); }
template <typename stce_9782> inline typename internal::trait_value<stce_9782>::stce_10029
  value(const stce_9782 &stce_9843) { return internal::trait_value<stce_9782>::value(stce_9843); }
namespace internal {
template <typename stce_9782, typename enable_if>
  struct trait_value {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef type& stce_9802;
  typedef const type& stce_10029;
  static inline stce_9802 value( type &value) { return value; }
  static inline stce_10029 value(const type &value) { return value; }
};
template <typename stce_9782>
  struct trait_value<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type && !is_gbcp_type<stce_9782>::value && !dco::mode<stce_9782>::is_intermediate_type>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::value_t& stce_9802;
  typedef typename dco::mode<type>::value_t const& stce_10029;
  static inline stce_9802 value( type &value) { return value._value(); }
  static inline stce_10029 value(const type &value) { return value._value(); }
};
template <typename stce_9782>
  struct trait_value<stce_9782, typename dco::enable_if<!is_gbcp_type<stce_9782>::value && dco::mode<stce_9782>::is_intermediate_type>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::value_t const& stce_9802;
  typedef typename dco::mode<type>::value_t const& stce_10029;
  static inline stce_9802 value( type &value) { return value._value(); }
  static inline stce_10029 value(const type &value) { return value._value(); }
};
template <typename stce_9782>
  struct trait_value<stce_9782, typename dco::enable_if<is_gbcp_type<stce_9782>::value>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::value_t& stce_9802;
  typedef typename dco::mode<type>::value_t const& stce_10029;
  static inline stce_9802 value( type &stce_10346) { return dco::value(stce_10346._value()); }
  static inline stce_10029 value(const type &stce_10346) { return dco::value(stce_10346._value()); }
};
template <typename stce_9782>
  struct trait_value<std::vector<stce_9782>, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type>::type> {
  typedef std::vector<stce_9782> stce_10283;
  typedef stce_10340<stce_9797, stce_9782, stce_10283> stce_9802;
  typedef const stce_10340<stce_9797, stce_9782, stce_10283> stce_10029;
  static inline stce_9802 value( stce_10283 &stce_9842) { return stce_9802(stce_9842); }
  static inline stce_10029 value(const stce_10283 &stce_9842) { return stce_10029(stce_9842); }
};
template<stce_9796 stce_10335, class stce_9782, typename stce_10347, typename stce_10348>
  struct trait_value<stce_10340<stce_10335, stce_9782, stce_10347, stce_10348>, void > {
  typedef stce_10340<stce_10335, stce_9782, stce_10347, stce_10348> stce_10283;
  typedef stce_10340<stce_9797, typename stce_10283::value_t, stce_10283> stce_9802;
  typedef stce_10340<stce_9797, typename stce_10283::value_t, stce_10283> stce_10029;
  static inline stce_9802 value(const stce_10283 &stce_9842) { return stce_9802(stce_9842); }
};
}
}
       
namespace dco {
template <typename stce_9782> inline typename internal::stce_9785<stce_9782>::stce_9802
  passive_value( stce_9782 &stce_9843) { return internal::stce_9785<stce_9782>::value(stce_9843); }
template <typename stce_9782> inline typename internal::stce_9785<stce_9782>::stce_10029
  passive_value(const stce_9782 &stce_9843) { return internal::stce_9785<stce_9782>::value(stce_9843); }
namespace internal {
template <typename stce_9782, typename enable_if>
  struct stce_9785 {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef type& stce_9802;
  typedef type const& stce_10029;
  static inline stce_9802 value( type &value) { return value; }
  static inline stce_10029 value(const type &value) { return value; }
};
template <typename stce_9782>
  struct stce_9785<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type && !dco::mode<stce_9782>::is_intermediate_type>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::active_value_t active_value_t;
  typedef typename dco::mode<type>::value_t value_t;
  typedef typename dco::mode<type>::passive_t passive_t;
  typedef passive_t& stce_9802;
  typedef passive_t const& stce_10029;
  static inline stce_9802 value( type &value) { return stce_9785<value_t>::value(dco::value(value)); }
  static inline stce_10029 value(const type &value) { return stce_9785<value_t>::value(dco::value(value)); }
};
template <typename stce_9782>
  struct stce_9785<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_intermediate_type>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::active_value_t active_value_t;
  typedef typename dco::mode<type>::value_t value_t;
  typedef typename dco::mode<type>::passive_t passive_t;
  typedef passive_t const& stce_9802;
  typedef passive_t const& stce_10029;
  static inline stce_9802 value( type &value) { return stce_9785<value_t>::value(dco::value(value)); }
  static inline stce_10029 value(const type &value) { return stce_9785<value_t>::value(dco::value(value)); }
};
template <typename stce_9782>
  struct stce_9785<std::vector<stce_9782>, typename dco::enable_if<!dco::mode<stce_9782>::is_dco_type>::type> {
  typedef std::vector<stce_9782> stce_10283;
  typedef stce_10283& stce_9802;
  typedef const stce_10283& stce_10029;
  static inline stce_9802 value( stce_10283 &value) { return value; }
  static inline stce_10029 value(const stce_10283 &value) { return value; }
};
template <typename stce_9782>
  struct stce_9785<std::vector<stce_9782>, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type>::type> {
  typedef std::vector<stce_9782> stce_10283;
  typedef stce_10340<stce_9800, stce_9782, stce_10283> stce_9802;
  typedef const stce_10340<stce_9800, stce_9782, stce_10283> stce_10029;
  static inline stce_9802 value( stce_10283 &stce_9842) { return stce_9802(stce_9842); }
  static inline stce_10029 value(const stce_10283 &stce_9842) { return stce_10029(stce_9842); }
};
template<stce_9796 stce_10335, class stce_9782, typename stce_10347, typename stce_10348>
  struct stce_9785<stce_10340<stce_10335, stce_9782, stce_10347, stce_10348>, void > {
  typedef stce_10340<stce_10335, stce_9782, stce_10347, stce_10348> stce_10283;
  typedef stce_10340<stce_9800, typename stce_10283::value_t, stce_10283> stce_9802;
  typedef stce_10340<stce_9800, typename stce_10283::value_t, stce_10283> stce_10029;
  static inline stce_9802 value(const stce_10283 &stce_9842) { return stce_9802(stce_9842); }
};
}
}
       
namespace dco {
template <typename stce_9782> inline typename internal::stce_9784<stce_9782>::stce_9802
  derivative(const stce_9782 &stce_9843) { return internal::stce_9784<stce_9782>::value(stce_9843); }
template <typename stce_9782, typename stce_9783> inline typename internal::stce_9784<stce_9782, stce_9783>::stce_10339
  derivative(const stce_9782 &stce_9843, stce_9783& adjoints) { return internal::stce_9784<stce_9782, stce_9783>::value(stce_9843, &adjoints); }
namespace internal {
template <typename stce_9782, typename stce_9783, typename enable_if>
  struct stce_9784 {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef type stce_9802;
  static inline stce_9802 value(const type &) { return type(); }
  static inline stce_9802 value(const type &, stce_9783*) { return type(); }
};
template <typename stce_9782>
  struct stce_9784<stce_9782, void, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type && !is_gbcp_type<stce_9782>::value>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::derivative_t& stce_9802;
  static inline stce_9802 value(type const& value) { return value.stce_9989(); }
};
template <typename stce_9782>
  struct stce_9784<stce_9782, void, typename dco::enable_if<(is_gbcp_type<stce_9782>::value && dco::mode<stce_9782>::order > 0)>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::value_t& stce_9802;
  static inline stce_9802 value(const type &value) { return dco::derivative(value._value()); }
};
template <typename stce_9782>
  struct stce_9784<stce_9782, void, typename dco::enable_if<(is_gbcp_type<stce_9782>::value && dco::mode<stce_9782>::order == 0)>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::value_t stce_9802;
  static inline stce_9802 value(const type &) { return stce_9802(); }
};
template <typename stce_9782, typename stce_10349, typename stce_10350, int stce_10351>
  struct stce_9784<stce_9782, adjoint_vector<stce_10349, stce_10350, stce_10351>, typename dco::enable_if<(stce_10351 == -1)>::type> {
  typedef adjoint_vector<stce_10349, stce_10350, stce_10351> stce_10348;
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef stce_10350& stce_10339;
  static inline stce_10339 value(const type &value, stce_10348* adjoints) { return adjoints->derivative(value)[0]; }
};
template <typename stce_9782, typename stce_10349, typename stce_10350, int stce_10351>
  struct stce_9784<stce_9782, adjoint_vector<stce_10349, stce_10350, stce_10351>, typename dco::enable_if<(stce_10351 > 0)>::type> {
  typedef adjoint_vector<stce_10349, stce_10350, stce_10351> stce_10348;
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename stce_10348::value_t& stce_10339;
  static inline stce_10339 value(const type &value, stce_10348* adjoints) { return adjoints->derivative(value); }
};
template <typename stce_9782>
  struct stce_9784<std::vector<stce_9782>, void, typename dco::enable_if<!dco::mode<stce_9782>::is_dco_type>::type> {
  typedef std::vector<stce_9782> stce_10283;
  typedef stce_10283 stce_9802;
  static inline stce_9802 value(const stce_10283 &value) { return stce_9802(value.size(), 0.0); }
};
template <typename stce_9782>
  struct stce_9784<std::vector<stce_9782>, void, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type>::type> {
  typedef std::vector<stce_9782> stce_10283;
  typedef stce_10340<stce_9798, stce_9782, stce_10283> stce_9802;
  static inline stce_9802 value(const stce_10283 &stce_9842) { return stce_9802(stce_9842); }
};
template <typename stce_9782, typename stce_10349, typename stce_10350, int stce_10351>
  struct stce_9784<std::vector<stce_9782>, adjoint_vector<stce_10349, stce_10350, stce_10351>, typename dco::enable_if<(stce_10351 == -1)>::type> {
  typedef adjoint_vector<stce_10349, stce_10350, stce_10351> stce_10348;
  typedef std::vector<stce_9782> stce_10283;
  typedef stce_10340<stce_9798, stce_9782, stce_10283> stce_9802;
  typedef stce_10340<stce_9799, stce_9782, stce_10283, stce_10348> stce_10339;
  static inline stce_9802 value(const stce_10283 &stce_9842) { return stce_9802(stce_9842); }
  static inline stce_10339 value(const stce_10283 &stce_9842, stce_10348* adjoints) { return stce_10339(stce_9842, adjoints); }
};
template<stce_9796 stce_10335, class stce_9782, typename stce_10347, typename stce_10348>
  struct stce_9784<stce_10340<stce_10335, stce_9782, stce_10347, stce_10348>, void > {
  typedef stce_10340<stce_10335, stce_9782, stce_10347, stce_10348> stce_10283;
  typedef stce_10340<stce_9798, typename stce_10283::value_t, stce_10283> stce_9802;
  static inline stce_9802 value(const stce_10283 &stce_9842) { return stce_9802(stce_9842); }
};
}
}
       
namespace dco {
template <typename stce_9782>
  typename internal::stce_9786<stce_9782>::stce_9802 tape(const stce_9782 &stce_9843) {
  return internal::stce_9786<stce_9782>::value(stce_9843);
}
template <typename stce_9782>
  typename internal::stce_9787<stce_9782>::stce_10029 tape_index(const stce_9782 &stce_9843) {
  return internal::stce_9787<stce_9782>::value(stce_9843);
}
template <typename stce_9782>
  typename internal::stce_9787<stce_9782>::stce_9802 tape_index(stce_9782 &stce_9843) {
  return internal::stce_9787<stce_9782>::value(stce_9843);
}
namespace internal {
template<typename stce_9782, typename enable_if> struct stce_9786 {
  typedef void* stce_9802;
  static stce_9802 value(const stce_9782&) { return static_cast<void *>(0); }
};
template<typename stce_9782> struct stce_9786<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_adjoint_type>::type> {
  typedef typename dco::mode<stce_9782>::tape_t* stce_9802;
  static stce_9802 value(const stce_9782& value) { return value.tape(); }
};
template<typename stce_9782> struct stce_9786<std::vector<stce_9782>, typename dco::enable_if<dco::mode<stce_9782>::is_adjoint_type>::type> {
  typedef typename dco::mode<stce_9782>::tape_t* stce_9802;
  static stce_9802 value(const std::vector<stce_9782>& value) {
    stce_9802 tape = 0;
    for (size_t stce_9832 = 0; stce_9832 < value.size(); stce_9832++) {
      if (stce_9786<stce_9782>::value(value[stce_9832]) != 0) {
        tape = stce_9786<stce_9782>::value(value[stce_9832]);
        break;
      }
    }
    return tape;
  }
};
template<typename stce_9782, typename enable_if> struct stce_9787 {
  typedef DCO_TAPE_INT stce_9802;
  typedef DCO_TAPE_INT stce_10338;
  typedef const DCO_TAPE_INT stce_10029;
  static DCO_TAPE_INT value(const stce_9782&) { return 0; }
};
template<typename stce_9782> struct stce_9787<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_adjoint_type>::type> {
  typedef DCO_TAPE_INT& stce_9802;
  typedef DCO_TAPE_INT stce_10338;
  typedef const DCO_TAPE_INT stce_10029;
  static DCO_TAPE_INT value(const stce_9782& value) { return value._tape_index(); }
  static DCO_TAPE_INT& value(stce_9782& value) { return value._tape_index(); }
};
template <typename stce_9782>
  struct stce_9787<std::vector<stce_9782>, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type>::type> {
  typedef std::vector<stce_9782> stce_10283;
  typedef stce_10340<stce_9801, stce_9782, stce_10283> stce_9802;
  typedef const stce_10340<stce_9801, stce_9782, stce_10283> stce_10029;
  static inline stce_9802 value( stce_10283 &stce_9842) { return stce_9802(stce_9842); }
  static inline stce_10029 value(const stce_10283 &stce_9842) { return stce_10029(stce_9842); }
};
}
}
namespace dco {
template <typename stce_9920, typename... stce_10352>
  int fprintf(stce_9920& stream, char const* stce_10353, stce_10352&&... stce_10354) {
  return ::fprintf(stream, stce_10353, dco::passive_value(stce_10354)...);
}
template <typename... stce_10352>
  int printf(char const * const stce_10353, stce_10352&&... stce_10354) {
  return ::printf(stce_10353, dco::passive_value(stce_10354)...);
}
template <typename... stce_10352>
  int sprintf(char* stream, char const* stce_10353, stce_10352&&... stce_10354) {
  return ::sprintf(stream, stce_10353, dco::passive_value(stce_10354)...);
}
namespace stce_10355 {
inline void init() {}
template <typename stce_9782, typename... stce_10352>
  typename dco::enable_if<dco::mode<typename std::remove_pointer<stce_9782>::type>::is_dco_type>::type init(stce_9782 &stce_9843, stce_10352&... stce_10354) {
  *stce_9843 = 0;
  init(stce_10354...);
}
template <typename stce_9782, typename... stce_10352>
  typename dco::enable_if<!dco::mode<typename std::remove_pointer<stce_9782>::type>::is_dco_type>::type init(stce_9782&, stce_10352&... stce_10354) {
  init(stce_10354...);
}
}
template <typename stce_9920, typename... stce_10352>
  int sscanf(stce_9920& stream, char const* stce_10353, stce_10352&&... stce_10354) {
  stce_10355::init(stce_10354...);
  return ::sscanf(stream, stce_10353, &dco::passive_value(*stce_10354)...);
}
template <typename stce_9920, typename... stce_10352>
  int fscanf(stce_9920& stream, char const* stce_10353, stce_10352&&... stce_10354) {
  stce_10355::init(stce_10354...);
  return ::fscanf(stream, stce_10353, &dco::passive_value(*stce_10354)...);
}
}
       
namespace dco {
template <typename tape_t>
struct stce_10356 {
  typedef std::map<std::string, mem_long_t> stce_10357;
  stce_10357 stce_10358;
  tape_t* stce_9952;
  stce_10356(tape_t* tape) : stce_9952(tape) {
    stce_10358["max_tape_memory"] = 0;
    stce_10358["max_adjoint_memory"] = 0;
    stce_10358["max_checkpoint_memory"] = 0;
    stce_10358["max_total_memory"] = 0;
  }
  mem_long_t operator[](std::string stce_10359) { return stce_10358[stce_10359]; }
  void stce_10360() {
    mem_long_t stce_10361 = dco::size_of(stce_9952, tape_t::size_of_stack);
    if (stce_10361 > stce_10358["max_tape_memory"]) {
      stce_10358["max_tape_memory"] = stce_10361;
    }
    stce_10361 = dco::size_of(stce_9952, tape_t::size_of_internal_adjoint_vector);
    if (stce_10361 > stce_10358["max_adjoint_memory"]) {
      stce_10358["max_adjoint_memory"] = stce_10361;
    }
    stce_10361 = dco::size_of(stce_9952, tape_t::size_of_checkpoints);
    if (stce_10361 > stce_10358["max_checkpoint_memory"]) {
      stce_10358["max_checkpoint_memory"] = stce_10361;
    }
    stce_10358["max_total_memory"] = 0;
    for (stce_10357::iterator stce_9832 = stce_10358.begin(); stce_9832 != stce_10358.end(); ++stce_9832) {
      stce_10358["max_total_memory"] += stce_9832->second;
    }
  }
  mem_long_t stce_10362() { return stce_10358["max_tape_memory"]; }
};
}
namespace dco {
namespace internal {
template<memory_model::TYPE MEMORY_MODEL,
         typename stce_10363,
         memory_model::TYPE stce_9789>
  class tape : public dco::stce_10079<typename helper::stce_9836<typename stce_10363::stce_9764>::stce_9838>,
               public stce_10149<MEMORY_MODEL,
                                     tape<MEMORY_MODEL, stce_10363, stce_9789>,
                                     stce_9789> {
public:
  typedef stce_10363 stce_9762;
  typedef stce_10149<MEMORY_MODEL, tape, stce_9789> stce_10283;
  typedef typename stce_10283::iterator iterator_t;
  typedef iterator_t position_t;
  typedef typename stce_9762::stce_9764 stce_9764;
  typedef typename stce_9762::stce_9763 stce_9763;
  typedef typename helper::stce_9836<typename stce_10363::stce_9764>::stce_9838 stce_10085;
  typedef tape<MEMORY_MODEL, stce_9762, stce_9789> tape_t;
  typedef stce_9793<typename stce_9762::stce_9764, stce_10283> vector_t;
  vector_t stce_10290;
  stce_10356<tape> stce_10364;
  stce_10356<tape>& stce_10365() { return stce_10364; }
  const stce_10356<tape>& stce_10365() const { return stce_10364; }
  tape& get_tape() { return *this; }
private:
  int stce_10366;
  bool stce_10367;
  bool stce_10368;
  DCO_TAPE_INT stce_10309;
  tape_options stce_10197;
  bool stce_10369;
  void stce_10313(const iterator_t &stce_9862, const iterator_t &stce_9863 = stce_10283::start()) {
    stce_10290.stce_10303(stce_9862.index(), stce_9863.index());
  }
  void stce_10370() {
    stce_10290.resize(this->stce_10185());
    stce_10290.stce_10306();
  }
  template <typename stce_10319, typename adjoint_vector>
    DCO_STRONG_INLINE void stce_10371(stce_10319 stce_9817, const stce_10319& stce_9863,
                                      adjoint_vector& stce_10321) {
    if (stce_10369) { dco::interpret< true >(stce_9817, stce_9863, stce_10321); }
    else { dco::interpret< false>(stce_9817, stce_9863, stce_10321); }
  }
  template <typename stce_10319, typename adjoint_vector>
    DCO_STRONG_INLINE void stce_10372(stce_10319 stce_9817, const stce_10319& stce_9863,
                                      adjoint_vector& stce_10321) {
    if (stce_10369) { dco::stce_10329< true >(stce_9817, stce_9863, stce_10321); }
    else { dco::stce_10329< false>(stce_9817, stce_9863, stce_10321); }
  }
public:
  const tape_options& stce_10170() const { return stce_10197; }
  bool stce_10373() { return this->start() == this->current(); }
  tape(tape_options stce_10170 = tape_options()) :
    stce_10283(stce_10170, *this), stce_10290(1, *this, stce_10170),
    stce_10364(this),
    stce_10366(dco::build_info::COMPATIBILITY_VERSION),
    stce_10367(false), stce_10368(true), stce_10309(0), stce_10197(stce_10170),
    stce_10369(false) { }
  ~tape() {
    for(typename std::vector<stce_10374>::iterator stce_10164 = stce_10375.begin();
        stce_10164 != stce_10375.end(); ++stce_10164) {
      stce_10164->stce_10376();
    }
  }
  bool is_compatible(const int stce_10377) {
    return stce_10377 == stce_10366;
  }
  bool is_active() { return stce_9758 ? stce_10368 : true; }
  typedef helper::stce_10084<tape> stce_10378;
  static void stce_10379(stce_10378* stce_10380) { stce_10380->get_tape()->stce_10369 = true; }
  static void stce_10381(stce_10378* stce_10380) { stce_10380->get_tape()->stce_10369 = false; }
  void set_sparse_interpret() {
    if (stce_10369 == false) {
      stce_10378* stce_10380 = this->create_callback_object<stce_10378>();
      this->insert_callback(stce_10381, stce_10380);
    }
    stce_10369 = true;
  }
  void unset_sparse_interpret() {
    if (stce_10369 == true) {
      stce_10378* stce_10380 = this->create_callback_object<stce_10378>();
      this->insert_callback(stce_10379, stce_10380);
    }
    stce_10369 = false;
  }
  bool& sparse_interpret() { return stce_10369; }
  void switch_to_active() {
    if (!stce_9758) {
      throw dco::exception::create<std::runtime_error>("activity check disabled but switched", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 160);
    }
    if (!stce_10368) { stce_10368 = true; }
  }
  void switch_to_passive() {
    if (!stce_9758) {
      throw dco::exception::create<std::runtime_error>("activity check disabled but switched", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 167);
    }
    if (stce_10368) { stce_10368 = false; }
  }
  template <bool stce_10382>
  DCO_TAPE_INT stce_10383() {
    if(!stce_10382 && this->empty()) { return stce_10283::stce_9928().stce_10168(); }
    else {
      typename stce_10283::stce_9924 stce_10384(stce_10283::stce_9928());
      stce_10384.insert(0, 0);
      return stce_10384.index();
    }
  }
  template <typename stce_9931> void register_variable(stce_9931& stce_9843) {
    stce_9931::register_variable(stce_9843, stce_10383< false>(), this);
  }
  template <typename stce_9931> void register_variable(std::vector<stce_9931>& stce_9843) {
    this->register_variable(stce_9843.data(), stce_9843.size());
  }
  template <typename stce_10385> void register_variable(stce_10385 first, stce_10385 stce_10386) {
    while(first != stce_10386) {
      this->register_variable(*first++);
    }
  }
  template <typename stce_10385, typename stce_10387>
    void register_variable(stce_10385 stce_10388, stce_10385 stce_10389, stce_10387 stce_10390) {
    while(stce_10388 != stce_10389) {
      *stce_10390 = *stce_10388++;
      this->register_variable(*stce_10390++);
    }
  }
  template <typename stce_9931> void register_variable(stce_9931* stce_9843, size_t stce_10102) {
    for (size_t stce_9832 = 0; stce_9832 < stce_10102; ++stce_9832) {
      stce_9931::register_variable(stce_9843[stce_9832], stce_10383< false>(), this);
    }
  }
  template<typename stce_9931> void register_output_variable(stce_9931 &stce_9843) { stce_9843 = 1*stce_9843; }
  template<typename stce_9931> void register_output_variable(std::vector<stce_9931>& stce_9843) {
    this->register_output_variable(stce_9843.data(), stce_9843.size());
  }
  template<typename stce_9931> void register_output_variable(stce_9931* stce_9843, size_t stce_10102) {
    for (size_t stce_9832 = 0; stce_9832 < stce_10102; ++stce_9832) {
      this->register_output_variable(stce_9843[stce_9832]);
    }
  }
  template <typename stce_10385> void register_output_variable(stce_10385 first, stce_10385 stce_10386) {
    while(first != stce_10386) {
      this->register_output_variable(*first++);
    }
  }
  template <stce_9796 stce_10391, typename stce_9782, typename stce_10331>
    void register_output_variable(stce_10340<stce_10391, stce_9782, stce_10331> stce_9843) {
    for (size_t stce_9832 = 0; stce_9832 < stce_9843.size(); ++stce_9832) {
      this->register_output_variable(stce_9843[stce_9832]);
    }
  }
  void reset_to(const iterator_t &stce_9863) {
    stce_10265(stce_9863);
    stce_10392(stce_9863);
  }
  void reset() { reset_to(stce_10283::start()); }
  inline iterator_t current() const { return get_position(); }
  inline iterator_t begin() const { return stce_10283::start(); }
  vector_t& adjoints() { return stce_10290; }
  void interpret_adjoint(interpretation_settings stce_10311 = interpretation_settings()) {
    iterator_t stce_9863 = stce_10283::start();
    stce_10312(get_position(), stce_9863, stce_10311);
  }
  void interpret_adjoint_to(const iterator_t &stce_9863) {
    interpretation_settings stce_10311;
    stce_10312(get_position(), stce_9863, stce_10311);
  }
  void interpret_adjoint_from(const iterator_t &stce_9862) {
    iterator_t stce_9863 = stce_10283::start();
    interpretation_settings stce_10311;
    stce_10312(stce_9862, stce_9863, stce_10311);
  }
  void interpret_adjoint_from_to(const iterator_t &stce_9862, const iterator_t &stce_9863) {
    interpretation_settings stce_10311;
    stce_10312(stce_9862, stce_9863, stce_10311);
  }
  void interpret_adjoint_and_reset_to(const iterator_t &stce_9863) {
    iterator_t stce_9862(get_position());
    interpretation_settings stce_10311;
    stce_10311.reset = true;
    stce_10311.stce_10065 = true;
    stce_10312(stce_9862, stce_9863, stce_10311);
    stce_10265(stce_9863);
    stce_10392(stce_9863);
  }
  void interpret_adjoint_and_zero_adjoints_to(const iterator_t &stce_9863) {
    iterator_t stce_9862(get_position());
    interpretation_settings stce_10311;
    stce_10311.reset = false;
    stce_10311.stce_10065 = true;
    stce_10312(stce_9862, stce_9863, stce_10311);
  }
  void interpret_adjoint_and_zero_adjoints_from_to(const iterator_t &stce_9862, const iterator_t &stce_9863) {
    interpretation_settings stce_10311;
    stce_10311.reset = false;
    stce_10311.stce_10065 = true;
    stce_10312(stce_9862, stce_9863, stce_10311);
  }
  void zero_adjoints() {
    iterator_t stce_9863 = stce_10283::start();
    stce_10313(get_position(), stce_9863);
  }
  void zero_adjoints_to(const iterator_t &stce_9863) {
    stce_10313(get_position(), stce_9863);
  }
  void zero_adjoints_from(const iterator_t &stce_9862) {
    iterator_t stce_9863 = stce_10283::start();
    stce_10313(stce_9862, stce_9863);
  }
  void zero_adjoints_from_to(const iterator_t &stce_9862, const iterator_t &stce_9863) {
    stce_10313(stce_9862, stce_9863);
  }
  void stce_10265(const iterator_t &stce_9863) {
    bool alloc = stce_10290.stce_10291(stce_9863.index());
    stce_10290.resize(stce_9863.index()+1);
    if(alloc) { stce_10290.stce_10306(); }
    stce_10283::erase(stce_9863);
  }
public:
  iterator_t get_position() const {
    return stce_10283::current();
  }
  stce_10085 &stce_10393(const DCO_TAPE_INT tape_index) {
    stce_9764 &stce_10082 = _adjoint(tape_index);
    stce_10085 *stce_10394 = reinterpret_cast<stce_10085*>(&stce_10082);
    return stce_10394[stce_10309];
  }
  stce_10085 stce_10393(const DCO_TAPE_INT tape_index) const {
    stce_9764 stce_10082(_adjoint(tape_index));
    stce_10085 *stce_10394 = reinterpret_cast<stce_10085 *>(&stce_10082);
    return stce_10394[stce_10309];
  }
  stce_9764 &_adjoint(DCO_TAPE_INT tape_index) {
    if(!stce_10290.stce_10291(get_position().index()+1))
      stce_10370();
    return stce_10290[tape_index];
  }
  stce_9764 _adjoint(const DCO_TAPE_INT tape_index) const {
    return stce_10290.stce_10291(tape_index)?stce_10290[tape_index]: 0;
  }
  enum stce_10395 {
    size_of_stack = 1,
    size_of_allocated_stack = 2,
    size_of_internal_adjoint_vector = 4,
    size_of_checkpoints = 8,
    stce_10076 = size_of_stack | size_of_internal_adjoint_vector
  };
  mem_long_t stce_10396() const {
    mem_long_t stce_10397 = 0;
    for (size_t stce_9832 = 0; stce_9832 < stce_10375.size(); stce_9832++) {
      stce_10397 += stce_10375[stce_9832].stce_10327()->size_in_byte();
    }
    return stce_10397;
  }
  mem_long_t stce_10077(const int stce_10398 = stce_10076) const {
    mem_long_t stce_10399 = 0;
    if (stce_10398 & size_of_allocated_stack) {
      stce_10399 += static_cast<size_t>(stce_10283::stce_10182());
    }
    else if (stce_10398 & size_of_stack) {
      stce_10399 += static_cast<size_t>(stce_10283::size_in_byte());
    }
    if (stce_10398 & size_of_internal_adjoint_vector) {
      stce_10399 += stce_10290.stce_10301();
    }
    if (stce_10398 & size_of_checkpoints) {
      stce_10399 += stce_10396();
    }
    return stce_10399;
  }
  static tape_t *create(tape_options stce_10170 = tape_options()) {
    tape_t* tape = new tape_t(stce_10170);
    return tape;
  }
  static void remove(tape_t *&tape) {
    if (tape == 0)
      return;
    delete tape;
    tape = 0;
  }
  virtual stce_10085 stce_10080(const DCO_TAPE_INT stce_9808) {
    return stce_10393(stce_9808);
  }
  virtual void stce_10081(const DCO_TAPE_INT stce_9808, const stce_10085 &stce_10082) {
    stce_10393(stce_9808) += stce_10082;
  }
  typedef dco::helper::stce_10084<tape_t> callback_object_t;
  template <typename stce_10400>
    class stce_10401 {
public:
    typedef void (*stce_10402)(tape_t &stce_10403, const interpretation_settings &stce_9898, stce_10400 *stce_10404);
    typedef void (*stce_10405)(tape_t &stce_10403, stce_10400 *stce_10404);
    typedef void (*stce_10406)(stce_10400 *stce_10404);
  };
  class stce_10407 {
public:
    virtual void stce_10408(tape_t &stce_10403, const interpretation_settings &stce_9898, callback_object_t *stce_10404) = 0;
    virtual ~stce_10407() {};
  };
  template <typename stce_10400>
    class stce_10409 : public stce_10407 {
private:
    union {
      typename stce_10401<stce_10400>::stce_10406 stce_10410;
      typename stce_10401<stce_10400>::stce_10405 stce_10411;
      typename stce_10401<stce_10400>::stce_10402 stce_10412;
    } stce_10410;
    int stce_10413;
public:
    stce_10409(typename stce_10401<stce_10400>::stce_10406 stce_10414) : stce_10413(0) {
      stce_10410.stce_10410 = stce_10414;
    }
    stce_10409(typename stce_10401<stce_10400>::stce_10405 stce_10414) : stce_10413(1) {
      stce_10410.stce_10411 = stce_10414;
    }
    stce_10409(typename stce_10401<stce_10400>::stce_10402 stce_10414) : stce_10413(2) {
      stce_10410.stce_10412 = stce_10414;
    }
    void stce_10408(tape_t &stce_10403, const interpretation_settings &stce_9898, callback_object_t *stce_10404) {
      stce_10400 *stce_10415 = static_cast<stce_10400 *>(stce_10404);
      switch (stce_10413) {
      case 0:
        stce_10410.stce_10410(stce_10415);
        break;
      case 1:
        stce_10410.stce_10411(stce_10403, stce_10415);
        break;
      case 2:
        stce_10410.stce_10412(stce_10403, stce_9898, stce_10415);
        break;
      default:
        if (helper::stce_9822) {
          throw dco::exception::create<std::runtime_error>("unknown callback to run.");
        }
        break;
      }
    }
    ~stce_10409() {}
  };
  class stce_10374 {
    callback_object_t *stce_10404;
    stce_10407 *stce_10416;
    iterator_t stce_10417;
public:
    stce_10374(const iterator_t& stce_10158) :
      stce_10404(0), stce_10416(0), stce_10417(stce_10158){
    }
    bool operator<(const stce_10374& stce_9933) const
    {
      return stce_10417 < stce_9933.stce_10417;
    }
    void stce_10376() {
      delete stce_10404;
      if (stce_10416)
        delete stce_10416;
    }
    template<typename stce_10400, typename stce_10418>
      void stce_10419(stce_10418 stce_10414) {
      if (stce_10416)
        throw dco::exception::create<std::runtime_error>(
                                                         "currently not supported to insert external_adjoint_object_bases twice.");
      stce_10416 = new stce_10409<stce_10400>(stce_10414);
    }
    callback_object_t *&stce_10327() {
      return stce_10404;
    }
    callback_object_t * stce_10327() const {
      return stce_10404;
    }
    iterator_t &stce_10326() {
      return stce_10417;
    }
    iterator_t stce_10326() const {
      return stce_10417;
    }
    void stce_10408(tape_t &stce_10403, const interpretation_settings &stce_9898) {
      if (stce_10416)
        stce_10416->stce_10408(stce_10403, stce_9898, stce_10404);
    }
  };
  std::vector<stce_10374> stce_10375;
  size_t number_of_callbacks() {
    return stce_10375.size();
  }
  template <class stce_10420>
   void stce_10421(stce_10420 *stce_10404) {
       stce_10374 stce_10422(stce_10283::current());
       stce_10375.push_back(stce_10422);
       stce_10375.back().stce_10327() = stce_10404;
       stce_10404->set_tape(this);
       typedef typename stce_10420::active_t stce_10423;
       typedef dco::stce_10137::stce_10140<tape_t,stce_10423> stce_10424;
       stce_10424* stce_10136 = stce_10404;
       insert_callback(&dco::stce_10137::stce_10135<tape_t,stce_10423>, stce_10136);
   }
  template<class stce_10420, typename stce_10418>
    void insert_callback(stce_10418 stce_10416, stce_10420 *stce_10380) {
    if (stce_10375.back().stce_10327() == stce_10380) {
      stce_10375.back().template stce_10419<stce_10420>(stce_10416);
      stce_10375.back().stce_10326() = stce_10283::current();
      stce_10383< true>();
    } else {
      throw dco::exception::create<std::runtime_error>("please always insert most recently created external function object.");
    }
  }
  template<class stce_10420, typename FCN_PARAMETERS>
    stce_10420 *create_callback_object(const FCN_PARAMETERS &parameters) {
    const stce_9844<FCN_PARAMETERS> stce_10425(parameters);
    stce_10374 stce_10422(stce_10283::current());
    stce_10375.push_back(stce_10422);
    stce_10420 *stce_10404 = stce_10425.template create<stce_10420>();
    stce_10375.back().stce_10327() = stce_10404;
    stce_10404->set_tape(this);
    return stce_10404;
  }
  template<class stce_10420>
    stce_10420 *create_callback_object() {
    void *stce_9998;
    return create_callback_object<stce_10420>(stce_9998);
  }
  template <class stce_10420>
    inline stce_10420* create_ACM_object() {
    stce_10420* stce_9840 = create_callback_object<stce_10420>();
    typedef typename stce_10420::active_t stce_10423;
    typedef dco::ACM::baseclass<tape_t,stce_10423> stce_10426;
    stce_10426* stce_10136 = stce_9840;
    insert_callback(&dco::ACM::stce_10135<tape_t,stce_10423>, stce_10136);
    return stce_9840;
  }
  template <class stce_10420, typename FCN_PARAMETERS>
    inline stce_10420 *create_ACM_object(const FCN_PARAMETERS &parameters) {
    const stce_9844<FCN_PARAMETERS> stce_10425(parameters);
    stce_10374 stce_10422(stce_10283::current());
    stce_10375.push_back(stce_10422);
    stce_10420 *stce_10404 = stce_10425.template create<stce_10420>();
    stce_10375.back().stce_10327() = stce_10404;
    stce_10404->set_tape(this);
    typedef typename stce_10420::active_t stce_10423;
    typedef dco::ACM::baseclass<tape_t,stce_10423> stce_10426;
    stce_10426* stce_10136 = stce_10404;
    insert_callback(&dco::ACM::stce_10135<tape_t,stce_10423>, stce_10136);
    return stce_10404;
  }
    template <class stce_10420>
    void insert_ACM_object(stce_10420 *stce_10404) {
        stce_10374 stce_10422(stce_10283::current());
        stce_10375.push_back(stce_10422);
        stce_10375.back().stce_10327() = stce_10404;
        stce_10404->set_tape(this);
        typedef typename stce_10420::active_t stce_10423;
        typedef dco::ACM::baseclass<tape_t,stce_10423> stce_10426;
        stce_10426* stce_10136 = stce_10404;
        insert_callback(&dco::ACM::stce_10135<tape_t,stce_10423>, stce_10136);
    }
  void stce_10392(const iterator_t &stce_9863) {
    typedef typename std::vector<stce_10374>::iterator stce_10325;
    stce_10325 stce_10054 = stce_10375.begin();
    for (; stce_10054 != stce_10375.end() && stce_10054->stce_10326().index() < stce_9863.index(); ++stce_10054) {}
    for (stce_10325 stce_10164 = stce_10054; stce_10164 != stce_10375.end(); ++stce_10164)
      stce_10164->stce_10376();
    stce_10375.erase(stce_10054, stce_10375.end());
  }
public:
  void write_to_dot(const std::string &filename = "tape.dot")
  {
    write_to_dot(stce_10283::current(), stce_10283::start(), filename);
  }
  void write_to_dot(iterator_t stce_9862)
  {
    write_to_dot(stce_9862, stce_10283::start());
  }
  void write_to_dot(iterator_t stce_9817, iterator_t stce_9863,
                   const std::string &filename = "tape.dot") {
    stce_9848 stce_9849(filename);
    stce_10320(stce_9817, stce_9863 , stce_10290, stce_10375, stce_9849);
  }
  void write_to_csv(const std::string &filename = "tape.csv")
  {
    write_to_csv(stce_10283::current(), stce_10283::start(), filename);
  }
  void write_to_csv(iterator_t stce_9862)
  {
    write_to_csv(stce_9862, stce_10283::start());
  }
  void write_to_csv(iterator_t stce_9817, iterator_t stce_9863,
                   const std::string &filename = "tape.csv") {
    if(stce_10375.size() != 0)
      throw dco::exception::create<std::runtime_error>
        ("write_to_csv not implemented for callbacks.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 679);
    stce_10328(stce_9817, stce_9863, filename);
  }
    void stce_10427(const iterator_t& stce_9817, stce_10374& stce_10428,
                          const interpretation_settings &stce_10311)
  {
    if(stce_10428.stce_10326() != stce_9817)
      throw dco::exception::create<std::runtime_error>("invoking callback at wrong position", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 688);
    if (stce_10311.reset) {
      stce_10265(stce_9817);
      stce_10383< true>();
    }
    int stce_9837 = dco::helper::stce_9836<typename stce_9762::stce_9764>::stce_9837;
    for (DCO_TAPE_INT stce_10429 = 0; stce_10429 < stce_9837; stce_10429++) {
      stce_10309 = stce_10429;
      stce_10428.stce_10327()->stce_10087 = this;
      stce_10428.stce_10408(*this, stce_10311);
    }
  }
  template<class stce_10430>
    void stce_10427(stce_10430 &stce_10394, const iterator_t& stce_9817, stce_10374& stce_10428,
                          const interpretation_settings &stce_10311)
  {
    if(stce_10428.stce_10326() != stce_9817)
      throw dco::exception::create<std::runtime_error>("invoking callback at wrong position", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 712);
    if (stce_10311.reset) {
      stce_10265(stce_9817);
      stce_10383< true>();
    }
    int stce_9837 = dco::helper::stce_9836<typename stce_9762::stce_9764>::stce_9837;
    for (DCO_TAPE_INT stce_10429 = 0; stce_10429 < stce_9837; stce_10429++) {
      stce_10309 = stce_10429;
      stce_10428.stce_10327()->stce_10087 = &stce_10394;
      stce_10428.stce_10408(*this, stce_10311);
    }
  }
  struct stce_10431 {
    tape* stce_9952;
    bool stce_10369;
    stce_10431(tape* stce_9905) : stce_9952(stce_9905), stce_10369(stce_9905->sparse_interpret()) {}
    ~stce_10431() { stce_9952->sparse_interpret() = stce_10369; }
  };
  void stce_10312(iterator_t stce_9817, const iterator_t &stce_9863,
                                   const interpretation_settings &stce_10311) {
    stce_10431 stce_10432(this);
    if (stce_9817 > stce_10283::current()) {
      throw dco::exception::create<std::runtime_error>
        ("you try to use a tape position outside of the current tape.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 747);
    }
    if (stce_9863 > stce_9817) {
      throw dco::exception::create<std::runtime_error>
        ("adjoint interpretation: from < to.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 752);
    }
    stce_10370();
    typedef typename std::vector<stce_10374>::reverse_iterator stce_10325;
    for (stce_10325 stce_10164 = stce_10375.rbegin();
         stce_10164 != stce_10375.rend() && stce_9863 <= stce_10164->stce_10326();
         ++stce_10164) {
      const iterator_t& stce_10433 = stce_10164->stce_10326();
      if(stce_9817 < stce_10433) { continue; }
      stce_10283::stce_10183(stce_9817);
      stce_10371(stce_9817, stce_10433, stce_10290);
      stce_9817 = stce_10433;
      std::size_t stce_10434 = std::distance(stce_10375.rbegin(), stce_10164);
      if(stce_10311.stce_10065)
        stce_10290.stce_10303(stce_9817.index(), stce_10433.index());
      stce_10283::stce_10184();
      stce_10427(stce_9817, *stce_10164, stce_10311);
      stce_10164 = stce_10375.rbegin() + stce_10434;
      if (stce_10311.reset) {
        stce_10265(stce_9817);
        stce_10164->stce_10376();
        stce_10375.erase( --(stce_10164.base()) );
      }
    }
    stce_10283::stce_10183(stce_9817);
    stce_10371(stce_9817, stce_9863, stce_10290);
    if(stce_10311.stce_10065)
      stce_10290.stce_10303(stce_9817.index(), stce_9863.index());
    stce_10283::stce_10184();
  }
  template<class stce_10430>
    void stce_10316(stce_10430 &stce_10394,
                                       iterator_t stce_9817, const iterator_t &stce_9863,
                                       const interpretation_settings &stce_10311)
  {
    if (stce_9863 > stce_9817) {
      throw dco::exception::create<std::runtime_error>("adjoint interpretation: from < to.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 817);
    }
    if (stce_9817 > stce_10283::current()) {
      throw dco::exception::create<std::runtime_error>("you try to use a tape position outside of the current tape.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 820);
    }
    typedef typename std::vector<stce_10374>::reverse_iterator stce_10325;
    for (stce_10325 stce_10164 = stce_10375.rbegin(); stce_10164 != stce_10375.rend() && stce_9863 <= stce_10164->stce_10326();
         ++stce_10164) {
      const iterator_t& stce_10433 = stce_10164->stce_10326();
      if(stce_9817 < stce_10433)
        continue;
      stce_10283::stce_10183(stce_9817);
      if(stce_10311.stce_10065)
        stce_10372(stce_9817, stce_10433, stce_10394);
      else {
        stce_10371(stce_9817, stce_10433, stce_10394);
      }
      stce_9817 = stce_10433;
      std::size_t stce_10434 = std::distance(stce_10375.rbegin(), stce_10164);
      stce_10427(stce_10394, stce_9817, *stce_10164, stce_10311);
      if (stce_9817 != stce_10283::start()) {
        --stce_9817;
      }
      stce_10164 = stce_10375.rbegin() + stce_10434;
      if (stce_10311.reset) {
        stce_10265(stce_9817);
        stce_10164->stce_10376();
        stce_10375.erase( --(stce_10164.base()) );
      }
    }
    stce_10283::stce_10183(stce_9817);
    if(stce_10311.stce_10065) { stce_10372(stce_9817, stce_9863, stce_10394); }
    else { stce_10371(stce_9817, stce_9863, stce_10394); }
    if (stce_10311.reset) { stce_10265(stce_9817); }
    stce_10283::stce_10184();
  }
};
}
}
       
namespace dco {
template <size_t stce_10435>
struct subrange_t {
  static const size_t stce_10436 = stce_10435;
  subrange_t(size_t stce_10102) : stce_10102(stce_10102) {}
  size_t stce_10102;
};
template <size_t stce_10435>
struct range_t {
  subrange_t<stce_10435> const& stce_9898;
  size_t global_index, index;
  range_t(subrange_t<stce_10435> const& stce_9898) : stce_9898(stce_9898), global_index(0), index(0) {}
  range_t& operator++() {
    global_index += stce_9898.stce_10436;
    ++index;
    return *this;
  }
  bool operator!=(range_t const& stce_10437) const {
    if (stce_10437.index != index) return true;
    return false;
  }
  range_t & operator*() { return *this; }
  range_t const& operator*() const { return *this; }
  void stce_10193() {
    global_index = stce_9898.stce_10102;
    index = static_cast<size_t>(std::ceil(static_cast<double>(stce_9898.stce_10102)/stce_10435));
  }
  bool is_last() const {
    return global_index + stce_10435 >= stce_9898.stce_10102;
  }
};
template <size_t stce_10435>
struct inner_range_t {
  range_t<stce_10435> const& stce_9820;
  size_t global_index, sub_index;
  inner_range_t(range_t<stce_10435> const& stce_9820) : stce_9820(stce_9820), global_index(stce_9820.global_index), sub_index(0) {}
  inner_range_t& operator++() {
    ++sub_index;
    ++global_index;
    return *this;
  }
  void stce_10193() {
    if ((stce_9820.index+1)*stce_9820.stce_9898.stce_10436 > stce_9820.stce_9898.stce_10102) {
      sub_index = stce_9820.stce_9898.stce_10102 % stce_9820.stce_9898.stce_10436;
    } else {
      sub_index = stce_9820.stce_9898.stce_10436;
    }
  }
  bool operator!=(inner_range_t const& stce_10437) const {
    if (stce_10437.sub_index != sub_index) return true;
    return false;
  }
  inner_range_t & operator*() { return *this; }
  inner_range_t const& operator*() const { return *this; }
};
template <size_t stce_10435>
  static inner_range_t<stce_10435> begin(range_t<stce_10435> const& stce_9820) { return inner_range_t<stce_10435>(stce_9820); }
template <size_t stce_10435>
  static range_t<stce_10435> begin(subrange_t<stce_10435> const& stce_9898) { return range_t<stce_10435>(stce_9898); }
template <size_t stce_10435>
  static inner_range_t<stce_10435> end(range_t<stce_10435> const& stce_9820) {
  inner_range_t<stce_10435> stce_10438(stce_9820);
  stce_10438.stce_10193();
  return stce_10438;
}
template <size_t stce_10435>
static range_t<stce_10435> end(subrange_t<stce_10435> const& stce_9898) {
  range_t<stce_10435> stce_10438(stce_9898);
  stce_10438.stce_10193();
  return stce_10438;
}
template <size_t stce_10435>
  static subrange_t<stce_10435> subranges(size_t stce_10102) {
  return subrange_t<stce_10435>(stce_10102);
}
}
namespace dco {
  namespace stce_10439 {
    template<class stce_9771>struct stce_10440{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { return static_cast<stce_9771>(stce_9912 + stce_9913); } template<class stce_10441,class stce_10030>static inline const stce_9771 stce_10039(const stce_9771 &_value , const stce_10441 &stce_9912 , const stce_10030 &stce_9913 ) { (void) _value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(static_cast<stce_9771>(1.0)); } template<class stce_10441,class stce_10030>static inline const stce_9771 stce_10040(const stce_9771 &_value , const stce_10441 &stce_9912 , const stce_10030 &stce_9913 ) { (void) _value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(static_cast<stce_9771>(1.0)); } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << stce_9957 << "+" << stce_9958; return stce_9840.str(); } static std::string stce_9914() { return "static_cast<DCO_TAPE_REAL>(1.0)"; } static std::string stce_9916() { return "static_cast<DCO_TAPE_REAL>(1.0)"; } };
    template<class stce_9771>struct stce_10442{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { return static_cast<stce_9771>(stce_9912 - stce_9913); } template<class stce_10441,class stce_10030>static inline const stce_9771 stce_10039(const stce_9771 &_value , const stce_10441 &stce_9912 , const stce_10030 &stce_9913 ) { (void) _value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(static_cast<stce_9771>(1.0)); } template<class stce_10441,class stce_10030>static inline const stce_9771 stce_10040(const stce_9771 &_value , const stce_10441 &stce_9912 , const stce_10030 &stce_9913 ) { (void) _value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(static_cast<stce_9771>(-1.0)); } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << stce_9957 << "-" << stce_9958; return stce_9840.str(); } static std::string stce_9914() { return "static_cast<DCO_TAPE_REAL>(1.0)"; } static std::string stce_9916() { return "static_cast<DCO_TAPE_REAL>(-1.0)"; } };
    template<class stce_9771>struct stce_10443{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { return static_cast<stce_9771>(stce_9912 * stce_9913); } template<class stce_10441,class stce_10030>static inline const stce_9771 stce_10039(const stce_9771 &_value , const stce_10441 &stce_9912 , const stce_10030 &stce_9913 ) { (void) _value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(stce_9913); } template<class stce_10441,class stce_10030>static inline const stce_9771 stce_10040(const stce_9771 &_value , const stce_10441 &stce_9912 , const stce_10030 &stce_9913 ) { (void) _value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(stce_9912); } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << stce_9957 << "*" << stce_9958; return stce_9840.str(); } static std::string stce_9914() { return "arg2"; } static std::string stce_9916() { return "arg1"; } };
    template<class stce_9771>struct stce_10444{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { return static_cast<stce_9771>(stce_9912 / stce_9913); } template<class stce_10441,class stce_10030>static inline const stce_9771 stce_10039(const stce_9771 &_value , const stce_10441 &stce_9912 , const stce_10030 &stce_9913 ) { (void) _value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) / stce_9913); } template<class stce_10441,class stce_10030>static inline const stce_9771 stce_10040(const stce_9771 &_value , const stce_10441 &stce_9912 , const stce_10030 &stce_9913 ) { (void) _value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(-_value / stce_9913); } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << stce_9957 << "/" << stce_9958; return stce_9840.str(); } static std::string stce_9914() { return "static_cast<DCO_TAPE_REAL>(1.0) / arg2"; } static std::string stce_9916() { return "-_value / arg2"; } };
    template<class stce_9771>struct stce_10445{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { return stce_9912 + stce_9913; } template<class stce_10441>static inline const stce_9771 stce_10039(const stce_9771 &_value , const stce_10441 &stce_9912 , const stce_9771 &stce_9913 ) { (void)_value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(1.0); } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << stce_9957 << "+" << stce_9958; return stce_9840.str(); } static std::string stce_9914() { return "static_cast<DCO_TAPE_REAL>(1.0)"; } }; template<class stce_9771>struct stce_10446{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { return stce_9912 + stce_9913; } template<class stce_10030>static inline const stce_9771 stce_10040(const stce_9771 &_value , const stce_9771 &stce_9912 , const stce_10030 &stce_9913 ) { (void) _value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(1.0); } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << stce_9957 << "+" << stce_9958; return stce_9840.str(); } static std::string stce_9916() { return "static_cast<DCO_TAPE_REAL>(1.0)"; } };
    template<class stce_9771>struct stce_10447{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { return stce_9912 - stce_9913; } template<class stce_10441>static inline const stce_9771 stce_10039(const stce_9771 &_value , const stce_10441 &stce_9912 , const stce_9771 &stce_9913 ) { (void)_value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(1.0); } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << stce_9957 << "-" << stce_9958; return stce_9840.str(); } static std::string stce_9914() { return "static_cast<DCO_TAPE_REAL>(1.0)"; } }; template<class stce_9771>struct stce_10448{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { return stce_9912 - stce_9913; } template<class stce_10030>static inline const stce_9771 stce_10040(const stce_9771 &_value , const stce_9771 &stce_9912 , const stce_10030 &stce_9913 ) { (void) _value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(-1.0); } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << stce_9957 << "-" << stce_9958; return stce_9840.str(); } static std::string stce_9916() { return "static_cast<DCO_TAPE_REAL>(-1.0)"; } };
    template<class stce_9771>struct stce_10449{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { return stce_9912 * stce_9913; } template<class stce_10441>static inline const stce_9771 stce_10039(const stce_9771 &_value , const stce_10441 &stce_9912 , const stce_9771 &stce_9913 ) { (void)_value; (void)stce_9912; (void)stce_9913; return stce_9913; } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << stce_9957 << "*" << stce_9958; return stce_9840.str(); } static std::string stce_9914() { return "arg2"; } }; template<class stce_9771>struct stce_10450{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { return stce_9912 * stce_9913; } template<class stce_10030>static inline const stce_9771 stce_10040(const stce_9771 &_value , const stce_9771 &stce_9912 , const stce_10030 &stce_9913 ) { (void) _value; (void)stce_9912; (void)stce_9913; return stce_9912; } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << stce_9957 << "*" << stce_9958; return stce_9840.str(); } static std::string stce_9916() { return "arg1"; } };
    template<class stce_9771>struct stce_10451{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { return stce_9912 / stce_9913; } template<class stce_10441>static inline const stce_9771 stce_10039(const stce_9771 &_value , const stce_10441 &stce_9912 , const stce_9771 &stce_9913 ) { (void)_value; (void)stce_9912; (void)stce_9913; return static_cast<stce_9771>(1.0) / stce_9913; } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << stce_9957 << "/" << stce_9958; return stce_9840.str(); } static std::string stce_9914() { return "static_cast<DCO_TAPE_REAL>(1.0) / arg2"; } }; template<class stce_9771>struct stce_10452{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { return stce_9912 / stce_9913; } template<class stce_10030>static inline const stce_9771 stce_10040(const stce_9771 &_value , const stce_9771 &stce_9912 , const stce_10030 &stce_9913 ) { (void) _value; (void)stce_9912; (void)stce_9913; return -(stce_9912 / stce_9913) / stce_9913; } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << stce_9957 << "/" << stce_9958; return stce_9840.str(); } static std::string stce_9916() { return "-(arg1 / arg2) / arg2"; } };
    template<class stce_9771>struct stce_10453{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using std::sin; return sin(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using std::cos; (void)_value; return static_cast<stce_9771>(cos(stce_9843)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "sin" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("cos(x)"); return stce_9840; } };
    template<class stce_9771>struct stce_10454{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using std::cos; return cos(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using std::sin; (void)_value; return static_cast<stce_9771>(-sin(stce_9843)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "cos" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("-sin(x)"); return stce_9840; } };
    template<class stce_9771>struct stce_10455{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using std::tan; return tan(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using std::tan; (void)_value; return static_cast<stce_9771>((static_cast<stce_9771>(1.0) + (tan(stce_9843)*tan(stce_9843)))); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "tan" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("(static_cast<DCO_TAPE_REAL>(1.0) + (tan(x)*tan(x)))"); return stce_9840; } };
    template<class stce_9771>struct stce_10456{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using std::cosh; return cosh(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using std::sinh; (void)_value; return static_cast<stce_9771>(sinh(stce_9843)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "cosh" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("sinh(x)"); return stce_9840; } };
    template<class stce_9771>struct stce_10457{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using std::sinh; return sinh(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using std::cosh; (void)_value; return static_cast<stce_9771>(cosh(stce_9843)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "sinh" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("cosh(x)"); return stce_9840; } };
    template<class stce_9771>struct stce_10458{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using std::tanh; return tanh(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using std::tanh; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) - tanh(stce_9843)*tanh(stce_9843)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "tanh" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("static_cast<DCO_TAPE_REAL>(1.0) - tanh(x)*tanh(x)"); return stce_9840; } };
    template<class stce_9771>struct stce_10459{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using std::asin; return asin(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using std::sqrt; (void)_value; return static_cast<stce_9771>(1 / sqrt(static_cast<stce_9771>(1.0) - stce_9843*stce_9843)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "asin" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("1 / sqrt(static_cast<DCO_TAPE_REAL>(1.0) - x*x)"); return stce_9840; } };
    template<class stce_9771>struct stce_10460{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using std::acos; return acos(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using std::sqrt; (void)_value; return static_cast<stce_9771>(-1 / sqrt(static_cast<stce_9771>(1.0) - stce_9843*stce_9843)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "acos" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("-1 / sqrt(static_cast<DCO_TAPE_REAL>(1.0) - x*x)"); return stce_9840; } };
    template<class stce_9771>struct stce_10461{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using std::atan; return atan(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using std::atan; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) / (static_cast<stce_9771>(1.0) + stce_9843*stce_9843)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "atan" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("static_cast<DCO_TAPE_REAL>(1.0) / (static_cast<DCO_TAPE_REAL>(1.0) + x*x)"); return stce_9840; } };
    template<class stce_9771>struct stce_10462{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using std::exp; return exp(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using std::exp; (void)_value; return static_cast<stce_9771>(exp(stce_9843)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "exp" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("exp(x)"); return stce_9840; } };
    template<class stce_9771>struct stce_10463{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using std::sqrt; return sqrt(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using std::sqrt; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) / (static_cast<stce_9771>(2.0) * sqrt(stce_9843 + ((HAS_SQRT_AVOID_DIVZERO && stce_9843==0) ? DCO_SQRT_EPS : 0.0)))); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "sqrt" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("static_cast<DCO_TAPE_REAL>(1.0) / (static_cast<DCO_TAPE_REAL>(2.0) * sqrt(x + ((HAS_SQRT_AVOID_DIVZERO && x==0) ? DCO_SQRT_EPS : 0.0)))"); return stce_9840; } };
    template<class stce_9771>struct stce_10464{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using std::log; return log(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using std::log; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) / stce_9843); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "log" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("static_cast<DCO_TAPE_REAL>(1.0) / x"); return stce_9840; } };
    template<class stce_9771>struct stce_10465{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using ::erf; return erf(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using ::exp; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(2.0) / sqrt(static_cast<stce_9771>(3.14159265358979323846264338327950288)) * exp(-stce_9843 * stce_9843)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "erf" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("static_cast<DCO_TAPE_REAL>(2.0) / sqrt(static_cast<DCO_TAPE_REAL>(3.14159265358979323846264338327950288)) * exp(-x * x)"); return stce_9840; } };
    template<class stce_9771>struct stce_10466{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using ::erfc; return erfc(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using ::exp; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(-2.0) / sqrt(static_cast<stce_9771>(3.14159265358979323846264338327950288)) * exp(-stce_9843 * stce_9843)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "erfc" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("static_cast<DCO_TAPE_REAL>(-2.0) / sqrt(static_cast<DCO_TAPE_REAL>(3.14159265358979323846264338327950288)) * exp(-x * x)"); return stce_9840; } };
    template<class stce_9771>struct stce_10467{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using ::asinh; return asinh(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using ::sqrt; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.) / sqrt(static_cast<stce_9771>(1.) + (stce_9843*stce_9843))); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "asinh" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("static_cast<DCO_TAPE_REAL>(1.) / sqrt(static_cast<DCO_TAPE_REAL>(1.) + (x*x))"); return stce_9840; } };
    template<class stce_9771>struct stce_10468{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using ::acosh; return acosh(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using ::sqrt; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.) / sqrt((stce_9843*stce_9843) - static_cast<stce_9771>(1.))); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "acosh" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("static_cast<DCO_TAPE_REAL>(1.) / sqrt((x*x) - static_cast<DCO_TAPE_REAL>(1.))"); return stce_9840; } };
    template<class stce_9771>struct stce_10469{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using ::atanh; return atanh(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using ::atanh; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.) / (static_cast<stce_9771>(1.) - (stce_9843*stce_9843))); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "atanh" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("static_cast<DCO_TAPE_REAL>(1.) / (static_cast<DCO_TAPE_REAL>(1.) - (x*x))"); return stce_9840; } };
    template<class stce_9771>struct stce_10470{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using ::expm1; return expm1(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using ::exp; (void)_value; return static_cast<stce_9771>(exp(stce_9843)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "expm1" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("exp(x)"); return stce_9840; } };
    template<class stce_9771>struct stce_10471{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using ::log1p; return log1p(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using ::log1p; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) / (stce_9843 + 1)); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "log1p" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("static_cast<DCO_TAPE_REAL>(1.0) / (x + 1)"); return stce_9840; } };
    template<class stce_9771>struct stce_10472{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &arg) { using ::log10; return log10(arg); } template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9843) { using ::log10; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) / (stce_9843*static_cast<stce_9771>(log(10)))); } static std::string stce_10037(const std::string& stce_9843) { std::stringstream stce_9840; stce_9840 << "log10" << "("<< stce_9843 << ")"; return stce_9840.str(); } static std::string stce_10038() { std::string stce_9840 = std::string("static_cast<DCO_TAPE_REAL>(1.0) / (x*static_cast<DCO_TAPE_REAL>(log(10)))"); return stce_9840; } };
    template<class stce_9771>struct stce_10473 {
      template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &stce_9912) {
        return -stce_9912;
      }
      template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9912 ) {
        (void)_value;
        (void)stce_9912;
        return static_cast<stce_9771>(-1.0);
      }
      static std::string stce_10037(const std::string& stce_9843) {
        return "-1 *" + stce_9843;
      }
      static std::string stce_10038() {
        return "-1";
      }
    };
    template<class stce_9771>struct stce_10474 {
      template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &stce_9912) {
        return stce_9912;
      }
      template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9912 ) {
        (void)_value;
        (void)stce_9912;
        return static_cast<stce_9771>(1.0);
      }
      static std::string stce_10037(const std::string& stce_9843) {
        return stce_9843;
      }
      static std::string stce_10038() {
        return "1";
      }
    };
    using ::fabs;
    template<class stce_9771>struct stce_10475 {
      template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &stce_9912) {
        return fabs(stce_9912);
      }
      template<class stce_9782>static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9912 ) {
        (void) _value;
        if (stce_9912 < 0) return static_cast<stce_9771>(-1.0);
        else return static_cast<stce_9771>(1.0);
      }
      static std::string stce_10037(const std::string& stce_9843) {
        return "fabs(" + stce_9843 + ")";
      }
      static std::string stce_10038() {
        return "((x >= 0) ? 1 : -1 )";
      }
    };
    using std::abs;
    template<class stce_9771> struct stce_10476 {
      template<class stce_9782> static inline const stce_9771 stce_10035(const stce_9782 &stce_9912) {
        return abs(stce_9912);
      }
      template<class stce_9782> static inline const stce_9771 stce_10036(const stce_9771 &_value , const stce_9782 &stce_9912 ) {
        (void) _value;
        if (stce_9912 < 0) return static_cast<stce_9771>(-1.0);
        else return static_cast<stce_9771>(1.0);
      }
      static std::string stce_10037(const std::string& stce_9843) {
        return "abs(" + stce_9843 + ")";
      }
      static std::string stce_10038() {
        return "((x >= 0) ? 1 : -1 )";
      }
    };
    using std::atan2; template<class stce_9771>struct stce_10477{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { (void) stce_9912; (void) stce_9913; return atan2(stce_9912,stce_9913); } template<class stce_10441,class stce_10030>static inline const stce_9771 stce_10039(stce_9771 const& _value , const stce_10441 &stce_9912, const stce_10030 &stce_9913) { (void) _value; (void) stce_9912; (void) stce_9913; return stce_9913 / (stce_9913 * stce_9913 + stce_9912 * stce_9912); } template<class stce_10441,class stce_10030>static inline const stce_9771 stce_10040(stce_9771 const& _value , const stce_10441 &stce_9912, const stce_10030 &stce_9913) { (void) _value; (void) stce_9912; (void) stce_9913; return -stce_9912 / (stce_9913 * stce_9913 + stce_9912 * stce_9912); } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << "atan2" << "(" << stce_9957 << "," << stce_9958 <<")"; return stce_9840.str(); } static std::string stce_9914() { return "arg2 / (arg2 * arg2 + arg1 * arg1)"; } static std::string stce_9916() { return "-arg1 / (arg2 * arg2 + arg1 * arg1)"; } }; template<class stce_9771>struct stce_10478{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &stce_9912, const stce_9771 &stce_9913) { (void) stce_9912; (void) stce_9913; return atan2(stce_9912,stce_9913); } template<class stce_9782>static inline const stce_9771 stce_10039(stce_9771 const& _value , const stce_9782 &stce_9912, const stce_9771 &stce_9913) { (void) _value; (void) stce_9912; (void) stce_9913; return stce_9913 / (stce_9913 * stce_9913 + stce_9912 * stce_9912); } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << "atan2" << "(" << stce_9957 << "," << stce_9958 <<")"; return stce_9840.str(); } static std::string stce_9914() { return "arg2 / (arg2 * arg2 + arg1 * arg1)"; } }; template<class stce_9771>struct stce_10479{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9771 &stce_9912, const stce_9782 &stce_9913) { (void) stce_9912; (void) stce_9913; return atan2(stce_9912,stce_9913); } template<class stce_9782>static inline const stce_9771 stce_10040(stce_9771 const& _value , const stce_9771 &stce_9912, const stce_9782 &stce_9913) { (void) _value; (void) stce_9912; (void) stce_9913; return -stce_9912 / (stce_9913 * stce_9913 + stce_9912 *stce_9912); } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << "atan2" << "(" << stce_9957 << "," << stce_9958 <<")"; return stce_9840.str(); } static std::string stce_9916() { return "-arg1 / (arg2 * arg2 + arg1 *arg1)"; } };
    template<class stce_9771>
    struct stce_10480 {
      template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) {
        using std::pow;
        return pow(stce_9912, stce_9913);
      }
      template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10039(stce_9771 const& _value , const stce_10441 &stce_9912, const stce_10030 &stce_9913) {
        using std::pow;
        (void)_value;
        return stce_9913 * pow(stce_9912, stce_9913 - static_cast<stce_9771>(1.0));
      }
      template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10040(stce_9771 const& _value, const stce_10441 &stce_9912, const stce_10030&) {
        using std::pow;
        using std::log;
        return log(stce_9912) * _value;
      }
      static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) {
        return "pow(" + stce_9957 + "," + stce_9958 + ")";
      }
      static std::string stce_9914() {
        return "arg2 * pow(arg1, arg2 - 1.0)";
      }
      static std::string stce_9916() {
        return "log(arg1) * pow(arg1, arg2)";
      }
    };
    template<class stce_9771>
    struct stce_10481 {
      template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &stce_9912, const stce_9771 &stce_9913) {
        using std::pow;
        return pow(stce_9912, stce_9913);
      }
      template<class stce_9782>static inline const stce_9771 stce_10039(stce_9771 const& _value , const stce_9782 &stce_9912, const stce_9771 &stce_9913) {
        using std::pow;
        (void) _value;
        return stce_9913 * pow(stce_9912, stce_9913 - static_cast<stce_9771>(1.0));
      }
      static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) {
        return "pow(" + stce_9957 + "," + stce_9958 + ")";
      }
      static std::string stce_9914() {
        return "arg2 * pow(arg1, arg2 - 1.0)";
      }
    };
    template<class stce_9771>
    struct stce_10482 {
      template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9771 &stce_9912, const stce_9782 &stce_9913) {
        using std::pow;
        return pow(stce_9912, stce_9913);
      }
      template<class stce_9782>static inline const stce_9771 stce_10040(stce_9771 const& _value , const stce_9771 &stce_9912, const stce_9782&) {
        using std::pow;
        using std::log;
        return log(stce_9912) * _value;
      }
      static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) {
        return "pow(" + stce_9957 + "," + stce_9958 + ")";
      }
      static std::string stce_9914() {
        return "log(arg1) * pow(arg1, arg2)";
      }
    };
  template<typename stce_9804, size_t stce_9805, typename stce_10483> struct stce_10484 {
    template<class stce_9782> static inline const vector_ns::vector_t<stce_9804, stce_9805>
      stce_10036(const vector_ns::vector_t<stce_9804, stce_9805> &_value, const stce_9782 &stce_9843) {
      vector_ns::vector_t<stce_9804, stce_9805> stce_9840;
      for (size_t stce_9832 = 0; stce_9832 < stce_9805; ++stce_9832) { stce_9840[stce_9832] = stce_10483::stce_10036(_value[stce_9832], stce_9843[stce_9832]); }
      return stce_9840;
    }
    static std::string stce_10037(const std::string&) { return "NOT IMPLEMENTED FOR VECTOR TYPES."; }
    static std::string stce_10038() { return "NOT IMPLEMENTED FOR VECTOR TYPES."; }
  };
  template<typename stce_9804, size_t stce_9805> struct stce_10476<vector_ns::vector_t<stce_9804, stce_9805>>
    : public stce_10484<stce_9804, stce_9805, stce_10476<stce_9804>> {
    template<class stce_9782> static inline const vector_ns::vector_t<stce_9804, stce_9805> stce_10035(const stce_9782 &stce_9912) { return abs(stce_9912); }
  };
  template<typename stce_9804, size_t stce_9805> struct stce_10475<vector_ns::vector_t<stce_9804, stce_9805>>
    : public stce_10484<stce_9804, stce_9805, stce_10475<stce_9804>> {
    template<class stce_9782> static inline const vector_ns::vector_t<stce_9804, stce_9805> stce_10035(const stce_9782 &stce_9912) { return fabs(stce_9912); }
  };
  template<typename stce_9804, size_t stce_9805> struct stce_10463<vector_ns::vector_t<stce_9804, stce_9805>>
    : public stce_10484<stce_9804, stce_9805, stce_10463<stce_9804>> {
    template<class stce_9782>static inline const vector_ns::vector_t<stce_9804, stce_9805> stce_10035(const stce_9782 &arg) { return sqrt(arg); }
  };
    using ::hypot; template<class stce_9771>struct stce_10485{ template<class stce_10441, class stce_10030>static inline const stce_9771 stce_10035(const stce_10441 &stce_9912, const stce_10030 &stce_9913) { (void) stce_9912; (void) stce_9913; return hypot(stce_9912,stce_9913); } template<class stce_10441,class stce_10030>static inline const stce_9771 stce_10039(stce_9771 const& _value , const stce_10441 &stce_9912, const stce_10030 &stce_9913) { (void) _value; (void) stce_9912; (void) stce_9913; return stce_9912 / _value; } template<class stce_10441,class stce_10030>static inline const stce_9771 stce_10040(stce_9771 const& _value , const stce_10441 &stce_9912, const stce_10030 &stce_9913) { (void) _value; (void) stce_9912; (void) stce_9913; return stce_9913 / _value; } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << "hypot" << "(" << stce_9957 << "," << stce_9958 <<")"; return stce_9840.str(); } static std::string stce_9914() { return "arg1 / _value"; } static std::string stce_9916() { return "arg2 / _value"; } }; template<class stce_9771>struct stce_10486{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9782 &stce_9912, const stce_9771 &stce_9913) { (void) stce_9912; (void) stce_9913; return hypot(stce_9912,stce_9913); } template<class stce_9782>static inline const stce_9771 stce_10039(stce_9771 const& _value , const stce_9782 &stce_9912, const stce_9771 &stce_9913) { (void) _value; (void) stce_9912; (void) stce_9913; return stce_9912 / _value; } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << "hypot" << "(" << stce_9957 << "," << stce_9958 <<")"; return stce_9840.str(); } static std::string stce_9914() { return "arg1 / _value"; } }; template<class stce_9771>struct stce_10487{ template<class stce_9782>static inline const stce_9771 stce_10035(const stce_9771 &stce_9912, const stce_9782 &stce_9913) { (void) stce_9912; (void) stce_9913; return hypot(stce_9912,stce_9913); } template<class stce_9782>static inline const stce_9771 stce_10040(stce_9771 const& _value , const stce_9771 &stce_9912, const stce_9782 &stce_9913) { (void) _value; (void) stce_9912; (void) stce_9913; return stce_9913 / _value; } static std::string stce_10037(const std::string& stce_9957,const std::string& stce_9958) { std::stringstream stce_9840; stce_9840 << "hypot" << "(" << stce_9957 << "," << stce_9958 <<")"; return stce_9840.str(); } static std::string stce_9916() { return "arg2 / _value"; } };
  }
}
namespace dco {
namespace internal {
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10473<stce_9771> > operator -( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10473<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10473<stce_9771> > operator -( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10473<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10473<stce_9771> > operator -( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10473<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10473<stce_9771> > operator -( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10473<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10473<stce_9771> > operator -( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10473<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10474<stce_9771> > operator +( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10474<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10474<stce_9771> > operator +( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10474<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10474<stce_9771> > operator +( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10474<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10474<stce_9771> > operator +( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10474<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10474<stce_9771> > operator +( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10474<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10453<stce_9771> > sin( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10453<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10453<stce_9771> > sin( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10453<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10453<stce_9771> > sin( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10453<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10453<stce_9771> > sin( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10453<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10453<stce_9771> > sin( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10453<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10454<stce_9771> > cos( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10454<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10454<stce_9771> > cos( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10454<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10454<stce_9771> > cos( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10454<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10454<stce_9771> > cos( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10454<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10454<stce_9771> > cos( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10454<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10455<stce_9771> > tan( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10455<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10455<stce_9771> > tan( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10455<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10455<stce_9771> > tan( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10455<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10455<stce_9771> > tan( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10455<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10455<stce_9771> > tan( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10455<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10456<stce_9771> > cosh( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10456<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10456<stce_9771> > cosh( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10456<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10456<stce_9771> > cosh( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10456<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10456<stce_9771> > cosh( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10456<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10456<stce_9771> > cosh( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10456<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10457<stce_9771> > sinh( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10457<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10457<stce_9771> > sinh( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10457<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10457<stce_9771> > sinh( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10457<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10457<stce_9771> > sinh( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10457<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10457<stce_9771> > sinh( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10457<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10459<stce_9771> > asin( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10459<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10459<stce_9771> > asin( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10459<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10459<stce_9771> > asin( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10459<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10459<stce_9771> > asin( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10459<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10459<stce_9771> > asin( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10459<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10460<stce_9771> > acos( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10460<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10460<stce_9771> > acos( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10460<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10460<stce_9771> > acos( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10460<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10460<stce_9771> > acos( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10460<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10460<stce_9771> > acos( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10460<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10462<stce_9771> > exp( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10462<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10462<stce_9771> > exp( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10462<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10462<stce_9771> > exp( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10462<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10462<stce_9771> > exp( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10462<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10462<stce_9771> > exp( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10462<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10461<stce_9771> > atan( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10461<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10461<stce_9771> > atan( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10461<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10461<stce_9771> > atan( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10461<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10461<stce_9771> > atan( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10461<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10461<stce_9771> > atan( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10461<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10458<stce_9771> > tanh( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10458<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10458<stce_9771> > tanh( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10458<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10458<stce_9771> > tanh( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10458<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10458<stce_9771> > tanh( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10458<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10458<stce_9771> > tanh( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10458<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10463<stce_9771> > sqrt( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10463<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10463<stce_9771> > sqrt( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10463<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10463<stce_9771> > sqrt( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10463<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10463<stce_9771> > sqrt( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10463<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10463<stce_9771> > sqrt( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10463<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10464<stce_9771> > log( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10464<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10464<stce_9771> > log( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10464<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10464<stce_9771> > log( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10464<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10464<stce_9771> > log( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10464<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10464<stce_9771> > log( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10464<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10465<stce_9771> > erf( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10465<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10465<stce_9771> > erf( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10465<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10465<stce_9771> > erf( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10465<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10465<stce_9771> > erf( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10465<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10465<stce_9771> > erf( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10465<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10466<stce_9771> > erfc( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10466<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10466<stce_9771> > erfc( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10466<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10466<stce_9771> > erfc( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10466<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10466<stce_9771> > erfc( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10466<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10466<stce_9771> > erfc( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10466<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10470<stce_9771> > expm1( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10470<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10470<stce_9771> > expm1( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10470<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10470<stce_9771> > expm1( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10470<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10470<stce_9771> > expm1( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10470<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10470<stce_9771> > expm1( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10470<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10467<stce_9771> > asinh( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10467<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10467<stce_9771> > asinh( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10467<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10467<stce_9771> > asinh( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10467<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10467<stce_9771> > asinh( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10467<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10467<stce_9771> > asinh( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10467<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10468<stce_9771> > acosh( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10468<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10468<stce_9771> > acosh( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10468<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10468<stce_9771> > acosh( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10468<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10468<stce_9771> > acosh( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10468<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10468<stce_9771> > acosh( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10468<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10469<stce_9771> > atanh( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10469<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10469<stce_9771> > atanh( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10469<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10469<stce_9771> > atanh( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10469<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10469<stce_9771> > atanh( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10469<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10469<stce_9771> > atanh( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10469<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10471<stce_9771> > log1p( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10471<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10471<stce_9771> > log1p( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10471<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10471<stce_9771> > log1p( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10471<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10471<stce_9771> > log1p( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10471<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10471<stce_9771> > log1p( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10471<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10472<stce_9771> > log10( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10472<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10472<stce_9771> > log10( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10472<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10472<stce_9771> > log10( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10472<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10472<stce_9771> > log10( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10472<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10472<stce_9771> > log10( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10472<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10475<stce_9771> > fabs( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10475<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10475<stce_9771> > fabs( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10475<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10475<stce_9771> > fabs( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10475<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10475<stce_9771> > fabs( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10475<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10475<stce_9771> > fabs( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10475<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10476<stce_9771> > abs( const dco::internal::active_type<stce_9771, stce_9939> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::stce_10439::stce_10476<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10476<stce_9771> > abs( const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::stce_10439::stce_10476<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10476<stce_9771> > abs( const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::stce_10439::stce_10476<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10476<stce_9771> > abs( const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::stce_10439::stce_10476<stce_9771> >(stce_9957); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10476<stce_9771> > abs( const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::stce_10439::stce_10476<stce_9771> >(stce_9957); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> > operator + (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10440<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10445<stce_9771> > operator + (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10445<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10446<stce_9771> > operator + (const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10446<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10445<stce_9771> > operator + (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10445<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10446<stce_9771> > operator + (const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10446<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10445<stce_9771> > operator + (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10445<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10446<stce_9771> > operator + (const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10446<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10445<stce_9771> > operator + (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10445<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10446<stce_9771> > operator + (const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10446<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10445<stce_9771> > operator + (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10445<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10446<stce_9771> > operator + (const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10446<stce_9771> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10440<typename dco::mode<stce_9824>::active_value_t> > >::type operator + (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10440<typename dco::mode<stce_9824>::active_value_t> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10440<typename dco::mode<stce_9782>::active_value_t> > >::type operator + (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10440<typename dco::mode<stce_9782>::active_value_t> >(stce_9957, stce_9958); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> > operator - (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10442<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10447<stce_9771> > operator - (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10447<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10448<stce_9771> > operator - (const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10448<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10447<stce_9771> > operator - (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10447<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10448<stce_9771> > operator - (const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10448<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10447<stce_9771> > operator - (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10447<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10448<stce_9771> > operator - (const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10448<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10447<stce_9771> > operator - (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10447<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10448<stce_9771> > operator - (const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10448<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10447<stce_9771> > operator - (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10447<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10448<stce_9771> > operator - (const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10448<stce_9771> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10442<typename dco::mode<stce_9824>::active_value_t> > >::type operator - (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10442<typename dco::mode<stce_9824>::active_value_t> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10442<typename dco::mode<stce_9782>::active_value_t> > >::type operator - (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10442<typename dco::mode<stce_9782>::active_value_t> >(stce_9957, stce_9958); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> > operator * (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10443<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10449<stce_9771> > operator * (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10449<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10450<stce_9771> > operator * (const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10450<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10449<stce_9771> > operator * (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10449<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10450<stce_9771> > operator * (const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10450<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10449<stce_9771> > operator * (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10449<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10450<stce_9771> > operator * (const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10450<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10449<stce_9771> > operator * (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10449<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10450<stce_9771> > operator * (const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10450<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10449<stce_9771> > operator * (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10449<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10450<stce_9771> > operator * (const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10450<stce_9771> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10443<typename dco::mode<stce_9824>::active_value_t> > >::type operator * (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10443<typename dco::mode<stce_9824>::active_value_t> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10443<typename dco::mode<stce_9782>::active_value_t> > >::type operator * (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10443<typename dco::mode<stce_9782>::active_value_t> >(stce_9957, stce_9958); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> > operator / (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10444<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10451<stce_9771> > operator / (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10451<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10452<stce_9771> > operator / (const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10452<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10451<stce_9771> > operator / (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10451<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10452<stce_9771> > operator / (const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10452<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10451<stce_9771> > operator / (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10451<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10452<stce_9771> > operator / (const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10452<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10451<stce_9771> > operator / (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10451<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10452<stce_9771> > operator / (const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10452<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10451<stce_9771> > operator / (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10451<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10452<stce_9771> > operator / (const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10452<stce_9771> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10444<typename dco::mode<stce_9824>::active_value_t> > >::type operator / (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10444<typename dco::mode<stce_9824>::active_value_t> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10444<typename dco::mode<stce_9782>::active_value_t> > >::type operator / (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10444<typename dco::mode<stce_9782>::active_value_t> >(stce_9957, stce_9958); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> > atan2 (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10477<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10478<stce_9771> > atan2 (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10478<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10479<stce_9771> > atan2 (const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10479<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10478<stce_9771> > atan2 (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10478<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10479<stce_9771> > atan2 (const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10479<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10478<stce_9771> > atan2 (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10478<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10479<stce_9771> > atan2 (const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10479<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10478<stce_9771> > atan2 (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10478<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10479<stce_9771> > atan2 (const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10479<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10478<stce_9771> > atan2 (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10478<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10479<stce_9771> > atan2 (const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10479<stce_9771> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10477<typename dco::mode<stce_9824>::active_value_t> > >::type atan2 (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10477<typename dco::mode<stce_9824>::active_value_t> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10477<typename dco::mode<stce_9782>::active_value_t> > >::type atan2 (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10477<typename dco::mode<stce_9782>::active_value_t> >(stce_9957, stce_9958); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> > pow (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10480<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10481<stce_9771> > pow (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10481<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10482<stce_9771> > pow (const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10482<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10481<stce_9771> > pow (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10481<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10482<stce_9771> > pow (const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10482<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10481<stce_9771> > pow (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10481<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10482<stce_9771> > pow (const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10482<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10481<stce_9771> > pow (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10481<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10482<stce_9771> > pow (const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10482<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10481<stce_9771> > pow (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10481<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10482<stce_9771> > pow (const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10482<stce_9771> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10480<typename dco::mode<stce_9824>::active_value_t> > >::type pow (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10480<typename dco::mode<stce_9824>::active_value_t> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10480<typename dco::mode<stce_9782>::active_value_t> > >::type pow (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10480<typename dco::mode<stce_9782>::active_value_t> >(stce_9957, stce_9958); }
template<class stce_9771, class stce_9939 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::active_type<stce_9771, stce_10492> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::active_type<stce_9771, stce_10492>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> > hypot (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489> &stce_9958) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>,dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>, dco::stce_10439::stce_10485<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10486<stce_9771> > hypot (const dco::internal::active_type<stce_9771, stce_9939> &stce_9957, const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10486<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10487<stce_9771> > hypot (const typename dco::internal::active_type<stce_9771, stce_9939>::VALUE_TYPE &stce_9957, const dco::internal::active_type<stce_9771, stce_9939> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9939>, dco::stce_10439::stce_10487<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10486<stce_9771> > hypot (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9957, const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10486<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10487<stce_9771> > hypot (const typename dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>, dco::stce_10439::stce_10487<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10486<stce_9771> > hypot (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10486<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10487<stce_9771> > hypot (const typename dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>, dco::stce_10439::stce_10487<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10486<stce_9771> > hypot (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9957, const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10486<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10487<stce_9771> > hypot (const typename dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>, dco::stce_10439::stce_10487<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10486<stce_9771> > hypot (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9957, const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9958) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10486<stce_9771> >(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10487<stce_9771> > hypot (const typename dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::VALUE_TYPE &stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9958) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>, dco::stce_10439::stce_10487<stce_9771> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10485<typename dco::mode<stce_9824>::active_value_t> > >::type hypot (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9824>::active_value_t, typename dco::mode<stce_9824>::active_t, stce_9824, dco::stce_10439::stce_10485<typename dco::mode<stce_9824>::active_value_t> >(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10485<typename dco::mode<stce_9782>::active_value_t> > >::type hypot (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::internal::stce_9779<typename dco::mode<stce_9782>::active_value_t, stce_9782, typename dco::mode<stce_9782>::active_t, dco::stce_10439::stce_10485<typename dco::mode<stce_9782>::active_value_t> >(stce_9957, stce_9958); }
template<class stce_9771, class stce_9939 > static inline auto operator == (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > static inline auto operator == (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator == (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > static inline auto operator == (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > static inline auto operator == (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > static inline auto operator == (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator == (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator == (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator == (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator == (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > static inline auto operator == (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator == (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator == (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator == (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator == (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > static inline auto operator == (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator == (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator == (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator == (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator == (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > static inline auto operator == (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator == (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator == (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator == (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator == (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > static inline auto operator == (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const typename dco::mode<dco::internal::active_type<stce_9771, stce_9939> >::passive_t& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9959(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9959(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > static inline auto operator == (const typename dco::mode<dco::internal::active_type<stce_9771, stce_9939> >::passive_t& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9960(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9960(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > static inline auto operator == (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9959(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9959(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > static inline auto operator == (const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9960(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9960(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline auto operator == (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9959(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9959(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline auto operator == (const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9960(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9960(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > static inline auto operator == (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9959(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9959(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > static inline auto operator == (const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9960(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9960(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > static inline auto operator == (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9959(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9959(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > static inline auto operator == (const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9960(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9960(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const stce_9782& stce_9958) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9959(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const stce_9782& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9960(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9959(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const stce_9782& stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9958) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9960(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9959(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const stce_9782& stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9958) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9960(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9959(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const stce_9782& stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9958) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9960(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9959(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const stce_9782& stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9958) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9960(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> static inline typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), bool>::type operator == (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::mode<stce_9824>::active_t::data_t::stce_9993::stce_9956(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> static inline typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), bool>::type operator == (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::mode<stce_9782>::active_t::data_t::stce_9993::stce_9956(stce_9957,stce_9958); }
template<class stce_9771, class stce_9939 > static inline auto operator != (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > static inline auto operator != (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator != (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > static inline auto operator != (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > static inline auto operator != (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > static inline auto operator != (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator != (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator != (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator != (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator != (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > static inline auto operator != (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator != (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator != (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator != (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator != (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > static inline auto operator != (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator != (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator != (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator != (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator != (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > static inline auto operator != (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator != (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator != (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator != (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator != (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > static inline auto operator != (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const typename dco::mode<dco::internal::active_type<stce_9771, stce_9939> >::passive_t& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9962(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9962(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > static inline auto operator != (const typename dco::mode<dco::internal::active_type<stce_9771, stce_9939> >::passive_t& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9963(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9963(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > static inline auto operator != (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9962(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9962(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > static inline auto operator != (const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9963(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9963(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline auto operator != (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9962(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9962(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline auto operator != (const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9963(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9963(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > static inline auto operator != (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9962(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9962(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > static inline auto operator != (const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9963(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9963(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > static inline auto operator != (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9962(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9962(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > static inline auto operator != (const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9963(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9963(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const stce_9782& stce_9958) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9962(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const stce_9782& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9963(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9962(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const stce_9782& stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9958) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9963(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9962(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const stce_9782& stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9958) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9963(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9962(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const stce_9782& stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9958) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9963(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9962(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const stce_9782& stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9958) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9963(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> static inline typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), bool>::type operator != (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::mode<stce_9824>::active_t::data_t::stce_9993::stce_9961(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> static inline typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), bool>::type operator != (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::mode<stce_9782>::active_t::data_t::stce_9993::stce_9961(stce_9957,stce_9958); }
template<class stce_9771, class stce_9939 > static inline auto operator < (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > static inline auto operator < (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator < (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > static inline auto operator < (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > static inline auto operator < (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > static inline auto operator < (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator < (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator < (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator < (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator < (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > static inline auto operator < (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator < (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator < (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator < (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator < (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > static inline auto operator < (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator < (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator < (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator < (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator < (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > static inline auto operator < (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator < (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator < (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator < (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator < (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > static inline auto operator < (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const typename dco::mode<dco::internal::active_type<stce_9771, stce_9939> >::passive_t& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9965(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9965(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > static inline auto operator < (const typename dco::mode<dco::internal::active_type<stce_9771, stce_9939> >::passive_t& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9966(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9966(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > static inline auto operator < (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9965(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9965(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > static inline auto operator < (const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9966(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9966(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline auto operator < (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9965(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9965(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline auto operator < (const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9966(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9966(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > static inline auto operator < (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9965(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9965(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > static inline auto operator < (const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9966(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9966(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > static inline auto operator < (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9965(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9965(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > static inline auto operator < (const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9966(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9966(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const stce_9782& stce_9958) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9965(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const stce_9782& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9966(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9965(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const stce_9782& stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9958) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9966(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9965(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const stce_9782& stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9958) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9966(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9965(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const stce_9782& stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9958) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9966(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9965(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const stce_9782& stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9958) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9966(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> static inline typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), bool>::type operator < (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::mode<stce_9824>::active_t::data_t::stce_9993::stce_9964(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> static inline typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), bool>::type operator < (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::mode<stce_9782>::active_t::data_t::stce_9993::stce_9964(stce_9957,stce_9958); }
template<class stce_9771, class stce_9939 > static inline auto operator <= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > static inline auto operator <= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator <= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > static inline auto operator <= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > static inline auto operator <= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > static inline auto operator <= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > static inline auto operator <= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > static inline auto operator <= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > static inline auto operator <= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator <= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > static inline auto operator <= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const typename dco::mode<dco::internal::active_type<stce_9771, stce_9939> >::passive_t& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9968(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9968(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > static inline auto operator <= (const typename dco::mode<dco::internal::active_type<stce_9771, stce_9939> >::passive_t& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9969(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9969(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > static inline auto operator <= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9968(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9968(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > static inline auto operator <= (const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9969(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9969(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline auto operator <= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9968(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9968(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline auto operator <= (const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9969(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9969(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > static inline auto operator <= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9968(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9968(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > static inline auto operator <= (const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9969(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9969(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > static inline auto operator <= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9968(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9968(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > static inline auto operator <= (const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9969(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9969(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const stce_9782& stce_9958) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9968(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const stce_9782& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9969(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9968(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const stce_9782& stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9958) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9969(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9968(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const stce_9782& stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9958) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9969(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9968(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const stce_9782& stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9958) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9969(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9968(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const stce_9782& stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9958) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9969(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> static inline typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), bool>::type operator <= (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::mode<stce_9824>::active_t::data_t::stce_9993::stce_9967(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> static inline typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), bool>::type operator <= (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::mode<stce_9782>::active_t::data_t::stce_9993::stce_9967(stce_9957,stce_9958); }
template<class stce_9771, class stce_9939 > static inline auto operator > (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > static inline auto operator > (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator > (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > static inline auto operator > (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > static inline auto operator > (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > static inline auto operator > (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator > (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator > (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator > (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator > (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > static inline auto operator > (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator > (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator > (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator > (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator > (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > static inline auto operator > (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator > (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator > (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator > (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator > (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > static inline auto operator > (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator > (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator > (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator > (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator > (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > static inline auto operator > (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const typename dco::mode<dco::internal::active_type<stce_9771, stce_9939> >::passive_t& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9971(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9971(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > static inline auto operator > (const typename dco::mode<dco::internal::active_type<stce_9771, stce_9939> >::passive_t& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9972(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9972(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > static inline auto operator > (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9971(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9971(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > static inline auto operator > (const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9972(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9972(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline auto operator > (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9971(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9971(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline auto operator > (const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9972(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9972(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > static inline auto operator > (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9971(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9971(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > static inline auto operator > (const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9972(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9972(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > static inline auto operator > (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9971(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9971(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > static inline auto operator > (const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9972(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9972(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const stce_9782& stce_9958) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9971(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const stce_9782& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9972(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9971(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const stce_9782& stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9958) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9972(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9971(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const stce_9782& stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9958) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9972(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9971(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const stce_9782& stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9958) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9972(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9971(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const stce_9782& stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9958) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9972(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> static inline typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), bool>::type operator > (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::mode<stce_9824>::active_t::data_t::stce_9993::stce_9970(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> static inline typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), bool>::type operator > (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::mode<stce_9782>::active_t::data_t::stce_9993::stce_9970(stce_9957,stce_9958); }
template<class stce_9771, class stce_9939 > static inline auto operator >= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10488, class stce_10489 > static inline auto operator >= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator >= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10490, class stce_10489 > static inline auto operator >= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, class stce_10491, class stce_10489 > static inline auto operator >= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10492 > static inline auto operator >= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10492 > static inline auto operator >= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10492 > static inline auto operator >= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10492 > static inline auto operator >= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::active_type<stce_9771, stce_10492>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10488, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9776<stce_9771, stce_10488, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10491, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9779<stce_9771, stce_10490, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10490, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9780<stce_9771, stce_10490, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, class stce_10491, class stce_10489 > static inline auto operator >= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const dco::internal::stce_9781<stce_9771, stce_10491, stce_10489>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > static inline auto operator >= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const typename dco::mode<dco::internal::active_type<stce_9771, stce_9939> >::passive_t& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9974(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9974(stce_9957,stce_9958); } template<class stce_9771, class stce_9939 > static inline auto operator >= (const typename dco::mode<dco::internal::active_type<stce_9771, stce_9939> >::passive_t& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) -> decltype(dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9975(stce_9957, stce_9958)) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9975(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > static inline auto operator >= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9974(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9974(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943 > static inline auto operator >= (const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9975(stce_9957, stce_9958)) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9975(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline auto operator >= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9974(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9974(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline auto operator >= (const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9975(stce_9957, stce_9958)) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9975(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > static inline auto operator >= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9974(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9974(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943 > static inline auto operator >= (const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9975(stce_9957, stce_9958)) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9975(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > static inline auto operator >= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> >::passive_t& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9974(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9974(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943 > static inline auto operator >= (const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> >::passive_t& stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9958) -> decltype(dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9975(stce_9957, stce_9958)) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9975(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const dco::internal::active_type<stce_9771, stce_9939>& stce_9957, const stce_9782& stce_9958) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9974(stce_9957,stce_9958); } template<class stce_9771, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const stce_9782& stce_9957, const dco::internal::active_type<stce_9771, stce_9939>& stce_9958) { return dco::internal::active_type<stce_9771, stce_9939>::data_t::stce_9993::stce_9975(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9974(stce_9957,stce_9958); } template<class stce_9771, class stce_9949, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const stce_9782& stce_9957, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9958) { return dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>::data_t::stce_9993::stce_9975(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9974(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const stce_9782& stce_9957, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9958) { return dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>::data_t::stce_9993::stce_9975(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9974(stce_9957,stce_9958); } template<class stce_9771, class stce_9941, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const stce_9782& stce_9957, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9958) { return dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>::data_t::stce_9993::stce_9975(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9957, const stce_9782& stce_9958) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9974(stce_9957,stce_9958); } template<class stce_9771, class stce_9942, class stce_9943, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const stce_9782& stce_9957, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9958) { return dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>::data_t::stce_9993::stce_9975(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> static inline typename enable_if<(dco::mode<stce_9824>::order > 0 && (dco::internal::is_gbcp_type<stce_9782>::value || (dco::mode<stce_9782>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9782::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9782>::active_t, typename dco::mode<stce_9824>::active_t>::value ), bool>::type operator >= (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::mode<stce_9824>::active_t::data_t::stce_9993::stce_9973(stce_9957,stce_9958); } template<typename stce_9782, typename stce_9824> static inline typename enable_if<(dco::mode<stce_9782>::order > 0 && (dco::internal::is_gbcp_type<stce_9824>::value || (dco::mode<stce_9824>::is_intermediate_type && dco::internal::is_gbcp_type<typename stce_9824::data_t>::value) ) && dco::helper::is_same<typename dco::mode<stce_9824>::active_t, typename dco::mode<stce_9782>::active_t>::value ), bool>::type operator >= (const stce_9782 &stce_9957, const stce_9824 &stce_9958) { return dco::mode<stce_9782>::active_t::data_t::stce_9993::stce_9973(stce_9957,stce_9958); }
template<class stce_9771, class stce_9939 >
static inline void stce_10493(dco::internal::active_type<stce_9771, stce_9939> &stce_9843) {
    double stce_10063 = 0;
    get(stce_9843, stce_10063);
    stce_9843 = stce_10063;
}
template<class stce_9771, class stce_9939 >
    static inline std::istream &operator >> (std::istream &in, dco::internal::active_type<stce_9771, stce_9939> &stce_9843) {
      stce_9771 &stce_10063 = stce_9843._value();
      in >> stce_10063;
      return in;
    }
    using std::ceil; using std::floor; using std::isfinite; using std::isnan; using std::isnormal; using std::isinf; template<class stce_9771, class stce_9939 > static inline bool isnan(const dco::internal::active_type<stce_9771, stce_9939>& stce_9843) { return isnan(stce_9843._value()); } template<class stce_9771, class stce_9939 > static inline bool isnormal(const dco::internal::active_type<stce_9771, stce_9939>& stce_9843) { return isnormal(stce_9843._value()); } template<class stce_9771, class stce_9939 > static inline bool isinf(const dco::internal::active_type<stce_9771, stce_9939>& stce_9843) { return isinf(stce_9843._value()); } template<class stce_9771, class stce_9939 > static inline double ceil(const dco::internal::active_type<stce_9771, stce_9939>& stce_9843) { return ceil(stce_9843._value()); } template<class stce_9771, class stce_9939 > static inline double floor(const dco::internal::active_type<stce_9771, stce_9939>& stce_9843) { return floor(stce_9843._value()); } template<class stce_9771, class stce_9939 > static inline bool isfinite(const dco::internal::active_type<stce_9771, stce_9939>& stce_9843) { return isfinite(stce_9843._value()); }
    using std::ceil; using std::floor; using std::isfinite; using std::isnan; using std::isnormal; using std::isinf; template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline bool isnan(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9843) { return isnan(stce_9843._value()); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline bool isnormal(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9843) { return isnormal(stce_9843._value()); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline bool isinf(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9843) { return isinf(stce_9843._value()); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline double ceil(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9843) { return ceil(stce_9843._value()); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline double floor(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9843) { return floor(stce_9843._value()); } template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline bool isfinite(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9843) { return isfinite(stce_9843._value()); }
    using std::ceil; using std::floor; using std::isfinite; using std::isnan; using std::isnormal; using std::isinf; template<class stce_9771, class stce_9941, class stce_9943 > static inline bool isnan(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9843) { return isnan(stce_9843._value()); } template<class stce_9771, class stce_9941, class stce_9943 > static inline bool isnormal(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9843) { return isnormal(stce_9843._value()); } template<class stce_9771, class stce_9941, class stce_9943 > static inline bool isinf(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9843) { return isinf(stce_9843._value()); } template<class stce_9771, class stce_9941, class stce_9943 > static inline double ceil(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9843) { return ceil(stce_9843._value()); } template<class stce_9771, class stce_9941, class stce_9943 > static inline double floor(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9843) { return floor(stce_9843._value()); } template<class stce_9771, class stce_9941, class stce_9943 > static inline bool isfinite(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9843) { return isfinite(stce_9843._value()); }
    using std::ceil; using std::floor; using std::isfinite; using std::isnan; using std::isnormal; using std::isinf; template<class stce_9771, class stce_9942, class stce_9943 > static inline bool isnan(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9843) { return isnan(stce_9843._value()); } template<class stce_9771, class stce_9942, class stce_9943 > static inline bool isnormal(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9843) { return isnormal(stce_9843._value()); } template<class stce_9771, class stce_9942, class stce_9943 > static inline bool isinf(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9843) { return isinf(stce_9843._value()); } template<class stce_9771, class stce_9942, class stce_9943 > static inline double ceil(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9843) { return ceil(stce_9843._value()); } template<class stce_9771, class stce_9942, class stce_9943 > static inline double floor(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9843) { return floor(stce_9843._value()); } template<class stce_9771, class stce_9942, class stce_9943 > static inline bool isfinite(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9843) { return isfinite(stce_9843._value()); }
    using std::ceil; using std::floor; using std::isfinite; using std::isnan; using std::isnormal; using std::isinf; template<class stce_9771, class stce_9949, class stce_9943 > static inline bool isnan(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9843) { return isnan(stce_9843._value()); } template<class stce_9771, class stce_9949, class stce_9943 > static inline bool isnormal(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9843) { return isnormal(stce_9843._value()); } template<class stce_9771, class stce_9949, class stce_9943 > static inline bool isinf(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9843) { return isinf(stce_9843._value()); } template<class stce_9771, class stce_9949, class stce_9943 > static inline double ceil(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9843) { return ceil(stce_9843._value()); } template<class stce_9771, class stce_9949, class stce_9943 > static inline double floor(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9843) { return floor(stce_9843._value()); } template<class stce_9771, class stce_9949, class stce_9943 > static inline bool isfinite(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9843) { return isfinite(stce_9843._value()); }
    template<class stce_9771, class stce_9939 > static inline std::ostream& operator << (std::ostream& out, const dco::internal::active_type<stce_9771, stce_9939>& stce_9843) { out << stce_9843._value(); return out; }
    template<class stce_9771, class stce_9941, class stce_9942, class stce_9943 > static inline std::ostream& operator << (std::ostream& out, const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943>& stce_9843) { out << stce_9843._value(); return out; }
    template<class stce_9771, class stce_9941, class stce_9943 > static inline std::ostream& operator << (std::ostream& out, const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943>& stce_9843) { out << stce_9843._value(); return out; }
    template<class stce_9771, class stce_9942, class stce_9943 > static inline std::ostream& operator << (std::ostream& out, const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943>& stce_9843) { out << stce_9843._value(); return out; }
    template<class stce_9771, class stce_9949, class stce_9943 > static inline std::ostream& operator << (std::ostream& out, const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943>& stce_9843) { out << stce_9843._value(); return out; }
template <typename stce_9782, typename stce_9824, typename stce_10494=decltype(passive_value(std::declval<stce_9782>()) < passive_value(std::declval<stce_9824>()))>
  typename enable_if<mode<stce_9782>::is_dco_type || mode<stce_9824>::is_dco_type,
                     typename stce_9867<(mode<stce_9782>::order > mode<stce_9824>::order), typename mode<stce_9782>::type, typename mode<stce_9824>::type>::type
                                       >::type min(stce_9782 const& stce_9905, stce_9824 const& stce_10495) {
  if (passive_value(stce_9905) < passive_value(stce_10495)) return stce_9905;
  else return stce_10495;
}
template <typename stce_9782, typename stce_9824, typename stce_10494=decltype(passive_value(std::declval<stce_9782>()) > passive_value(std::declval<stce_9824>()))>
  typename enable_if<mode<stce_9782>::is_dco_type || mode<stce_9824>::is_dco_type,
                     typename stce_9867<(mode<stce_9782>::order > mode<stce_9824>::order), typename mode<stce_9782>::type, typename mode<stce_9824>::type>::type
                                             >::type max(stce_9782 const& stce_9905, stce_9824 const& stce_10495) {
  if (passive_value(stce_9905) > passive_value(stce_10495)) return stce_9905;
  else return stce_10495;
}
}
}
namespace dco {
  namespace instrument {
    extern bool dco_instrument;
    struct stce_10496 {
      struct iterator_t {
        DCO_TAPE_INT _progvarcounter;
        DCO_TAPE_INT _stackcounter;
        iterator_t(DCO_TAPE_INT stce_10497, DCO_TAPE_INT stce_10498) : _progvarcounter(stce_10497), _stackcounter(stce_10498) {}
        DCO_TAPE_INT stce_10499() {
          return _progvarcounter + _stackcounter * 2;
        }
        DCO_TAPE_INT &stce_10500() {
          return _progvarcounter;
        }
      };
      struct stce_10501 {
        DCO_TAPE_INT id;
        std::string stce_10502;
        DCO_TAPE_INT stce_10503;
        iterator_t stce_10164;
        DCO_TAPE_INT stce_10301;
        DCO_TAPE_INT stce_10504;
        DCO_TAPE_INT stce_10133;
        DCO_TAPE_INT stce_10134;
        DCO_TAPE_INT stce_10505;
        DCO_TAPE_INT stce_10506;
        stce_10501 *stce_10403;
        stce_10501(stce_10496 *tape, const std::string &name, stce_10501 *stce_10507) : id(++tape->stce_9857),
          stce_10502(name),
          stce_10503(-1),
          stce_10164(tape->get_position()),
          stce_10301(0),
          stce_10504(0),
          stce_10133(0),
          stce_10134(0),
          stce_10505(0),
          stce_10506(0),
          stce_10403(stce_10507) {
        }
      };
    public:
      struct callback_object_t {
        virtual ~callback_object_t() {
        }
        virtual void stce_10088(stce_9848 &stce_9849) {
          (void)stce_9849;
        }
      };
      struct interpretation_settings {
      };
      typedef void (*stce_10508)(stce_10496 &stce_10403, const interpretation_settings &stce_10311, callback_object_t *stce_10404);
      void stce_10509(stce_10508 stce_10510, callback_object_t *stce_10404) {
        (void) stce_10510;
        (void) stce_10404;
      }
      static stce_10496 *get_tape() {
        return NULL;
      }
      DCO_TAPE_INT tape_index() const {
        return 0;
      }
      DCO_TAPE_INT stce_9857;
      std::ofstream stce_10511;
      std::ofstream stce_10512;
      std::ofstream stce_10513;
      std::ofstream stce_10514;
      std::ofstream stce_10515;
      std::ofstream stce_10516;
      DCO_TAPE_INT _progvarcounter;
      DCO_TAPE_INT _stackcounter;
      stce_10501 *stce_10169;
      std::vector<stce_10501 *> stce_10517;
      std::stack<stce_10501 *> stce_10518;
      DCO_TAPE_INT stce_10519;
      stce_10496() : stce_9857(0),
        _progvarcounter(0),
        _stackcounter(0),
        stce_10169(NULL),
        stce_10519(0) { }
      ~stce_10496() {
        _finish_instrument();
      }
      bool is_active() {
        return dco_instrument;
      }
      void switch_to_active() {
        dco_instrument = true;
      }
      void switch_to_passive() {
        dco_instrument = false;
      }
      static void remove(stce_10496 *stce_9782) {
        delete stce_9782;
      }
      static stce_10496 *create(const DCO_TAPE_INT size = 0, const DCO_TAPE_INT stce_10520 = 0) {
        (void) size;
        (void) stce_10520;
        stce_10496 *stce_9840 = new stce_10496();
        stce_9840->_init_instrument();
        return stce_9840;
      }
      static stce_10496 *create(std::string stce_10521) {
        stce_10496 *stce_9840 = new stce_10496();
        stce_9840->_init_instrument(stce_10521);
        return stce_9840;
      }
      template<class stce_9771, class stce_9939>
      void register_variable(dco::internal::active_type<stce_9771, stce_9939> &stce_9843) {
        _progvarcounter++;
        stce_9843.tape_index = _progvarcounter;
        stce_9843.stce_10522 = stce_9857;
      }
      template<class stce_9771, class stce_9939>
      void register_variable(dco::internal::active_type<stce_9771, stce_9939> *stce_9843, DCO_TAPE_INT stce_10102) {
        for (DCO_TAPE_INT stce_9832 = 0; stce_9832 < stce_10102; stce_9832++)
          this->register_variable(stce_9843[stce_9832]);
      }
      template<class stce_9771, class stce_9939>
      void register_variable(std::vector<dco::internal::active_type<stce_9771, stce_9939> > &stce_9843) {
        register_variable(&(stce_9843[0]), stce_9843.size());
      }
      template<class stce_9771, class stce_9939>
      void register_output_variable(dco::internal::active_type<stce_9771, stce_9939> &stce_9843) {
        stce_9843 = 1 * stce_9843;
      }
      template<class stce_9771, class stce_9939>
      void register_output_variable(dco::internal::active_type<stce_9771, stce_9939> *stce_9843, const DCO_TAPE_INT stce_10102) {
        if (stce_10102 == 0) return;
        for (DCO_TAPE_INT stce_9832 = 0; stce_9832 < stce_10102; ++stce_9832)
          stce_9843[stce_9832] = 1 * stce_9843[stce_9832];
      }
      template<class stce_9771, class stce_9939>
      void register_output_variable(std::vector<dco::internal::active_type<stce_9771, stce_9939> > &stce_9843) {
        register_output_variable(&(stce_9843[0]), stce_9843.size());
      }
      void _init_instrument(const std::string &stce_10521 = "") {
        stce_9857 = 0;
        std::string stce_10523 = stce_10521;
        if (stce_10521.size() != 0) {
          stce_10523 += "_";
        }
        stce_10511.open((stce_10523 + "nodes.dot").c_str());
        stce_10512.open((stce_10523 + "edges.dot").c_str());
        stce_10511 << "digraph {" << std::endl;
        stce_10513.open((stce_10523 + "nodes.csv").c_str());
        stce_10514.open((stce_10523 + "edges.csv").c_str());
        stce_10515.open((stce_10523 + "ext_nodes.csv").c_str());
        stce_10516.open((stce_10523 + "ext_edges.csv").c_str());
        dco_instrument = true;
        std::string name = "ROOT";
        std::string stce_10403 = "caller";
        _dco_instrument_enter_function(name, stce_10403);
      }
      void _finish_instrument() {
        std::string name = "ROOT";
        std::string stce_10403 = "caller";
        _dco_instrument_leave_function(name, stce_10403);
        for (size_t stce_9832 = 0; stce_9832 < stce_10517.size(); ++stce_9832) {
          stce_10501 *stce_9817 = stce_10517[stce_9832];
          stce_10511 << stce_9817->id << "[shape=box,label=\"" << stce_9817->stce_10502 << "("
                    << stce_9817->stce_10301 << "," << stce_9817->stce_10133 << "," << stce_9817->stce_10134 << ")\"]" << std::endl;
          stce_10513 << stce_9817->stce_10505 << " "
                    << stce_9817->stce_10301 << " "
                    << stce_9817->stce_10506 << " "
                    << stce_9817->stce_10133 << " "
                    << stce_9817->stce_10134 << " "
                    << stce_9817->stce_10502 << std::endl;
          delete stce_9817;
        }
        stce_10517.clear();
        stce_10512 << "}" << std::endl;
        stce_10511.close();
        stce_10512.close();
        stce_10513.close();
        stce_10514.close();
        stce_10515.close();
        stce_10516.close();
        dco_instrument = false;
      }
      iterator_t get_position() {
        return iterator_t(_progvarcounter, _stackcounter);
      }
      void interpret_adjoint() { }
      void interpret_adjoint_to(iterator_t &stce_10164) {
        (void) stce_10164;
      }
      void reset() {
        std::cerr << "dco instrument -- possible error due to reset() call" << std::endl;
      }
      void reset_to(iterator_t &stce_10164) {
        std::cerr << "dco instrument -- possible error due to reset() call" << std::endl;
        (void) stce_10164;
      }
      void _dco_instrument_enter_function(const std::string &stce_10502 , const std::string &stce_10403 ) {
        (void) stce_10502;
        (void) stce_10403;
        stce_10501 *stce_9817 = new stce_10501(this, stce_10502, stce_10169);
        stce_10517.push_back(stce_9817);
        stce_10169 = stce_9817;
        if (! stce_10518.empty()) {
          stce_10501 *stce_10386 = stce_10518.top();
          DCO_TAPE_INT size = get_position().stce_10499() - stce_10386->stce_10164.stce_10499();
          stce_10512 << stce_10386->id << " -> " << stce_9817->id << "[label=\"" << size << "\"]" << std::endl;
          stce_10514 << stce_10386->id << " " << stce_9817->id << " " << size << std::endl;
          stce_10386->stce_10505++;
          stce_10386->stce_10301 += size;
          stce_10386->stce_10504 += size;
          stce_10386->stce_10164 = get_position();
        }
        stce_10518.push(stce_9817);
        stce_10519 = stce_9817->stce_10164.stce_10500();
        stce_9817->stce_10503 = stce_10519;
      }
      void _dco_instrument_leave_function(const std::string &stce_10502 , const std::string &stce_10403 ) {
        (void) stce_10502;
        (void) stce_10403;
        stce_10501 *stce_9817 = stce_10518.top();
        DCO_TAPE_INT size = get_position().stce_10499() - stce_9817->stce_10164.stce_10499();
        stce_9817->stce_10301 += size;
        stce_9817->stce_10504 += size;
        stce_10518.pop();
        if (!stce_10518.empty()) {
          stce_10501 *stce_10524 = stce_10518.top();
          stce_10524->stce_10164 = get_position();
          stce_10524->stce_10504 += stce_9817->stce_10504;
          stce_10519 = stce_10524->stce_10503;
          stce_10169 = stce_10524;
        }
      }
    };
  template<class stce_9771, stce_10496 *&stce_10525>
    struct stce_10526 {
    typedef internal::stce_9953 stce_9993;
    typedef stce_9803<stce_9771> mode_t;
  stce_9771 stce_10527;
  DCO_TAPE_INT tape_index;
  DCO_TAPE_INT _tape_index;
  mutable DCO_TAPE_INT stce_10522;
  mutable stce_10496::stce_10501 *stce_10403;
  stce_10526 &operator =(const stce_10526 &stce_9894) {
    tape_index = stce_9894.tape_index;
    stce_10522 = stce_9894.stce_10522;
    stce_10403 = stce_9894.stce_10403;
    _tape_index = stce_9894._tape_index;
    return *this;
  }
  stce_10526() :
      stce_10527(0.0), tape_index(0), _tape_index(0), stce_10522(0), stce_10403(0) {
    if (stce_10525)
      stce_10522 = stce_10525->stce_9857;
  }
  const stce_9771 &get_derivative() const {
    return stce_10527;
  }
  stce_9771 &get_derivative() {
    return stce_10527;
  }
  void clear() {
    tape_index = 0;
    stce_10522 = 0;
  }
  DCO_TAPE_INT stce_10528() const {
    return stce_9755 && tape_index == 0 ? 0 : 1;
  }
  template<class stce_9939 >
  static DCO_TAPE_INT stce_10529(const dco::internal::active_type<stce_9771, stce_9939> &stce_9843) {
    const stce_10526 &data = stce_9843;
    if (data.tape_index != 0) {
      if (data.stce_10403 != 0) {
        if (stce_10525->stce_10169 != data.stce_10403) {
          data.stce_10403->stce_10134++;
          data.stce_10403 = 0;
        }
      }
      if (data.tape_index <= stce_10525->stce_10519) {
        stce_10496::stce_10501 *stce_10530 = stce_10525->stce_10518.top();
        while (stce_10530) {
          if (data.stce_10522 < stce_10530->id) {
            stce_10530->stce_10133++;
            stce_10530 = stce_10530->stce_10403;
          } else {
            stce_10530 = NULL;
          }
        }
        data.stce_10522 = stce_10525->stce_10518.top()->id;
      }
    }
    return stce_9843.stce_10528();
  }
  template<class stce_9941, class stce_9942, class stce_9943 >
  static DCO_TAPE_INT stce_10529(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843) {
    return stce_10529(stce_9843.stce_9945) + stce_10529(stce_9843.stce_9947);
  }
  template<class stce_9949, class stce_9943 >
  static DCO_TAPE_INT stce_10529(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843) {
    return stce_10529(stce_9843.stce_9950);
  }
  template<class stce_9941, class stce_9943 >
  static DCO_TAPE_INT stce_10529(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843) {
    return stce_10529(stce_9843.stce_9945);
  }
  template<class stce_9942, class stce_9943 >
  static DCO_TAPE_INT stce_10529(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843) {
    return stce_10529(stce_9843.stce_9947);
  }
  template<class stce_9939 >
  static DCO_TAPE_INT stce_10531(const dco::internal::active_type<stce_9771, stce_9939> &stce_9843) {
    (void) stce_9843;
    return 0;
  }
  template<class stce_9941, class stce_9942, class stce_9943 >
  static DCO_TAPE_INT stce_10531(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843) {
    return stce_10531(stce_9843.stce_9945) + stce_10531(stce_9843.stce_9947) + 1;
  }
  template<class stce_9949, class stce_9943 >
  static DCO_TAPE_INT stce_10531(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843) {
    return stce_10531(stce_9843.stce_9950) + 1;
  }
  template<class stce_9941, class stce_9943 >
  static DCO_TAPE_INT stce_10531(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843) {
    return stce_10531(stce_9843.stce_9945) + 1;
  }
  template<class stce_9942, class stce_9943 >
  static DCO_TAPE_INT stce_10531(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843) {
    return stce_10531(stce_9843.stce_9947) + 1;
  }
  template<class stce_9892>
  DCO_ALWAYS_INLINE void stce_9893(const stce_9892 &stce_9894) {
    DCO_TAPE_INT stce_10532 = stce_10529(stce_9894);
    if (dco_instrument && (stce_10525->stce_10169 != NULL)) {
      DCO_TAPE_INT stce_10533 = stce_10531(stce_9894);
      stce_10525->stce_10169->stce_10506 += (stce_10533);
    }
    if (NULL != stce_10525 && !stce_10525->is_active()) {
      clear();
    }
    if (stce_10532 <= 0) {
      clear();
      return;
    }
    stce_10525->_stackcounter += (stce_10532 + 1);
    stce_10525->_progvarcounter++;
    tape_index = stce_10525->_progvarcounter;
    stce_10403 = stce_10525->stce_10169;
  }
};
  }
}
namespace dco {
  namespace internal {
    template<class stce_9771>
    struct stce_10534 {
    private:
      struct stce_10535 {
        DCO_TAPE_INT stce_10532;
        DCO_TAPE_INT stce_10536;
      };
      stce_10535 *stce_10537;
      DCO_TAPE_INT stce_10538;
      int *stce_10539;
      DCO_TAPE_INT stce_10540;
      stce_10534(const stce_10534<stce_9771> &) {
        throw dco::exception::create<std::runtime_error>("dco_error: Not implemented yet! Please report, why needed!", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape_blob_pattern.hpp", 30);
      }
      stce_10534(DCO_TAPE_INT stce_10532, DCO_TAPE_INT stce_10497) {
        stce_10537 = new stce_10535[stce_9752(stce_10497)];
        stce_10538 = 0;
        stce_10539 = new int[size_t(stce_10532)];
        stce_10540 = 0;
      }
    public:
      struct stce_10541 {
        friend struct stce_10534<stce_9771>;
      private:
        DCO_TAPE_INT stce_10542;
        DCO_TAPE_INT stce_10543;
        stce_10541(DCO_TAPE_INT stce_10544, DCO_TAPE_INT stce_10545) : stce_10542(stce_10544), stce_10543(stce_10545) {}
      public:
        stce_10541() : stce_10542(0), stce_10543(0) {};
        DCO_TAPE_INT stce_10546() {
          return stce_10542;
        }
        DCO_TAPE_INT stce_10547() {
          if (stce_10543 < 0) return -stce_10543;
          return stce_10543;
        }
        bool stce_10548() {
          return stce_10543 < 0;
        }
      };
      struct iterator_t {
      private:
        DCO_TAPE_INT _progvarcounter;
      public:
        DCO_TAPE_INT stce_10497() const {
          return _progvarcounter;
        }
        iterator_t(const DCO_TAPE_INT stce_10549) : _progvarcounter(stce_10549) {}
      };
      iterator_t get_position() {
        return iterator_t(stce_10538);
      }
      typedef std::vector<stce_10541> sparse_jacobi;
      template<class stce_9782>
      sparse_jacobi *build_pattern(stce_9782 *stce_10546, DCO_TAPE_INT stce_10102, const iterator_t &stce_10550) {
        sparse_jacobi *stce_9840 = new sparse_jacobi();
        for (DCO_TAPE_INT stce_9832 = 0; stce_9832 < stce_10102; ++stce_9832) {
          stce_10551(stce_9832 + 1, stce_10546[stce_9832].tape_index, stce_10550, stce_9840);
        }
        return stce_9840;
      }
      void create_pattern_file(sparse_jacobi *stce_10552, std::string filename) {
        std::ofstream out(filename.c_str());
        for (size_t stce_9832 = 0; stce_9832 < stce_10552->size(); ++stce_9832) {
          stce_10541 &stce_10553 = (*stce_10552)[stce_9832];
          out << stce_10553.stce_10546() << " " << stce_10553.stce_10547() << " 1" << std::endl;
        }
        out.close();
      }
      void stce_10554(std::string filename) {
        std::ofstream out(filename.c_str());
        out << "digraph {" << std::endl;
        for (DCO_TAPE_INT stce_9832 = 1; stce_9832 <= stce_10538; ++stce_9832) {
          out << stce_9832 << std::endl;
        }
        for (DCO_TAPE_INT stce_9832 = 1; stce_9832 <= stce_10538; ++stce_9832) {
          int stce_10532 = stce_10537[stce_9832].stce_10532;
          int stce_10555 = stce_10537[stce_9832].stce_10536;
          for (int stce_10556 = 0; stce_10556 < stce_10532; ++stce_10556) {
            int stce_10054 = stce_10539[stce_10556 + stce_10555];
            out << stce_9832 << " -> " << stce_10054 << std::endl;
          }
        }
        out << "}" << std::endl;
        out.close();
      }
      void stce_10551(DCO_TAPE_INT stce_10544, DCO_TAPE_INT stce_9857, const iterator_t &stce_10550, sparse_jacobi *stce_10557) {
        if (stce_9857 == 0) return;
        if (stce_9857 <= stce_10550.stce_10497()) {
          stce_10557->push_back(stce_10541(stce_10544, stce_9857));
        } else {
          DCO_TAPE_INT stce_10532 = stce_10537[stce_9857].stce_10532;
          DCO_TAPE_INT stce_10555 = stce_10537[stce_9857].stce_10536;
          for (int stce_9832 = 0; stce_9832 < stce_10532; ++stce_9832) {
            stce_10551(stce_10544, stce_10539[stce_9832 + stce_10555], stce_10550, stce_10557);
          }
        }
      }
      static stce_10534 *create(DCO_TAPE_INT size, DCO_TAPE_INT stce_10520 = 0) {
        if (stce_10520 == 0) stce_10520 = size / 2;
        return new stce_10534(size, stce_10520);
      }
      template<class stce_9939 >
      void register_variable(dco::internal::active_type<stce_9771, stce_9939> &stce_10558) {
        stce_9939 &data = stce_10558.data();
        stce_10538++;
        data.tape_index = stce_10538;
      }
      template<stce_10534 *&global_tape>
      struct stce_10559 {
        typedef stce_9803<stce_9771> mode_t;
        DCO_TAPE_INT tape_index;
        stce_10559() : tape_index(0) {}
        void clear() {
          tape_index = 0;
        }
        template<class stce_9939 >
        static void interpret(const dco::internal::active_type<stce_9771, stce_9939> &stce_9843, bool stce_10548) {
          (void) stce_10548;
          global_tape->stce_10539[global_tape->stce_10540] = stce_9843.tape_index;
          global_tape->stce_10540++;
        }
        template<class stce_9941, class stce_9942, class stce_9943 >
        static void interpret(const dco::internal::stce_9779<stce_9771, stce_9941, stce_9942, stce_9943> &stce_9843, bool stce_10548) {
          interpret(stce_9843.stce_9945, stce_10548);
          interpret(stce_9843.stce_9947, stce_10548);
        }
        template<class stce_9949, class stce_9943 >
        static void interpret(const dco::internal::stce_9776<stce_9771, stce_9949, stce_9943> &stce_9843, bool stce_10548) {
          interpret(stce_9843.stce_9950, stce_10548);
        }
        template<class stce_9941, class stce_9943 >
        static void interpret(const dco::internal::stce_9780<stce_9771, stce_9941, stce_9943> &stce_9843, bool stce_10548) {
          interpret(stce_9843.stce_9945, stce_10548);
        }
        template<class stce_9942, class stce_9943 >
        static void interpret(const dco::internal::stce_9781<stce_9771, stce_9942, stce_9943> &stce_9843, bool stce_10548) {
          interpret(stce_9843.stce_9947, stce_10548);
        }
        template<class stce_9892, class stce_9794>
        static void stce_9893(const stce_9892 &stce_9894, stce_9794 &stce_10560) {
          const int stce_10536 = global_tape->stce_10540;
          interpret(stce_9894, false);
          stce_10559 &data = const_cast<stce_10559 &>(stce_10560);
          if (stce_10536 == global_tape->stce_10540) {
            data.clear();
          } else {
            global_tape->stce_10538++;
            data.tape_index = global_tape->stce_10538;
            global_tape->stce_10537[global_tape->stce_10538].stce_10532 = global_tape->stce_10540 - stce_10536;
            global_tape->stce_10537[global_tape->stce_10538].stce_10536 = stce_10536;
          }
        }
        template<class stce_9892>
        void stce_9893(const stce_9892 &stce_9894) {
          stce_9893(stce_9894, *this);
        }
      };
    };
  }
}
namespace dco {
  namespace instrument {
    typedef stce_10496 tape_t;
    extern tape_t *global_tape;
    typedef stce_10526<double, global_tape> data;
    typedef internal::active_type<double, data> type;
    typedef helper::stce_10108<type, tape_t> external_adjoint_object_t;
  }
  namespace p1t {
    typedef internal::stce_10534<double > tape_t;
    extern tape_t *global_tape;
    typedef tape_t::stce_10559<global_tape> data;
    typedef internal::active_type<double, data> type;
  }
  namespace p1f {
    typedef internal::stce_10060<32, double> data;
    typedef internal::active_type<double, data> type;
  }
  namespace p1f {
    inline void set(type &stce_9843, const bool stce_9753, const DCO_TAPE_INT stce_9808)
    {
      data &stce_10063 = stce_9843;
      stce_10063.stce_10062[stce_9808] = stce_9753;
    }
     inline void get(const type &stce_9843, bool &stce_9753, const DCO_TAPE_INT stce_9808)
    {
      stce_9753 = stce_9843.stce_10062[stce_9808];
    }
     inline void get(type &stce_9843, double &stce_9753)
    {
      stce_9753 = stce_9843._value();
    }
     inline std::string stce_10561(const type &stce_9843)
    {
      return stce_9843.stce_10062.to_string<char, std::char_traits<char>, std::allocator<char> >();
    }
  }
  namespace p1t {
     inline void set(type &stce_9843, const double stce_9753)
    {
      double &value = const_cast<double &>(stce_9843._value());
      value = stce_9753;
    }
     inline void get(type &stce_9843, double &stce_9753)
    {
      stce_9753 = stce_9843._value();
    }
  }
  namespace instrument {
  inline void set(type &stce_9843, const double stce_9753, const DCO_TAPE_INT what = 0)
    {
      stce_1000002();;
      if (what == 0) {
        double &value = const_cast<double &>(stce_9843._value());
        value = stce_9753;
      } else if (what == -1) {
      } else {
        throw exception::create<std::invalid_argument>("invalid kind of what (musst be 0 or -1)");
      }
    }
  inline void get(const type &stce_9843, double &stce_9753, const DCO_TAPE_INT what = 0)
    {
      stce_1000002();;
      if (what == 0)
        stce_9753 = stce_9843._value();
      else if (what == -1) {
        stce_9753 = 0;
      } else {
        throw exception::create<std::invalid_argument>("invalid kind of what (musst be 0 or -1)");
      }
    }
  }
}
       
       
namespace dco {
namespace internal {
template<class stce_10423>
  struct jacobian_preaccumulator_t: stce_10079<typename dco::mode<stce_10423>::stce_10020> {
  typedef typename dco::mode<stce_10423>::stce_10020 stce_9764;
  typedef typename dco::mode<stce_10423>::tape_t tape_t;
  typedef typename tape_t::iterator_t stce_10562;
  typedef typename tape_t::stce_9924 stce_10563;
  tape_t* const stce_9952;
  stce_10562 stce_10564;
  DCO_TAPE_INT stce_10565;
  typedef stce_9793<stce_9764, tape_t> vector_t;
  vector_t stce_10290;
  bool stce_10566;
  void init() {}
  typedef stce_10149<memory_model::BLOB_TAPE, tape_t> stce_10283;
  typedef typename stce_10283::stce_9924 stce_10567;
  typedef typename stce_10283::iterator stce_10568;
  stce_10283* stce_10145;
  stce_10567* stce_9927;
  typedef typename stce_10423::data_t stce_10569;
  std::vector<stce_10569*> stce_10570;
  DCO_TAPE_INT stce_10309;
  jacobian_preaccumulator_t(tape_t * const tape) :
    stce_9952(tape),
    stce_10290(0, *tape, tape->stce_10170()),
    stce_10145(0),
    stce_9927(0),
    stce_10309(0)
  {
    if (stce_9952) {
      stce_10564 = stce_9952->get_position();
      stce_10565 = stce_10564.index();
      stce_10566 = true;
    } else {
      stce_10566 = false;
    }
  }
  ~jacobian_preaccumulator_t() {
    if (stce_10566)
      this->finish();
    if (stce_10145)
      delete stce_10145;
  }
  void start() {
    if (stce_9952) {
      stce_10564 = stce_9952->get_position();
      stce_10565 = stce_10564.index();
      stce_10566 = true;
    } else {
      stce_10566 = false;
    }
  }
  void stce_10571(stce_10423 &stce_9843) {
    if (dco::tape_index(stce_9843) && dco::tape_index(stce_9843) <= stce_10565) {
      throw dco::exception::create<std::runtime_error>("jacobian_preaccumulator_t: registering output, which is indeed an input (was alive already before starting point).");
    }
    stce_10570.push_back(&static_cast<stce_10569&>(stce_9843));
  }
  void register_output(stce_10423 &stce_9843)
  {
    if (stce_9952) { stce_10571(stce_9843); }
  }
  void register_output(std::vector<stce_10423> &stce_9843)
  {
    if (stce_9952) {
      stce_10570.reserve(stce_10570.size() + stce_9843.size());
      for (size_t stce_9832 = 0; stce_9832 < stce_9843.size(); ++stce_9832) {
        stce_10571(stce_9843[stce_9832]);
      }
    }
  }
  void finish() {
    if (!stce_10566) return;
    stce_10566 = false;
    if (!stce_9952) return;
    if (!stce_10570.size()) {
      stce_9952->reset_to(stce_10564);
      return;
    }
    const stce_10562 stce_10572 = stce_9952->get_position();
    const DCO_TAPE_INT stce_10573 = stce_10572.index() - stce_10564.index();
    stce_10290.resize(stce_10573 + 1);
    const DCO_TAPE_INT stce_10574 = stce_9952->size(stce_10572, stce_10564);
    const DCO_TAPE_INT stce_10575 = stce_10574 * (1 + static_cast<DCO_TAPE_INT>(stce_10570.size()));
    if (stce_10145) delete stce_10145;
    tape_options stce_10437;
    stce_10437.blob_size_in_byte() = stce_10575 * stce_10283::entry::size_of();
    stce_10145 = new stce_10283(stce_10437, *stce_9952);
    stce_10145->stce_10200().index() = stce_10572.index();
    stce_10290.stce_10306();
    for(typename std::vector<stce_10569*>::iterator stce_10164 = stce_10570.begin(); stce_10164 != stce_10570.end(); ++stce_10164) {
      stce_10569& stce_9817 = **stce_10164;
      if (stce_9817._tape_index() == 0) continue;
      stce_10290.stce_10303(stce_10573-1,0);
      stce_10049(stce_9817._tape_index()) = 1.0;
      stce_9927 = new stce_10567(stce_10145->stce_9928());
      stce_9952->stce_10316(*this,
                                           stce_10572, stce_10564,
                                           interpretation_settings(false, true));
      delete stce_9927;
    }
    stce_9952->reset_to(stce_10564);
    stce_10568 stce_10234 = stce_10145->current();
    for (typename std::vector<stce_10569*>::reverse_iterator stce_10164 = stce_10570.rbegin();
         stce_10164 != stce_10570.rend(); ++stce_10164) {
      stce_10563 stce_10576(stce_9952->stce_9928());
      for (;; --stce_10234) {
        stce_10576.insert(stce_10234->arg(), stce_10234->pval());
        if(stce_10234->stce_10157())
          break;
      }
      --stce_10234;
      (*stce_10164)->_tape_index() = stce_10576.index();
    }
    stce_10570.clear();
  }
  stce_9764& stce_10049(const DCO_TAPE_INT stce_9808) {
    return stce_10290[stce_9808 - stce_10565];
  }
  template<class stce_10315>
    void stce_10081(const DCO_TAPE_INT &stce_9808, const stce_10315 &stce_9985,
                            const stce_9764 &stce_10082) {
    stce_10081(stce_9808, stce_9985 * stce_10082);
  }
  void stce_10303(DCO_TAPE_INT stce_9862, DCO_TAPE_INT stce_9863) {
    stce_10290.stce_10303(stce_9862,stce_9863);
  }
  void stce_10295(DCO_TAPE_INT stce_9808) {
    stce_10290.stce_10295(stce_9808 - stce_10565);
  }
  stce_9764& operator[](const DCO_TAPE_INT stce_9808) {
    return stce_10290[stce_9808 - stce_10565];
  }
  template <typename iterator_t>
  void stce_10296(iterator_t& stce_9817) {
    stce_10081(stce_9817->arg(), stce_9817->pval()*(*this)[stce_9817.index()]);
  }
  virtual stce_9764 stce_10080(const DCO_TAPE_INT stce_9808) {
    return stce_10290[stce_9808 - stce_10565];
  }
  virtual void stce_10081(const DCO_TAPE_INT stce_9808, const stce_9764 &stce_10082) {
    if (stce_10082 == 0)
      return;
    if (stce_9808 <= stce_10565) {
      stce_9927->insert(stce_9808, stce_10082);
    } else {
      stce_10049(stce_9808) += stce_10082;
    }
  }
};
template<>
struct jacobian_preaccumulator_t<void> {
  jacobian_preaccumulator_t(void*) {}
  void start() {}
  template<typename stce_9782> void register_output(const stce_9782&) {}
  void finish() {}
};
}
}
namespace dco {
template <typename stce_9782,
          typename stce_9760 = stce_9782,
          typename stce_9761 = stce_9760,
          memory_model::TYPE MEMORY_MODEL = DCO_TAPE_MEMORY_MODEL>
  class ga1s {
  typedef stce_9762<stce_9782, stce_9760, stce_9761> stce_10018;
public:
  static const memory_model::TYPE TAPE_MODEL = MEMORY_MODEL;
  typedef dco::internal::single_tape_data<ga1s> data_t;
  typedef typename stce_10018::value_t value_t;
  typedef typename stce_10018::value_t active_value_t;
  typedef typename stce_10018::stce_9764 derivative_t;
  typedef derivative_t stce_10020;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef internal::tape<TAPE_MODEL, stce_10018> tape_t;
  static DCOi_THREAD_LOCAL tape_t *global_tape;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10084<tape_t> callback_object_t;
  typedef dco::helper::stce_10090<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10108<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
template <class stce_9782, class stce_9760, class stce_9761, memory_model::TYPE MEMORY_MODEL>
  DCOi_THREAD_LOCAL typename ga1s<stce_9782,stce_9760,stce_9761,MEMORY_MODEL>::tape_t *ga1s<stce_9782,stce_9760,stce_9761,MEMORY_MODEL>::global_tape = 0;
}
       
namespace dco {
template <class stce_9782, class stce_9760 = stce_9782, class stce_9761 = stce_9760>
  class ga1s_mod {
  typedef stce_9762<stce_9782, stce_9760, stce_9761> stce_10018;
public:
  typedef dco::internal::single_tape_data<ga1s_mod> data_t;
  static const memory_model::TYPE TAPE_MODEL = memory_model::stce_9770;
  static const memory_model::TYPE stce_10577 = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10018::value_t value_t;
  typedef typename stce_10018::value_t active_value_t;
  typedef typename stce_10018::stce_9764 derivative_t;
  typedef derivative_t stce_10020;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef dco::internal::tape<TAPE_MODEL, stce_10018, stce_10577> tape_t;
  static DCOi_THREAD_LOCAL tape_t *global_tape;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10084<tape_t> callback_object_t;
  typedef dco::helper::stce_10090<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10108<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
template <class stce_9782, class stce_9760, class stce_9761>
  DCOi_THREAD_LOCAL typename ga1s_mod<stce_9782,stce_9760,stce_9761>::tape_t *ga1s_mod<stce_9782,stce_9760,stce_9761>::global_tape = 0;
}
       
namespace dco {
template <typename stce_9782,
          typename stce_9760 = stce_9782,
          typename stce_9761 = stce_9760>
class gtas {
  typedef stce_9762<stce_9782, stce_9760, stce_9761> stce_10018;
public:
  typedef dco::internal::stce_10001<gtas> data_t;
  static const memory_model::TYPE TAPE_MODEL = memory_model::stce_9769;
  typedef typename stce_10018::value_t value_t;
  typedef typename stce_10018::value_t active_value_t;
  typedef typename stce_10018::stce_9764 derivative_t;
  typedef derivative_t stce_10020;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = true;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef internal::tape<TAPE_MODEL, stce_10018> tape_t;
  static tape_t *global_tape;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10084<tape_t> callback_object_t;
  typedef dco::helper::stce_10090<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10108<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
  typedef typename tape_t::logger logger;
  enum stce_10578 {
    check_branch = 1,
    check_fd = 2,
    check_identity = 4,
    stce_10579 = check_branch | check_fd | check_identity
  };
  static double& fd_eps() {
    assert(global_tape);
    return global_tape->stce_10252;
  }
  static double& ta_eps() {
    assert(global_tape);
    return global_tape->stce_10253;
  }
  static void init(const int stce_10398 = stce_10579) {
    assert(global_tape);
    global_tape->template stce_10268<data_t>();
    global_tape->stce_10249 = (stce_10398 & check_branch) != 0;
    global_tape->stce_10250 = (stce_10398 & check_fd) != 0;
    global_tape->stce_10010 = (stce_10398 & check_identity) != 0;
    global_tape->stce_10255 = 0;
    global_tape->stce_10256 = -1;
  }
  static DCO_TAPE_INT& throw_on_event() {
    assert(global_tape);
    return global_tape->stce_10256;
  }
  static bool& throw_on_warning() {
    assert(global_tape);
    return global_tape->stce_10251;
  }
};
template<typename stce_9782, typename stce_10030, typename stce_10580> typename gtas<stce_9782, stce_10030, stce_10580>::tape_t *gtas<stce_9782, stce_10030, stce_10580>::global_tape = 0;
}
       
namespace dco {
template <class stce_9782, class stce_9760 = stce_9782, class stce_9761 = stce_9760>
  class ga1sm {
  typedef stce_9762<stce_9782, stce_9760, stce_9761> stce_10018;
public:
  typedef dco::internal::stce_9996<ga1sm> data_t;
  static const memory_model::TYPE TAPE_MODEL = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10018::value_t value_t;
  typedef typename stce_10018::value_t active_value_t;
  typedef typename stce_10018::stce_9764 derivative_t;
  typedef derivative_t stce_10020;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef internal::tape<TAPE_MODEL, stce_10018> tape_t;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10084<tape_t> callback_object_t;
  typedef dco::helper::stce_10090<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10108<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
}
       
namespace dco {
template <class stce_9782, class stce_9760 = stce_9782, class stce_9761 = stce_9760>
  class ga1sm_mod {
  typedef stce_9762<stce_9782, stce_9760, stce_9761> stce_10018;
public:
  typedef dco::internal::stce_9996<ga1sm_mod> data_t;
  static const memory_model::TYPE TAPE_MODEL = memory_model::stce_9770;
  static const memory_model::TYPE stce_10577 = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10018::value_t value_t;
  typedef typename stce_10018::value_t active_value_t;
  typedef typename stce_10018::stce_9764 derivative_t;
  typedef derivative_t stce_10020;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef dco::internal::tape<TAPE_MODEL, stce_10018, stce_10577> tape_t;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10084<tape_t> callback_object_t;
  typedef dco::helper::stce_10090<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10108<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
}
namespace dco {
template <class stce_9782, const int stce_10581=1, class stce_9760=stce_9782, class stce_9761=stce_9760>
  class ga1vm {
  typedef stce_9762<stce_9782, stce_9760, dco::helper::stce_9827<stce_9761, stce_10581> > stce_10018;
public:
  typedef dco::internal::stce_9996<ga1vm> data_t;
  static const memory_model::TYPE TAPE_MODEL = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10018::value_t value_t;
  typedef typename stce_10018::value_t active_value_t;
  typedef typename stce_10018::stce_9764 derivative_t;
  typedef stce_9761 stce_10020;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef internal::tape<TAPE_MODEL, stce_10018> tape_t;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10084<tape_t> callback_object_t;
  typedef dco::helper::stce_10090<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10108<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
}
namespace dco {
template <class stce_9782, const int stce_10581 = 1, class stce_9760 = stce_9782, class stce_9761 = stce_9760>
  class ga1vm_mod {
  typedef stce_9762<stce_9782, stce_9760, dco::helper::stce_9827<stce_9761, stce_10581> > stce_10018;
public:
  typedef dco::internal::stce_9996<ga1vm_mod> data_t;
  static const memory_model::TYPE TAPE_MODEL = memory_model::stce_9770;
  static const memory_model::TYPE stce_10577 = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10018::value_t value_t;
  typedef typename stce_10018::value_t active_value_t;
  typedef typename stce_10018::stce_9764 derivative_t;
  typedef typename derivative_t::value_t stce_10020;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef dco::internal::tape<TAPE_MODEL, stce_10018, stce_10577> tape_t;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10084<tape_t> callback_object_t;
  typedef dco::helper::stce_10090<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10108<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
}
namespace dco {
template <class stce_9782, const int stce_10581 = 1, class stce_9760 = stce_9782, class stce_9761 = stce_9760>
  class ga1v {
  typedef stce_9762<stce_9782, stce_9760, dco::helper::stce_9827<stce_9761, stce_10581> > stce_10018;
public:
  typedef dco::internal::single_tape_data<ga1v> data_t;
  static const memory_model::TYPE TAPE_MODEL = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10018::value_t value_t;
  typedef typename stce_10018::value_t active_value_t;
  typedef typename stce_10018::stce_9764 derivative_t;
  typedef stce_9761 stce_10020;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef internal::tape<TAPE_MODEL, stce_10018> tape_t;
  static DCOi_THREAD_LOCAL tape_t *global_tape;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10084<tape_t> callback_object_t;
  typedef dco::helper::stce_10090<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10108<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
template <class stce_9782, const int stce_10581, class stce_9760, class stce_9761>
  DCOi_THREAD_LOCAL typename ga1v<stce_9782,stce_10581,stce_9760,stce_9761>::tape_t *ga1v<stce_9782,stce_10581,stce_9760,stce_9761>::global_tape;
}
namespace dco {
template <class stce_9782, const int stce_10581 = 1, class stce_9760 = stce_9782, class stce_9761 = stce_9760>
  class ga1v_mod {
  typedef stce_9762<stce_9782, stce_9760, dco::helper::stce_9827<stce_9761, stce_10581> > stce_10018;
public:
  typedef dco::internal::single_tape_data<ga1v_mod> data_t;
  static const memory_model::TYPE TAPE_MODEL = memory_model::stce_9770;
  static const memory_model::TYPE stce_10577 = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10018::value_t value_t;
  typedef typename stce_10018::value_t active_value_t;
  typedef typename stce_10018::stce_9764 derivative_t;
  typedef typename derivative_t::value_t stce_10020;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef dco::internal::tape<TAPE_MODEL, stce_10018, stce_10577> tape_t;
  static DCOi_THREAD_LOCAL tape_t *global_tape;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10084<tape_t> callback_object_t;
  typedef dco::helper::stce_10090<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10108<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
template <class stce_9782, const int stce_10581, class stce_9760, class stce_9761>
  DCOi_THREAD_LOCAL typename ga1v_mod<stce_9782,stce_10581,stce_9760,stce_9761>::tape_t *ga1v_mod<stce_9782,stce_10581,stce_9760,stce_9761>::global_tape;
}
namespace dco
{
  template <typename stce_10582>
  struct ginstrument {
    typedef instrument::type type;
    typedef type active_t;
    typedef instrument::tape_t tape_t;
    typedef instrument::data data;
    typedef instrument::external_adjoint_object_t external_adjoint_object_t;
    typedef double active_value_t;
    typedef double value_t;
    typedef double passive_t;
    static tape_t* &global_tape;
  };
  template <typename stce_10582> typename ginstrument<stce_10582>::tape_t *&ginstrument<stce_10582>::global_tape = instrument::global_tape;
}
       
namespace dco {
template <class stce_9782, typename stce_10583 = stce_9782>
  class gt1s {
public:
  typedef dco::internal::ts_data<gt1s> data_t;
  typedef stce_9782 value_t;
  typedef stce_9782 active_value_t;
  typedef stce_10583 derivative_t;
  typedef stce_10583 stce_10020;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<stce_9782, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = false;
  static const bool is_tangent_type = true;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef void tape_t;
  typedef void local_gradient_t;
  typedef void external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<void> jacobian_preaccumulator_t;
};
}
       
namespace dco {
template <class stce_9782,
          int stce_10584=1>
  class gt1v {
public:
    typedef dco::internal::stce_10056<gt1v> data_t;
    typedef stce_9782 value_t;
    typedef stce_9782 active_value_t;
    typedef dco::helper::stce_9827<stce_9782, stce_10584> derivative_t;
    typedef stce_9782 stce_10020;
    typedef typename mode<value_t>::passive_t passive_t;
    typedef dco::internal::active_type<stce_9782, data_t> type;
    typedef type active_t;
    static const bool is_dco_type = true;
    static const bool is_adjoint_type = false;
    static const bool is_tangent_type = true;
    static const int order = dco::mode<stce_9782>::order + 1;
    typedef void tape_t;
    typedef void local_gradient_t;
    typedef void external_adjoint_object_t;
    typedef internal::jacobian_preaccumulator_t<void> jacobian_preaccumulator_t;
  };
}
namespace dco {
template <typename active_t, typename stce_10585 = void>
  struct gbcp;
template <typename stce_10423>
  struct gbcp<stce_10423, typename dco::enable_if<(dco::mode<stce_10423>::is_dco_type && !internal::is_gbcp_type<stce_10423>::value)>::type> {
  typedef stce_10423 binary_t;
  typedef typename dco::mode<binary_t>::value_t active_t;
  typedef active_t active_value_t;
  typedef typename dco::mode<active_value_t>::value_t value_t;
  typedef typename mode<active_value_t>::passive_t passive_t;
  typedef internal::jacobian_preaccumulator_t<void> jacobian_preaccumulator_t;
  typedef dco::internal::stce_10017<gbcp> data_t;
  typedef dco::internal::active_type<active_value_t, data_t> type;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = dco::mode<active_value_t>::is_adjoint_type;
  static const bool is_tangent_type = dco::mode<active_value_t>::is_tangent_type;
  static const int gbcp_order = 1;
  static const int order = dco::mode<binary_t>::order - 1;
};
template <typename stce_10586>
  struct gbcp<stce_10586, typename dco::enable_if<(dco::mode<stce_10586>::is_dco_type && internal::is_gbcp_type<stce_10586>::value)>::type> {
  typedef stce_10586 binary_t;
  typedef typename dco::mode<binary_t>::active_t active_t;
  typedef typename dco::mode<binary_t>::active_value_t active_value_t;
  typedef typename dco::mode<binary_t>::value_t value_t;
  typedef typename dco::mode<binary_t>::passive_t passive_t;
  typedef dco::internal::stce_10017<gbcp> data_t;
  typedef dco::internal::active_type<active_value_t, data_t> type;
  typedef internal::jacobian_preaccumulator_t<void> jacobian_preaccumulator_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = dco::mode<binary_t>::is_adjoint_type;
  static const bool is_tangent_type = dco::mode<binary_t>::is_tangent_type;
  static const int gbcp_order = dco::mode<binary_t>::gbcp_order+1;
  static const int order = dco::mode<binary_t>::order - 1;
};
template <typename stce_10352>
  struct gbcp<stce_10352, typename dco::enable_if<(!dco::mode<stce_10352>::is_dco_type )>::type> {
  typedef internal::jacobian_preaccumulator_t<void> jacobian_preaccumulator_t;
  typedef stce_10352 binary_t;
  typedef binary_t active_t;
  typedef binary_t active_value_t;
  typedef binary_t value_t;
  typedef binary_t passive_t;
  typedef binary_t type;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = false;
  static const bool is_tangent_type = false;
  static const int gbcp_order = 1;
  static const int order = 0;
};
}
       
namespace dco {
  typedef ga1s<double> a1w;
  typedef double w_rtype;
}
namespace dco {
  namespace instrument {
    class instrumentation_data {
      const std::string stce_10203;
      const std::string stce_10587;
      bool stce_10588;
    public:
      instrumentation_data(std::string stce_10589) : stce_10203(stce_10589), stce_10587("caller"), stce_10588(false) {
        if (dco::instrument::dco_instrument && dco::instrument::global_tape) {
          dco::instrument::global_tape->_dco_instrument_enter_function(stce_10203, stce_10587);
          stce_10588 = true;
        }
      }
      ~instrumentation_data() {
        if (stce_10588) {
          dco::instrument::global_tape->_dco_instrument_leave_function(stce_10203, stce_10587);
        }
      }
    };
    inline void stce_10590(std::string stce_10589) {
      std::string stce_10587 = "caller";
      dco::instrument::global_tape->_dco_instrument_enter_function(stce_10589, stce_10587);
    }
    inline void stce_10591(std::string stce_10589) {
      std::string stce_10587 = "caller";
      dco::instrument::global_tape->_dco_instrument_leave_function(stce_10589, stce_10587);
    }
  }
}

