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
}
extern "C" {
  extern int stce_1000002(int stce_9804=false);
  extern void* stce_1000001(size_t *stce_9805, int);
  extern void* stce_999999(size_t, int);
  extern int stce_1000003(int *stce_9806);
  extern long stce_1000004(long *stce_9806);
  extern void stce_1000000(void* stce_9807);
}
template <typename stce_9782> inline stce_9782 stce_9808(stce_9782 &stce_9806);
template<> inline int stce_9808<int>(int &stce_9806) {
  return stce_1000003(&stce_9806);
}
template<>
  inline long stce_9808<long>(long &stce_9806) {
  return stce_1000004(&stce_9806);
}
namespace dco {
template<memory_model::TYPE MEMORY_MODEL, class stce_9788, memory_model::TYPE stce_9789>
  bool is_null(internal::tape<MEMORY_MODEL, stce_9788, stce_9789> *tape) {
  return (tape == NULL);
}
 static unsigned int stce_9809=0xABCDE;
 static unsigned int stce_9810(unsigned int *stce_9811, unsigned int stce_9812) {
    int stce_9813;
    stce_9813 = *stce_9811 & 1;
    *stce_9811 >>=1;
    if(stce_9813==1)
        *stce_9811=stce_9812;
    return *stce_9811;
 }
 static inline double stce_9814() {
      stce_9810(&stce_9809, 0xB4BCD35C);
      unsigned int stce_9815 = stce_9810(&stce_9809, 0xB4BCD35C);
      return double(stce_9815)/RAND_MAX;
 }
  template <class stce_9782>
  std::string
  type_name()
  {
    return typeid(stce_9782).name();
  }
namespace internal {
template <typename stce_9782, bool>
struct stce_9816 { typedef stce_9782 type; };
template <typename stce_9782>
  struct stce_9816<stce_9782, true> { typedef stce_9782& type; };
}
namespace helper {
static bool stce_9817 = true;
template<typename stce_9782>
  struct stce_9818 {
  static const bool stce_9753 = true;
};
template<class stce_9782, class stce_9819>
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
template<class stce_9820, const int stce_9821=1>
  struct stce_9822 {
  static const int stce_9823 = stce_9821;
  typedef stce_9820 value_t;
  stce_9820 stce_9824[stce_9821];
  void stce_9825(const stce_9820 &stce_9826) {
    for (DCO_TAPE_INT stce_9827 = 0; stce_9827 < stce_9821; ++stce_9827)
      stce_9824[stce_9827] = stce_9826;
  }
  stce_9822(const stce_9820 &stce_9826) {
    stce_9825(stce_9826);
  }
  stce_9822() {
    stce_9825(0);
  }
  stce_9820 &operator[](const DCO_TAPE_INT stce_9806) {
    return stce_9824[stce_9806];
  }
  const stce_9820 &operator[](const DCO_TAPE_INT stce_9806) const {
    return stce_9824[stce_9806];
  }
  stce_9822 &operator= (const stce_9820 &stce_9828) {
    for (int stce_9827 = 0; stce_9827 < stce_9821; ++stce_9827)
      this->stce_9824[stce_9827] = stce_9828;
    return *this;
  }
  stce_9822 &operator += (const stce_9822 &stce_9828) {
    for (int stce_9827 = 0; stce_9827 < stce_9821; ++stce_9827)
      this->stce_9824[stce_9827] += stce_9828.stce_9824[stce_9827];
    return *this;
  }
};
template <typename stce_9820, const int stce_9821>
  bool operator==(stce_9822<stce_9820, stce_9821> const& stce_9829, int const& stce_9830) {
  for (int stce_9827 = 0; stce_9827 < stce_9821; ++stce_9827) {
    if (stce_9829[stce_9827] != stce_9830) {
      return false;
    }
  }
  return true;
}
template <typename stce_9820, typename stce_9831, const int stce_9821>
  bool operator==(stce_9822<stce_9820, stce_9821> const& stce_9829, stce_9822<stce_9831, stce_9821> const& stce_9830) {
  for (int stce_9827 = 0; stce_9827 < stce_9821; ++stce_9827) {
    if (stce_9829[stce_9827] != stce_9830[stce_9827]) {
      return false;
    }
  }
  return true;
}
template<class stce_9782>
struct stce_9832 {
  static const int stce_9833;
  typedef stce_9782 stce_9834;
};
template<class stce_9782>
const int stce_9832<stce_9782>::stce_9833 = 1;
template<class stce_9820, int stce_9821>
struct stce_9832<stce_9822<stce_9820, stce_9821> > {
  static const int stce_9833;
  typedef stce_9820 stce_9834;
};
template<class stce_9820, int stce_9821>
const int stce_9832<stce_9822<stce_9820, stce_9821> > ::stce_9833 = stce_9821;
template<class stce_9820, int stce_9821>
  static inline
  stce_9822<stce_9820, stce_9821>
  operator + (const stce_9822<stce_9820, stce_9821> &stce_9835, const stce_9822<stce_9820, stce_9821> &stce_9828) {
  stce_9822<stce_9820, stce_9821> stce_9836;
  for (int stce_9827 = 0; stce_9827 < stce_9821; ++stce_9827)
    stce_9836.stce_9824[stce_9827] = stce_9835.stce_9824[stce_9827] + stce_9828.stce_9824[stce_9827];
  return stce_9836;
}
template<class stce_9820, int stce_9821>
static inline stce_9822<stce_9820, stce_9821> operator *(const stce_9822<stce_9820, stce_9821> &stce_9835,
        const stce_9822<stce_9820, stce_9821> &stce_9828)
{
    stce_9822<stce_9820, stce_9821> stce_9836;
    for (DCO_TAPE_INT stce_9827 = 0; stce_9827 < stce_9821; ++stce_9827)
        stce_9836.stce_9824[stce_9827] = stce_9835.stce_9824[stce_9827] * stce_9828.stce_9824[stce_9827];
    return stce_9836;
}
template<class stce_9820, class stce_9831, int stce_9821>
static inline stce_9822<stce_9820, stce_9821> operator*(const stce_9831 &stce_9837,
        const stce_9822<stce_9820, stce_9821> &stce_9838)
{
    stce_9822<stce_9820, stce_9821> stce_9836;
    for (DCO_TAPE_INT stce_9827 = 0; stce_9827 < stce_9821; ++stce_9827)
        stce_9836.stce_9824[stce_9827] = stce_9838.stce_9824[stce_9827] * stce_9837;
    return stce_9836;
}
template<class stce_9820, class stce_9831, int stce_9821>
static inline stce_9822<stce_9820, stce_9821> operator *(const stce_9822<stce_9820, stce_9821> &stce_9835,
        const stce_9831 &stce_9828)
{
    stce_9822<stce_9820, stce_9821> stce_9836;
    for (DCO_TAPE_INT stce_9827 = 0; stce_9827 < stce_9821; ++stce_9827)
        stce_9836.stce_9824[stce_9827] = stce_9835.stce_9824[stce_9827] * stce_9828;
    return stce_9836;
}
template<class stce_9820, int stce_9821>
  static inline std::ostream &operator << (std::ostream &out, const stce_9822<stce_9820, stce_9821> &stce_9839) {
  out << "(";
  for (DCO_TAPE_INT stce_9827 = 0; stce_9827 < stce_9821 - 1; ++stce_9827) {
    out << stce_9839.stce_9824[stce_9827] << ";";
  }
  out << stce_9839.stce_9824[stce_9821 - 1];
  out << ")";
  return out;
}
}
template<typename TYPE>
  class stce_9840 {
  const TYPE &stce_9782;
public:
  stce_9840(const TYPE &stce_9841) : stce_9782(stce_9841) {}
  template <typename stce_9842>
    stce_9842 *create() const {
    return new stce_9842(stce_9782);
  }
};
template<>
  class stce_9840<void *> {
public:
  stce_9840(void *stce_9843) {
    (void) stce_9843;
  }
  template <typename stce_9842>
    stce_9842 *create() const {
    return new stce_9842();
  }
};
class stce_9844 {
  std::ofstream stce_9845;
  std::string stce_9846;
  std::string stce_9847;
  std::stringstream stce_9848;
  std::stringstream stce_9849;
public:
  stce_9844(const std::string &filename) : stce_9846("digraph {\n"), stce_9847("\n}\n") {
    stce_9848.clear();
    stce_9849.clear();
    stce_9845.open(filename.c_str());
  }
  void stce_9850(std::vector<DCO_TAPE_INT> &stce_9851, std::vector<DCO_TAPE_INT> &stce_9852, DCO_TAPE_INT stce_9853) {
    stce_9848 << stce_9853 << "[shape=box,label=\"(" << stce_9853 << ", ext. function)\",color=\"red\"]\n";
    for (std::vector<DCO_TAPE_INT>::iterator stce_9827 = stce_9851.begin(); stce_9827 != stce_9851.end(); ++stce_9827) {
      stce_9849 << *stce_9827 << " -> " << stce_9853 << "[label=\"unknown" << "\",color=\"red\"]\n";
    }
    for (std::vector<DCO_TAPE_INT>::iterator stce_9827 = stce_9852.begin(); stce_9827 != stce_9852.end(); ++stce_9827) {
      stce_9849 << stce_9853 << " -> " << *stce_9827 << "[label=\"unknown" << "\",color=\"red\"]\n";
    }
  }
  void stce_9854(DCO_TAPE_INT stce_9853) {
    if(stce_9853 == 0){
      return;
    }
    stce_9848 << stce_9853 << "[shape=box,label=\"" << stce_9853 << "\", color=\"green\"]\n";
  }
  template <typename stce_9855>
  void stce_9854(DCO_TAPE_INT stce_9853, stce_9855 stce_9856) {
    if(stce_9853 == 0){
      return;
    }
    stce_9848 << stce_9853 << "[shape=box,label=\"(" << stce_9853 << ", " << stce_9856 << ")\", color=\"green\"]\n";
  }
  template <typename stce_9855>
  void stce_9857(DCO_TAPE_INT stce_9858, DCO_TAPE_INT stce_9859, stce_9855 stce_9860) {
    if(stce_9858 == 0 || stce_9859 == 0){
      return;
    }
    stce_9849 << stce_9858 << " -> " << stce_9859 << "[label=\"" << stce_9860 << "\"]\n";
  }
  ~stce_9844() {
    stce_9845 << stce_9846;
    stce_9845 << stce_9848.str();
    stce_9845 << stce_9849.str();
    stce_9845 << stce_9847;
    stce_9845.close();
  }
};
template<bool stce_9861, class stce_9782, class stce_9862> struct stce_9863 { typedef stce_9782 type; };
template< class stce_9782, class stce_9862> struct stce_9863<false, stce_9782, stce_9862> { typedef stce_9862 type; };
template<bool stce_9861, typename stce_9782 = void> struct enable_if {};
template< typename stce_9782 > struct enable_if<true, stce_9782> { typedef stce_9782 type; };
template <typename stce_9782> struct type_identity { typedef stce_9782 type; };
template<class stce_9782> struct remove_const { typedef stce_9782 type; };
template<class stce_9782> struct remove_const<const stce_9782> { typedef stce_9782 type; };
template<class stce_9782> struct remove_const<const stce_9782&> { typedef stce_9782& type; };
template<class stce_9782> struct stce_9864 { typedef stce_9782 type; };
template<class stce_9782> struct stce_9864<stce_9782&> { typedef stce_9782 type; };
template<typename stce_9782> struct remove_all { typedef stce_9782 type; };
template<typename stce_9782> struct remove_all<const stce_9782> { typedef typename remove_all<stce_9782>::type type; };
template<typename stce_9782> struct remove_all<stce_9782 const&> { typedef typename remove_all<stce_9782>::type type; };
template<typename stce_9782> struct remove_all<stce_9782&> { typedef typename remove_all<stce_9782>::type type; };
template<typename stce_9782> struct remove_all<stce_9782 const*> { typedef typename remove_all<stce_9782>::type type; };
template<typename stce_9782> struct remove_all<stce_9782*> { typedef typename remove_all<stce_9782>::type type; };
template<typename stce_9782> struct stce_9865 { enum { value = false }; };
template<typename stce_9782> struct stce_9865<stce_9782*> { enum { value = true }; };
template <typename stce_9866, typename stce_9867>
  struct stce_9868 {
  static const bool value = (stce_9866::order < stce_9867::order);
};
template <>
  struct stce_9868<void, void> {
  static const bool value = false;
};
static inline void *alloc(size_t size, bool stce_9869 = false) {
  return stce_999999(size, static_cast<int>(stce_9869));
}
static inline void dealloc(void* stce_9807) {
  stce_1000000(stce_9807);
}
}
       
namespace dco {
  namespace internal {
    template<class stce_9771, class stce_9870 = stce_9771>
    struct stce_9871 {
      DCO_TAPE_INT stce_9872;
      DCO_TAPE_INT stce_9873;
      DCO_TAPE_INT _progvarcounter;
      std::string stce_9874;
      std::ofstream stce_9875;
      std::ofstream stce_9876;
    private:
      stce_9871(const std::string stce_9877) : stce_9872(0), stce_9873(0), _progvarcounter(0), stce_9874(stce_9877) {
      }
    public:
      static stce_9871 *create(const std::string stce_9877, const std::string filename) {
        const char *stce_9846 = "//=================================================//\n"
                             "// This file is generated by dco/c++.              //\n"
                             "// Further inquiries info@stce.rwth-aachen.de      //\n"
                             "//=================================================//\n\n"
                             "#include <cmath>\n\n"
                             "typedef double DCO_TAPE_REAL;\n\n";
        stce_9871 *stce_9836 = new stce_9871(stce_9877);
        stce_9836->stce_9875.open(filename.c_str());
        stce_9836->stce_9875 << stce_9846 << "void " << stce_9877 << "(const double *x, double *x_a1s, double *y, const double *y_a1s) {\n double va0_a1s=0; \n";
        stce_9836->stce_9876.open("ad_tmp");
        return stce_9836;
      }
      static void remove(stce_9871 *tape) {
        tape->stce_9876.close();
        std::ifstream stce_9878("ad_tmp");
        std::vector<std::string> stce_9879;
        std::string stce_9880;
        while (std::getline(stce_9878, stce_9880)) {
          stce_9879.insert(stce_9879.begin(), stce_9880);
        }
        tape->stce_9875 << "// BEGIN OF ADJOINT SECTION\n";
        for (size_t stce_9827 = 0; stce_9827 < stce_9879.size(); ++stce_9827) {
          tape->stce_9875 << stce_9879[stce_9827] << "\n";
        }
        tape->stce_9875 << "}\n";
        tape->stce_9875.close();
      }
      template<class stce_9881>
      void register_variable(stce_9881 &stce_9882) {
        _progvarcounter++;
        stce_9882._tape_index() = _progvarcounter;
        stce_9875 << "const double &va" << _progvarcounter << "=x[" << stce_9872 << "];\n";
        stce_9875 << "double &va" << _progvarcounter << "_a1s=x_a1s[" << stce_9872 << "];\n";
        stce_9872++;
      }
      template<class stce_9881>
      void register_output_variable(stce_9881 &stce_9882) {
        stce_9875 << "y[" << stce_9873 << "]=va" << stce_9882._tape_index() << ";\n";
        stce_9875 << "va" << stce_9882._tape_index() << "_a1s+=y_a1s[" << stce_9873 << "];\n";
        stce_9873++;
      }
    };
    class stce_9883
    {
      template<typename>
      friend class stce_9884;
      DCO_TAPE_INT stce_9885;
      stce_9883(): stce_9885(0) {}
    public:
      std::string stce_9886() const {
        std::stringstream stce_9807;
        stce_9807 << "va" << stce_9885;
        return stce_9807.str();
      }
      DCO_TAPE_INT& _tape_index() {
        return stce_9885;
      }
      DCO_TAPE_INT _tape_index() const {
        return stce_9885;
      }
      void clear() {
        _tape_index() = 0;
      }
    };
    template<class stce_9887>
    class stce_9884: public stce_9883{
    public:
      typedef stce_9887 mode_t;
      stce_9884(){}
      template<class stce_9888>
      DCO_ALWAYS_INLINE void stce_9889(const stce_9888 &stce_9890) {
        stce_9890.stce_9891(stce_9887::global_tape->stce_9875);
        _tape_index() = stce_9887::global_tape->_progvarcounter;
        std::stringstream stce_9892;
        stce_9892 << "va" << ++stce_9887::global_tape->_progvarcounter;
       stce_9887::global_tape->stce_9876 << stce_9890.stce_9886() << "_a1s+=" <<
            stce_9892.str() << "_a1s;\t";
        stce_9890.stce_9893(stce_9887::global_tape->stce_9876);
        stce_9887::global_tape->stce_9876 << "//EOL\n";
        _tape_index()++;
        stce_9887::global_tape->stce_9875 << "double " << stce_9886()
            << "=" << stce_9890.stce_9886() << ";\n";
        stce_9887::global_tape->stce_9875 << "double "
            << stce_9886() << "_a1s=0.0;\n";
      }
    };
  }
}
namespace dco {
namespace internal {
static inline bool is_character(const std::string &stce_9894, size_t stce_9806) {
  if(stce_9806 >=stce_9894.size()) return false;
  return std::isalpha(stce_9894[stce_9806]) != 0 ;
}
static inline std::string replace_all_in_string(std::string& str, const std::string &stce_9858, const std::string &stce_9859) {
  size_t stce_9895=0;
  while ((stce_9895 = str.find(stce_9858, stce_9895)) != std::string::npos) {
    if((stce_9895 > 0 && is_character(str, stce_9895-1))
      || is_character(str, stce_9895 + stce_9858.length())) {
        stce_9895++;
        continue;
    }
    str.replace(stce_9895, stce_9858.length(), stce_9859);
    stce_9895 += stce_9859.length();
  }
  return str;
}
template<class stce_9771>
class stce_9896 {
  static DCO_TAPE_INT stce_9897;
  DCO_TAPE_INT stce_9898;
public:
  stce_9896() :
      stce_9898(stce_9897++) {
  }
   std::string stce_9886() const {
    std::stringstream stce_9807;
    stce_9807 << "im" << stce_9898;
    return stce_9807.str();
  }
protected:
  void stce_9899(std::ofstream &out, std::string stce_9900) const {
    out << "const double " << stce_9886() << "=" << stce_9900 << ";\n";
    out << "double " << stce_9886() << "_a1s=0.0;\n";
  }
  template<class stce_9782>
  std::string str(stce_9782 arg) const {
    return arg.stce_9886();
  }
  std::string str(stce_9771 stce_9901) const {
    std::stringstream str;
    str << std::scientific << std::setprecision(17) << stce_9901;
    return str.str();
  }
  template <class stce_9902>
  void stce_9903(std::string& stce_9900, const stce_9902& arg, const std::string stce_9904) const {
    replace_all_in_string(stce_9900, stce_9904, str(arg));
  }
  void stce_9903(std::string& stce_9900, const stce_9771& arg, const std::string stce_9904) const {
    replace_all_in_string(stce_9900, stce_9904, str(arg));
  }
  template <class stce_9775, class stce_9905, class stce_9906>
  void stce_9907(std::ofstream &out, const stce_9905& stce_9908, const stce_9906& stce_9909) const {
    std::string stce_9900 = stce_9775::stce_9910();
    this->stce_9903(stce_9900, stce_9908, "arg1");
    this->stce_9903(stce_9900, stce_9909, "arg2");
    out << this->str(stce_9908) << "_a1s +=(" << stce_9900 << ")*" << this->stce_9886() << "_a1s;\t";
    stce_9908.stce_9893(out);
  }
  template <class stce_9775, class stce_9906>
  void stce_9907(std::ofstream&, const stce_9771&, const stce_9906&) const {}
  template <class stce_9775, class stce_9905, class stce_9906>
  void stce_9911(std::ofstream &out, const stce_9905& stce_9908, const stce_9906& stce_9909) const {
    std::string stce_9900 = stce_9775::stce_9912();
    this->stce_9903(stce_9900, stce_9908, "arg1");
    this->stce_9903(stce_9900, stce_9909, "arg2");
    out << this->str(stce_9909) << "_a1s +=(" << stce_9900 << ")*" << this->stce_9886() << "_a1s;\t";
    stce_9909.stce_9893(out);
  }
  template <class stce_9775, class stce_9905>
  void stce_9911(std::ofstream&, const stce_9905&, const stce_9771&) const {}
};
template<class stce_9771>
DCO_TAPE_INT stce_9896<stce_9771>::stce_9897 = 0;
}
}
       
       
       
namespace dco {
  class exception {
  public:
    template <typename stce_9913>
    static stce_9913 create(std::string stce_9914, std::string stce_9915 = "", int stce_9880 = 0) {
      std::stringstream stce_9916;
      stce_9916 << "--- dco/c++ --- " << stce_9914;
      if (stce_9915 != "") stce_9916 << " --- " << stce_9915 << ":" << stce_9880 << ".";
      DCO_LOG(dco::logERROR) << "EXCEPTION thrown: " << stce_9916.str();
      return stce_9913(stce_9916.str());
    }
    template <typename stce_9913>
    static stce_9913 stce_9917(std::string stce_9914, std::string stce_9915 = "", int stce_9880 = 0) {
      std::stringstream stce_9916;
      stce_9916 << "--- dco/c++ --- " << stce_9914;
      if (stce_9915 != "") stce_9916 << " --- " << stce_9915 << ":" << stce_9880 << ".";
      DCO_LOG(dco::logERROR) << "EXCEPTION thrown: " << stce_9916.str();
      return stce_9913();
    }
};
}
       
namespace dco {
namespace folding {
template <typename stce_9782, typename enable_if = void>
  struct is_zero_trait {
  static bool get(const stce_9782& stce_9839) {
    return stce_9839 == 0;
  }
};
template <typename stce_9782>
  bool is_zero(const stce_9782& stce_9839) {
  return is_zero_trait<stce_9782>::get(stce_9839);
}
template <typename stce_9782>
  struct is_zero_trait<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_tangent_type>::type> {
  static bool get(const stce_9782& stce_9839) {
    return is_zero(value(stce_9839)) && is_zero(derivative(stce_9839));
  }
};
template <typename stce_9782>
  struct is_zero_trait<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_adjoint_type>::type> {
  static bool get(const stce_9782& stce_9839) {
    return is_zero(value(stce_9839)) && !tape_index(stce_9839);
  }
};
}
}
namespace dco {
namespace folding {
template<class stce_9918>
  class stce_9919 {
  typedef typename stce_9918::stce_9764 stce_9764;
  typedef typename stce_9918::stce_9763 stce_9763;
  typedef typename stce_9918::stce_9920 stce_9921;
public:
  stce_9919(stce_9918* tape) : stce_9922(1), stce_9923(tape->stce_9924()) {}
  template<class stce_9925>
    DCO_ALWAYS_INLINE stce_9925& stce_9926(stce_9763 const& pval) {
    stce_9922 *= pval;
    return static_cast<stce_9925&>(*this);
  }
  DCO_ALWAYS_INLINE stce_9763& pval() { return stce_9922; }
  DCO_ALWAYS_INLINE DCO_TAPE_INT index() const { return stce_9923.index(); }
  template <class stce_9927>
    DCO_ALWAYS_INLINE void insert(stce_9927 const& stce_9839) {
    const DCO_TAPE_INT stce_9928 = stce_9839._tape_index();
    if (!stce_9755 || stce_9928 != 0) {
      if (!stce_9756 || !is_zero(stce_9922)) {
        stce_9923.insert(stce_9928, stce_9922);
      }
    }
  }
protected:
  stce_9763 stce_9922;
private:
  stce_9919();
  stce_9919(stce_9919& stce_9929);
  stce_9921 stce_9923;
};
template<class stce_9918>
  class stce_9930 : public stce_9919<stce_9918> {
  typedef typename stce_9918::stce_9764 stce_9764;
  typedef typename stce_9918::stce_9920 stce_9921;
public:
  stce_9930(stce_9918* tape) : stce_9919<stce_9918>(tape),
    stce_9931(0) {
  }
  template <class stce_9927>
    DCO_ALWAYS_INLINE void insert(stce_9927 stce_9839) {
    stce_9919<stce_9918>::insert(stce_9839);
    stce_9931 += this->stce_9922 * stce_9839.stce_9932;
  }
  stce_9764 stce_9933() const {
    return stce_9931;
  }
private:
  stce_9930();
  stce_9930(stce_9930& stce_9929);
  stce_9764 stce_9931;
};
template<class stce_9934, class stce_9759, class stce_9935>
  DCO_ALWAYS_INLINE void interpret(const dco::internal::active_type<stce_9759, stce_9935> &stce_9839,
                                   stce_9934& stce_9936) {
  stce_9936.insert(stce_9839);
}
template<class stce_9934, class stce_9759, class stce_9937, class stce_9938, class stce_9939>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9779<stce_9759, stce_9937, stce_9938, stce_9939> &stce_9839,
                                   stce_9934& stce_9936) {
  stce_9759 stce_9940 = stce_9936.pval();
  interpret<stce_9934>(stce_9839.stce_9941, stce_9936.template stce_9926<stce_9934>(stce_9839.stce_9942()));
  stce_9936.pval() = stce_9940;
  interpret<stce_9934>(stce_9839.stce_9943, stce_9936.template stce_9926<stce_9934>(stce_9839.stce_9944()));
  stce_9936.pval() = stce_9940;
}
template<class stce_9934, class stce_9759, class stce_9945, class stce_9939>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9776<stce_9759, stce_9945, stce_9939> &stce_9839,
                                   stce_9934& stce_9936) {
  stce_9759 stce_9940 = stce_9936.pval();
  interpret<stce_9934>(stce_9839.stce_9946, stce_9936.template stce_9926<stce_9934>(stce_9839.pval()));
  stce_9936.pval() = stce_9940;
}
template<class stce_9934, class stce_9759, class stce_9937, class stce_9939>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9780<stce_9759, stce_9937, stce_9939> &stce_9839,
                                   stce_9934& stce_9936) {
  stce_9759 stce_9940 = stce_9936.pval();
  interpret<stce_9934>(stce_9839.stce_9941, stce_9936.template stce_9926<stce_9934>(stce_9839.stce_9942()));
  stce_9936.pval() = stce_9940;
}
template<class stce_9934, class stce_9759, class stce_9938, class stce_9939>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9781<stce_9759, stce_9938, stce_9939> &stce_9839,
                                   stce_9934& stce_9936) {
  stce_9759 stce_9940 = stce_9936.pval();
  interpret<stce_9934>(stce_9839.stce_9943, stce_9936.template stce_9926<stce_9934>(stce_9839.stce_9944()));
  stce_9936.pval() = stce_9940;
}
template<class stce_9918>
  class stce_9947 {
public:
  stce_9947() :
    stce_9948(0) {
  }
  template <class stce_9927>
    DCO_ALWAYS_INLINE void insert(stce_9927 stce_9839) {
    stce_9948 = stce_9839.tape();
  }
  stce_9918* tape() const {
    return stce_9948;
  }
private:
  stce_9947(stce_9947& stce_9929);
  stce_9918* stce_9948;
};
template<class stce_9759, class stce_9937, class stce_9938, class stce_9939, class stce_9918>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9779<stce_9759, stce_9937, stce_9938, stce_9939> &stce_9839,
                                   stce_9947<stce_9918>& stce_9936) {
  interpret(stce_9839.stce_9941, stce_9936);
  if(stce_9936.tape())
    {
      if(IS_DEBUG)
        {
          stce_9918 * tape = stce_9936.tape();
          interpret(stce_9839.stce_9943, stce_9936);
          if (stce_9936.tape() && tape != stce_9936.tape())
            throw dco::exception::create<std::runtime_error>(
                                                             "multiple tape error: different tapes on right-hand side of statement",
                                                             "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_folding.hpp", 159);
        }
      return;
    }
  interpret(stce_9839.stce_9943, stce_9936);
}
template<class stce_9759, class stce_9945, class stce_9939, class stce_9918>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9776<stce_9759, stce_9945, stce_9939> &stce_9839,
                                   stce_9947<stce_9918>& stce_9936) {
  interpret(stce_9839.stce_9946, stce_9936);
}
template<class stce_9759, class stce_9937, class stce_9939, class stce_9918>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9780<stce_9759, stce_9937, stce_9939> &stce_9839,
                                   stce_9947<stce_9918>& stce_9936) {
  interpret(stce_9839.stce_9941, stce_9936);
}
template<class stce_9759, class stce_9938, class stce_9939, class stce_9918>
  DCO_ALWAYS_INLINE void interpret(
                                   const dco::internal::stce_9781<stce_9759, stce_9938, stce_9939> &stce_9839,
                                   stce_9947<stce_9918>& stce_9936) {
  interpret(stce_9839.stce_9943, stce_9936);
}
}
}
namespace dco {
namespace internal {
struct stce_9949 {
  template<class stce_9950, class stce_9951> static bool stce_9952(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953._value() == stce_9954._value(); } template<class stce_9950, class stce_9951> static bool stce_9955(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953._value() == stce_9954; } template<class stce_9950, class stce_9951> static bool stce_9956(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953 == stce_9954._value(); }
  template<class stce_9950, class stce_9951> static bool stce_9957(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953._value() != stce_9954._value(); } template<class stce_9950, class stce_9951> static bool stce_9958(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953._value() != stce_9954; } template<class stce_9950, class stce_9951> static bool stce_9959(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953 != stce_9954._value(); }
  template<class stce_9950, class stce_9951> static bool stce_9960(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953._value() < stce_9954._value(); } template<class stce_9950, class stce_9951> static bool stce_9961(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953._value() < stce_9954; } template<class stce_9950, class stce_9951> static bool stce_9962(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953 < stce_9954._value(); }
  template<class stce_9950, class stce_9951> static bool stce_9963(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953._value() <= stce_9954._value(); } template<class stce_9950, class stce_9951> static bool stce_9964(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953._value() <= stce_9954; } template<class stce_9950, class stce_9951> static bool stce_9965(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953 <= stce_9954._value(); }
  template<class stce_9950, class stce_9951> static bool stce_9966(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953._value() > stce_9954._value(); } template<class stce_9950, class stce_9951> static bool stce_9967(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953._value() > stce_9954; } template<class stce_9950, class stce_9951> static bool stce_9968(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953 > stce_9954._value(); }
  template<class stce_9950, class stce_9951> static bool stce_9969(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953._value() >= stce_9954._value(); } template<class stce_9950, class stce_9951> static bool stce_9970(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953._value() >= stce_9954; } template<class stce_9950, class stce_9951> static bool stce_9971(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { return stce_9953 >= stce_9954._value(); }
};
template<class stce_9972, stce_9972 *&stce_9948>
struct stce_9973 {
  template<class stce_9950, class stce_9951> static bool stce_9952(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953._value() == stce_9954._value(); const bool stce_9974 = stce_9953.stce_9975(0) == stce_9954.stce_9975(0); const bool stce_9976 = stce_9953.stce_9975(1) == stce_9954.stce_9975(1); stce_9948->stce_9977("compare_AA", stce_9836, stce_9974,stce_9976); return stce_9836; } template<class stce_9950, class stce_9951> static bool stce_9955(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953._value() == stce_9954; const bool stce_9974 = stce_9953.stce_9975(0) == stce_9954; const bool stce_9976 = stce_9953.stce_9975(1) == stce_9954; stce_9948->stce_9977("compare_AP", stce_9836, stce_9974, stce_9976); return stce_9836; } template<class stce_9950, class stce_9951> static bool stce_9956(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953 == stce_9954._value(); const bool stce_9974 = stce_9953 == stce_9954.stce_9975(0); const bool stce_9976 = stce_9953 == stce_9954.stce_9975(1); stce_9948->stce_9977("compare_PA", stce_9836, stce_9974, stce_9976); return stce_9836; }
  template<class stce_9950, class stce_9951> static bool stce_9957(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953._value() != stce_9954._value(); const bool stce_9974 = stce_9953.stce_9975(0) != stce_9954.stce_9975(0); const bool stce_9976 = stce_9953.stce_9975(1) != stce_9954.stce_9975(1); stce_9948->stce_9977("compare_AA", stce_9836, stce_9974,stce_9976); return stce_9836; } template<class stce_9950, class stce_9951> static bool stce_9958(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953._value() != stce_9954; const bool stce_9974 = stce_9953.stce_9975(0) != stce_9954; const bool stce_9976 = stce_9953.stce_9975(1) != stce_9954; stce_9948->stce_9977("compare_AP", stce_9836, stce_9974, stce_9976); return stce_9836; } template<class stce_9950, class stce_9951> static bool stce_9959(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953 != stce_9954._value(); const bool stce_9974 = stce_9953 != stce_9954.stce_9975(0); const bool stce_9976 = stce_9953 != stce_9954.stce_9975(1); stce_9948->stce_9977("compare_PA", stce_9836, stce_9974, stce_9976); return stce_9836; }
  template<class stce_9950, class stce_9951> static bool stce_9960(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953._value() < stce_9954._value(); const bool stce_9974 = stce_9953.stce_9975(0) < stce_9954.stce_9975(0); const bool stce_9976 = stce_9953.stce_9975(1) < stce_9954.stce_9975(1); stce_9948->stce_9977("compare_AA", stce_9836, stce_9974,stce_9976); return stce_9836; } template<class stce_9950, class stce_9951> static bool stce_9961(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953._value() < stce_9954; const bool stce_9974 = stce_9953.stce_9975(0) < stce_9954; const bool stce_9976 = stce_9953.stce_9975(1) < stce_9954; stce_9948->stce_9977("compare_AP", stce_9836, stce_9974, stce_9976); return stce_9836; } template<class stce_9950, class stce_9951> static bool stce_9962(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953 < stce_9954._value(); const bool stce_9974 = stce_9953 < stce_9954.stce_9975(0); const bool stce_9976 = stce_9953 < stce_9954.stce_9975(1); stce_9948->stce_9977("compare_PA", stce_9836, stce_9974, stce_9976); return stce_9836; }
  template<class stce_9950, class stce_9951> static bool stce_9963(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953._value() <= stce_9954._value(); const bool stce_9974 = stce_9953.stce_9975(0) <= stce_9954.stce_9975(0); const bool stce_9976 = stce_9953.stce_9975(1) <= stce_9954.stce_9975(1); stce_9948->stce_9977("compare_AA", stce_9836, stce_9974,stce_9976); return stce_9836; } template<class stce_9950, class stce_9951> static bool stce_9964(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953._value() <= stce_9954; const bool stce_9974 = stce_9953.stce_9975(0) <= stce_9954; const bool stce_9976 = stce_9953.stce_9975(1) <= stce_9954; stce_9948->stce_9977("compare_AP", stce_9836, stce_9974, stce_9976); return stce_9836; } template<class stce_9950, class stce_9951> static bool stce_9965(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953 <= stce_9954._value(); const bool stce_9974 = stce_9953 <= stce_9954.stce_9975(0); const bool stce_9976 = stce_9953 <= stce_9954.stce_9975(1); stce_9948->stce_9977("compare_PA", stce_9836, stce_9974, stce_9976); return stce_9836; }
  template<class stce_9950, class stce_9951> static bool stce_9966(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953._value() > stce_9954._value(); const bool stce_9974 = stce_9953.stce_9975(0) > stce_9954.stce_9975(0); const bool stce_9976 = stce_9953.stce_9975(1) > stce_9954.stce_9975(1); stce_9948->stce_9977("compare_AA", stce_9836, stce_9974,stce_9976); return stce_9836; } template<class stce_9950, class stce_9951> static bool stce_9967(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953._value() > stce_9954; const bool stce_9974 = stce_9953.stce_9975(0) > stce_9954; const bool stce_9976 = stce_9953.stce_9975(1) > stce_9954; stce_9948->stce_9977("compare_AP", stce_9836, stce_9974, stce_9976); return stce_9836; } template<class stce_9950, class stce_9951> static bool stce_9968(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953 > stce_9954._value(); const bool stce_9974 = stce_9953 > stce_9954.stce_9975(0); const bool stce_9976 = stce_9953 > stce_9954.stce_9975(1); stce_9948->stce_9977("compare_PA", stce_9836, stce_9974, stce_9976); return stce_9836; }
  template<class stce_9950, class stce_9951> static bool stce_9969(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953._value() >= stce_9954._value(); const bool stce_9974 = stce_9953.stce_9975(0) >= stce_9954.stce_9975(0); const bool stce_9976 = stce_9953.stce_9975(1) >= stce_9954.stce_9975(1); stce_9948->stce_9977("compare_AA", stce_9836, stce_9974,stce_9976); return stce_9836; } template<class stce_9950, class stce_9951> static bool stce_9970(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953._value() >= stce_9954; const bool stce_9974 = stce_9953.stce_9975(0) >= stce_9954; const bool stce_9976 = stce_9953.stce_9975(1) >= stce_9954; stce_9948->stce_9977("compare_AP", stce_9836, stce_9974, stce_9976); return stce_9836; } template<class stce_9950, class stce_9951> static bool stce_9971(const stce_9950 &stce_9953, const stce_9951 &stce_9954) { const bool stce_9836 = stce_9953 >= stce_9954._value(); const bool stce_9974 = stce_9953 >= stce_9954.stce_9975(0); const bool stce_9976 = stce_9953 >= stce_9954.stce_9975(1); stce_9948->stce_9977("compare_PA", stce_9836, stce_9974, stce_9976); return stce_9836; }
};
}
}
namespace dco {
namespace internal {
template<class stce_9978>
class local_gradient_t {
private:
  typedef typename dco::mode<stce_9978>::tape_t tape_t;
  typedef typename tape_t::stce_9763 stce_9763;
  typedef typename tape_t::stce_9920 stce_9920;
  stce_9978& stce_9979;
  stce_9920 stce_9923;
  local_gradient_t(const local_gradient_t&);
public:
  local_gradient_t(stce_9978 &stce_9980)
    : stce_9979(stce_9980)
    , stce_9923(stce_9980.tape()->stce_9924())
  {}
  ~local_gradient_t() {
    stce_9979._tape_index() = stce_9923.index();
  }
  template<typename stce_9794>
    void put(const stce_9794 &stce_9839, const stce_9763& stce_9981) {
    if (!stce_9755 || stce_9839._tape_index() != 0) {
      bool stce_9982 = (stce_9923.index() == 0);
      stce_9923.insert(stce_9839._tape_index(), stce_9981);
      if (stce_9982)
        stce_9979.stce_9983(stce_9979._value(), stce_9814());
    }
  }
};
template<class stce_9887>
struct stce_9984 {
  typedef typename stce_9887::tape_t::stce_9764 derivative_t;
protected:
  mutable DCO_TAPE_INT stce_9885;
public:
  stce_9984() : stce_9885(0){}
  DCO_STRONG_INLINE stce_9984 &operator =(const stce_9984 &stce_9839) {
    stce_9885 = stce_9839.stce_9885;
    return *this;
  }
  void clear() {
    stce_9885 = 0;
  }
  derivative_t& stce_9985() const {
    return _adjoint();
  }
  derivative_t& _adjoint() const {
    if (IS_DEBUG && !stce_9887::global_tape)
      throw dco::exception::create<std::runtime_error>("Tape not created yet. Please use tape_t::create() before accessing derivative components of adjoint types.");
    return stce_9887::global_tape->_adjoint(stce_9885);
  }
  void stce_9986(const derivative_t &stce_9843) const {
    if (IS_DEBUG && !stce_9887::global_tape)
      throw dco::exception::create<std::runtime_error>("Tape not created yet. Please use tape_t::create() before accessing derivative components of adjoint types.");
    stce_9887::global_tape->_adjoint(stce_9885) = stce_9843;
  }
  bool _is_registered() const {
    return stce_9885 == 0 ? false : true;
  }
  DCO_TAPE_INT & _tape_index() { return stce_9885; }
  DCO_TAPE_INT const& _tape_index() const { return stce_9885; }
  template<class stce_9987>
  DCO_STRONG_INLINE void stce_9983(const stce_9987 &, const double &) { }
  template <typename stce_9927>
    static DCO_STRONG_INLINE void register_variable(stce_9927& stce_9839, DCO_TAPE_INT index, typename stce_9887::tape_t*) {
    if (IS_DEBUG && !stce_9887::global_tape) {
      throw dco::exception::create<std::runtime_error>("Global tape not yet created (is NULL).");
    }
    stce_9839.stce_9885 = index;
  }
  static typename stce_9887::tape_t *tape() {
    return stce_9887::global_tape;
  }
  void set_tape(typename stce_9887::tape_t *) {}
};
template<typename stce_9887>
  class single_tape_data : public stce_9984<stce_9887> {
  typedef stce_9984<stce_9887> stce_9988;
public:
  typedef stce_9887 mode_t;
  typedef stce_9949 stce_9989;
  typedef typename stce_9887::tape_t tape_t;
  typedef typename tape_t::stce_9764 derivative_t;
  single_tape_data() {}
  DCO_STRONG_INLINE single_tape_data &operator =(const single_tape_data &stce_9839) {
    stce_9988::operator= (static_cast<const stce_9988&>(stce_9839));
    return *this;
  }
  template<class stce_9888>
  DCO_ALWAYS_INLINE void stce_9889(const stce_9888 &stce_9990) {
    if (!stce_9887::global_tape || !stce_9887::global_tape->is_active()) {
      this->clear();
      return;
    }
    using namespace folding;
    stce_9919<tape_t> stce_9991(stce_9887::global_tape);
    interpret(stce_9990, stce_9991);
    this->stce_9885 = stce_9991.index();
  }
};
template<typename stce_9887>
  class stce_9992 {
  typedef typename stce_9887::tape_t tape_t;
  mutable DCO_TAPE_INT stce_9885;
  mutable tape_t *stce_9993;
public:
  typedef stce_9887 mode_t;
  typedef stce_9949 stce_9989;
  typedef typename tape_t::stce_9764 derivative_t;
  stce_9992() : stce_9885(0), stce_9993(0) {}
  DCO_STRONG_INLINE stce_9992 &operator =(const stce_9992 &stce_9839) {
   stce_9885 = stce_9839.stce_9885;
   stce_9993 = stce_9839.stce_9993;
    return *this;
  }
  void clear() {
    stce_9885 = 0;
    stce_9993 = 0;
  }
  derivative_t& stce_9985() const {
    return _adjoint();
  }
  derivative_t &_adjoint() const {
    static derivative_t stce_9994 = 0;
    if (stce_9993 == 0) {
      stce_9994 = 0;
      return stce_9994;
    } else {
      return stce_9993->_adjoint(stce_9885);
    }
  }
  void stce_9986(const derivative_t &stce_9843) const {
    if (stce_9993 != 0)
      stce_9993->_adjoint(stce_9885) = stce_9843;
  }
  bool _is_registered() const {
    return stce_9885 == 0 ? false : true;
  }
  DCO_TAPE_INT &_tape_index() {
    return stce_9885;
  }
  const DCO_TAPE_INT &_tape_index() const {
    return stce_9885;
  }
  template<typename stce_9987>
  DCO_STRONG_INLINE void stce_9983(const stce_9987 &, const double &) { }
  template <typename stce_9927>
    static DCO_STRONG_INLINE void register_variable(stce_9927& stce_9839, DCO_TAPE_INT index, tape_t* tape) {
    stce_9839.stce_9993 = tape;
    stce_9839.stce_9885 = index;
  }
  tape_t *tape() const {
    return stce_9993;
  }
  void set_tape(tape_t *stce_9782) {
    stce_9993 = stce_9782;
  }
  template<class stce_9888>
  DCO_ALWAYS_INLINE void stce_9889(const stce_9888 &stce_9990) {
    using namespace folding;
    stce_9947<tape_t> stce_9947;
    interpret(stce_9990, stce_9947);
    tape_t *stce_9995 = stce_9947.tape();
    if (!stce_9995 || !stce_9995->is_active()) {
       clear();
       return;
     }
    stce_9919<tape_t> stce_9991(stce_9995);
    interpret(stce_9990, stce_9991);
    stce_9885 = stce_9991.index();
    stce_9993 = stce_9995;
  }
};
template<typename stce_9887, size_t stce_9996=2>
class stce_9997: public stce_9984<stce_9887> {
  typedef stce_9984<stce_9887> stce_9988;
public:
  typedef stce_9887 mode_t;
  typedef stce_9973<typename stce_9887::tape_t, stce_9887::global_tape> stce_9989;
  typedef typename stce_9887::tape_t tape_t;
  typedef typename tape_t::stce_9764 derivative_t;
  mutable derivative_t stce_9932;
  mutable derivative_t stce_9998[stce_9996];
  void stce_9999(const derivative_t &value, const double &stce_9933) const {
    if(stce_9887::global_tape->stce_10000 == 0)
      this->stce_9998[0] = this->stce_9998[1] = value;
    stce_9887::global_tape->stce_10001(this, stce_9933);
    stce_9932 = stce_9933;
  }
public:
  stce_9997(): stce_9932(0) {}
  DCO_STRONG_INLINE stce_9997 &operator =(const stce_9997 &stce_9839) {
    stce_9988::operator= (static_cast<const stce_9988&>(stce_9839));
    stce_9932 = stce_9839.stce_9932;
    return *this;
  }
  inline void clear() {
    stce_9988::clear();
    stce_9932 = 0;
  }
  template<class stce_9987>
  DCO_STRONG_INLINE void stce_9983(const stce_9987 &value, const double &stce_9933) {
    stce_9999(value, stce_9933);
  }
  template <typename stce_9927>
    static DCO_STRONG_INLINE void register_variable(stce_9927& stce_9839, DCO_TAPE_INT index, tape_t* tape) {
    stce_9988::register_variable(stce_9839, index, tape);
    stce_9839.stce_9999(trait_value<stce_9927>::value(stce_9839), stce_9814());
  }
  const derivative_t& stce_9975(const derivative_t &value, const int stce_9806) const {
      if(this->stce_9885==0) return value;
      return stce_9998[stce_9806];
  }
  template<class stce_9888>
    DCO_ALWAYS_INLINE void stce_9889(const stce_9888 &stce_9990) {
    if (!stce_9887::global_tape || !stce_9887::global_tape->is_active()) {
      clear();
      return;
    }
    {
      using namespace folding;
      stce_9930<tape_t> stce_9991(stce_9887::global_tape);
      interpret(stce_9990, stce_9991);
      this->stce_9885 = stce_9991.index();
      if (this->stce_9885 == 0)
        return;
      stce_9932 = stce_9991.stce_9933();
    }
    stce_9998[0] = stce_9990.stce_9975(0);
    stce_9998[1] = stce_9990.stce_9975(1);
    derivative_t stce_10002 = stce_9998[1] - stce_9998[0];
    double stce_10003 = stce_10002/(2*stce_9887::global_tape->stce_10000);
    stce_9887::global_tape->stce_10004(stce_9932, stce_10003);
    stce_9887::global_tape->stce_10005();
    if(stce_9887::global_tape->stce_10006) {
      const double stce_10007=3.1415692;
      stce_9887::global_tape->_adjoint(this->stce_9885) = stce_10007;
      stce_9887::global_tape->interpret_adjoint();
      double stce_10008=0;
      for(size_t stce_9827=0;stce_9827<stce_9887::global_tape->stce_10009.size();++stce_9827)
        stce_10008 += stce_9887::global_tape->_adjoint(stce_9887::global_tape->stce_10009[stce_9827]) * stce_9887::global_tape->stce_10010[stce_9827];
      stce_9887::global_tape->zero_adjoints();
      stce_9887::global_tape->stce_10011(stce_10008, stce_9932 * stce_10007, stce_10003 * stce_10007);
    }
    stce_9887::global_tape->stce_10005();
  }
};
template <typename stce_10012>
  struct stce_10013 : public dco::mode<typename stce_10012::binary_t>::data_t {
  typedef stce_10012 mode_t;
  static const int gbcp_order = mode_t::gbcp_order;
  static const int order = mode_t::order;
  template<class stce_9888>
    DCO_ALWAYS_INLINE void stce_9889(stce_9888 const&) {
    }
  void clear() {}
  stce_10013& operator=(stce_10013 const&) {
    return *this;
  }
};
}
}
namespace dco
{
  template <class stce_9782>
  class ga1s_codegen {
    typedef stce_9762<stce_9782> stce_10014;
public:
    static const memory_model::TYPE TAPE_MODEL = memory_model::stce_9767;
    typedef dco::internal::stce_9884<ga1s_codegen> data;
    typedef typename stce_10014::value_t value_t;
    typedef typename stce_10014::value_t active_value_t;
    typedef typename stce_10014::stce_9764 derivative_t;
    typedef typename stce_10014::stce_9763 stce_9763;
    typedef typename mode<value_t>::passive_t passive_t;
    typedef dco::internal::active_type<derivative_t, data> type;
    typedef type active_t;
    static const bool is_dco_type = true;
    static const bool is_adjoint_type = true;
    static const bool is_tangent_type = false;
    static const int order = dco::mode<stce_9782>::order + 1;
    typedef typename internal::stce_9871<derivative_t, stce_9763> tape_t;
    static DCOi_THREAD_LOCAL tape_t *global_tape;
    typedef typename internal::local_gradient_t<type> local_gradient_t;
  };
template <class stce_9782>
  DCOi_THREAD_LOCAL typename ga1s_codegen<stce_9782>::tape_t *ga1s_codegen<stce_9782>::global_tape;
}
       
       
       
namespace dco {
  template <class stce_9759>
  class stce_10015 {
  public:
    typedef stce_9759 active_value_t;
    typedef stce_9759 type;
    typedef type active_t;
    typedef stce_9759 value_t;
    typedef stce_9759 passive_t;
    typedef void derivative_t;
    typedef void stce_10016;
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
  typename stce_10015<stce_9759>::tape_t *stce_10015<stce_9759>::global_tape = 0;
}
namespace dco {
template <class stce_9795>
  struct mode : public dco::stce_10015<stce_9795> {
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
template <typename stce_9782, typename stce_10017, typename stce_10018>
  struct mode<internal::stce_9776<stce_9782, stce_10017, stce_10018> > :
    public internal::stce_9776<stce_9782, stce_10017, stce_10018>::data_t::mode_t {
  static const bool is_intermediate_type = true;
};
template<class stce_9782, class stce_10019, class stce_10020, class stce_10018>
  struct mode<internal::stce_9779<stce_9782, stce_10019, stce_10020, stce_10018> > :
    public internal::stce_9779<stce_9782, stce_10019, stce_10020, stce_10018>::data_t::mode_t {
  static const bool is_intermediate_type = true;
};
template<class stce_9782, class stce_10017, class stce_10018>
  struct mode<internal::stce_9780<stce_9782, stce_10017, stce_10018> > :
    public internal::stce_9780<stce_9782, stce_10017, stce_10018>::data_t::mode_t {
  static const bool is_intermediate_type = true;
};
template<class stce_9782, class stce_10017, class stce_10018>
  struct mode<internal::stce_9781<stce_9782, stce_10017, stce_10018> > :
    public internal::stce_9781<stce_9782, stce_10017, stce_10018>::data_t::mode_t {
  static const bool is_intermediate_type = true;
};
template<class stce_9771, class DATA_TYPE>
  struct mode<internal::stce_9773<stce_9771, DATA_TYPE> > :
    public DATA_TYPE::mode_t {
  static const bool is_intermediate_type = true;
};
}
       
namespace dco {
template <typename stce_9782, typename stce_9819>
  stce_9782& create_gbcp(stce_9819& stce_10021, typename enable_if<!stce_9865<stce_9782>::value, stce_9782>::type * = 0) {
      return reinterpret_cast<stce_9782&>(stce_10021);
  }
template <typename stce_9782, typename stce_9819>
  stce_9782 create_gbcp(stce_9819& stce_10021, typename enable_if<stce_9865<stce_9782>::value, stce_9782>::type * = 0) {
      return reinterpret_cast<stce_9782>(&stce_10021);
  }
namespace internal {
template<typename stce_9782>
  struct stce_10022 {
  static const bool value = false;
};
template <typename mode_t>
  struct stce_10022<stce_10013<mode_t> > {
  static const bool value = true;
};
template <typename stce_10023, typename data_t>
  struct stce_10022<active_type<stce_10023, data_t> > {
  static const bool value = stce_10022<data_t>::value;
};
template<typename stce_9782>
  struct is_gbcp_type {
  static const bool value = stce_10022<typename dco::remove_all<stce_9782>::type>::value;
};
template<typename stce_9782, int stce_9839> struct stce_10024;
template<typename stce_9782> struct stce_10024<stce_9782, 1> {
    typedef typename mode<stce_9782>::active_value_t& stce_9802;
    typedef const typename mode<stce_9782>::active_value_t& stce_10025;
    static stce_9802 value(stce_9782& at) { return at._value(); }
    static stce_10025 value(const stce_9782& at) { return at._value(); }
};
template<typename stce_9782, int stce_9839> struct stce_10024 {
    typedef typename stce_10024<typename mode<stce_9782>::active_value_t, stce_9839-1>::stce_9802 stce_9802;
    typedef typename stce_10024<typename mode<stce_9782>::active_value_t, stce_9839-1>::stce_10025 stce_10025;
    static stce_9802 value(stce_9782& at) { return stce_10024<typename mode<stce_9782>::active_value_t, stce_9839-1>::value(at._value()); }
    static stce_10025 value(const stce_9782& at) { return stce_10024<typename mode<stce_9782>::active_value_t, stce_9839-1>::value(at._value()); }
};
}
}
namespace dco {
template <typename stce_9782,typename stce_10026=stce_9782>
  struct dco_type_constructable_from {};
template <typename stce_9782>
  struct dco_explicit_type_cast_to {
  const static bool value = false;
};
}
namespace dco {
namespace internal {
template<class stce_9771>
struct stce_10027 {
  const stce_9771 stce_10028;
  stce_10027(stce_9771 const& value) :
      stce_10028(value) {
  }
  inline const stce_9771 &_value() const { return stce_10028; }
};
template<class stce_9771, class DATA_TYPE>
struct stce_9773 : public stce_10027<stce_9771> {
  stce_9773(const stce_9771& value) :
    stce_10027<stce_9771>(value) {
  }
};
template<class stce_9771>
  struct stce_9773<stce_9771, typename ga1s_codegen<stce_9771>::data>
: public stce_10027<stce_9771>, public stce_9896<stce_9771> {
  stce_9773(const stce_9771& value) :
    stce_10027<stce_9771>(value) {
  }
};
template<class stce_9771, class stce_9774, class stce_9775>
struct stce_9776: public stce_9773<stce_9771, typename stce_9774::data_t>
  {
      const stce_9774 DCO_TEMPORARY_REFORCOPY stce_9946;
      typedef typename stce_9774::data_t data_t;
      typedef typename data_t::mode_t::value_t VALUE_TYPE;
      typedef typename data_t::mode_t::passive_t stce_10029;
    stce_9776(const stce_9774 &arg)
      :stce_9773<stce_9771, data_t>(stce_9775::stce_10030(arg._value())),
       stce_9946(arg) {
      }
    inline const stce_9771 stce_9975(const int stce_9806) const {
        return stce_9775::stce_10030(stce_9946.stce_9975(stce_9806));
    }
    inline const stce_9771 pval() const {
      return stce_9775::stce_10031(this->_value(), stce_9946._value());
    }
    void stce_9891(std::ofstream &out) const {
      stce_9946.stce_9891(out);
      this->stce_9899(out, stce_9775::stce_10032(this->str(stce_9946)));
    }
    void stce_9893(std::ofstream &out) const {
      std::string stce_9900 = stce_9775::stce_10033();
      this->stce_9903(stce_9900, stce_9946, "x");
      out << this->str(stce_9946) << "_a1s +=(" << stce_9900 << ")*" << this->stce_9886() << "_a1s;\t";
      stce_9946.stce_9893(out);
    }
};
    template<class stce_9771, class stce_9777, class stce_9778, class stce_9775>
    struct stce_9779: public stce_9773<stce_9771, typename stce_9777::data_t> {
      const stce_9777 DCO_TEMPORARY_REFORCOPY stce_9941;
      const stce_9778 DCO_TEMPORARY_REFORCOPY stce_9943;
      typedef typename stce_9777::data_t data_t;
      typedef typename data_t::mode_t::value_t VALUE_TYPE;
      typedef typename data_t::mode_t::passive_t stce_10029;
      stce_9779(const stce_9777 &stce_9908, const stce_9778 &stce_9909)
        :stce_9773<stce_9771, data_t>(stce_9775::stce_10030(stce_9908._value(), stce_9909._value())),
        stce_9941(stce_9908),
        stce_9943(stce_9909) {}
      inline const stce_9771 stce_9942() const {
        return stce_9775::stce_10034(this->_value(), stce_9941._value(), stce_9943._value());
      }
      inline const stce_9771 stce_9944() const {
        return stce_9775::stce_10035(this->_value(), stce_9941._value(), stce_9943._value());
      }
      inline const stce_9771 stce_9975(const int stce_9806) const {
        return stce_9775::stce_10030(stce_9941.stce_9975(stce_9806),stce_9943.stce_9975(stce_9806));
      }
      void stce_9891(std::ofstream &out) const {
        stce_9941.stce_9891(out);
        stce_9943.stce_9891(out);
        this->stce_9899(out, stce_9775::stce_10032(this->str(stce_9941),
                                                                                   this->str(stce_9943)));
      }
      void stce_9893(std::ofstream &out) const {
        this->template stce_9907<stce_9775>(out, stce_9941, stce_9943);
        this->template stce_9911<stce_9775>(out, stce_9941, stce_9943);
      }
    };
    template<class stce_9771, class stce_9777, class stce_9775>
    struct stce_9780 : public stce_9773<stce_9771, typename stce_9777::data_t> {
    const stce_9777 DCO_TEMPORARY_REFORCOPY stce_9941;
    const stce_9771 stce_9943;
    typedef typename stce_9777::data_t data_t;
    typedef typename data_t::mode_t::value_t VALUE_TYPE;
    typedef typename data_t::mode_t::passive_t stce_10029;
      stce_9780(const stce_9777 &stce_9908, const stce_9771 &stce_9909)
      :stce_9773<stce_9771, data_t>(stce_9775::stce_10030(stce_9908._value(), stce_9909)),
        stce_9941(stce_9908), stce_9943(stce_9909) {
      }
    inline const stce_9771 stce_9975(const int stce_9806) const {
      return stce_9775::stce_10030(stce_9941.stce_9975(stce_9806), stce_9943);
    }
    inline const stce_9771 stce_9942() const {
        return stce_9775::stce_10034(this->stce_10028, stce_9941._value(), stce_9943);
    }
    void stce_9891(std::ofstream &out) const {
      stce_9941.stce_9891(out);
      this->stce_9899(out, stce_9775::stce_10032(this->str(stce_9941),
          this->str(stce_9943)));
    }
    void stce_9893(std::ofstream &out) const {
      this->template stce_9907<stce_9775>(out, stce_9941, stce_9943);
      this->template stce_9911<stce_9775>(out, stce_9941, stce_9943);
    }
};
    template<class stce_9771, class stce_9778, class stce_9775>
      struct stce_9781: public stce_9773<stce_9771, typename stce_9778::data_t> {
      const stce_9771 stce_9941;
      const stce_9778 DCO_TEMPORARY_REFORCOPY stce_9943;
      typedef typename stce_9778::data_t data_t;
      typedef typename data_t::mode_t::value_t VALUE_TYPE;
      typedef typename data_t::mode_t::passive_t stce_10029;
      stce_9781(const stce_9771 &stce_9908, const stce_9778 &stce_9909)
        :stce_9773<stce_9771, data_t>(stce_9775::stce_10030(stce_9908, stce_9909._value())),
        stce_9941(stce_9908), stce_9943(stce_9909) {
      }
      inline const stce_9771 stce_9975(const int stce_9806) const {
        return stce_9775::stce_10030(stce_9941, stce_9943.stce_9975(stce_9806));
      }
      inline const stce_9771 stce_9944() const {
        return stce_9775::stce_10035(this->stce_10028, stce_9941, stce_9943._value());
      }
      void stce_9891(std::ofstream &out) const {
        stce_9943.stce_9891(out);
        this->stce_9899(out, stce_9775::stce_10032(this->str(stce_9941),
                                                                                   this->str(stce_9943)));
      }
      void stce_9893(std::ofstream &out) const {
        this->template stce_9907<stce_9775>(out, stce_9941, stce_9943);
        this->template stce_9911<stce_9775>(out, stce_9941, stce_9943);
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
  class stce_10036 {
protected:
  stce_9782 stce_10028;
public:
  stce_10036() : stce_10028(0) {}
  template <typename stce_9819>
    stce_10036(const stce_9819& stce_9839) : stce_10028(stce_9839) {}
};
template<class stce_9771, class stce_9772>
struct active_type : private stce_10036<stce_9771>, public stce_9772 {
  typedef stce_10036<stce_9771> stce_10037;
  using stce_10037::stce_10028;
public:
  typedef stce_9772 data_t;
  typedef typename data_t::mode_t::value_t VALUE_TYPE;
  typedef typename data_t::mode_t::passive_t stce_10029;
  DCO_STRONG_INLINE const stce_9771& _value() const {
    return stce_10028;
  }
  stce_9772& data() {
    return *static_cast<stce_9772*>(this);
  }
  DCO_STRONG_INLINE stce_9771 &_value() {
    return stce_10028;
  }
  DCO_STRONG_INLINE const stce_9771 stce_9975(const int stce_9806) const {
        return stce_9772::stce_9975(stce_10028, stce_9806);
    }
  DCO_STRONG_INLINE active_type() :
    stce_10037(static_cast<stce_9771>(0.0)) {
  }
  template<typename TYPE>
    DCO_STRONG_INLINE active_type(const TYPE &stce_9753, typename dco_type_constructable_from<TYPE>::type* = 0) :
      stce_10037(stce_9753) { }
  template<typename TYPE>
    DCO_STRONG_INLINE typename dco::enable_if<dco_type_constructable_from<TYPE>::value, active_type>::type& operator=(const TYPE &stce_9753) {
    stce_10028 = stce_9753;
    return *this;
  }
  template<typename TYPE>
  DCO_STRONG_INLINE active_type(const TYPE& stce_9753, typename dco::enable_if<stce_9868<typename TYPE::data_t::mode_t, typename data_t::mode_t>::value, void*>::type = NULL) :
     stce_10037(stce_9753) { }
  DCO_STRONG_INLINE active_type(const stce_10029 &stce_9753) : stce_10037(stce_9753) { }
  DCO_STRONG_INLINE active_type &operator =(const active_type &stce_9839) {
    stce_9772::operator= (static_cast<const stce_9772&>(stce_9839));
    stce_10028 = stce_9839.stce_10028;
    return *this;
  }
  DCO_STRONG_INLINE active_type &operator =(const stce_10029 &stce_9990) {
    this->stce_10028 = stce_9990;
    this->clear();
    return *this;
  }
  template<class stce_10038, class stce_10039>
  DCO_STRONG_INLINE typename dco::enable_if<stce_9868<typename stce_10039::mode_t, typename data_t::mode_t>::value, active_type>::type& operator =(const active_type<stce_10038, stce_10039> &stce_9990) {
    stce_10028 = stce_9990;
    this->clear();
    return *this;
  }
  template<typename TYPE>
  DCO_STRONG_INLINE typename dco_type_constructable_from<TYPE, active_type>::type& operator=(const TYPE &stce_9753){
    stce_10028 = stce_9753;
    this->clear();
    return *this;
  }
    template<class stce_9937, class stce_9938, class stce_9939 > DCO_ALWAYS_INLINE active_type(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839) { this->stce_9889(stce_9839); this->stce_10028 = stce_9839.stce_10028; } template<class stce_9937, class stce_9938, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator=(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839) { this->stce_9889(stce_9839); this->stce_10028 = stce_9839.stce_10028; return *this; }
    template<class stce_9937, class stce_9939 > DCO_ALWAYS_INLINE active_type(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839) { this->stce_9889(stce_9839); this->stce_10028 = stce_9839.stce_10028; } template<class stce_9937, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator=(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839) { this->stce_9889(stce_9839); this->stce_10028 = stce_9839.stce_10028; return *this; }
    template<class stce_9938, class stce_9939 > DCO_ALWAYS_INLINE active_type(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839) { this->stce_9889(stce_9839); this->stce_10028 = stce_9839.stce_10028; } template<class stce_9938, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator=(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839) { this->stce_9889(stce_9839); this->stce_10028 = stce_9839.stce_10028; return *this; }
    template<class stce_9945, class stce_9939 > DCO_ALWAYS_INLINE active_type(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839) { this->stce_9889(stce_9839); this->stce_10028 = stce_9839.stce_10028; } template<class stce_9945, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator=(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839) { this->stce_9889(stce_9839); this->stce_10028 = stce_9839.stce_10028; return *this; }
    template<class stce_10039> DCO_ALWAYS_INLINE active_type& operator += (const active_type<stce_9771, stce_10039> &stce_9839){ *this = *this + stce_9839; return *this; } template<class stce_9937, class stce_9938, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator += (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839){ *this = *this + stce_9839; return *this; } template<class stce_9937, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator += (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839){ *this = *this + stce_9839; return *this; } template<class stce_9938, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator += (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839){ *this = *this + stce_9839; return *this; } template<class stce_9945, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator += (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839){ *this = *this + stce_9839; return *this; } DCO_ALWAYS_INLINE active_type& operator += (const VALUE_TYPE &stce_9839) { this->_value() += stce_9839; return *this; }
    template<class stce_10039> DCO_ALWAYS_INLINE active_type& operator -= (const active_type<stce_9771, stce_10039> &stce_9839){ *this = *this - stce_9839; return *this; } template<class stce_9937, class stce_9938, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator -= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839){ *this = *this - stce_9839; return *this; } template<class stce_9937, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator -= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839){ *this = *this - stce_9839; return *this; } template<class stce_9938, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator -= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839){ *this = *this - stce_9839; return *this; } template<class stce_9945, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator -= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839){ *this = *this - stce_9839; return *this; } DCO_ALWAYS_INLINE active_type& operator -= (const VALUE_TYPE &stce_9839) { this->_value() -= stce_9839; return *this; }
    template<class stce_10039> DCO_ALWAYS_INLINE active_type& operator *= (const active_type<stce_9771, stce_10039> &stce_9839){ *this = *this * stce_9839; return *this; } template<class stce_9937, class stce_9938, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator *= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839){ *this = *this * stce_9839; return *this; } template<class stce_9937, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator *= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839){ *this = *this * stce_9839; return *this; } template<class stce_9938, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator *= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839){ *this = *this * stce_9839; return *this; } template<class stce_9945, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator *= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839){ *this = *this * stce_9839; return *this; } DCO_STRONG_INLINE active_type& operator *= (const VALUE_TYPE &stce_9839) { *this = *this * stce_9839; return *this; }
    template<class stce_10039> DCO_ALWAYS_INLINE active_type& operator /= (const active_type<stce_9771, stce_10039> &stce_9839){ *this = *this / stce_9839; return *this; } template<class stce_9937, class stce_9938, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator /= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839){ *this = *this / stce_9839; return *this; } template<class stce_9937, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator /= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839){ *this = *this / stce_9839; return *this; } template<class stce_9938, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator /= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839){ *this = *this / stce_9839; return *this; } template<class stce_9945, class stce_9939 > DCO_ALWAYS_INLINE active_type& operator /= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839){ *this = *this / stce_9839; return *this; } DCO_STRONG_INLINE active_type& operator /= (const VALUE_TYPE &stce_9839) { *this = *this / stce_9839; return *this; }
  DCO_STRONG_INLINE active_type &operator ++() {
    ++this->stce_10028;
    return *this;
  }
  DCO_STRONG_INLINE active_type &operator --() {
    --this->stce_10028;
    return *this;
  }
  DCO_STRONG_INLINE active_type operator ++(int) {
    active_type stce_9836(*this);
    ++this->stce_10028;
    return stce_9836;
  }
  DCO_STRONG_INLINE active_type operator --(int) {
    active_type stce_9836(*this);
    --this->stce_10028;
    return stce_9836;
  }
  std::string stce_9886() const {
    if (stce_9772::_tape_index() == 0) {
      std::stringstream stce_10040;
      return stce_10040.str();
    } else
      return stce_9772::stce_9886();
  }
  void stce_9891(const std::ofstream &) const {
  }
  void stce_9893(const std::ofstream &) const {
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
  template<typename stce_9887>
    struct ts_data {
    typedef stce_9887 mode_t;
    typedef stce_9949 stce_9989;
    typedef typename mode_t::derivative_t derivative_t;
    mutable derivative_t stce_10041;
    ts_data() : stce_10041(0) {}
    ts_data &operator = (const ts_data &stce_10042) {
      stce_10041 = stce_10042.stce_10041;
      return *this;
    }
    void stce_9986(const derivative_t &stce_9753) {
        stce_10041 = stce_9753;
      }
      void stce_10043(derivative_t &stce_9753) const {
        stce_9753 = stce_10041;
      }
      void clear() {
        stce_10041 = 0;
      }
      derivative_t &stce_9985() const
      {
        stce_1000002();;
        return stce_10041;
      }
      template<class stce_9888>
      void stce_9889(const stce_9888 &stce_9990) {
        stce_10041 = stce_10044(stce_9990, 1.0);
      }
      template<class stce_9771, class stce_9935 >
      static derivative_t stce_10044(const dco::internal::active_type<stce_9771, stce_9935> &stce_9839, const derivative_t &pval) {
        return stce_9839.stce_10041 * pval;
      }
      template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 >
      static derivative_t stce_10044(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839, const derivative_t &pval) {
        return stce_10044(stce_9839.stce_9941, stce_9839.stce_9942() * pval) + stce_10044(stce_9839.stce_9943, stce_9839.stce_9944() * pval);
      }
      template<class stce_9771, class stce_9945, class stce_9939 >
      static derivative_t stce_10044(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839, const derivative_t &pval) {
        return stce_10044(stce_9839.stce_9946, stce_9839.pval() * pval);
      }
      template<class stce_9771, class stce_9937, class stce_9939 >
      static derivative_t stce_10044(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839, const derivative_t &pval) {
        return stce_10044(stce_9839.stce_9941, stce_9839.stce_9942() * pval);
      }
      template<class stce_9771, class stce_9938, class stce_9939 >
      static derivative_t stce_10044(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839, const derivative_t &pval) {
        return stce_10044(stce_9839.stce_9943, stce_9839.stce_9944() * pval);
      }
      typedef void tape_t;
      typedef DCO_TAPE_INT stce_10045;
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
  template<typename stce_9887>
    struct stce_10046 {
    typedef stce_9887 mode_t;
    typedef stce_9949 stce_9989;
    typedef typename mode_t::derivative_t derivative_t;
    typedef typename derivative_t::value_t stce_10016;
private:
    typedef typename mode_t::value_t stce_9771;
public:
      static const int stce_9833;
      mutable derivative_t tlms;
      stce_10046() {
        for (int stce_9827 = 0; stce_9827 < stce_9833; ++stce_9827) tlms[stce_9827] = 0;
      }
    derivative_t &stce_9985() const {
        return tlms;
      }
      void stce_9986(const stce_10016 &stce_9753, const int stce_10047) {
        stce_1000002();;
        tlms[stce_10047] = stce_9753;
      }
      void stce_10043(stce_10016 &stce_9753, const int stce_10047) const {
        stce_1000002();;
        stce_9753 = tlms[stce_10047];
      }
      void clear() {
        for (int stce_9827 = 0; stce_9827 < stce_9833; ++stce_9827) {
          tlms[stce_9827] = 0.0;
        }
      }
      template<class stce_9888, class stce_9794>
      static void stce_9889(const stce_9888 &stce_9990, stce_9794 &stce_10048) {
        stce_10046 &data = const_cast<stce_10046 &>(stce_10048);
        for (int stce_9827 = 0; stce_9827 < stce_9833; ++stce_9827) {
          data.tlms[stce_9827] = stce_10044(stce_9990, 1.0, stce_9827);
        }
      }
      template<class stce_9935 >
      static stce_10016 stce_10044(const dco::internal::active_type<stce_9771, stce_9935> &stce_9839, const stce_10016 &pval, const int stce_10047) {
        return stce_9839.tlms[stce_10047] * pval;
      }
      template<class stce_9937, class stce_9938, class stce_9939 >
      static stce_10016 stce_10044(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839, const stce_10016 &pval, const int stce_10047) {
        return stce_10044(stce_9839.stce_9941, stce_9839.stce_9942() * pval, stce_10047) + stce_10044(stce_9839.stce_9943, stce_9839.stce_9944() * pval, stce_10047);
      }
      template<class stce_9945, class stce_9939 >
      static stce_10016 stce_10044(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839, const stce_10016 &pval, const int stce_10047) {
        return stce_10044(stce_9839.stce_9946, stce_9839.pval() * pval, stce_10047);
      }
      template<class stce_9937, class stce_9939 >
      static stce_10016 stce_10044(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839, const stce_10016 &pval, const int stce_10047) {
        return stce_10044(stce_9839.stce_9941, stce_9839.stce_9942() * pval, stce_10047);
      }
      template<class stce_9938, class stce_9939 >
      static stce_10016 stce_10044(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839, const stce_10016 &pval, const int stce_10047) {
        return stce_10044(stce_9839.stce_9943, stce_9839.stce_9944() * pval, stce_10047);
      }
      typedef void stce_10049;
      typedef DCO_TAPE_INT stce_10045;
      void *stce_9948() const {
        return NULL;
      }
      DCO_TAPE_INT _tape_index() const {
        return 0;
      }
    };
  template<typename stce_9887>
    struct stce_10050 {
      typedef stce_9887 mode_t;
      typedef stce_9949 stce_9989;
      typedef typename mode_t::derivative_t derivative_t;
      typedef typename derivative_t::value_t stce_10016;
protected:
      typedef typename mode_t::value_t stce_9771;
public:
      static const int stce_9833 = derivative_t::stce_9823;
      mutable bool stce_10051;
      mutable derivative_t tlms;
      stce_10050() : stce_10051(false) {
        for (int stce_9827 = 0; stce_9827 < stce_9833; ++stce_9827) tlms[stce_9827] = 0;
      }
      void stce_9986(const stce_10016 &stce_9753, const int stce_10047) {
        if (!stce_10051) {
          stce_10051 = true;
        }
        tlms[stce_10047] = stce_9753;
      }
      derivative_t &stce_9985() const {
        stce_1000002();;
        stce_10051 = true;
        return tlms;
      }
      void stce_10043(stce_10016 &stce_9753, const int stce_10047) const {
        stce_1000002();;
        if (stce_10051) {
          stce_9753 = tlms[stce_10047];
        } else {
          stce_9753 = 0;
        }
      }
      void clear() {
        for (int stce_9827 = 0; stce_9827 < stce_9833; ++stce_9827) {
          tlms[stce_9827] = 0.0;
        }
      }
      template<class stce_9888>
      DCO_ALWAYS_INLINE void stce_9889(const stce_9888 &stce_9990) {
        bool stce_10052 = this->stce_10053(stce_9990);
        if (stce_10052) {
          for (int stce_9827 = 0; stce_9827 < stce_9833; ++stce_9827) {
            tlms[stce_9827] = stce_10044(stce_9990, 1.0, stce_9827);
          }
          stce_10051 = true;
        } else
          stce_10051 = false;
      }
      template<class stce_9935 >
      static stce_10016 stce_10044(const dco::internal::active_type<stce_9771, stce_9935> &stce_9839, const stce_10016 &pval, const int stce_10047) {
        if (stce_9839.stce_10051) {
          return stce_9839.tlms[stce_10047] * pval;
        } else
          return 0;
      }
      template<class stce_9937, class stce_9938, class stce_9939 >
      static stce_10016 stce_10044(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839, const stce_10016 &pval, const int stce_10047) {
        return stce_10044(stce_9839.stce_9941, stce_9839.stce_9942() * pval, stce_10047) + stce_10044(stce_9839.stce_9943, stce_9839.stce_9944() * pval, stce_10047);
      }
      template<class stce_9945, class stce_9939 >
      static stce_10016 stce_10044(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839, const stce_10016 &pval, const int stce_10047) {
        return stce_10044(stce_9839.stce_9946, stce_9839.pval() * pval, stce_10047);
      }
      template<class stce_9937, class stce_9939 >
      static stce_10016 stce_10044(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839, const stce_10016 &pval, const int stce_10047) {
        return stce_10044(stce_9839.stce_9941, stce_9839.stce_9942() * pval, stce_10047);
      }
      template<class stce_9938, class stce_9939 >
      static stce_10016 stce_10044(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839, const stce_10016 &pval, const int stce_10047) {
        return stce_10044(stce_9839.stce_9943, stce_9839.stce_9944() * pval, stce_10047);
      }
      template<class stce_9935 >
      static bool stce_10053(const dco::internal::active_type<stce_9771, stce_9935> &stce_9839) {
        return stce_9839.stce_10051;
      }
      template<class stce_9937, class stce_9938, class stce_9939 >
      static bool stce_10053(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839) {
        return stce_10053(stce_9839.stce_9941) || stce_10053(stce_9839.stce_9943);
      }
      template<class stce_9945, class stce_9939 >
      static bool stce_10053(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839) {
        return stce_10053(stce_9839.stce_9946);
      }
      template<class stce_9937, class stce_9939 >
      static bool stce_10053(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839) {
        return stce_10053(stce_9839.stce_9941);
      }
      template<class stce_9938, class stce_9939 >
      static bool stce_10053(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839) {
        return stce_10053(stce_9839.stce_9943);
      }
      typedef void stce_10049;
      typedef void tape_t;
      typedef DCO_TAPE_INT stce_10045;
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
    struct stce_10054 {
    typedef stce_9949 stce_9989;
    typedef std::bitset<size> stce_10055;
    typedef stce_9803<stce_9771> mode_t;
      stce_10055 stce_10056;
      stce_10054() {}
      DCO_STRONG_INLINE stce_10054 &operator = (const stce_10054 &stce_10042) {
        stce_10056 = stce_10042.stce_10056;
        return *this;
      }
       void stce_9986(const size_t &stce_9806) {
        stce_10056.set(stce_9806);
      }
       bool stce_10043(const size_t &stce_9806) const {
        return stce_10056[stce_9806];
      }
       void clear() {
        stce_10056.reset();
      }
      template<class stce_9888>
      DCO_ALWAYS_INLINE void stce_9889(const stce_9888 &stce_9990) {
        stce_10054 stce_10057;
        stce_10058(stce_9990, stce_10057.stce_10056);
        stce_10056 = stce_10057.stce_10056;
      }
      template<class stce_9935 >
      static void stce_10058(const dco::internal::active_type<stce_9771, stce_9935> &stce_9839, stce_10055 &stce_10048) {
        stce_10048 |= stce_9839.stce_10056;
      }
      template<class stce_9937, class stce_9938, class stce_9939 >
      static void stce_10058(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839, stce_10055 &stce_10048) {
        stce_10058(stce_9839.stce_9941, stce_10048);
        stce_10058(stce_9839.stce_9943, stce_10048);
      }
      template<class stce_9945, class stce_9939 >
      static void stce_10058(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839, stce_10055 &stce_10048) {
        stce_10058(stce_9839.stce_9946, stce_10048);
      }
      template<class stce_9937, class stce_9939 >
      static void stce_10058(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839, stce_10055 &stce_10048) {
        stce_10058(stce_9839.stce_9941, stce_10048);
      }
      template<class stce_9938, class stce_9939 >
      static void stce_10058(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839, stce_10055 &stce_10048) {
        stce_10058(stce_9839.stce_9943, stce_10048);
      }
    };
  }
}
       
       
namespace dco {
  namespace internal {
struct interpretation_settings {
  bool reset;
  bool stce_10059;
  interpretation_settings(bool stce_10060=false, bool stce_10061=false)
    : reset(stce_10060), stce_10059(stce_10061) {
  }
};
}
}
namespace dco {
    class tape_options {
    private:
      size_t stce_10062, stce_10063;
      bool stce_10064;
      std::string stce_10065;
      bool stce_10066;
      bool stce_10067;
    public:
      tape_options()
 : stce_10062(128*1024*1024),
          stce_10063(0),
   stce_10064(false),
   stce_10065("dco_tape_data.bin"),
   stce_10066(false),
          stce_10067(false) { }
      bool operator==(const tape_options& stce_9929) {
        return
          stce_10062 == stce_9929.stce_10062 &&
          stce_10063 == stce_9929.stce_10063 &&
          stce_10064 == stce_9929.stce_10064 &&
          stce_10066 == stce_9929.stce_10066 &&
          stce_10065 == stce_9929.stce_10065 &&
          stce_10067 == stce_9929.stce_10067;
      }
      bool operator!=(const tape_options& stce_9929) { return !operator==(stce_9929); }
      DCO_TAPE_INT chunk_size() const {
        if (helper::stce_9817) {
          throw dco::exception::create<std::runtime_error>("chunk_size removed; use chunk_size_in_byte() instead");
        }
        return -1;
      }
      size_t chunk_size_in_byte() const { return stce_10062; }
      size_t& chunk_size_in_byte() { return stce_10062; }
      void set_chunk_size_in_byte (mem_long_t size) {
        if (static_cast<double>(size) / static_cast<double>(std::numeric_limits<size_t>::max()) > 1.0) {
          DCO_LOG(logWARNING) << "set_chunk_size: requesting too much memory for your system (size_t too small). Resetting to std::numeric_limits<size_t>::max()";
          stce_10062 = std::numeric_limits<size_t>::max();
        } else {
          stce_10062 = static_cast<size_t>(size);
        }
      }
      void set_chunk_size_in_kbyte(double size) { set_chunk_size_in_byte(1024*static_cast<mem_long_t>(size)); }
      void set_chunk_size_in_mbyte(double size) { set_chunk_size_in_kbyte(1024.*size); }
      void set_chunk_size_in_gbyte(double size) { set_chunk_size_in_mbyte(1024.*size); }
      size_t blob_size_in_byte() const { return stce_10063; }
      size_t& blob_size_in_byte() { return stce_10063; }
      void set_blob_size_in_byte (mem_long_t size) {
        if (static_cast<double>(size) / static_cast<double>(std::numeric_limits<size_t>::max()) > 1.0) {
          DCO_LOG(logWARNING) << "set_blob_size: requesting too much memory for your system (size_t too small). Resetting to std::numeric_limits<size_t>::max()";
          stce_10063 = std::numeric_limits<size_t>::max();
        } else {
          stce_10063 = static_cast<size_t>(size);
        }
      }
      void set_blob_size_in_kbyte(double size) { set_blob_size_in_byte(1024*static_cast<mem_long_t>(size)); }
      void set_blob_size_in_mbyte(double size) { set_blob_size_in_kbyte(1024.*size); }
      void set_blob_size_in_gbyte(double size) { set_blob_size_in_mbyte(1024.*size); }
      bool &deallocation_on_reset() { return stce_10064; }
      const bool &deallocation_on_reset() const { return stce_10064; }
      std::string &filename() { return stce_10065; }
      const std::string &filename() const { return stce_10065; }
      bool& write_to_file() { return stce_10066; }
      bool write_to_file() const { return stce_10066; }
      bool& alloc_use_huge_pages() { return stce_10067; }
      bool alloc_use_huge_pages() const { return stce_10067; }
};
}
       
       
namespace dco {
template <typename stce_9782, typename stce_10068 = void> struct trait_size_of;
template<typename stce_9782> mem_long_t size_of(const stce_9782& stce_9839) {
  return trait_size_of<stce_9782>::get(stce_9839);
}
template<typename stce_9782> mem_long_t size_of(const stce_9782& stce_9839, int stce_10069) {
  return trait_size_of<stce_9782>::get(stce_9839, stce_10069);
}
template <typename stce_9782, typename stce_10068> struct trait_size_of {
  static mem_long_t get(const stce_9782&) {
    if (helper::stce_9817) {
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
                        const int stce_10069 = tape_t::stce_10070) {
  return tape->stce_10071(stce_10069);
  }
};
template<typename stce_9782> struct trait_size_of<std::vector<stce_9782> > {
  static mem_long_t get(const std::vector<stce_9782>& stce_9839) { return stce_9839.size() * sizeof(stce_9782) + sizeof(stce_9839); }
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
template<> struct trait_size_of<float> { static mem_long_t get(const float&) { return sizeof(float); } };
template<> struct trait_size_of<double> { static mem_long_t get(const double&) { return sizeof(double); } };
template<> struct trait_size_of<long double> { static mem_long_t get(const long double&) { return sizeof(long double); } };
template<> struct trait_size_of<bool> { static mem_long_t get(const bool&) { return sizeof(bool); } };
}
       
namespace dco {
namespace helper {
template <typename, typename stce_10017>
  struct stce_10072 {
  stce_10072(stce_10017*) {}
  template <typename stce_9782> void register_input(const stce_9782&) { }
  template <typename stce_9782> void register_output(const stce_9782&) { }
};
}
}
       
       
namespace dco {
template<class stce_9761>
  struct stce_10073 {
  virtual ~stce_10073(){};
  virtual stce_9761 stce_10074(const DCO_TAPE_INT stce_9806)=0;
  virtual void stce_10075(const DCO_TAPE_INT stce_9806, const stce_9761 &stce_10076)=0;
};
}
namespace dco {
namespace helper {
template<typename stce_10077>
  class stce_10078 {
  template<memory_model::TYPE, typename, memory_model::TYPE>
    friend class internal::tape;
  typedef typename stce_10077::stce_10079 stce_9764;
protected:
  stce_10077 *stce_10080;
public:
  virtual ~stce_10078() {
  }
  dco::stce_10073<stce_9764> *stce_10081;
public:
  void set_tape(stce_10077 *stce_9901) {
    if (stce_10080 != 0) { return; }
    else { stce_10080 = stce_9901; }
  }
  inline stce_10077 *get_tape() { return stce_10080; }
  stce_10078() : stce_10080(0) {}
  virtual void stce_10082(stce_9844 &stce_9845, DCO_TAPE_INT stce_9853) {
    (void)stce_9845;
    (void)stce_9853;
  }
  virtual mem_long_t size_in_byte() { return sizeof(stce_10077); }
};
}
}
       
namespace dco {
namespace helper {
template<class stce_10083, class stce_10077>
  class stce_10084 : public stce_10078<stce_10077> {
private:
  class stce_10085 {
public:
    virtual ~stce_10085() {};
    virtual mem_long_t size_in_byte() = 0;
    virtual void stce_10086(const std::string&, const std::string&) {};
  };
  template <typename stce_10087>
    class stce_10088 : public stce_10085 {
public:
    const stce_10087 _data;
    stce_10088(const stce_10087 &stce_10089) : _data(stce_10089) {}
    const stce_10087 &stce_10090() const { return _data; }
  };
  template <typename stce_10087>
    class stce_10091 : public stce_10088<stce_10087> {
public:
    stce_10091(stce_10087 data) : stce_10088<stce_10087>(data) {}
    virtual ~stce_10091() { }
    virtual mem_long_t size_in_byte() {
      return size_of(stce_10088<stce_10087>::_data);
    }
    virtual void stce_10086(const std::string& stce_10092, const std::string& stce_10093) {
      if (false && IS_DEBUG && stce_10093 != stce_10092) {
        std::stringstream stce_10057;
        stce_10057 << "read_data error: tried to read " << stce_10092 << " instead of " << stce_10093 << "!";
        throw dco::exception::create<std::runtime_error>(stce_10057.str());
      }
    }
  };
  template <typename stce_10087>
    class stce_10094 : public stce_10088<stce_10087*> {
public:
    const DCO_TAPE_INT stce_10095;
    stce_10094(const stce_10087 *data, DCO_TAPE_INT stce_10096) : stce_10088<stce_10087*>(new stce_10087[size_t(stce_10096)]), stce_10095(stce_10096) {
      for (DCO_TAPE_INT stce_9827 = 0; stce_9827 < stce_10096; ++stce_9827)
        this->_data[stce_9827] = data[stce_9827];
    }
    stce_10094(const stce_10087 *data, const DCO_TAPE_INT stce_10097, const DCO_TAPE_INT stce_10096) :
      stce_10088<stce_10087 *>(new stce_10087[stce_9752(stce_10096)]), stce_10095(stce_10096) {
      for (DCO_TAPE_INT stce_9827 = 0, stce_9806 = 0; stce_9827 < stce_10096; ++stce_9827, stce_9806 += stce_10097)
        this->_data[stce_9827] = data[stce_9806];
    }
    virtual ~stce_10094() {
      delete [] this->_data;
    }
    virtual mem_long_t size_in_byte() {
      return stce_10095 * static_cast<DCO_TAPE_INT>(sizeof(stce_10087));
    }
    virtual void stce_10086(const std::string& stce_10092, const std::string& stce_10093) {
      (void) stce_10092; (void) stce_10093;
      if (false && IS_DEBUG && stce_10093+"*" != stce_10092) {
        std::stringstream stce_10057;
        stce_10057 << "read_data error: tried to read " << stce_10092 << " instead of " << stce_10093 << "!";
        throw dco::exception::create<std::runtime_error>(stce_10057.str());
      }
    }
  };
  DCO_TAPE_INT stce_10098;
  std::vector<stce_10085 *> stce_10099;
  std::vector<std::string> stce_10100;
protected:
  virtual ~stce_10084() {
    for (size_t stce_9827 = 0; stce_9827 < stce_10099.size(); stce_9827++)
      delete stce_10099[stce_9827];
    stce_10099.clear();
  }
public:
  stce_10084(): stce_10078<stce_10077>(), stce_10098(0) {}
  virtual mem_long_t size_in_byte() {
    mem_long_t stce_9916 = stce_10078<stce_10077>::size_in_byte();
    for (size_t stce_9827 = 0; stce_9827 < stce_10099.size(); stce_9827++) {
      stce_9916 += stce_10099[stce_9827]->size_in_byte();
    }
    stce_9916 += stce_10099.size() * sizeof(stce_10085*);
    stce_9916 += stce_10100.size() * sizeof(std::string);
    return stce_9916;
  }
  template<typename stce_10087>
    void write_data(const stce_10087 &stce_10101) {
    stce_10099.push_back(new stce_10091<stce_10087>(stce_10101));
    if (IS_DEBUG)
      stce_10100.push_back(type_name<stce_10087>());
  }
  template<typename stce_10087>
    void write_data(const stce_10087 *const stce_10101, const DCO_TAPE_INT stce_10096) {
    stce_10099.push_back(new stce_10094<stce_10087>(stce_10101, stce_10096));
    if (IS_DEBUG)
      stce_10100.push_back(type_name<stce_10087>());
  }
  template<typename stce_10087>
    void write_data(const stce_10087 *const &stce_10101, const DCO_TAPE_INT stce_10097, const DCO_TAPE_INT stce_10096) {
    stce_10099.push_back(new stce_10094<stce_10087>(stce_10101, stce_10097, stce_10096));
    if (IS_DEBUG) {
      stce_10100.push_back(type_name<stce_10087>());
    }
  }
  template<typename stce_10087>
    const stce_10087 &read_data() {
    const stce_10087 &stce_10101 = static_cast<stce_10088<stce_10087>* >(stce_10099[stce_10098])->stce_10090();
    if (IS_DEBUG)
      stce_10099[stce_10098]->stce_10086(type_name<stce_10087>(), stce_10100.at(stce_10098));
    ++stce_10098;
    if (stce_10098 == static_cast<DCO_TAPE_INT>(stce_10099.size())) stce_10098 = 0;
    return stce_10101;
  }
};
}
}
       
namespace dco {
namespace helper {
template<class stce_10083, class stce_10077>
  class stce_10102 :
           public stce_10084<stce_10083, stce_10077>,
           public stce_10072<stce_10077, stce_10102<stce_10083, stce_10077> > {
  typedef stce_10072<stce_10077, stce_10102> stce_10103;
protected:
  std::vector<DCO_TAPE_INT> stce_9851;
  std::vector<DCO_TAPE_INT> stce_9852;
  DCO_TAPE_INT stce_10104;
  DCO_TAPE_INT stce_10105;
public:
  stce_10077* &tape() { return this->stce_10080; }
  virtual mem_long_t size_in_byte() {
    mem_long_t stce_9916 = stce_10084<stce_10083, stce_10077>::size_in_byte();
    stce_9916 += stce_9851.size() * sizeof(DCO_TAPE_INT);
    stce_9916 += stce_9852.size() * sizeof(DCO_TAPE_INT);
    stce_9916 += 2 * sizeof(DCO_TAPE_INT);
    return stce_9916;
  }
  size_t stce_10106() {
    return stce_9851.size();
  }
  size_t stce_10107() {
    return stce_9852.size();
  }
  void stce_10082(stce_9844 &stce_9845, DCO_TAPE_INT stce_9853) {
    stce_9845.stce_9850(stce_9851, stce_9852, stce_9853 - static_cast<DCO_TAPE_INT>(stce_9852.size()) - 1);
  }
public:
  void stce_10108(const stce_10083 &stce_9839) {
    if ((stce_9839.tape() != 0) && (this->stce_10080 != stce_9839.tape()))
      throw dco::exception::create<std::runtime_error>("impossible binding tape - wrong tape in variable!");
  }
protected:
  ~stce_10102() {
    if (IS_DEBUG && stce_10104 != 0) {
      DCO_LOG(dco::logERROR) << "not all input adjoints incremented in external adjoint object.";
    }
    if (IS_DEBUG && stce_10105 != 0) {
      DCO_LOG(dco::logERROR) << "not all output adjoints read in external adjoint object.";
    }
  }
public:
  stce_10102(const std::pair<DCO_TAPE_INT, DCO_TAPE_INT> &stce_10109):
    stce_10084<stce_10083, stce_10077>(),
    stce_10103(this),
    stce_10104(0),
    stce_10105(0) {
    stce_9851.reserve(stce_10109.first);
    stce_9852.reserve(stce_10109.second);
  }
  stce_10102(): stce_10084<stce_10083, stce_10077>(),
    stce_10103(this),
    stce_10104(0), stce_10105(0) { }
  typename stce_10083::VALUE_TYPE register_input(const stce_10083 &stce_9839) {
    stce_10108(stce_9839);
    const DCO_TAPE_INT stce_9806 = stce_9839._tape_index();
    stce_9851.push_back(stce_9806);
    stce_10103::register_input(stce_9839);
    return stce_9839._value();
  }
  void register_input(const stce_10083 *const stce_9839, typename stce_10083::VALUE_TYPE *stce_10110, const DCO_TAPE_INT stce_10096) {
    stce_9851.reserve(stce_9851.size() + stce_10096);
    for (int stce_9827 = 0; stce_9827 < stce_10096; ++stce_9827) {
      stce_10108(stce_9839[stce_9827]);
      const DCO_TAPE_INT stce_9806 = stce_9839[stce_9827]._tape_index();
      stce_9851.push_back(stce_9806);
      stce_10110[stce_9827] = stce_9839[stce_9827]._value();
      stce_10103::register_input(stce_9839[stce_9827]);
    }
  }
  void register_input(const std::vector<stce_10083> &stce_9839, std::vector<typename stce_10083::VALUE_TYPE> &stce_10110) {
    assert(stce_9839.size() == stce_10110.size());
    register_input(&(stce_9839[0]), &(stce_10110[0]), stce_9839.size());
  }
  std::vector<typename stce_10083::VALUE_TYPE> register_input(const std::vector<stce_10083> &stce_9839) {
    std::vector<typename stce_10083::VALUE_TYPE> stce_10110(stce_9839.size());
    register_input(stce_9839, stce_10110);
    return stce_10110;
  }
  void register_output(stce_10083 *stce_10111, const size_t stce_10096) {
    if (this->stce_10080 == NULL) {
      throw dco::exception::create<std::runtime_error>("impossible binding output - no tape available");
    } else {
      register_output(NULL, stce_10111, stce_10096);
    }
  }
  void register_output(const typename stce_10083::VALUE_TYPE *const stce_10112, stce_10083 *stce_10111, const size_t stce_10096) {
    if (this->stce_10080 == NULL) {
      throw dco::exception::create<std::runtime_error>("impossible binding output - no tape available");
    } else {
      stce_9852.reserve(stce_9852.size() + stce_10096);
      for (size_t stce_9827 = 0; stce_9827 < stce_10096; ++stce_9827) {
        if (stce_10112) {
          stce_10111[stce_9827] = stce_10112[stce_9827];
        }
        this->stce_10080->register_variable(stce_10111[stce_9827]);
        stce_9852.push_back(stce_10111[stce_9827]._tape_index());
        stce_10103::register_output(stce_10111[stce_9827]);
      }
    }
  }
  void register_output(const std::vector<typename stce_10083::VALUE_TYPE> &stce_10112, std::vector<stce_10083> &stce_10111) {
    assert(stce_10112.size() == stce_10111.size());
    register_output(&(stce_10112[0]), &(stce_10111[0]), stce_10112.size());
  }
  std::vector<stce_10083> register_output(const std::vector<typename stce_10083::VALUE_TYPE> &stce_10112) {
    std::vector<stce_10083> stce_10111(stce_10112.size());
    register_output(stce_10112, stce_10111);
    return stce_10111;
  }
  void register_output(std::vector<stce_10083> &stce_10111) {
    register_output(&(stce_10111[0]), stce_10111.size());
  }
  stce_10083 register_output(const typename stce_10083::VALUE_TYPE &stce_10113, stce_10077 *tape = NULL) {
    stce_10083 stce_10114;
    if (tape != NULL) {
      if (this->stce_10080 != NULL && this->stce_10080 != tape) {
        throw dco::exception::create<std::runtime_error>("impossible binding output in external function (register_output) - tape of inputs and outputs differ!");
      }
      this->stce_10080 = tape;
    }
    if (this->stce_10080 != NULL) {
      stce_10114 = stce_10113;
      this->stce_10080->register_variable(stce_10114);
    } else
      throw dco::exception::create<std::runtime_error>("impossible binding output in external function - no tape available");
    stce_9852.push_back(stce_10114._tape_index());
    stce_10103::register_output(stce_10114);
    return stce_10114;
  }
  inline typename dco::mode<stce_10083>::stce_10016 get_output_adjoint() {
    DCO_TAPE_INT stce_9806 = stce_10105;
    stce_10105++;
    if (stce_10105 == static_cast<DCO_TAPE_INT>(stce_9852.size()))
      stce_10105 = 0;
    typename dco::mode<stce_10083>::stce_10016 back = this->stce_10081->stce_10074(stce_9852[static_cast<size_t>(stce_9806)]);
    return back;
  }
  void get_output_adjoint(typename stce_10083::VALUE_TYPE *stce_9807, const size_t stce_10096) {
    DCO_TAPE_INT stce_9806 = stce_10105;
    for (size_t stce_9827 = 0; stce_9827 < stce_10096; ++stce_9827) {
      stce_9807[stce_9827] = get_output_adjoint();
      stce_9806++;
    }
    if (stce_10105 == static_cast<DCO_TAPE_INT>(stce_9852.size())) {
      stce_10105 = 0;
    }
  }
  void get_output_adjoint(std::vector<typename stce_10083::VALUE_TYPE> &stce_9807) {
    assert(stce_9807.size());
    get_output_adjoint(&(stce_9807[0]), stce_9807.size());
  }
  void increment_input_adjoint(const typename stce_10083::VALUE_TYPE *const stce_10076, const DCO_TAPE_INT stce_10096) {
    for (DCO_TAPE_INT stce_9827 = 0; stce_9827 < stce_10096; ++stce_9827) {
      this->stce_10081->stce_10075(stce_9851[stce_10104 + stce_9827], stce_10076[stce_9827]);
    }
    stce_10104 += stce_10096;
    if (stce_10104 == static_cast<DCO_TAPE_INT>(stce_9851.size()))
      stce_10104 = 0;
  }
  void increment_input_adjoint(const std::vector<typename stce_10083::VALUE_TYPE> &stce_10076) {
    assert(stce_10076.size() != 0);
    increment_input_adjoint(&(stce_10076[0]), static_cast<DCO_TAPE_INT>(stce_10076.size()));
  }
  bool stce_10115() {
    if (stce_10104 == 0) return true;
    else return false;
  }
  bool stce_10116() {
    if (stce_10105 == 0) return true;
    else return false;
  }
  void increment_input_adjoint(const typename dco::mode<stce_10083>::stce_10016 &stce_10076) {
    DCO_TAPE_INT stce_9806 = stce_10104;
    stce_10104++;
    if (stce_10104 == static_cast<DCO_TAPE_INT>(stce_9851.size()))
      stce_10104 = 0;
    this->stce_10081->stce_10075(stce_9851[static_cast<size_t>(stce_9806)], stce_10076);
  }
};
}
}
       
namespace dco {
namespace ACM {
template<class stce_10077, class stce_9978>
  struct baseclass : dco::helper::stce_10078<stce_10077> {
private:
  std::vector<stce_9978*> stce_10117;
  std::vector<stce_9978*> stce_10118;
  std::vector<DCO_TAPE_INT> stce_10119;
  std::vector<DCO_TAPE_INT> stce_10120;
  std::vector<void*> _data;
  int stce_10121;
public:
  typedef stce_9978 active_t;
  typedef stce_10077 tape_t;
  baseclass(tape_t *tape) : stce_10121(0) {
    this->set_tape(tape);
  }
  size_t register_input(stce_9978 &stce_9839) {
    size_t stce_9836 = stce_10117.size();
    stce_10117.push_back(&stce_9839);
    stce_10120.push_back(stce_9839._tape_index());
    return stce_9836;
  }
  size_t register_input(std::vector<stce_9978> &stce_10122) {
    size_t stce_9836 = stce_10117.size();
    for(int stce_9827=0;stce_9827<stce_10122.size();++stce_9827) {
      const stce_9978 &stce_10123 = stce_10122[stce_9827];
      this->register_input(stce_10123);
    }
    return stce_9836;
  }
  size_t register_output(stce_9978 &stce_9839) {
    if(stce_9839._tape_index()==0) {
      this->get_tape()->register_variable(stce_9839);
    }
    size_t stce_9836 = stce_10118.size();
    stce_10118.push_back(&stce_9839);
    stce_10119.push_back(stce_9839._tape_index());
    return stce_9836;
  }
  int register_output(std::vector<stce_9978> &stce_10122) {
    int stce_9836 = stce_10118.size();
    for(int stce_9827=0;stce_9827<stce_10122.size();++stce_9827) {
      if(stce_10122[stce_9827]._tape_index()==0) {
        this->get_tape()->register_variable(stce_10122[stce_9827]);
      }
      stce_10119.push_back(stce_10122[stce_9827]._tape_index());
      stce_10118.push_back(&stce_10122[stce_9827]);
    }
    return stce_9836;
  }
  template<typename stce_9782> int write_data(const stce_9782 &data) {
    stce_9782* stce_10124 = new stce_9782(data);
    int stce_9836 = _data.size();
    _data.push_back(stce_10124);
    return stce_9836;
  }
  template<typename stce_9782> const stce_9782& read_data(const int stce_10125=-1) {
    int stce_9806=0;
    if(stce_10125==-1) {
      stce_9806 = stce_10121%_data.size(); stce_10121++;
    }
    void* stce_10126 = _data[stce_9806];
    return *static_cast<stce_9782*>(stce_10126);
  }
  size_t stce_10127() {
    return stce_10117.size();
  }
  size_t stce_10128() {
    return stce_10117.size();
  }
  typename stce_9978::VALUE_TYPE& input_value(const int stce_9806) {
    return stce_10117[stce_9806]->_value();
  }
  typename stce_9978::VALUE_TYPE& output_value(const int stce_9806) {
    return stce_10118[stce_9806]->_value();
  }
  typename dco::mode<stce_9978>::stce_10016& input_adjoint(const int stce_9806) {
    return this->get_tape()->_adjoint( stce_10120[stce_9806] );
  }
  typename dco::mode<stce_9978>::stce_10016& output_adjoint(const int stce_9806) {
    return this->get_tape()->_adjoint( stce_10119[stce_9806] );
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
template<class stce_10077, class stce_9978>
  static void stce_10129( baseclass<stce_10077, stce_9978> *stce_10130 ) {
  stce_10130->evaluate_adjoint();
}
}
}
namespace dco {
namespace stce_10131 {
    template<class stce_9987>
    struct stce_10132 {
        std::vector<stce_9987*> stce_10117;
        std::vector<stce_9987*> stce_10118;
        std::vector<void*> _data;
        int stce_10121;
        stce_10132() : stce_10121(0) {}
        virtual ~stce_10132() {}
        template<typename stce_9782> int write_data(const stce_9782 &data) {
            stce_9782* stce_10124 = new stce_9782(data);
            int stce_9836 = _data.size();
            _data.push_back(stce_10124);
            return stce_9836;
        }
        template<typename stce_9782> const stce_9782& read_data(const int stce_10125=-1) {
            int stce_9806=0;
            if(stce_10125==-1) {
                stce_9806 = stce_10121%_data.size(); stce_10121++;
            }
            void* stce_10126 = _data[stce_9806];
            return *static_cast<stce_9782*>(stce_10126);
        }
        size_t stce_10127() {
            return stce_10117.size();
        }
        size_t stce_10128() {
            return stce_10118.size();
        }
        stce_9987& input_value(const int stce_9806) {
            return *stce_10117[stce_9806];
        }
        stce_9987& output_value(const int stce_9806) {
            return *stce_10118[stce_9806];
        }
        virtual void evaluate_primal() {
          std::cout << "WARNING: you need to implement evaluate_primal!!!!\n";
        }
    };
    template<class stce_9987>
    struct stce_10133 : stce_10132<stce_9987> {
        int register_input(stce_9987 &stce_9839) {
            int stce_9836 = this->stce_10117.size();
            this->stce_10117.push_back(&stce_9839);
            return stce_9836;
        }
        int register_input(std::vector<stce_9987> &stce_10122) {
            int stce_9836 = this->stce_10117.size();
            for(int stce_9827=0;stce_9827<stce_10122.size();++stce_9827) {
                const stce_9987 &stce_10123 = stce_10122[stce_9827];
                this->register_input(stce_10123);
            }
            return stce_9836;
        }
        int register_output(stce_9987 &stce_9839) {
            int stce_9836 = this->stce_10118.size();
            this->stce_10118.push_back(&stce_9839);
            return stce_9836;
        }
        int register_output(std::vector<stce_9987> &stce_10122) {
            int stce_9836 = this->stce_10118.size();
            for(int stce_9827=0;stce_9827<stce_10122.size();++stce_9827) {
                register_output(stce_10122[stce_9827]);
            }
            return stce_9836;
        }
    };
    template<class stce_10077, class stce_9978>
    struct stce_10134 : public dco::helper::stce_10078<stce_10077>, public stce_10132<stce_9978> {
    private:
        std::vector<DCO_TAPE_INT> stce_10119;
        std::vector<DCO_TAPE_INT> stce_10120;
    public:
        typedef stce_9978 active_t;
        typedef stce_10077 tape_t;
        stce_10134(tape_t *tape) {
            this->set_tape(tape);
        }
        stce_9978& stce_10135(const int stce_9806) {
            return *(this->stce_10117[stce_9806]);
        }
        stce_9978& stce_10136(const int stce_9806) {
            return *(this->stce_10118[stce_9806]);
        }
        typename stce_9978::VALUE_TYPE& input_value(const int stce_9806) {
            return this->stce_10117[stce_9806]->_value();
        }
        typename stce_9978::VALUE_TYPE& output_value(const int stce_9806) {
            return this->stce_10118[stce_9806]->_value();
        }
        typename dco::mode<stce_9978>::stce_10016& input_adjoint(const int stce_9806) {
            return this->get_tape()->_adjoint( stce_10120[stce_9806] );
        }
        typename dco::mode<stce_9978>::stce_10016& output_adjoint(const int stce_9806) {
            return this->get_tape()->_adjoint( stce_10119[stce_9806] );
        }
        int register_input(stce_9978 &stce_9839) {
            int stce_9836 = this->stce_10117.size();
            this->stce_10117.push_back(&stce_9839);
            this->stce_10120.push_back(stce_9839._tape_index());
            return stce_9836;
        }
        int register_input(std::vector<stce_9978> &stce_10122) {
            int stce_9836 = this->stce_10117.size();
            for(int stce_9827=0;stce_9827<stce_10122.size();++stce_9827) {
                const stce_9978 &stce_10123 = stce_10122[stce_9827];
                this->register_input(stce_10123);
            }
            return stce_9836;
        }
        int register_output(stce_9978 &stce_9839) {
            if(stce_9839._tape_index()==0) {
                this->get_tape()->register_variable(stce_9839);
            }
            int stce_9836 = this->stce_10118.size();
            this->stce_10118.push_back(&stce_9839);
            this->stce_10119.push_back(stce_9839._tape_index());
            return stce_9836;
        }
        int register_output(std::vector<stce_9978> &stce_10122) {
            int stce_9836 = this->stce_10118.size();
            for(int stce_9827=0;stce_9827<stce_10122.size();++stce_9827) {
                if(stce_10122[stce_9827]._tape_index()==0) {
                    this->get_tape()->register_variable(stce_10122[stce_9827]);
                }
                this->stce_10119.push_back(stce_10122[stce_9827]._tape_index());
                this->stce_10118.push_back(&stce_10122[stce_9827]);
            }
            return stce_9836;
        }
        virtual void evaluate_augmented_primal() {
            std::cout << "WARNING: you need to implement evaluate_augmented_primal!!!!\n";
        }
        virtual void evaluate_adjoint() {
            std::cout << "WARNING: you need to implement evaluate_adjoint!!!!\n";
        }
    };
    template<class stce_10077, class stce_9978>
    static void stce_10129( stce_10134<stce_10077, stce_9978> *stce_10130 ) {
        stce_10130->evaluate_adjoint();
    }
}
}
       
       
       
namespace dco {
namespace internal {
template <typename stce_10137>
  class stce_10138 {
protected:
  stce_10137& stce_10139;
public:
  stce_10138(stce_10137& stce_10140) : stce_10139(stce_10140) { }
  stce_10137& stce_10140() { return stce_10139; }
};
}
}
namespace dco {
namespace internal {
template<memory_model::TYPE memory_model,
         typename stce_10141,
         memory_model::TYPE stce_10142=memory_model::stce_9766> class stce_10143 {};
}
}
namespace dco {
namespace internal {
template<class stce_10141, memory_model::TYPE stce_10144>
  class stce_10143<memory_model::BLOB_TAPE_SPLINT, stce_10141, stce_10144>
  : public stce_10138<stce_10141>, object_logger
{
  typedef stce_10138<stce_10141> stce_10145;
  typedef typename stce_9790<stce_10141>::type stce_10146;
public:
  struct entry {
    stce_10146* stce_10147;
    DCO_TAPE_INT* stce_10148;
    DCO_ALWAYS_INLINE entry(stce_10146* stce_10149,
                            DCO_TAPE_INT* stce_10150)
      : stce_10147(stce_10149),
        stce_10148(stce_10150) {};
    DCO_ALWAYS_INLINE DCO_TAPE_INT arg() {
      return (*stce_10148 < 0) ? std::abs(stce_9752(*stce_10148))-1 : stce_9752(*stce_10148);
    }
    DCO_ALWAYS_INLINE stce_10146& pval() { return *stce_10147; }
    DCO_ALWAYS_INLINE stce_10146 pval() const { return *stce_10147; }
    DCO_ALWAYS_INLINE bool stce_10151() { return *stce_10148 < 0; }
    DCO_ALWAYS_INLINE static bool stce_10151(DCO_TAPE_INT* stce_10152) { return *stce_10152 < 0; }
    DCO_ALWAYS_INLINE static DCO_TAPE_INT stce_10153(DCO_TAPE_INT stce_10048) { return -stce_10048-1; }
    static size_t size_of() { return sizeof(stce_10146) + sizeof(DCO_TAPE_INT); }
  };
  class iterator {
public:
    DCO_TAPE_INT stce_10154;
    entry stce_10155;
    template<class stce_10156>
      void stce_10157(stce_10156 &stce_10158) const {
      assert(sizeof(iterator) <= sizeof(stce_10156)) ;
      stce_10159(&stce_10158, this, sizeof(iterator));
    }
    template<class stce_10156>
      void stce_10160(const stce_10156 &stce_10158) {
      assert(sizeof(iterator) <= sizeof(stce_10156)) ;
      stce_10159(static_cast<void*>(this), &stce_10158, sizeof(iterator));
    }
    DCO_ALWAYS_INLINE const entry* operator->() const { return &stce_10155; }
    DCO_ALWAYS_INLINE entry* operator->() { return &stce_10155; }
    DCO_ALWAYS_INLINE void stce_10161() {
      stce_10155.stce_10147++;
      stce_10155.stce_10148--;
      if (stce_9754 && static_cast<void*>(stce_10155.stce_10147) >= static_cast<void*>(stce_10155.stce_10148)) {
        throw dco::exception::create<std::runtime_error>("Blob tape (splint) container overflow. Allocate bigger blob or use chunk tape container instead.");
      }
    }
    DCO_ALWAYS_INLINE DCO_TAPE_INT stce_10162() { return ++stce_10154; }
    DCO_ALWAYS_INLINE explicit iterator() : stce_10154(-1), stce_10155(NULL, NULL) { }
    DCO_ALWAYS_INLINE explicit iterator(stce_10146* stce_10149, DCO_TAPE_INT* stce_10150, DCO_TAPE_INT index = DCO_TAPE_INT(-1)) :
      stce_10154(index), stce_10155(stce_10149, stce_10150) { }
    DCO_ALWAYS_INLINE DCO_TAPE_INT index() const { return stce_10154; }
    DCO_ALWAYS_INLINE DCO_TAPE_INT& index() { return stce_10154; }
    bool stce_10151() const {
      return entry::stce_10151(stce_10155.stce_10148);
    }
    DCO_ALWAYS_INLINE void operator--() {
      if(entry::stce_10151(stce_10155.stce_10148))
        stce_10154--;
      stce_10155.stce_10147--;
      stce_10155.stce_10148++;
    }
    DCO_ALWAYS_INLINE iterator operator--(int) {
      iterator stce_10057(*this);
      operator--();
      return stce_10057;
    }
    DCO_ALWAYS_INLINE bool operator==(const iterator& stce_9929) const {
      return stce_10155.stce_10147 == stce_9929.stce_10155.stce_10147;
    }
    DCO_ALWAYS_INLINE bool operator!=(const iterator& stce_9929) const {
      return !operator==(stce_9929);
    }
    DCO_ALWAYS_INLINE bool operator<(const iterator& stce_9929) const {
      return stce_10155.stce_10147 < stce_9929.stce_10155.stce_10147;
    }
    DCO_ALWAYS_INLINE bool operator>=(const iterator& stce_9929) const {
      return !operator<(stce_9929);
    }
    DCO_ALWAYS_INLINE bool operator>(const iterator& stce_9929) const {
      return stce_10155.stce_10147 > stce_9929.stce_10155.stce_10147;
    }
    DCO_ALWAYS_INLINE bool operator<=(const iterator& stce_9929) const {
      return !operator>(stce_9929);
    }
  };
  class stce_9920
  {
    iterator& stce_10163;
    DCO_TAPE_INT stce_10154;
public:
    stce_9920(iterator& current) :
      stce_10163(current), stce_10154(0){
    }
    iterator& current() { return stce_10163; }
    stce_9920& operator=(stce_9920 stce_9929) {
      swap(stce_9929);
      return *this;
    }
    DCO_TAPE_INT index() const {
      return stce_10154;
    }
    void swap(stce_9920& stce_9929)
    {
      std::swap(stce_10163, stce_9929.stce_10163);
      std::swap(stce_10154, stce_9929.stce_10154);
    }
    DCO_ALWAYS_INLINE_IF_AGRESSIVE void insert(const DCO_TAPE_INT& stce_10048, const stce_10146& pval) {
      if (IS_DEBUG && stce_10163.index() < stce_10048) {
        throw dco::exception::create<std::runtime_error>("DCO_DEBUG: tape index of read variable bigger than current progvarcounter.");
      }
      if(stce_10154 == 0) {
        stce_10163.stce_10161();
        new (stce_10163.stce_10155.stce_10147) stce_10146(pval);
        new (stce_10163.stce_10155.stce_10148) DCO_TAPE_INT(entry::stce_10153(stce_10048));
        stce_10154 = stce_10163.stce_10162();
        DCO_LOG(logDEBUG2) << "inserting first entry, index now = " << stce_10154 << "; pval = " << pval << ", target = " << stce_10048;
        return;
      }
      if(stce_9757 && stce_10048 == *(stce_10163.stce_10155.stce_10148))
        {
          *(stce_10163.stce_10155.stce_10147) += pval;
          return;
        }
      stce_10163.stce_10161();
      new (stce_10163.stce_10155.stce_10147) stce_10146(pval);
      new (stce_10163.stce_10155.stce_10148) DCO_TAPE_INT(stce_10048);
    }
private:
    stce_9920(const stce_9920& stce_9929);
  };
protected:
  void init(tape_options const& stce_10164, DCO_TAPE_INT stce_10165 = 0) {
    size_t stce_10166;
    if(stce_10164.write_to_file())
      stce_10166 = stce_10164.blob_size_in_byte();
    else
      stce_10166 = get_allocation_size(stce_10164.blob_size_in_byte());
    size_t stce_10167 = stce_10166 / entry::size_of();
    DCO_LOG(logDEBUG1) << "BLOB SPLINT container: user requested allocation size = " << stce_10166
                       << "; size of element = " << entry::size_of()
                       << "; max elements in memory = " << stce_10167;
    if (static_cast<double>(stce_10167) / static_cast<double>(std::numeric_limits<DCO_INTEGRAL_TAPE_INT>::max()) > 1.0) {
      stce_10167 = std::numeric_limits<DCO_TAPE_INT>::max();
      DCO_LOG(logDEBUG1) << "BLOB SPLINT container: too many elements requested for currently used DCO_TAPE_INT. Setting to: " << stce_10167 * 1.e-6 << "M elements";
    }
    DCO_TAPE_INT stce_10168 = static_cast<DCO_TAPE_INT>(stce_10167);
    size_t stce_10169 = stce_10168 * entry::size_of();
    DCO_LOG(logDEBUG) << "BLOB SPLINT container: trying allocation for = " << stce_10169 / 1024 << "K ~=~ " << stce_9752(stce_10168) * 1.e-6 << "M elements";
    stce_10139 = stce_1000001(&stce_10169, stce_10164.alloc_use_huge_pages());
    if (!stce_10139) {
      std::string stce_10170 = "";
      if (stce_10164.alloc_use_huge_pages()) {
        stce_10170 = " You are using huge pages; have you prepared your system correspondingly (see documentation as well)?";
      }
      throw dco::exception::create<std::runtime_error>
        ("Could not allocate memory. Use environment variables (DCO_MEM_RATIO, DCO_MAX_ALLOCATION) and see documentation." + stce_10170);
    }
    stce_10168 = static_cast<DCO_TAPE_INT>(stce_10169 / entry::size_of());
    DCO_LOG(logINFO) << "BLOB SPLINT container: actually allocated size = " << stce_10169 / 1024 << "K ~=~ " << stce_9752(stce_10168) * 1.e-6 << "M elements";
    stce_10171 = static_cast<stce_10146*>(stce_10139);
    stce_10172 = static_cast<DCO_TAPE_INT*>(stce_10139) + stce_10169 / sizeof(DCO_TAPE_INT) - 1;
    stce_10163 = iterator(stce_10171, stce_10172, stce_10165);
    stce_10173 = iterator(stce_10171, stce_10172, stce_10165);
    stce_10174 = iterator(stce_10171 + stce_9752(stce_10168) - 1, stce_10172 - stce_9752(stce_10168) + 1, std::numeric_limits<DCO_TAPE_INT>::max());
    new (stce_10163.stce_10155.stce_10147) stce_10146(0.0);
    new (stce_10163.stce_10155.stce_10148) DCO_TAPE_INT(entry::stce_10153(0));
  }
public:
  explicit stce_10143(tape_options const& stce_10164, stce_10141 &tape) :
    stce_10145(tape),
    object_logger("BLOB SPLINT container"),
    stce_10139(0), stce_10171(0), stce_10172(0),
    stce_10163(), stce_10173(), stce_10174() {
    assert(stce_10144!=memory_model::stce_9770);
    init(stce_10164);
  }
  explicit stce_10143(iterator& stce_9815, stce_10141 &tape) :
    stce_10145(tape),
    object_logger("BLOB SPLINT container"),
    stce_10139(0), stce_10171(0), stce_10172(0),
    stce_10163(stce_9815), stce_10173(), stce_10174() {
    assert(stce_10144==memory_model::stce_9770);
  }
  ~stce_10143() {
    dco::dealloc(stce_10139);
  }
  DCO_ALWAYS_INLINE iterator& stce_9924() { return stce_10163; }
  iterator start() { return stce_10173; }
  iterator start() const { return stce_10173; }
  iterator rend() { return start(); }
  iterator current() const { return stce_10163; }
  iterator current() { return stce_10163; }
  iterator rbegin() const { return current(); }
  iterator capacity() const { return stce_10174; }
  bool empty() const { return stce_10163 == start(); }
  DCO_TAPE_INT size(const iterator& stce_9858, const iterator& stce_9859) const {
    if(stce_9858 < stce_9859)
      throw dco::exception::create<std::runtime_error>("size from < to");
    return static_cast<DCO_TAPE_INT>(stce_9858->stce_10147 - stce_9859->stce_10147);
  }
  DCO_TAPE_INT size(iterator stce_9858) const {
    return size(stce_9858, start());
  }
  DCO_TAPE_INT size() const {
    return size(current(), start());
  }
  DCO_TAPE_INT stce_10175() const {
    return size(stce_10174, start());
  }
  size_t size_in_byte() const { return size(current(), start())*entry::size_of(); }
  size_t stce_10176() const { return size(stce_10174 , start())*entry::size_of(); }
  void erase(const iterator& stce_10152) {
    if(!stce_10171)
      return;
    if(stce_10152 > stce_10174)
      throw dco::exception::create<std::runtime_error>("pos behind stack!");
    if(stce_10152 > stce_10163)
      throw dco::exception::create<std::runtime_error>("reset to position out of tape!");
    stce_10163 = stce_10152;
  }
  void stce_10177(const iterator&) { }
  void stce_10178 ( ) { }
  DCO_TAPE_INT stce_10179() const { return stce_10163.index()+1; }
private:
  stce_10143(const stce_10143& stce_9929);
  void *stce_10139;
  stce_10146 * stce_10171;
  DCO_TAPE_INT * stce_10172;
protected:
  typename stce_9816
    <iterator, stce_10144==memory_model::stce_9770>::type stce_10163;
  iterator stce_10173;
  iterator stce_10174;
};
}
}
       
       
namespace dco {
namespace internal {
template<typename stce_10141, memory_model::TYPE stce_10144>
  class stce_10143<memory_model::BLOB_TAPE, stce_10141, stce_10144> :
    public stce_10138<stce_10141>, public object_logger {
  typedef stce_10138<stce_10141> stce_10145;
  typedef stce_10141 tape_t;
  typedef typename stce_9790<stce_10141>::type stce_10146;
public:
  enum stce_10180 {
    stce_10181,
    stce_10182
  };
      struct entry {
        DCO_TAPE_INT stce_9946;
        stce_10146 stce_9922;
      public:
        static size_t size_of() {
          return sizeof(entry);
        }
        entry(DCO_TAPE_INT arg = 0, stce_10146 const& pval = 0.) :
          stce_9946(arg), stce_9922(pval) {
        }
        DCO_STRONG_INLINE DCO_TAPE_INT arg() const {
          return (stce_9946 < 0)?std::abs(stce_9752(stce_9946)) -1:stce_9752(stce_9946);
        }
        DCO_STRONG_INLINE static DCO_TAPE_INT stce_10153(const DCO_TAPE_INT& arg) {
          return -arg-1;
        }
        DCO_STRONG_INLINE bool stce_10151() const {
          return stce_9946 < 0;
        }
        DCO_STRONG_INLINE void stce_10183(const stce_10146& stce_10097) {
          stce_9922 += stce_10097;
        }
        DCO_STRONG_INLINE stce_10146 pval() const
        {
          return stce_9922;
        }
      }
      ;
      class iterator {
        entry* stce_10184;
        DCO_TAPE_INT stce_10154;
        entry const* stce_10185;
      public:
        template<class stce_10156>
          void stce_10157(stce_10156 &stce_10158) const {
          assert(sizeof(iterator) <= sizeof(stce_10156)) ;
          stce_10159(&stce_10158, this, sizeof(iterator));
        }
        template<class stce_10156>
          void stce_10160(const stce_10156 &stce_10158) {
          assert(sizeof(iterator) <= sizeof(stce_10156)) ;
          stce_10159(static_cast<void*>(this), &stce_10158, sizeof(iterator));
        }
        void stce_10186(DCO_TAPE_INT stce_9806)
        {
          stce_10154 = stce_9806;
        }
        void stce_10161() {
          stce_10184++;
          if (stce_10184 > stce_10185) {
            throw dco::exception::create<std::runtime_error>("Blob tape container overflow. Allocate bigger blob or use chunk tape container instead.");
          }
        }
        DCO_TAPE_INT index() const { return stce_10154; }
        DCO_TAPE_INT& index() { return stce_10154; }
        DCO_TAPE_INT stce_10162() { return ++stce_10154; }
        void stce_10187(entry const* end) { stce_10185 = end; }
        explicit iterator(entry* stce_10152=0, DCO_TAPE_INT index = DCO_TAPE_INT(-1)) : stce_10184(stce_10152), stce_10154(index) { }
        void operator--() {
          if(stce_10184->stce_10151())
            stce_10154--;
          stce_10184--;
        }
        bool stce_10151() const {
          return stce_10184->stce_10151();
        }
        iterator operator--(int) {
          iterator stce_10057(*this);
          operator--();
          return stce_10057;
        }
        entry* operator*() { return stce_10184; }
        entry* operator->() { return stce_10184; }
        entry* operator*() const { return stce_10184; }
        entry* operator->() const { return stce_10184; }
        bool operator==(const iterator& stce_9929) const { return stce_10184 == stce_9929.stce_10184; }
        bool operator!=(const iterator& stce_9929) const { return !operator==(stce_9929); }
        bool operator< (const iterator& stce_9929) const { return stce_10184 < stce_9929.stce_10184; }
        bool operator>=(const iterator& stce_9929) const { return !operator<(stce_9929); }
        bool operator> (const iterator& stce_9929) const { return stce_10184 > stce_9929.stce_10184; }
        bool operator<=(const iterator& stce_9929) const { return !operator>(stce_9929); }
        friend std::ostream& operator<< ( std::ostream &stce_10188, const iterator &stce_10158 ) {
          stce_10188 << "(" << stce_10158.stce_10184 << ":" << stce_10158.stce_10154 << ")";
          return stce_10188;
        }
      };
      class stce_9920 {
        iterator& stce_10163;
        DCO_TAPE_INT stce_10154;
      public:
        stce_9920(iterator& current) : stce_10163(current), stce_10154(0) {}
        iterator& current() { return stce_10163; }
        stce_9920& operator=(stce_9920 stce_9929) {
          swap(stce_9929);
          return *this;
        }
        DCO_TAPE_INT index() const { return stce_10154; }
        void swap(stce_9920& stce_9929) {
          std::swap(stce_10163, stce_9929.stce_10163);
          std::swap(stce_10154, stce_9929.stce_10154);
        }
        DCO_ALWAYS_INLINE_IF_AGRESSIVE void insert(const DCO_TAPE_INT& stce_10048, const stce_10146& pval) {
          if (IS_DEBUG && stce_10163.index() < stce_10048) {
            throw dco::exception::create<std::runtime_error>("DCO_DEBUG: tape index of read variable bigger than current progvarcounter.");
          }
          if(stce_10154 == 0) {
              stce_10163.stce_10161();
              new (*stce_10163) entry(entry::stce_10153(stce_10048),pval);
              stce_10154 = stce_10163.stce_10162();
              return;
            }
          if(stce_9757 && stce_10048 == stce_10163->arg()) {
              stce_10163->stce_10183(pval);
              return;
            }
          stce_10163.stce_10161();
          new (*stce_10163) entry(stce_10048,pval);
        }
      private:
        stce_9920(const stce_9920& stce_9929);
      };
      explicit stce_10143(stce_10143 &stce_9858,
                              stce_10180 stce_10189)
        : stce_10145(stce_9858),
          object_logger("BLOB container"),
          stce_10139(0),
          stce_10163(), stce_10173(), stce_10174(),
          stce_10190(stce_10189),
          stce_10191(stce_9858.stce_10191) {
        if (stce_10189 == stce_10181) {
            stce_10139 = stce_9858.stce_10139;
            stce_10163 = iterator(stce_10139, stce_9858.current().index());
            stce_10173 = stce_10163;
            stce_10174 = stce_9858.stce_10174;
            stce_10163.stce_10187(*stce_10174);
          } else if (stce_10189 == stce_10182) {
          DCO_TAPE_INT stce_10167 = allocate(stce_10191);
          stce_10163 = iterator(stce_10139, (stce_9858.current().index()));
          stce_10173 = iterator(stce_10139 , (stce_9858.current().index()));
          stce_10174 = iterator(stce_10139 + static_cast<DCO_INTEGRAL_TAPE_INT>(stce_10167) - 1, std::numeric_limits<DCO_INTEGRAL_TAPE_INT>::max());
          stce_10163.stce_10187(*stce_10174);
          new (*stce_10163) entry(entry::stce_10153(0));
        }
      }
  DCO_TAPE_INT allocate(const tape_options& stce_10164) {
    size_t stce_10169 = get_allocation_size(stce_10164.blob_size_in_byte());
    stce_10139 = static_cast<entry*>(stce_1000001(&stce_10169, stce_10164.alloc_use_huge_pages()));
    if (!stce_10139) {
      throw dco::exception::stce_9917<std::bad_alloc>
        ("Could not allocate memory. Use environment variables (DCO_MEM_RATIO, DCO_MAX_ALLOCATION) and see documentation.");
    }
    size_t stce_10192 = stce_10169 / entry::size_of();
    if (static_cast<double>(stce_10192) / static_cast<double>(std::numeric_limits<DCO_INTEGRAL_TAPE_INT>::max()) > 1.0) {
      stce_10192 = std::numeric_limits<DCO_TAPE_INT>::max();
    }
    DCO_TAPE_INT stce_10167 = static_cast<DCO_TAPE_INT>(stce_10192);
    DCO_LOG(logINFO) << "BLOB container; size = " << stce_10164.blob_size_in_byte() << "b ~=~ "
                     << stce_10164.blob_size_in_byte() / entry::size_of() << "elements";
    DCO_LOG(logINFO) << "BLOB container; tried allocation size = " << stce_10169
                     << "b ~=~ " << stce_10169 / entry::size_of() << "elements";
    DCO_LOG(logDEBUG) << "sizeof(entry) = " << sizeof(entry) << " == " << entry::size_of();
    return stce_10167;
  }
  void init(tape_options const& stce_10164) {
    DCO_TAPE_INT stce_10167 = allocate(stce_10164);
    stce_10163 = iterator(stce_10139, 0);
    stce_10173 = iterator(stce_10139, 0);
    stce_10174 = iterator(stce_10139 + static_cast<DCO_INTEGRAL_TAPE_INT>(stce_10167) - 1, std::numeric_limits<DCO_INTEGRAL_TAPE_INT>::max());
    stce_10163.stce_10187(*stce_10174);
    new (*stce_10163) entry(entry::stce_10153(0));
  }
  explicit stce_10143(tape_options const& stce_10164,
                          tape_t& tape)
    : stce_10145(tape),
      object_logger("BLOB container"),
      stce_10139(0),
      stce_10163(), stce_10173(), stce_10174(),
      stce_10190(stce_10182),
      stce_10191(stce_10164) {
    assert(stce_10144!=memory_model::stce_9770);
    init(stce_10164);
  }
  explicit stce_10143(iterator& stce_9815,
                          tape_t& tape)
    : stce_10145(tape),
      object_logger("BLOB container"),
      stce_10139(0),
      stce_10163(stce_9815), stce_10173(), stce_10174(),
      stce_10190(stce_10182),
      stce_10191(tape.stce_10164()) {
    assert(stce_10144==memory_model::stce_9770);
  }
      ~stce_10143() {
        if (stce_10190 == stce_10182)
          dco::dealloc(stce_10139);
      }
      DCO_ALWAYS_INLINE iterator& stce_9924() { return stce_10163; }
      entry* &stce_10193() { return stce_10139; }
      iterator current() const { return stce_10163; }
      iterator current() { return stce_10163; }
      iterator& stce_10194() { return stce_10163; }
      iterator start() { return stce_10173; }
      iterator start() const { return stce_10173; }
      iterator rend() { return start(); }
      iterator rbegin() const { return current(); }
      iterator capacity() const { return stce_10174; }
      bool empty() const { return stce_10163 == start(); }
      DCO_TAPE_INT size(const iterator& stce_9858, const iterator& stce_9859) const {
        if(stce_9858 < stce_9859)
          throw dco::exception::create<std::runtime_error>("size from < to");
        return static_cast<DCO_TAPE_INT>(*stce_9858 - *stce_9859);
      }
  tape_options stce_10164() const { return stce_10191; }
  DCO_TAPE_INT size(iterator stce_9858) const { return size(stce_9858, start()); }
  DCO_TAPE_INT size() const { return size(current(), start()); }
  DCO_TAPE_INT stce_10175() const { return size(stce_10174, start()); }
      size_t size_in_byte() const { return size(current(), start())*sizeof(entry); }
      size_t stce_10176() const { return size(stce_10174 , start())*sizeof(entry); }
      void erase(const iterator& stce_10152) {
        if(!stce_10139)
          return;
        if(IS_DEBUG && stce_10152 > stce_10174)
          throw dco::exception::create<std::runtime_error>("pos behind stack!");
        if(stce_10152 > stce_10163)
          throw dco::exception::create<std::runtime_error>("reset to position out of tape!");
        stce_10163 = stce_10152;
        stce_10163.stce_10187(*stce_10174);
      }
        void erase() {
        if(!stce_10139)
          return;
        stce_10163 = stce_10173;
        stce_10163.stce_10187(*stce_10174);
      }
      void stce_10177(const iterator&) { }
      void stce_10178 ( ) { }
      DCO_TAPE_INT stce_10179() { return stce_10163.index()+1; }
private:
  stce_10143(const stce_10143& stce_9929);
  entry * stce_10139;
  typename stce_9816
    <iterator, stce_10144==memory_model::stce_9770>::type stce_10163;
  iterator stce_10173;
  iterator stce_10174;
  stce_10180 stce_10190;
  tape_options stce_10191;
};
}
}
namespace dco {
namespace internal {
template<typename stce_10141, memory_model::TYPE stce_10144>
  class stce_10143<memory_model::CHUNK_TAPE, stce_10141, stce_10144>
  : public stce_10138<stce_10141>, object_logger
{
  typedef stce_10138<stce_10141> stce_10145;
private:
  typedef typename stce_9790<stce_10141>::type stce_10146;
protected:
  typedef stce_10143<memory_model::BLOB_TAPE, stce_10141> stce_10195;
  static void stce_10196(const std::string& stce_10197, char* data, long size)
  {
    DCO_LOG(logDEBUG1) << "save_data_to_file: bytes to write = " << size << ", fname = " << stce_10197;
    std::ofstream out(stce_10197.c_str(), std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
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
  static DCO_INTEGRAL_TAPE_INT stce_10198(const std::string& stce_10197, char* data)
  {
    struct stat stce_10199;
    if (stat(stce_10197.c_str(), &stce_10199) == -1)
      throw dco::exception::create<std::runtime_error>(
          "could not check file size, perhaps index corruption?");
    DCO_INTEGRAL_TAPE_INT size = static_cast<DCO_INTEGRAL_TAPE_INT>(stce_10199.st_size);
    std::ifstream in(stce_10197.c_str(), std::ios_base::in | std::ios_base::binary);
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
  class stce_10200 : public stce_10195 {
    DCO_TAPE_INT stce_10201;
    tape_options stce_10191;
    bool stce_10202;
    stce_10200 *stce_10203, *stce_10204;
    DCO_TAPE_INT stce_10205;
public:
    typedef typename stce_10195::entry entry;
    typedef typename stce_10195::iterator iterator;
    DCO_TAPE_INT stce_10206() const { return stce_10201; }
    static std::string stce_10207(const stce_10200& stce_10195) {
      std::stringstream stce_10197;
      stce_10197 << stce_10195.stce_10208().filename() << '.' << stce_10195.stce_10206();
      return stce_10197.str();
    }
    stce_10200(stce_10200& stce_9858, typename stce_10195::stce_10180 stce_10209)
      : stce_10195(stce_9858, stce_10209), stce_10201(stce_9858.stce_10206()+1), stce_10191(stce_9858.stce_10208()), stce_10202(true),
        stce_10203(0), stce_10204(0), stce_10205(-1) { }
    stce_10200(tape_options stce_10164, stce_10141 &tape)
      : stce_10195(stce_10164, tape), stce_10201(0), stce_10191(stce_10164), stce_10202(true),
        stce_10203(0), stce_10204(0), stce_10205(-1) { }
    ~stce_10200() {
      if(stce_10208().write_to_file()) std::remove(stce_10207(*this).c_str());
    }
    DCO_TAPE_INT size(const iterator& stce_9858, const iterator& stce_9859) const {
      return stce_10195::size(stce_9858, stce_9859);
    }
    DCO_TAPE_INT size(iterator stce_9858) const {
      if(!*this->start())
        throw dco::exception::create<std::runtime_error>("given iterator does not point "
            "into given tape!");
      return stce_10195::size(stce_9858, this->start());
    }
    DCO_TAPE_INT size() const {
        return stce_10195::size(this->current(), this->start());
    }
    DCO_TAPE_INT stce_10175() const {
      if(!stce_10210())
        return 0;
      return stce_10195::stce_10175();
    }
    size_t size_in_byte() const {
      return stce_10195::size(this->current(), this->start()) * stce_10195::entry::size_of();
    }
    size_t stce_10176() const {
      if(!stce_10210()) { return 0; }
      return stce_10195::stce_10175() * stce_10195::entry::size_of();
    }
    tape_options stce_10208() const {
      return stce_10191;
    }
    bool stce_10210() const { return stce_10202; }
    void stce_10211(bool stce_10212) { stce_10202 = stce_10212; }
    void stce_10213() { stce_10205 = size(); }
    void stce_10214() { stce_10205 = -1; }
    bool stce_10215() { return stce_10205 != size(); }
    stce_10200* &stce_10216() { return stce_10203; }
    stce_10200* &stce_10217() { return stce_10204; }
  };
  typedef typename stce_10195::entry entry;
  typedef stce_10200* stce_10218;
  typedef typename stce_10195::iterator stce_10219;
  static void stce_10220(stce_10200& stce_10195) {
    if(stce_10195.stce_10208().write_to_file()) {
      assert(!stce_10195.stce_10210());
      DCO_LOG(logDEBUG1) << "wake_up: waking up chunk " << stce_10195.stce_10206() << " | " << stce_10195.stce_10193() << " from sleep.";
      DCO_INTEGRAL_TAPE_INT size = stce_10198(stce_10200::stce_10207(stce_10195),
                                                       reinterpret_cast<char*>(stce_10195.stce_10193()));
      (void)size;
      stce_10195.stce_10211(true);
    } else {
      DCO_LOG(logDEBUG1) << "wake_up: waking already awake chunk " << stce_10195.stce_10206() << " | " << stce_10195.stce_10193();
    }
  }
  static void sleep(stce_10200& stce_10195) {
    DCO_LOG(logDEBUG1) << "sleep: chunk " << stce_10195.stce_10206() << " | " << stce_10195.stce_10193() << " dirty = " << stce_10195.stce_10215();
    if(stce_10195.stce_10208().write_to_file() && stce_10195.stce_10215()) {
      assert(stce_10195.stce_10210());
      DCO_LOG(logDEBUG1) << "sleep: writing chunk " << stce_10195.stce_10206() << " | " << stce_10195.stce_10193();
      const long stce_10221 = stce_9752(stce_10195.size()+1) * sizeof(entry);
      DCO_LOG(logDEBUG1) << "sleep: chunk.size() = " << stce_10195.size() << ", bytes to write = " << stce_10221 << ", ptr diff = " << stce_10195.current() << ", " << stce_10195.start() << " = " << *stce_10195.current() - *stce_10195.start() << " " << sizeof(entry);
      entry* data = stce_10195.stce_10193();
      stce_10196(stce_10200::stce_10207(stce_10195),
                        reinterpret_cast<char*>(data), stce_10221);
      stce_10195.stce_10213();
    } else {
      if (!stce_10195.stce_10215()) {
        DCO_LOG(logDEBUG1) << "sleep: not writing data for chunk: " << stce_10195.stce_10206();
      }
    }
  }
  static void stce_10222(stce_10200& stce_10195, stce_10200& stce_10223) {
    if (!stce_10195.stce_10210()) {
      sleep(stce_10223);
      stce_10223.stce_10211(false);
      stce_10220(stce_10195);
    }
  }
  void stce_10222(stce_10200& stce_10195) {
    if (!stce_10195.stce_10210()) {
      DCO_LOG(logDEBUG1) << "activate: chunk " << stce_10195.stce_10206() << " not awake yet";
      stce_10200* stce_10223 = stce_10163.stce_10224();
      while (stce_10223 && !stce_10223->stce_10210()) { stce_10223 = stce_10223->stce_10216(); }
      if (stce_10223) {
        DCO_LOG(logDEBUG1) << "activate: currently awake " << stce_10223->stce_10206();
        stce_10222(stce_10195, *stce_10223);
      } else {
        stce_10220(stce_10195);
      }
    }
  }
  class iterator {
private:
    stce_10218 stce_10225;
    stce_10219 stce_10226;
public:
    template<class stce_10156>
      void stce_10160(const stce_10156 &stce_10158) {
      assert(sizeof(iterator) <= sizeof(stce_10156)) ;
      stce_10159(static_cast<void*>(this), &stce_10158, sizeof(iterator));
    }
    template<class stce_10156>
      void stce_10157(stce_10156 &stce_10158) const {
      assert(sizeof(iterator) <= sizeof(stce_10156)) ;
      stce_10159(&stce_10158, this, sizeof(iterator));
    }
    stce_10218 stce_10224() const {
      return stce_10225;
    }
    stce_10219 stce_10227() const {
      return stce_10226;
    }
    explicit iterator(const stce_10219& stce_10228, stce_10218 stce_10224) :
        stce_10225(stce_10224), stce_10226(stce_10228) {
    }
    explicit iterator(stce_10218 stce_10224) :
        stce_10225(stce_10224), stce_10226(stce_10224->current()) {
    }
    iterator(const iterator& stce_9929) :
        stce_10225(stce_9929.stce_10225), stce_10226(stce_9929.stce_10226) {
    }
    iterator() : stce_10225(), stce_10226() {}
    bool stce_10151() const {
      return stce_10226.stce_10151();
    }
    void operator--() {
      DCO_LOG(logDEBUG4) << stce_10226 << " =?= " << stce_10225->start();
      --stce_10226;
      if (stce_10226 == stce_10225->start()) {
          if(stce_10225->stce_10206()==0) return;
        stce_10225 = stce_10225->stce_10216();
        stce_10222(*stce_10225, *(stce_10225->stce_10217()));
        stce_10226 = stce_10225->current();
      }
    }
    DCO_TAPE_INT index() const { return stce_10226.index(); }
    DCO_TAPE_INT& index() { return stce_10226.index(); }
    DCO_TAPE_INT stce_10162() { stce_10225->stce_10194().stce_10162(); return stce_10226.stce_10162(); }
    stce_10218& stce_10195() { return stce_10225; }
    stce_10219& stce_10229() { return stce_10226; }
    iterator operator--(int) {
      iterator stce_10057(*this);
      operator--();
      return stce_10057;
    }
    entry* operator* () { return *stce_10226; }
    entry* operator* () const { return *stce_10226; }
    entry* operator->() { return *stce_10226; }
    entry* operator->() const { return *stce_10226; }
    bool operator==(const iterator& stce_9929) const {
      return stce_10226 == stce_9929.stce_10226 && stce_10225 == stce_9929.stce_10225;
    }
    bool operator!=(const iterator& stce_9929) const {
      return !operator==(stce_9929);
    }
    bool operator<(const iterator& stce_9929) const {
      return
          (stce_10225 == stce_9929.stce_10225) ?
              stce_10226 < stce_9929.stce_10226 : stce_10226.index() < stce_9929.stce_10226.index();
    }
    bool operator>=(const iterator& stce_9929) const {
      return !operator<(stce_9929);
    }
    bool operator>(const iterator& stce_9929) const {
      return
          (stce_10225 == stce_9929.stce_10225) ?
              stce_10226 > stce_9929.stce_10226 : stce_10226.index() > stce_9929.stce_10226.index();
    }
    bool operator<=(const iterator& stce_9929) const {
      return !operator>(stce_9929);
    }
  };
  class stce_9920
  {
    stce_10218& stce_10230;
    stce_10219& stce_10231;
    stce_10219 stce_10226;
    stce_10219 stce_10232;
    DCO_TAPE_INT stce_10154;
public:
    stce_9920(iterator& stce_9815) :
      stce_10230(stce_9815.stce_10195()),
      stce_10231(stce_9815.stce_10229()),
      stce_10226(stce_10230->current()),
      stce_10232(stce_10230->capacity()),
      stce_10154(0) {
    }
    stce_10219& current() { return stce_10226; }
    ~stce_9920() {
      if(static_cast<DCO_INTEGRAL_TAPE_INT>(stce_10154)) {
        stce_10231 = stce_10230->stce_10194() = stce_10226;
      }
    }
    DCO_TAPE_INT index() const { return stce_10154; }
    DCO_TAPE_INT& index() { return stce_10154; }
    DCO_ALWAYS_INLINE_IF_AGRESSIVE void insert(const DCO_TAPE_INT& stce_10048, const stce_10146& pval) {
      if (IS_DEBUG && stce_10226.index() < stce_10048) {
        throw dco::exception::create<std::runtime_error>("DCO_DEBUG: tape index of read variable bigger than current progvarcounter.");
      }
      if (stce_10154 == 0) {
        stce_10226.stce_10161();
        new (*stce_10226) entry(entry::stce_10153(stce_10048), pval);
        stce_10154 = stce_10226.stce_10162();
        DCO_LOG(logDEBUG2) << "inserting first entry, index now = " << stce_10154 << "; pval = " << pval << ", target = " << stce_10048 << " in chunk " << stce_10230 << " into memory at = " << stce_10226;
        return;
      }
      if(stce_9757 && stce_10226->arg() == stce_10048){
        DCO_LOG(logDEBUG3) << "merging a parallel edge";
        stce_10226->stce_10183(pval);
        return;
      }
      if (*stce_10226 >= *stce_10232) {
        DCO_LOG(logDEBUG1) << "insert: chunk full, updating current and getting next";
        stce_10231 = stce_10230->stce_10194() = stce_10226;
        stce_10233(stce_10230);
        stce_10226 = stce_10230->current();
        stce_10232 = stce_10230->capacity();
      }
      stce_10226.stce_10161();
      new (*stce_10226) entry(stce_10048, pval);
      DCO_LOG(logDEBUG2) << "inserting next entry, index now = " << stce_10154 << "; pval = " << pval << ", target = " << stce_10048 << " in chunk " << stce_10230 << " into memory at = " << stce_10226;
    }
  private:
    stce_9920(const stce_9920& stce_9929);
  };
public:
  void init(tape_options const& stce_10164) {
    tape_options stce_10234 = stce_10164;
    stce_10234.blob_size_in_byte() = stce_10164.chunk_size_in_byte();
    stce_10235 = new stce_10200(stce_10234, stce_10145::stce_10140());
    stce_10163 = iterator(stce_10235);
  }
  explicit stce_10143(tape_options const& stce_10164, stce_10141& tape) :
    stce_10145(tape),
    object_logger("CHUNK container (w/ options)"),
    stce_10235(0),
    stce_10163() {
    assert(stce_10144!=memory_model::stce_9770);
    init(stce_10164);
  }
  explicit stce_10143(iterator &stce_9815, stce_10141& tape) :
    stce_10145(tape),
    object_logger("CHUNK container (w/ options)"),
    stce_10235(0), stce_10163(stce_9815) {
    assert(stce_10144==memory_model::stce_9770);
  }
  ~stce_10143() {
    stce_10200* stce_10158 = stce_10235;
    while (stce_10158) {
      stce_10200* stce_10236 = stce_10158;
      stce_10158 = stce_10158->stce_10217();
      delete stce_10236;
    }
  }
  static void stce_10233(stce_10218& stce_9815)
  {
    DCO_LOG(logDEBUG1) << "current chunk: id = " << stce_9815->stce_10206() << ", index = " << stce_9815->current().index()
                       << ", data pointer = " << stce_9815->stce_10193();
    delete stce_9815->stce_10217();
    if (stce_9815->stce_10208().write_to_file()) {
      sleep(*stce_9815);
      stce_9815->stce_10211(false);
      stce_9815->stce_10217() = new stce_10200(*stce_9815, stce_10200::stce_10181);
    } else {
      stce_9815->stce_10217() = new stce_10200(*stce_9815, stce_10200::stce_10182);
    }
    stce_9815->stce_10217()->stce_10216() = stce_9815;
    stce_9815 = stce_9815->stce_10217();
    DCO_LOG(logDEBUG1) << "next chunk: id = " << stce_9815->stce_10206() << ", index = " << stce_9815->current().index()
                       << ", data pointer = " << stce_9815->stce_10193();
  }
  DCO_ALWAYS_INLINE iterator& stce_9924()
  {
    if (stce_10163.stce_10195()->current() >= stce_10163.stce_10195()->capacity()) {
      DCO_LOG(logDEBUG1) << "insert_ref: chunk full, getting next";
      stce_10233(stce_10163.stce_10195());
    }
    return stce_10163;
  }
  iterator start() { return iterator(stce_10235->start(), stce_10235); }
  iterator start() const { return iterator(stce_10235->start(), stce_10235); }
  iterator rend() { return start(); }
  iterator current() const { return stce_10163; }
  iterator rbegin() const { return current(); }
  bool empty() const { return current() == start(); }
  DCO_TAPE_INT size(const iterator& stce_9858, const iterator& stce_9859) const {
    const stce_10200& stce_10237 = *(stce_9858.stce_10224());
    if(stce_9858.stce_10224() == stce_9859.stce_10224())
      return stce_10237.size(stce_9858.stce_10227(), stce_9859.stce_10227());
    DCO_TAPE_INT size = 0;
    stce_10218 stce_10238 = stce_9858.stce_10224();
    size += stce_10238->size(stce_9858.stce_10227());
    DCO_LOG(logDEBUG2) << "size of last chunk = " << size;
    for (stce_10238 = stce_10238->stce_10216(); stce_10238 != stce_9859.stce_10224(); stce_10238 = stce_10238->stce_10216()) {
      if(stce_10238) {
        size += stce_10238->size();
        DCO_LOG(logDEBUG2) << "adding " << stce_10238->size() << " => size = " << size;
      }
    }
    const stce_10200& stce_10239 = *(stce_9859.stce_10224());
    size += stce_10239.size(stce_10239.current(), stce_9859.stce_10227());
    DCO_LOG(logDEBUG2) << "very first bit: adding " << stce_10239.size(stce_10239.current(), stce_9859.stce_10227()) << " => size = " << size;
    return size;
  }
  DCO_TAPE_INT size(iterator stce_9858) const { return size(stce_9858, start()); }
  DCO_TAPE_INT size() const { return size(current(), start()); }
  mem_long_t size_in_byte() const {
    return size(current(), start())*entry::size_of();
  }
  DCO_TAPE_INT stce_10175() const {
    DCO_TAPE_INT size = 0;
    for(stce_10218 stce_10240 = stce_10235; stce_10240 != NULL; stce_10240 = stce_10240->stce_10217())
      if(stce_10240)
        size += stce_10240->stce_10175();
    return size;
  }
  mem_long_t stce_10176() const {
    mem_long_t size = 0;
    for(stce_10218 stce_10240 = stce_10235; stce_10240 != NULL; stce_10240 = stce_10240->stce_10217())
      if(stce_10240) size += stce_10240->stce_10176();
    return size;
  }
  void erase(const iterator& stce_10152) {
    if (stce_10152 > current())
      throw dco::exception::create<std::invalid_argument>("erase pos out of tape!");
    stce_10218 stce_10241(stce_10235);
    while (stce_10241 != stce_10152.stce_10224())
      stce_10241 = stce_10241->stce_10217();
    stce_10163.stce_10195() = stce_10241;
    stce_10222(*stce_10163.stce_10195());
    stce_10163.stce_10195()->stce_10214();
    stce_10218 stce_10158 = stce_10241->stce_10217();
    stce_10241->stce_10217() = 0;
    while (stce_10158) {
      stce_10218 stce_10236 = stce_10158;
      stce_10158 = stce_10158->stce_10217();
      delete stce_10236;
    }
    stce_10163.stce_10195()->erase(stce_10152.stce_10227());
    stce_10163.stce_10229() = stce_10163.stce_10195()->current();
  }
  void stce_10177(const iterator& stce_9815) {
    DCO_LOG(logDEBUG1) << "chunk tape: prepare_interpretation; activate: " << stce_9815.stce_10224()->stce_10206() << ", is awake = " << stce_9815.stce_10224()->stce_10210();
    stce_10222(*stce_9815.stce_10224());
  }
  void stce_10178() {
    DCO_LOG(logDEBUG1) << "chunk tape: prepare_recording";
    stce_10222(*stce_10163.stce_10195());
  }
  DCO_TAPE_INT stce_10179() { return current().index()+1; }
private:
  stce_10200* stce_10235;
  typename stce_9816
    <iterator, stce_10144==memory_model::stce_9770>::type stce_10163;
};
}
}
       
namespace dco {
namespace internal {
template<class stce_10141>
  class stce_10143<memory_model::stce_9769, stce_10141>
  : public stce_10143<memory_model::CHUNK_TAPE, stce_10141>
{
  typedef stce_10141 tape_t;
  typedef typename stce_9790<stce_10141>::type stce_10146;
public:
  typedef stce_10143<memory_model::CHUNK_TAPE, stce_10141> stce_10242;
  typedef typename stce_10242::entry entry;
  typedef typename stce_10242::stce_9920 stce_9920;
  typedef typename stce_10242::iterator iterator;
  bool stce_10243, stce_10006, stce_10244, stce_10245;
  double stce_10246, stce_10247, stce_10000;
  std::vector<DCO_TAPE_INT> stce_10009;
  std::vector<stce_10146> stce_10010;
  std::vector<const void*> stce_10248;
  DCO_TAPE_INT stce_10249, stce_10250;
  struct stce_10251 {
    static std::string id() { return " gtas "; }
  };
  typedef logging<stce_10251> logger;
  explicit stce_10143(tape_options stce_10164, stce_10141 &tape) :
    stce_10242(stce_10164, tape),
    stce_10243(false),
    stce_10006(false),
    stce_10244(false),
    stce_10245(false),
    stce_10246(1e-6),
    stce_10247(1e-6),
    stce_10000(0.0),
    stce_10249(0),
    stce_10250(-1) {}
  void stce_10252(log_level_enum level) {
    if (stce_10245 && level <= dco::logWARNING) {
      throw dco::exception::create<std::runtime_error>("Warning occured.");
    }
  }
  void stce_10005() {
    if (stce_10250 == stce_10249) {
      stce_10249++;
      throw dco::exception::create<std::runtime_error>("Event counter reached.");
    }
    stce_10249++;
  };
  void stce_10004(stce_10146 stce_9933, double stce_10003)
  {
    if(stce_10244) {
      double stce_10253 = ::fabs(stce_9933 - stce_10003);
      log_level_enum stce_10254 = dco::logWARNING;
      if (stce_10253 < stce_10246) {stce_10254 = dco::logINFO;}
      logger(stce_10254).get_stream() << stce_10249 << ": fd vs. tangent check";
      logger(stce_10254).get_stream() << " delta = " << stce_10253 << " tangent = " << stce_9933
          << " finite differences = " << stce_10003;
      stce_10252(stce_10254);
    }
  }
  void stce_10011(double stce_10008, double stce_10255, double stce_10256)
  {
    const double stce_10257 = ::fabs(stce_10008-stce_10255);
    const double stce_10258 = ::fabs(stce_10008-stce_10256);
    log_level_enum stce_10254 = dco::logWARNING;
    if (stce_10257 < stce_10247) {stce_10254 = dco::logINFO;}
    logger(stce_10254).get_stream() << stce_10249 << ": Tangent/adjoint identity check:";
    logger(stce_10254).get_stream() << " tangent vs. adjoint = " << stce_10257;
    logger(stce_10254).get_stream() << "  <x_(1),x^(1)> = " << stce_10008;
    logger(stce_10254).get_stream() << "  <y_(1),y^(1)> = " << stce_10255;
    stce_10252(stce_10254);
    stce_10254 = dco::logWARNING;
    if (stce_10258 < stce_10246) {stce_10254 = dco::logINFO;}
    logger(stce_10254).get_stream() << " finite differences vs. adjoint = " << stce_10258;
    logger(stce_10254).get_stream() << "  <y_(1),y_{FD}^(1)> = " << stce_10256;
    stce_10252(stce_10254);
  }
  void stce_10259(const iterator &stce_9859) {
    assert(stce_10242::start() == stce_9859);
    stce_10242::stce_10259(stce_9859);
    stce_10009.clear();
    stce_10010.clear();
    stce_10248.clear();
    stce_10000=0;
  }
  void stce_9977(std::string stce_10260, const bool stce_9836, const bool stce_9974, const bool stce_9976)
  {
    if (!stce_10243) {
      stce_10005();
      return;
    }
    log_level_enum stce_10254 = dco::logWARNING;
    if (stce_9836 == stce_9974 && stce_9836 == stce_9976) {
      stce_10254 = dco::logINFO;
    }
    logger(stce_10254).get_stream() << stce_10249 << ": Branch check:" << stce_10260
        << std::boolalpha << " at current value returns " << stce_9836 << ", at minus h " << stce_9974
        << ", and at plus h " << stce_9976;
    stce_10252(stce_10254);
    stce_10005();
  }
  template <class data_t>
  void stce_10001(const data_t* data, const double &stce_9933) {
    if(stce_10000 == 0) {
        stce_10248.push_back(data);
        stce_10009.push_back(data->_tape_index());
        stce_10010.push_back(stce_9933);
    }
  }
  template <typename stce_10261>
  void stce_10262() {
      double stce_10263=0;
      for(size_t stce_9827=0;stce_9827<stce_10248.size();++stce_9827) {
        const stce_10261 &stce_9815 = *static_cast<const stce_10261*>(stce_10248[stce_9827]);
        stce_10263 = std::max(stce_10263, std::fabs(stce_9815.stce_9998[0]));
      }
      const double stce_10264 = 1e-14;
      stce_10000 = std::sqrt(stce_10264) * stce_10263;
      for(size_t stce_9827=0;stce_9827<stce_10248.size();++stce_9827) {
        const stce_10261 &stce_9815 = *static_cast<const stce_10261*>(stce_10248[stce_9827]);
        stce_9815.stce_9998[0] -= stce_10000 * stce_9815.stce_9932;
        stce_9815.stce_9998[1] += stce_10000 * stce_9815.stce_9932;
      }
  }
};
}
}
       
namespace dco {
namespace internal {
template<class stce_10141, memory_model::TYPE stce_10142>
  class stce_10143<memory_model::stce_9770, stce_10141, stce_10142>
  : public stce_10143<stce_10142, stce_10141, memory_model::stce_9770> {
  typedef typename stce_9790<stce_10141>::type stce_10146;
public:
  typedef stce_10143<stce_10142, stce_10141, memory_model::stce_9770> stce_10242;
  explicit stce_10143(tape_options stce_10164, stce_10141 &tape) :
    stce_10242(stce_10163, tape), stce_10163() {
    stce_10242::init(stce_10164);
  }
  struct iterator : public stce_10242::iterator {
    iterator() : stce_10265(1), stce_10266(0) {}
    iterator(const typename stce_10242::iterator& stce_10267) : stce_10242::iterator(stce_10267), stce_10265(-1), stce_10266(-1) {}
    DCO_TAPE_INT stce_10265;
    DCO_TAPE_INT stce_10266;
    const DCO_TAPE_INT& max_dist() const { return stce_10265; }
    void stce_10268(const DCO_TAPE_INT &stce_10269) {
      stce_10265 = stce_10269+1;
      DCO_LOG(logDEBUG3) << " updating max dist: " << stce_10265;
    }
    const DCO_TAPE_INT &stce_10270() const { return stce_10266; }
    DCO_TAPE_INT stce_10162() {
      stce_10266++;
      return stce_10242::iterator::stce_10162();
    }
  };
  class stce_9920 : public stce_10242::stce_9920 {
    iterator &stce_10163;
public:
    stce_9920(iterator& current) : stce_10242::stce_9920(current), stce_10163(current) { }
    iterator& current() { return stce_10163; }
    DCO_ALWAYS_INLINE void insert(const DCO_TAPE_INT& stce_10048, const stce_10146& pval) {
      DCO_LOG(logDEBUG4) << "target = " << stce_10048 << ", parametercounter = " << current().stce_10270() << ", max dist = " << current().max_dist();
      DCO_TAPE_INT stce_10271 = current().index();
      stce_10242::stce_9920::insert(stce_10048, pval);
      if(stce_10048 > current().stce_10270()) {
        DCO_TAPE_INT stce_10272 = stce_10271 - stce_10048+1;
        DCO_LOG(logDEBUG4) << "curdist = " << stce_10272;
        if(stce_10272 - current().max_dist() > 0) {
          current().stce_10268(stce_10272);
        }
      }
    }
private:
    stce_9920(const stce_9920& stce_9929);
  };
  DCO_ALWAYS_INLINE iterator& stce_9924() {
    return static_cast<iterator&>(stce_10242::stce_9924());
  }
  iterator current() const { return stce_10163; }
  iterator current() { return stce_10163; }
  iterator& stce_10194() { return stce_10163; }
  iterator rbegin() const { return stce_10163; }
  DCO_TAPE_INT stce_10179() {
    using ::log;
    using ::pow;
    DCO_TAPE_INT max_dist = stce_10163.max_dist();
    bool stce_10273 = (max_dist & (max_dist - 1)) == 0;
    if (HAS_BITWISE_MODULO && !stce_10273) {
      DCO_TAPE_INT stce_10089 = static_cast<DCO_TAPE_INT>(floor(log(static_cast<double>(max_dist))/log(2.0)))+1;
      stce_10163.stce_10268( (1 << stce_9752(stce_10089)) - 1 );
      DCO_LOG(logDEBUG2) << "rounding adjoint vector size to " << stce_10163.max_dist();
    }
    return stce_10163.stce_10270() + stce_10163.max_dist() + 1;
  }
private:
  iterator stce_10163;
};
}
namespace helper {
template<typename stce_9788, memory_model::TYPE stce_9789, typename external_adjoint_object_t>
  struct stce_10072<internal::tape<memory_model::stce_9770, stce_9788, stce_9789>, external_adjoint_object_t> {
  typedef internal::tape<memory_model::stce_9770, stce_9788, stce_9789> tape_t;
  DCO_TAPE_INT stce_10274;
  DCO_TAPE_INT stce_10275;
  tape_t* &stce_9948;
  stce_10072(external_adjoint_object_t* stce_10276) :
    stce_10274(std::numeric_limits<DCO_TAPE_INT>::max()),
    stce_10275(0),
    stce_9948(stce_10276->tape()) {}
  template <typename stce_9782> void register_input(const stce_9782& stce_9839) {
    DCO_TAPE_INT stce_9806 = stce_9839._tape_index();
    if(stce_9806 > stce_9948->current().stce_10270()) {
      stce_10274 = std::min(stce_10274, stce_9806);
    }
  }
  template <typename stce_9782> void register_output(const stce_9782& stce_10114) {
    stce_10275 = std::max(stce_10275, stce_10114._tape_index());
    if(stce_10274 != std::numeric_limits<DCO_TAPE_INT>::max()) {
        DCO_TAPE_INT max_dist = stce_10275 - stce_10274;
        if(max_dist > stce_9948->stce_10194().max_dist()) {
          stce_9948->stce_10194().stce_10268(max_dist);
        }
    }
  }
};
}
namespace internal {
template<typename stce_9791, typename stce_10141, memory_model::TYPE stce_10142>
  class stce_9793<stce_9791, internal::stce_10143<memory_model::stce_9770, stce_10141, stce_10142> >
  : public stce_9793<stce_9791, internal::stce_10143<stce_10142, stce_10141, memory_model::stce_9770> >
{
  typedef internal::stce_10143<memory_model::stce_9770, stce_10141, stce_10142> stce_10277;
  typedef typename stce_10277::stce_10242 stce_10242;
  typedef stce_9793<stce_9791, stce_10242> stce_10278;
  typedef typename stce_10278::stce_9764 stce_9764;
public:
  stce_10277 &stce_9948;
  stce_9793(DCO_TAPE_INT size, stce_10277 &tape, tape_options stce_10164) :
    stce_10278(size, static_cast<stce_10242&>(tape), stce_10164), stce_9948(tape),stce_10279(true), stce_10280(0) { }
  void init() {
    stce_10281 = stce_9948.current().stce_10270();
    if (HAS_BITWISE_MODULO) {
      stce_10282 = stce_9948.current().max_dist()-1;
    } else {
      stce_10282 = stce_9948.current().max_dist();
    }
    stce_10283 = this->stce_10284 + stce_9752(stce_10281) + 1;
  }
  bool stce_10285(DCO_TAPE_INT stce_10286) const {
    if(stce_10286> stce_9948.current().stce_10270()) {
      stce_10286 -= stce_9948.current().stce_10270();
      stce_10286 %= stce_9948.current().max_dist();
      stce_10286 += stce_9948.current().stce_10270() + 1;
    }
    return static_cast<DCO_TAPE_INT>(this->stce_10287) > stce_10286;
  }
  DCO_TAPE_INT stce_10288(DCO_TAPE_INT stce_9806) const {
    if(stce_9806 > stce_9948.current().stce_10270()) {
      stce_9806 -= stce_9948.current().stce_10270();
      stce_9806 %= stce_9948.current().max_dist();
      stce_9806 += stce_9948.current().stce_10270() + 1;
    }
    return stce_9806;
  }
  typename stce_10278::stce_9764& operator[](DCO_TAPE_INT stce_9806) {
    DCO_LOG(logDEBUG4) << "getting " << stce_9806;
    stce_9806 = stce_10288(stce_9806);
    DCO_LOG(logDEBUG4) << "matched to " << stce_9806;
    return stce_10278::operator[](stce_9806);
  }
  typename stce_10278::stce_9764 operator[](DCO_TAPE_INT stce_9806) const {
    DCO_LOG(logDEBUG4) << "getting " << stce_9806;
    stce_9806 = stce_10288(stce_9806);
    DCO_LOG(logDEBUG4) << "matched to " << stce_9806;
    return stce_10278::operator[](stce_9806);
  }
  stce_9764 at(DCO_TAPE_INT stce_9806) const { return stce_10285(stce_9806)?(*this)[stce_9806]:0; }
  stce_9764 &at(const DCO_TAPE_INT stce_9806) {
    if(!stce_10285(stce_9806))
      this->resize(stce_9948.stce_10179());
    return (*this)[stce_9806];
  }
  void stce_10289(DCO_TAPE_INT stce_9806) {
    if(!stce_10285(stce_9806)) {
      DCO_LOG(logDEBUG2) << "Resizing internal adjoint vector to " << stce_9806;
      this->resize(stce_9948.stce_10179());
    }
  }
  template <typename iterator_t>
  void stce_10290(iterator_t& stce_9815) {
    DCO_TAPE_INT stce_9859 = stce_9815->arg();
    stce_9764 pval = stce_9815->pval();
    DCO_TAPE_INT stce_9858 = stce_9815.index();
    if (stce_10279) {
      DCO_LOG(logDEBUG3) << "setting " << stce_9858 << " to 0 ";
      stce_9764& stce_10291 = HAS_BITWISE_MODULO ?
                           stce_10283[ (stce_9858 - stce_10281) & stce_10282 ] :
                           stce_10283[ (stce_9858 - stce_10281) % stce_10282 ];
      this->stce_10280 = stce_10291;
      stce_10291 = 0;
    }
    if(stce_9859 > stce_10281) {
      stce_9764& stce_10291 = HAS_BITWISE_MODULO ?
                           stce_10283[ (stce_9859 - stce_10281) & stce_10282 ] :
                           stce_10283[ (stce_9859 - stce_10281) % stce_10282 ];
      stce_10291 += pval*this->stce_10280;
    } else {
      (*this)[stce_9859] += pval*this->stce_10280;
    }
    stce_10279 = stce_9815.stce_10151();
  }
private:
  stce_9764* stce_10283;
  DCO_TAPE_INT stce_10281;
  DCO_TAPE_INT stce_10282;
  bool stce_10279;
  stce_9764 stce_10280;
};
}
}
namespace dco {
namespace internal {
template<class stce_9761, typename stce_10137>
  class stce_9793 : public stce_10138<stce_10137>
{
  typedef stce_10138<stce_10137> stce_10145;
public:
  typedef stce_9761 stce_9764;
protected:
  stce_9764 *stce_10284;
  DCO_TAPE_INT stce_10287;
  DCO_TAPE_INT stce_10292;
  DCO_TAPE_INT stce_10293;
  int stce_10294;
  tape_options stce_10191;
public:
  typedef stce_9764 value_t;
  void init() {}
  stce_9793(DCO_TAPE_INT size, stce_10137 &tape, tape_options stce_10164) :
    stce_10145(tape), stce_10284(0), stce_10287(0), stce_10292(0),
    stce_10293(size), stce_10294(0), stce_10191(stce_10164) { }
  ~stce_9793() {
    if (stce_10284) dco::dealloc(stce_10284);
  }
  DCO_TAPE_INT size() const { return stce_10293; }
  bool stce_10285(DCO_TAPE_INT stce_10286) const {
    return stce_10287 > stce_10286;
  }
  mem_long_t stce_10295() const {
    return static_cast<mem_long_t>(stce_10292) * sizeof(stce_9764);
  }
  stce_9764& operator[](DCO_TAPE_INT stce_9806) { return stce_10284[stce_9806]; }
  stce_9764 operator[](DCO_TAPE_INT stce_9806) const { return stce_10284[stce_9806]; }
  void stce_10290(typename stce_10137::iterator& stce_9815) {
    (*this)[stce_9815->arg()] += stce_9815->pval()*(*this)[stce_9815.index()];
  }
  void resize(DCO_TAPE_INT stce_10296) {
    DCO_LOG(logDEBUG2) << "Resizing internal adjoint vector to " << stce_10296;
    if (stce_10296 < stce_10293) {
      stce_10297(stce_10293, stce_10296-1);
    }
    stce_10293 = stce_10296;
  }
  void stce_10298(DCO_TAPE_INT stce_10096) {
    if (++stce_10294 % 100 == 0) {
      DCO_LOG(logWARNING) << "dco_adjoint_vector: resizing adjoint vector quite often (#=" << stce_10294
                          << ") -- might have severe performance impact (reallocation / initialization / copy)";
    }
    if (stce_10292 >= stce_10096) {
      DCO_LOG(logDEBUG1) << "dco_adjoint_vector: resizing adjoint vector to " << stce_10096 << "; no reallocation required.";
      for (DCO_TAPE_INT stce_9827 = stce_10287; stce_9827 < stce_10096; ++stce_9827) {
        stce_10284[stce_9827] = stce_9764();
      }
      stce_10287 = stce_10096;
      return;
    } else {
      DCO_LOG(logDEBUG1) << "dco_adjoint_vector: resizing adjoint vector to " << stce_10096 << "; reallocation required.";
      stce_9764* stce_10299 = static_cast<stce_9764*>(dco::alloc(stce_10096 * sizeof(stce_9764), stce_10191.alloc_use_huge_pages()));
      if (stce_10299 == 0) {
        throw dco::exception::create<std::runtime_error>("Unable to allocate memory for adjoint vector.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_lazy_adjoint_vector.hpp", 101);
      }
      if (stce_10284) {
        for (DCO_TAPE_INT stce_9827 = 0; stce_9827 < stce_10287; ++stce_9827) { stce_10299[stce_9827] = stce_10284[stce_9827]; }
        for (DCO_TAPE_INT stce_9827 = stce_10287; stce_9827 < stce_10096 ; ++stce_9827) { stce_10299[stce_9827] = stce_9764(); }
        dco::dealloc(stce_10284);
      } else {
        for (DCO_TAPE_INT stce_9827 = 0; stce_9827 < stce_10096 ; ++stce_9827) { stce_10299[stce_9827] = stce_9764(); }
      }
      stce_10287 = stce_10096;
      stce_10292 = stce_10096;
      stce_10284 = stce_10299;
    }
  }
  void stce_10289(DCO_TAPE_INT stce_9806) {
    if(!stce_10285(stce_9806)) { stce_10298(stce_9806); }
  }
  stce_9764 at(DCO_TAPE_INT stce_9806) const {
    return stce_10285(stce_9806)?stce_10284[stce_9806]:0;
  }
  stce_9764 &at(const DCO_TAPE_INT stce_9806) {
    if(!stce_10285(stce_9806)) { stce_10298(stce_9806+1); }
    return stce_10284[stce_9806];
  }
  void stce_10297(DCO_TAPE_INT stce_9858, DCO_TAPE_INT stce_9859) {
    if(stce_9858 < stce_9859) { std::swap(stce_9858, stce_9859); }
    if(!stce_10285(stce_9859)) { return; }
    stce_10300();
    for (DCO_TAPE_INT stce_9827 = std::min(stce_9858, stce_10293-1); stce_9827 > stce_9859; --stce_9827) {
      stce_10284[stce_9827] = 0;
    }
  }
  void stce_10300() { stce_10298(stce_10293); }
};
}
}
namespace dco {
template<class stce_10137, class stce_9761, int stce_10301=-1>
  struct adjoint_vector : dco::stce_10073<typename stce_10137::tape_t::stce_9764>, public internal::stce_10138<stce_10137> {
    static const int stce_9833 = stce_10301 < 0 ? -stce_10301 : stce_10301;
  typedef internal::stce_10138<stce_10137> stce_10145;
public:
  typedef dco::helper::stce_9822<stce_9761, stce_9833> value_t;
  typedef stce_9761 adjoint_real_t;
  typedef typename stce_10137::stce_9764 stce_9764;
private:
  value_t* stce_10284;
  DCO_TAPE_INT stce_10302;
public:
  int stce_10303;
  void init() {}
  adjoint_vector(stce_10137 *tape) : stce_10145(*tape), stce_10284(0), stce_10302(0), stce_10303(0) { }
  virtual ~adjoint_vector() {
    if (stce_10284) delete [] stce_10284;
  }
  template<class stce_9794>
    value_t& derivative(const stce_9794 &stce_9882) {
    stce_10304();
    const DCO_TAPE_INT stce_9806 = stce_9882._tape_index();
    return stce_10284[stce_9806];
  }
  void interpret_adjoint() {
    typename stce_10137::iterator stce_9859 = stce_10145::stce_10140().start();
    internal::interpretation_settings stce_10305;
    stce_10306(stce_10145::stce_10140().current(), stce_9859, stce_10305);
  }
  void interpret_adjoint_to(const typename stce_10137::iterator &stce_9859) {
    internal::interpretation_settings stce_10305;
    if (stce_9859 > stce_10145::stce_10140().current())
      throw dco::exception::create<std::runtime_error>("adjoint interpretation: from < to.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_adjoint_vector.hpp", 59);
    else
      stce_10306(stce_10145::stce_10140().current(), stce_9859, stce_10305);
  }
  void interpret_adjoint_from(const typename stce_10137::iterator &stce_9858) {
    typename stce_10137::iterator stce_9859 = stce_10145::stce_10140().start();
    assert(!(stce_9858 < stce_9859));
    typename internal::interpretation_settings stce_10305;
    stce_10306(stce_9858, stce_9859, stce_10305);
  }
  void interpret_adjoint_from_to(const typename stce_10137::iterator &stce_9858, const typename stce_10137::iterator &stce_9859) {
    typename internal::interpretation_settings stce_10305;
    if (stce_9859 > stce_9858)
      throw dco::exception::create<std::runtime_error>("adjoint interpretation: from < to.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_adjoint_vector.hpp", 75);
    else
      stce_10306(stce_9858, stce_9859, stce_10305);
  }
  void interpret_adjoint_and_zero_adjoints_to(const typename stce_10137::iterator &stce_9859) {
    typename stce_10137::iterator stce_9858(stce_10145::stce_10140().current());
    typename internal::interpretation_settings stce_10305;
    stce_10305.reset = false;
    stce_10305.stce_10059 = true;
    stce_10306(stce_9858, stce_9859, stce_10305);
  }
  void interpret_adjoint_and_zero_adjoints_from_to(const typename stce_10137::iterator &stce_9858, const typename stce_10137::iterator &stce_9859) {
    typename internal::interpretation_settings stce_10305;
    stce_10305.reset = false;
    stce_10305.stce_10059 = true;
    stce_10306(stce_9858, stce_9859, stce_10305);
  }
  void zero_adjoints() {
    typename stce_10137::iterator stce_9859 = stce_10145::stce_10140().start();
    stce_10307(stce_10145::stce_10140().current(), stce_9859);
  }
  void zero_adjoints_to(const typename stce_10137::iterator &stce_9859) {
    stce_10307(stce_10145::stce_10140().current(), stce_9859);
  }
  void zero_adjoints_from(const typename stce_10137::iterator &stce_9858) {
    typename stce_10137::iterator stce_9859 = stce_10145::stce_10140().start();
    stce_10307(stce_9858, stce_9859);
  }
  void zero_adjoints_from_to(const typename stce_10137::iterator &stce_9858, const typename stce_10137::iterator &stce_9859) {
    stce_10307(stce_9858, stce_9859);
  }
  value_t& operator[](DCO_TAPE_INT stce_9806) {
    return stce_10284[stce_9806];
  }
  value_t at(DCO_TAPE_INT stce_9806) const { return (stce_9806 < stce_10302)?stce_10284[stce_9806]:0; }
  value_t& at(DCO_TAPE_INT stce_9806) { stce_10289(stce_9806); return stce_10284[stce_9806]; }
  void stce_10289(DCO_TAPE_INT stce_9806) {
    if(stce_10302 < stce_9806) {
      stce_10302 = stce_9806;
    }
    stce_10304();
  }
  void stce_10290(typename stce_10137::iterator& stce_9815) {
    (*this)[stce_9815->arg()] += stce_9815->pval()*(*this)[stce_9815.index()];
  }
  value_t& stce_10043(const DCO_TAPE_INT stce_9806) {
    return stce_10284[stce_9806];
  }
  stce_9761& stce_10043(const DCO_TAPE_INT stce_9806, int stce_10308) {
    return stce_10284[stce_9806][stce_10308];
  }
  template<class stce_10309>
    void stce_10075(const DCO_TAPE_INT &stce_9806, const stce_10309 &stce_9981, const value_t &stce_10076) {
    stce_10284[stce_9806] += (stce_9981 * stce_10076);
  }
  virtual typename stce_10137::stce_9764 stce_10074(const DCO_TAPE_INT stce_9806) {
    return static_cast<typename stce_10137::stce_9764>(stce_10043(stce_9806)[stce_10303]);
  }
  virtual void stce_10075(const DCO_TAPE_INT stce_9806, const typename stce_10137::stce_9764 &stce_10076) {
    stce_10284[stce_9806][stce_10303] += stce_10076;
  }
  void stce_10307(const typename stce_10137::iterator &stce_9858, const typename stce_10137::iterator &stce_9859) {
    stce_10297(stce_9858.index(), stce_9859.index());
  }
  void stce_10297(DCO_TAPE_INT stce_9858, DCO_TAPE_INT stce_9859) {
    if(stce_9858 < stce_9859) { std::swap(stce_9858, stce_9859); }
    if(stce_10284==0) { return; }
    for(DCO_TAPE_INT stce_9827=stce_9859+1;stce_9827<=stce_9858;stce_9827++) {
      stce_10284[stce_9827] = 0;
    }
  }
  void stce_10306(const typename stce_10137::iterator &stce_9858,
                                   const typename stce_10137::iterator &stce_9859,
                                   typename internal::interpretation_settings stce_10305)
  {
    stce_10304();
    stce_10145::stce_10140().stce_10310(*this, stce_9858, stce_9859, stce_10305);
  }
  void stce_10304() {
    if(stce_10284==0) {
      stce_10302 = stce_10145::stce_10140().current().index()+1;
      stce_10284 = new value_t[size_t(stce_10302)]();
    }
    else {
      DCO_TAPE_INT stce_10311 = stce_10302 = stce_10145::stce_10140().current().index()+1;
      if(stce_10311 > stce_10302) {
        value_t *stce_10312 = new value_t[size_t(stce_10311)]();
        for(DCO_TAPE_INT stce_9827=0;stce_9827<stce_10302;++stce_9827)
          stce_10312[stce_9827] = stce_10284[stce_9827];
        delete [] stce_10284;
        stce_10284 = stce_10312;
      }
    }
  }
};
}
       
namespace dco {
template< typename stce_10313, typename adjoint_vector >
  void stce_10314(stce_10313 stce_9815, stce_10313 stce_9859, const adjoint_vector& stce_10315,
                 const std::string &filename = "tape.dot") {
  stce_10314(stce_9815, stce_9859 , stce_10315, stce_9844(filename));
}
template< typename stce_10313, typename adjoint_vector >
  void stce_10314(stce_10313 stce_9815, stce_10313 stce_9859, const adjoint_vector& stce_10315,
                 stce_9844 &stce_9845) {
  if (stce_9815.index() < stce_9859.index()) {
    std::swap(stce_9815, stce_9859);
  }
  DCO_TAPE_INT stce_10316 = stce_9815.index();
  stce_9845.stce_9854(stce_10316,stce_10315.at(stce_10316));
  for (; stce_9815 != stce_9859; --stce_9815) {
    stce_9845.stce_9857(stce_9815->arg(), stce_10316, stce_9815->pval());
    if(stce_9815->stce_10151()) {
      stce_10316--;
      stce_9845.stce_9854(stce_10316,stce_10315.at(stce_10316));
    }
  }
}
template< typename stce_10313, typename adjoint_vector, typename stce_10317 >
  void stce_10314(stce_10313 stce_9815, stce_10313 stce_9859, const adjoint_vector& stce_10315,
                 const stce_10317& stce_10318, stce_9844 &stce_9845) {
  if (stce_9815.index() < stce_9859.index()) {
    std::swap(stce_9815, stce_9859);
  }
  stce_10314(stce_9815, stce_9859, stce_10315, stce_9845);
  typedef typename stce_10317::const_reverse_iterator stce_10319;
  for (stce_10319 stce_10158 = stce_10318.rbegin(); stce_10158 != stce_10318.rend() && stce_9859 <= stce_10158->stce_10320(); ++stce_10158) {
    if (stce_9815 < stce_10158->stce_10320())
      continue;
    stce_10158->stce_10321()->stce_10082(stce_9845, stce_10158->stce_10320().index() + 1);
  }
}
template <typename stce_10313>
  void stce_10322(stce_10313 stce_9815, stce_10313 stce_9859,
                 std::string const& filename) {
  if (stce_9815.index() < stce_9859.index()) { std::swap(stce_9815, stce_9859); }
  std::ofstream stce_9915(filename.c_str());
  stce_9915 << stce_9815.index() << "\n";
  for (; stce_9815 != stce_9859; --stce_9815) {
    stce_9915 << stce_9815->arg() << ", " << stce_9815->pval() << "; " ;
    if(stce_9815->stce_10151())
      stce_9915 << "\n";
  }
}
template <bool USE_SPARSE_INTERPRET, typename stce_10313, typename adjoint_vector>
  DCO_STRONG_INLINE void interpret(stce_10313 stce_9815, const stce_10313& stce_9859,
                                            adjoint_vector& stce_10315) {
  stce_9815.index()--;
  stce_1000002();;
  stce_10315.stce_10289(stce_9808(stce_9752(stce_9815.index())));
  stce_10315.init();
  for (; stce_9815 != stce_9859; --stce_9815) {
    if(USE_SPARSE_INTERPRET && folding::is_zero(stce_10315[stce_9815.index()])) { continue; }
    DCO_LOG(logDEBUG3) << stce_10315[stce_9815->arg()] << " <=> adj[" << stce_9815->arg() << "] += " << stce_9815->pval() * stce_10315[stce_9815.index()] << " from idx " << stce_9815.index();
    stce_10315.stce_10290(stce_9815);
  }
}
template <bool USE_SPARSE_INTERPRET, typename stce_10313, typename adjoint_vector>
  DCO_STRONG_INLINE void stce_10323(stce_10313 stce_9815,
                                           const stce_10313& stce_9859,
                                           adjoint_vector& stce_10315) {
  stce_9815.index()--;
  stce_1000002();;
  stce_10315.stce_10289(stce_9808(stce_9752(stce_9815.index())));
  for (; stce_9815 != stce_9859; --stce_9815) {
    if(USE_SPARSE_INTERPRET && folding::is_zero(stce_10315[stce_9815.index()])) { continue; }
    DCO_LOG(logDEBUG3) << stce_10315[stce_9815->arg()] << " <=> adj[" << stce_9815->arg() << "] += " << stce_9815->pval() * stce_10315[stce_9815.index()] << " from idx " << stce_9815.index();
    stce_10315.stce_10290(stce_9815);
    if(stce_9815->stce_10151()) {
      stce_10315[stce_9815.index()] = 0;
    }
  }
}
template< typename stce_10313, typename adjoint_vector >
  DCO_STRONG_INLINE void stce_10324(stce_10313 stce_9815, const stce_10313& stce_9859, adjoint_vector& stce_10315) {
  stce_10315.stce_10289(stce_9815.index());
  for (; stce_9815 != stce_9859; --stce_9815) {
    std::cout << stce_9815.index() << ": " << stce_9815->arg() << ", " << stce_9815->pval() << ", " << stce_10315[stce_9815.index()] << std::endl;
  }
}
template< typename tape_t >
  DCO_STRONG_INLINE void stce_10324(const tape_t *tape) {
  typename tape_t::iterator stce_9815 = tape->get_position();
  for (; stce_9815 != tape->begin(); --stce_9815) {
    std::cout << stce_9815.index() << ": " << stce_9815->arg() << ", " << stce_9815->pval() << std::endl;
  }
}
}
       
       
       
namespace dco {
namespace internal {
template <typename stce_10325> struct stce_10326 {
  typedef stce_10325 type;
  typedef stce_10325& stce_10327;
  typedef const stce_10325 stce_10328;
};
template <typename stce_9782 > struct stce_10326<std::vector<stce_9782> > {
  typedef std::vector<stce_9782> & type;
  typedef std::vector<stce_9782> & stce_10327;
  typedef std::vector<stce_9782> const& stce_10328;
};
template <stce_9796 stce_10329, typename stce_9782, typename stce_10330 = void>
  struct stce_10331 { };
template <typename stce_9782> struct stce_10331<stce_9797, stce_9782> {
  typedef typename dco::mode<stce_9782>::value_t value_t;
  typedef typename trait_value<stce_9782>::stce_9802 stce_9802;
  typedef typename trait_value<stce_9782>::stce_10025 stce_10025;
  stce_9802 get( stce_9782& stce_9839) { return trait_value<stce_9782>::value(stce_9839); };
  stce_10025 get(const stce_9782& stce_9839) const { return trait_value<stce_9782>::value(stce_9839); };
};
template <typename stce_9782> struct stce_10331<stce_9800, stce_9782> {
  typedef typename dco::mode<stce_9782>::passive_t value_t;
  typedef typename stce_9785<stce_9782>::stce_9802 stce_9802;
  typedef typename stce_9785<stce_9782>::stce_10025 stce_10025;
  stce_9802 get( stce_9782& stce_9839) { return stce_9785<stce_9782>::value(stce_9839); };
  stce_10025 get(const stce_9782& stce_9839) const { return stce_9785<stce_9782>::value(stce_9839); };
};
template <typename stce_9782> struct stce_10331<stce_9798, stce_9782> {
  typedef typename dco::mode<stce_9782>::derivative_t value_t;
  typedef typename stce_9784<stce_9782>::stce_9802 stce_9802;
  typedef typename stce_9784<stce_9782>::stce_9802 stce_10025;
  stce_10025 get(const stce_9782& stce_9839) const { return stce_9784<stce_9782>::value(stce_9839); };
};
template <typename stce_9782> struct stce_10331<stce_9801, stce_9782> {
  typedef typename stce_9787<stce_9782>::stce_10332 stce_9802;
  typedef typename stce_9787<stce_9782>::stce_10025 stce_10025;
  typedef stce_9802 value_t;
  stce_9802 get( stce_9782& stce_9839) { return stce_9787<stce_9782>::value(stce_9839); };
  stce_10025 get(const stce_9782& stce_9839) const { return stce_9787<stce_9782>::value(stce_9839); };
};
template <typename stce_9782, typename stce_10330> struct stce_10331<stce_9799, stce_9782, stce_10330> {
  typedef typename dco::mode<stce_9782>::derivative_t value_t;
  typedef typename stce_9784<stce_9782, stce_10330>::stce_10333 stce_9802;
  typedef typename stce_9784<stce_9782, stce_10330>::stce_10333 stce_10025;
  stce_10330* stce_10284;
  stce_10331(stce_10330* adjoints) : stce_10284(adjoints) {}
  stce_9802 get( stce_9782& stce_9839) { return stce_9784<stce_9782, stce_10330>::value(stce_9839, stce_10284); };
  stce_10025 get(const stce_9782& stce_9839) const { return stce_9784<stce_9782, stce_10330>::value(stce_9839, stce_10284); };
};
template<stce_9796 stce_10329, class stce_9782, typename stce_10325, typename stce_9783 = void>
  struct stce_10334 : stce_10331<stce_10329, stce_9782, stce_9783> {
  typedef typename stce_10326<stce_10325>::type stce_10277;
  typedef typename stce_10326<stce_10325>::stce_10327 stce_10335;
  typedef typename stce_10326<stce_10325>::stce_10328 stce_10336;
  typedef stce_10331<stce_10329, stce_9782, stce_9783> stce_10337;
  typedef typename stce_10337::value_t value_t;
  stce_10277 stce_10338;
  stce_10334(stce_10336 stce_9838) : stce_10338(const_cast<stce_10335>(stce_9838)) {}
  stce_10334(stce_10336 stce_9838, stce_9783* stce_10339) : stce_10337(stce_10339), stce_10338(const_cast<stce_10335>(stce_9838)) {}
  stce_10334 &operator=(const std::vector<value_t> &stce_9890) {
    for (size_t stce_9827 = 0; stce_9827 < stce_10338.size(); ++stce_9827) {
      stce_10337::get(stce_10338[stce_9827]) = stce_9890[stce_9827];
    }
    return *this;
  }
  operator std::vector<value_t>() const {
    std::vector<value_t> stce_9836(stce_10338.size());
    for (size_t stce_9827 = 0; stce_9827 < stce_10338.size(); ++stce_9827) {
      stce_9836[stce_9827] = stce_10337::get(stce_10338[stce_9827]);
    }
    return stce_9836;
  }
  typename stce_10337::stce_9802 operator[](size_t stce_9827) { return stce_10337::get(stce_10338[stce_9827]); }
  typename stce_10337::stce_10025 operator[](size_t stce_9827) const {
    return stce_10337::get(stce_10338[stce_9827]);
  }
  size_t size() const { return stce_10338.size(); }
};
}
}
namespace dco {
template <typename stce_9782> inline typename internal::trait_value<stce_9782>::stce_9802
  value( stce_9782 &stce_9839) { return internal::trait_value<stce_9782>::value(stce_9839); }
template <typename stce_9782> inline typename internal::trait_value<stce_9782>::stce_10025
  value(const stce_9782 &stce_9839) { return internal::trait_value<stce_9782>::value(stce_9839); }
namespace internal {
template <typename stce_9782, typename enable_if>
  struct trait_value {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef type& stce_9802;
  typedef const type& stce_10025;
  static inline stce_9802 value( type &value) { return value; }
  static inline stce_10025 value(const type &value) { return value; }
};
template <typename stce_9782>
  struct trait_value<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type && !is_gbcp_type<stce_9782>::value && !dco::mode<stce_9782>::is_intermediate_type>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::value_t& stce_9802;
  typedef typename dco::mode<type>::value_t const& stce_10025;
  static inline stce_9802 value( type &value) { return value._value(); }
  static inline stce_10025 value(const type &value) { return value._value(); }
};
template <typename stce_9782>
  struct trait_value<stce_9782, typename dco::enable_if<!is_gbcp_type<stce_9782>::value && dco::mode<stce_9782>::is_intermediate_type>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::value_t const& stce_9802;
  typedef typename dco::mode<type>::value_t const& stce_10025;
  static inline stce_9802 value( type &value) { return value._value(); }
  static inline stce_10025 value(const type &value) { return value._value(); }
};
template <typename stce_9782>
  struct trait_value<stce_9782, typename dco::enable_if<is_gbcp_type<stce_9782>::value>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::value_t& stce_9802;
  typedef typename dco::mode<type>::value_t const& stce_10025;
  static inline stce_9802 value( type &stce_10340) { return dco::value(stce_10340._value()); }
  static inline stce_10025 value(const type &stce_10340) { return dco::value(stce_10340._value()); }
};
template <typename stce_9782>
  struct trait_value<std::vector<stce_9782>, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type>::type> {
  typedef std::vector<stce_9782> stce_10277;
  typedef stce_10334<stce_9797, stce_9782, stce_10277> stce_9802;
  typedef const stce_10334<stce_9797, stce_9782, stce_10277> stce_10025;
  static inline stce_9802 value( stce_10277 &stce_9838) { return stce_9802(stce_9838); }
  static inline stce_10025 value(const stce_10277 &stce_9838) { return stce_10025(stce_9838); }
};
template<stce_9796 stce_10329, class stce_9782, typename stce_10341, typename stce_10342>
  struct trait_value<stce_10334<stce_10329, stce_9782, stce_10341, stce_10342>, void > {
  typedef stce_10334<stce_10329, stce_9782, stce_10341, stce_10342> stce_10277;
  typedef stce_10334<stce_9797, typename stce_10277::value_t, stce_10277> stce_9802;
  typedef stce_10334<stce_9797, typename stce_10277::value_t, stce_10277> stce_10025;
  static inline stce_9802 value(const stce_10277 &stce_9838) { return stce_9802(stce_9838); }
};
}
}
       
namespace dco {
template <typename stce_9782> inline typename internal::stce_9785<stce_9782>::stce_9802
  passive_value( stce_9782 &stce_9839) { return internal::stce_9785<stce_9782>::value(stce_9839); }
template <typename stce_9782> inline typename internal::stce_9785<stce_9782>::stce_10025
  passive_value(const stce_9782 &stce_9839) { return internal::stce_9785<stce_9782>::value(stce_9839); }
namespace internal {
template <typename stce_9782, typename enable_if>
  struct stce_9785 {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef type& stce_9802;
  typedef type const& stce_10025;
  static inline stce_9802 value( type &value) { return value; }
  static inline stce_10025 value(const type &value) { return value; }
};
template <typename stce_9782>
  struct stce_9785<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type && !dco::mode<stce_9782>::is_intermediate_type>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::active_value_t active_value_t;
  typedef typename dco::mode<type>::value_t value_t;
  typedef typename dco::mode<type>::passive_t passive_t;
  typedef passive_t& stce_9802;
  typedef passive_t const& stce_10025;
  static inline stce_9802 value( type &value) { return stce_9785<value_t>::value(dco::value(value)); }
  static inline stce_10025 value(const type &value) { return stce_9785<value_t>::value(dco::value(value)); }
};
template <typename stce_9782>
  struct stce_9785<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_intermediate_type>::type> {
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename dco::mode<type>::active_value_t active_value_t;
  typedef typename dco::mode<type>::value_t value_t;
  typedef typename dco::mode<type>::passive_t passive_t;
  typedef passive_t const& stce_9802;
  typedef passive_t const& stce_10025;
  static inline stce_9802 value( type &value) { return stce_9785<value_t>::value(dco::value(value)); }
  static inline stce_10025 value(const type &value) { return stce_9785<value_t>::value(dco::value(value)); }
};
template <typename stce_9782>
  struct stce_9785<std::vector<stce_9782>, typename dco::enable_if<!dco::mode<stce_9782>::is_dco_type>::type> {
  typedef std::vector<stce_9782> stce_10277;
  typedef stce_10277& stce_9802;
  typedef const stce_10277& stce_10025;
  static inline stce_9802 value( stce_10277 &value) { return value; }
  static inline stce_10025 value(const stce_10277 &value) { return value; }
};
template <typename stce_9782>
  struct stce_9785<std::vector<stce_9782>, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type>::type> {
  typedef std::vector<stce_9782> stce_10277;
  typedef stce_10334<stce_9800, stce_9782, stce_10277> stce_9802;
  typedef const stce_10334<stce_9800, stce_9782, stce_10277> stce_10025;
  static inline stce_9802 value( stce_10277 &stce_9838) { return stce_9802(stce_9838); }
  static inline stce_10025 value(const stce_10277 &stce_9838) { return stce_10025(stce_9838); }
};
template<stce_9796 stce_10329, class stce_9782, typename stce_10341, typename stce_10342>
  struct stce_9785<stce_10334<stce_10329, stce_9782, stce_10341, stce_10342>, void > {
  typedef stce_10334<stce_10329, stce_9782, stce_10341, stce_10342> stce_10277;
  typedef stce_10334<stce_9800, typename stce_10277::value_t, stce_10277> stce_9802;
  typedef stce_10334<stce_9800, typename stce_10277::value_t, stce_10277> stce_10025;
  static inline stce_9802 value(const stce_10277 &stce_9838) { return stce_9802(stce_9838); }
};
}
}
       
namespace dco {
template <typename stce_9782> inline typename internal::stce_9784<stce_9782>::stce_9802
  derivative(const stce_9782 &stce_9839) { return internal::stce_9784<stce_9782>::value(stce_9839); }
template <typename stce_9782, typename stce_9783> inline typename internal::stce_9784<stce_9782, stce_9783>::stce_10333
  derivative(const stce_9782 &stce_9839, stce_9783& adjoints) { return internal::stce_9784<stce_9782, stce_9783>::value(stce_9839, &adjoints); }
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
  static inline stce_9802 value(type const& value) { return value.stce_9985(); }
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
template <typename stce_9782, typename stce_10343, typename stce_10344, int stce_10345>
  struct stce_9784<stce_9782, adjoint_vector<stce_10343, stce_10344, stce_10345>, typename dco::enable_if<(stce_10345 == -1)>::type> {
  typedef adjoint_vector<stce_10343, stce_10344, stce_10345> stce_10342;
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef stce_10344& stce_10333;
  static inline stce_10333 value(const type &value, stce_10342* adjoints) { return adjoints->derivative(value)[0]; }
};
template <typename stce_9782, typename stce_10343, typename stce_10344, int stce_10345>
  struct stce_9784<stce_9782, adjoint_vector<stce_10343, stce_10344, stce_10345>, typename dco::enable_if<(stce_10345 > 0)>::type> {
  typedef adjoint_vector<stce_10343, stce_10344, stce_10345> stce_10342;
  typedef typename dco::remove_const<stce_9782>::type type;
  typedef typename stce_10342::value_t& stce_10333;
  static inline stce_10333 value(const type &value, stce_10342* adjoints) { return adjoints->derivative(value); }
};
template <typename stce_9782>
  struct stce_9784<std::vector<stce_9782>, void, typename dco::enable_if<!dco::mode<stce_9782>::is_dco_type>::type> {
  typedef std::vector<stce_9782> stce_10277;
  typedef stce_10277 stce_9802;
  static inline stce_9802 value(const stce_10277 &value) { return stce_9802(value.size(), 0.0); }
};
template <typename stce_9782>
  struct stce_9784<std::vector<stce_9782>, void, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type>::type> {
  typedef std::vector<stce_9782> stce_10277;
  typedef stce_10334<stce_9798, stce_9782, stce_10277> stce_9802;
  static inline stce_9802 value(const stce_10277 &stce_9838) { return stce_9802(stce_9838); }
};
template <typename stce_9782, typename stce_10343, typename stce_10344, int stce_10345>
  struct stce_9784<std::vector<stce_9782>, adjoint_vector<stce_10343, stce_10344, stce_10345>, typename dco::enable_if<(stce_10345 == -1)>::type> {
  typedef adjoint_vector<stce_10343, stce_10344, stce_10345> stce_10342;
  typedef std::vector<stce_9782> stce_10277;
  typedef stce_10334<stce_9798, stce_9782, stce_10277> stce_9802;
  typedef stce_10334<stce_9799, stce_9782, stce_10277, stce_10342> stce_10333;
  static inline stce_9802 value(const stce_10277 &stce_9838) { return stce_9802(stce_9838); }
  static inline stce_10333 value(const stce_10277 &stce_9838, stce_10342* adjoints) { return stce_10333(stce_9838, adjoints); }
};
template<stce_9796 stce_10329, class stce_9782, typename stce_10341, typename stce_10342>
  struct stce_9784<stce_10334<stce_10329, stce_9782, stce_10341, stce_10342>, void > {
  typedef stce_10334<stce_10329, stce_9782, stce_10341, stce_10342> stce_10277;
  typedef stce_10334<stce_9798, typename stce_10277::value_t, stce_10277> stce_9802;
  static inline stce_9802 value(const stce_10277 &stce_9838) { return stce_9802(stce_9838); }
};
}
}
       
namespace dco {
template <typename stce_9782>
  typename internal::stce_9786<stce_9782>::stce_9802 tape(const stce_9782 &stce_9839) {
  return internal::stce_9786<stce_9782>::value(stce_9839);
}
template <typename stce_9782>
  typename internal::stce_9787<stce_9782>::stce_10025 tape_index(const stce_9782 &stce_9839) {
  return internal::stce_9787<stce_9782>::value(stce_9839);
}
template <typename stce_9782>
  typename internal::stce_9787<stce_9782>::stce_9802 tape_index(stce_9782 &stce_9839) {
  return internal::stce_9787<stce_9782>::value(stce_9839);
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
    for (size_t stce_9827 = 0; stce_9827 < value.size(); stce_9827++) {
      if (stce_9786<stce_9782>::value(value[stce_9827]) != 0) {
        tape = stce_9786<stce_9782>::value(value[stce_9827]);
        break;
      }
    }
    return tape;
  }
};
template<typename stce_9782, typename enable_if> struct stce_9787 {
  typedef DCO_TAPE_INT stce_9802;
  typedef DCO_TAPE_INT stce_10332;
  typedef const DCO_TAPE_INT stce_10025;
  static DCO_TAPE_INT value(const stce_9782&) { return 0; }
};
template<typename stce_9782> struct stce_9787<stce_9782, typename dco::enable_if<dco::mode<stce_9782>::is_adjoint_type>::type> {
  typedef DCO_TAPE_INT& stce_9802;
  typedef DCO_TAPE_INT stce_10332;
  typedef const DCO_TAPE_INT stce_10025;
  static DCO_TAPE_INT value(const stce_9782& value) { return value._tape_index(); }
  static DCO_TAPE_INT& value(stce_9782& value) { return value._tape_index(); }
};
template <typename stce_9782>
  struct stce_9787<std::vector<stce_9782>, typename dco::enable_if<dco::mode<stce_9782>::is_dco_type>::type> {
  typedef std::vector<stce_9782> stce_10277;
  typedef stce_10334<stce_9801, stce_9782, stce_10277> stce_9802;
  typedef const stce_10334<stce_9801, stce_9782, stce_10277> stce_10025;
  static inline stce_9802 value( stce_10277 &stce_9838) { return stce_9802(stce_9838); }
  static inline stce_10025 value(const stce_10277 &stce_9838) { return stce_10025(stce_9838); }
};
}
}
namespace dco {
}
       
namespace dco {
template <typename tape_t>
struct stce_10346 {
  typedef std::map<std::string, mem_long_t> stce_10347;
  stce_10347 stce_10348;
  tape_t* stce_9948;
  stce_10346(tape_t* tape) : stce_9948(tape) {
    stce_10348["max_tape_memory"] = 0;
    stce_10348["max_adjoint_memory"] = 0;
    stce_10348["max_checkpoint_memory"] = 0;
    stce_10348["max_total_memory"] = 0;
  }
  mem_long_t operator[](std::string stce_10349) { return stce_10348[stce_10349]; }
  void stce_10350() {
    mem_long_t stce_10351 = dco::size_of(stce_9948, tape_t::size_of_stack);
    if (stce_10351 > stce_10348["max_tape_memory"]) {
      stce_10348["max_tape_memory"] = stce_10351;
    }
    stce_10351 = dco::size_of(stce_9948, tape_t::size_of_internal_adjoint_vector);
    if (stce_10351 > stce_10348["max_adjoint_memory"]) {
      stce_10348["max_adjoint_memory"] = stce_10351;
    }
    stce_10351 = dco::size_of(stce_9948, tape_t::size_of_checkpoints);
    if (stce_10351 > stce_10348["max_checkpoint_memory"]) {
      stce_10348["max_checkpoint_memory"] = stce_10351;
    }
    stce_10348["max_total_memory"] = 0;
    for (stce_10347::iterator stce_9827 = stce_10348.begin(); stce_9827 != stce_10348.end(); ++stce_9827) {
      stce_10348["max_total_memory"] += stce_9827->second;
    }
  }
  mem_long_t stce_10352() { return stce_10348["max_tape_memory"]; }
};
}
namespace dco {
namespace internal {
template<memory_model::TYPE MEMORY_MODEL,
         typename stce_10353,
         memory_model::TYPE stce_9789>
  class tape : public dco::stce_10073<typename helper::stce_9832<typename stce_10353::stce_9764>::stce_9834>,
               public stce_10143<MEMORY_MODEL,
                                     tape<MEMORY_MODEL, stce_10353, stce_9789>,
                                     stce_9789> {
public:
  typedef stce_10353 stce_9762;
  typedef stce_10143<MEMORY_MODEL, tape, stce_9789> stce_10277;
  typedef typename stce_10277::iterator iterator_t;
  typedef iterator_t position_t;
  typedef typename stce_9762::stce_9764 stce_9764;
  typedef typename stce_9762::stce_9763 stce_9763;
  typedef typename helper::stce_9832<typename stce_10353::stce_9764>::stce_9834 stce_10079;
  typedef tape<MEMORY_MODEL, stce_9762, stce_9789> tape_t;
  typedef stce_9793<typename stce_9762::stce_9764, stce_10277> vector_t;
  vector_t stce_10284;
  stce_10346<tape> stce_10354;
  stce_10346<tape>& stce_10355() { return stce_10354; }
  const stce_10346<tape>& stce_10355() const { return stce_10354; }
  tape& get_tape() { return *this; }
private:
  int stce_10356;
  bool stce_10357;
  bool stce_10358;
  DCO_TAPE_INT stce_10303;
  tape_options stce_10191;
  bool stce_10359;
  void stce_10307(const iterator_t &stce_9858, const iterator_t &stce_9859 = stce_10277::start()) {
    stce_10284.stce_10297(stce_9858.index(), stce_9859.index());
  }
  void stce_10360() {
    stce_10284.resize(this->stce_10179());
    stce_10284.stce_10300();
  }
  template <typename stce_10313, typename adjoint_vector>
    DCO_STRONG_INLINE void stce_10361(stce_10313 stce_9815, const stce_10313& stce_9859,
                                      adjoint_vector& stce_10315) {
    if (stce_10359) { dco::interpret< true >(stce_9815, stce_9859, stce_10315); }
    else { dco::interpret< false>(stce_9815, stce_9859, stce_10315); }
  }
  template <typename stce_10313, typename adjoint_vector>
    DCO_STRONG_INLINE void stce_10362(stce_10313 stce_9815, const stce_10313& stce_9859,
                                      adjoint_vector& stce_10315) {
    if (stce_10359) { dco::stce_10323< true >(stce_9815, stce_9859, stce_10315); }
    else { dco::stce_10323< false>(stce_9815, stce_9859, stce_10315); }
  }
public:
  const tape_options& stce_10164() const { return stce_10191; }
  bool stce_10363() { return this->start() == this->current(); }
  tape(tape_options stce_10164 = tape_options()) :
    stce_10277(stce_10164, *this), stce_10284(1, *this, stce_10164),
    stce_10354(this),
    stce_10356(dco::build_info::COMPATIBILITY_VERSION),
    stce_10357(false), stce_10358(true), stce_10303(0), stce_10191(stce_10164),
    stce_10359(false) { }
  ~tape() {
    for(typename std::vector<stce_10364>::iterator stce_10158 = stce_10365.begin();
        stce_10158 != stce_10365.end(); ++stce_10158) {
      stce_10158->stce_10366();
    }
  }
  bool is_compatible(const int stce_10367) {
    return stce_10367 == stce_10356;
  }
  bool is_active() { return stce_9758 ? stce_10358 : true; }
  typedef helper::stce_10078<tape> stce_10368;
  static void stce_10369(stce_10368* stce_10370) { stce_10370->get_tape()->stce_10359 = true; }
  static void stce_10371(stce_10368* stce_10370) { stce_10370->get_tape()->stce_10359 = false; }
  void set_sparse_interpret() {
    if (stce_10359 == false) {
      stce_10368* stce_10370 = this->create_callback_object<stce_10368>();
      this->insert_callback(stce_10371, stce_10370);
    }
    stce_10359 = true;
  }
  void unset_sparse_interpret() {
    if (stce_10359 == true) {
      stce_10368* stce_10370 = this->create_callback_object<stce_10368>();
      this->insert_callback(stce_10369, stce_10370);
    }
    stce_10359 = false;
  }
  bool& sparse_interpret() { return stce_10359; }
  void switch_to_active() {
    if (!stce_9758) {
      throw dco::exception::create<std::runtime_error>("activity check disabled but switched", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 160);
    }
    if (!stce_10358) { stce_10358 = true; }
  }
  void switch_to_passive() {
    if (!stce_9758) {
      throw dco::exception::create<std::runtime_error>("activity check disabled but switched", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 167);
    }
    if (stce_10358) { stce_10358 = false; }
  }
  template <bool stce_10372>
  DCO_TAPE_INT stce_10373() {
    if(!stce_10372 && this->empty()) { return stce_10277::stce_9924().stce_10162(); }
    else {
      typename stce_10277::stce_9920 stce_10374(stce_10277::stce_9924());
      stce_10374.insert(0, 0);
      return stce_10374.index();
    }
  }
  template <typename stce_9927> void register_variable(stce_9927& stce_9839) {
    stce_9927::register_variable(stce_9839, stce_10373< false>(), this);
  }
  template <typename stce_9927> void register_variable(std::vector<stce_9927>& stce_9839) {
    this->register_variable(stce_9839.data(), stce_9839.size());
  }
  template <typename stce_10375> void register_variable(stce_10375 first, stce_10375 stce_10376) {
    while(first != stce_10376) {
      this->register_variable(*first++);
    }
  }
  template <typename stce_10375, typename stce_10377>
    void register_variable(stce_10375 stce_10378, stce_10375 stce_10379, stce_10377 stce_10380) {
    while(stce_10378 != stce_10379) {
      *stce_10380 = *stce_10378++;
      this->register_variable(*stce_10380++);
    }
  }
  template <typename stce_9927> void register_variable(stce_9927* stce_9839, size_t stce_10096) {
    for (size_t stce_9827 = 0; stce_9827 < stce_10096; ++stce_9827) {
      stce_9927::register_variable(stce_9839[stce_9827], stce_10373< false>(), this);
    }
  }
  template<typename stce_9927> void register_output_variable(stce_9927 &stce_9839) { stce_9839 = 1*stce_9839; }
  template<typename stce_9927> void register_output_variable(std::vector<stce_9927>& stce_9839) {
    this->register_output_variable(stce_9839.data(), stce_9839.size());
  }
  template<typename stce_9927> void register_output_variable(stce_9927* stce_9839, size_t stce_10096) {
    for (size_t stce_9827 = 0; stce_9827 < stce_10096; ++stce_9827) {
      this->register_output_variable(stce_9839[stce_9827]);
    }
  }
  template <typename stce_10375> void register_output_variable(stce_10375 first, stce_10375 stce_10376) {
    while(first != stce_10376) {
      this->register_output_variable(*first++);
    }
  }
  template <stce_9796 stce_10381, typename stce_9782, typename stce_10325>
    void register_output_variable(stce_10334<stce_10381, stce_9782, stce_10325> stce_9839) {
    for (size_t stce_9827 = 0; stce_9827 < stce_9839.size(); ++stce_9827) {
      this->register_output_variable(stce_9839[stce_9827]);
    }
  }
  void reset_to(const iterator_t &stce_9859) {
    stce_10259(stce_9859);
    stce_10382(stce_9859);
  }
  void reset() { reset_to(stce_10277::start()); }
  inline iterator_t current() const { return get_position(); }
  inline iterator_t begin() const { return stce_10277::start(); }
  vector_t& adjoints() { return stce_10284; }
  void interpret_adjoint(interpretation_settings stce_10305 = interpretation_settings()) {
    iterator_t stce_9859 = stce_10277::start();
    stce_10306(get_position(), stce_9859, stce_10305);
  }
  void interpret_adjoint_to(const iterator_t &stce_9859) {
    interpretation_settings stce_10305;
    stce_10306(get_position(), stce_9859, stce_10305);
  }
  void interpret_adjoint_from(const iterator_t &stce_9858) {
    iterator_t stce_9859 = stce_10277::start();
    interpretation_settings stce_10305;
    stce_10306(stce_9858, stce_9859, stce_10305);
  }
  void interpret_adjoint_from_to(const iterator_t &stce_9858, const iterator_t &stce_9859) {
    interpretation_settings stce_10305;
    stce_10306(stce_9858, stce_9859, stce_10305);
  }
  void interpret_adjoint_and_reset_to(const iterator_t &stce_9859) {
    iterator_t stce_9858(get_position());
    interpretation_settings stce_10305;
    stce_10305.reset = true;
    stce_10305.stce_10059 = true;
    stce_10306(stce_9858, stce_9859, stce_10305);
    stce_10259(stce_9859);
    stce_10382(stce_9859);
  }
  void interpret_adjoint_and_zero_adjoints_to(const iterator_t &stce_9859) {
    iterator_t stce_9858(get_position());
    interpretation_settings stce_10305;
    stce_10305.reset = false;
    stce_10305.stce_10059 = true;
    stce_10306(stce_9858, stce_9859, stce_10305);
  }
  void interpret_adjoint_and_zero_adjoints_from_to(const iterator_t &stce_9858, const iterator_t &stce_9859) {
    interpretation_settings stce_10305;
    stce_10305.reset = false;
    stce_10305.stce_10059 = true;
    stce_10306(stce_9858, stce_9859, stce_10305);
  }
  void zero_adjoints() {
    iterator_t stce_9859 = stce_10277::start();
    stce_10307(get_position(), stce_9859);
  }
  void zero_adjoints_to(const iterator_t &stce_9859) {
    stce_10307(get_position(), stce_9859);
  }
  void zero_adjoints_from(const iterator_t &stce_9858) {
    iterator_t stce_9859 = stce_10277::start();
    stce_10307(stce_9858, stce_9859);
  }
  void zero_adjoints_from_to(const iterator_t &stce_9858, const iterator_t &stce_9859) {
    stce_10307(stce_9858, stce_9859);
  }
  void stce_10259(const iterator_t &stce_9859) {
    bool alloc = stce_10284.stce_10285(stce_9859.index());
    stce_10284.resize(stce_9859.index()+1);
    if(alloc) { stce_10284.stce_10300(); }
    stce_10277::erase(stce_9859);
  }
public:
  iterator_t get_position() const {
    return stce_10277::current();
  }
  stce_10079 &stce_10383(const DCO_TAPE_INT tape_index) {
    stce_9764 &stce_10076 = _adjoint(tape_index);
    stce_10079 *stce_10384 = reinterpret_cast<stce_10079*>(&stce_10076);
    return stce_10384[stce_10303];
  }
  stce_10079 stce_10383(const DCO_TAPE_INT tape_index) const {
    stce_9764 stce_10076(_adjoint(tape_index));
    stce_10079 *stce_10384 = reinterpret_cast<stce_10079 *>(&stce_10076);
    return stce_10384[stce_10303];
  }
  stce_9764 &_adjoint(DCO_TAPE_INT tape_index) {
    if(!stce_10284.stce_10285(get_position().index()+1))
      stce_10360();
    return stce_10284[tape_index];
  }
  stce_9764 _adjoint(const DCO_TAPE_INT tape_index) const {
    return stce_10284.stce_10285(tape_index)?stce_10284[tape_index]: 0;
  }
  enum stce_10385 {
    size_of_stack = 1,
    size_of_allocated_stack = 2,
    size_of_internal_adjoint_vector = 4,
    size_of_checkpoints = 8,
    stce_10070 = size_of_stack | size_of_internal_adjoint_vector
  };
  mem_long_t stce_10386() const {
    mem_long_t stce_10387 = 0;
    for (size_t stce_9827 = 0; stce_9827 < stce_10365.size(); stce_9827++) {
      stce_10387 += stce_10365[stce_9827].stce_10321()->size_in_byte();
    }
    return stce_10387;
  }
  mem_long_t stce_10071(const int stce_10388 = stce_10070) const {
    mem_long_t stce_10389 = 0;
    if (stce_10388 & size_of_allocated_stack) {
      stce_10389 += static_cast<size_t>(stce_10277::stce_10176());
    }
    else if (stce_10388 & size_of_stack) {
      stce_10389 += static_cast<size_t>(stce_10277::size_in_byte());
    }
    if (stce_10388 & size_of_internal_adjoint_vector) {
      stce_10389 += stce_10284.stce_10295();
    }
    if (stce_10388 & size_of_checkpoints) {
      stce_10389 += stce_10386();
    }
    return stce_10389;
  }
  static tape_t *create(tape_options stce_10164 = tape_options()) {
    tape_t* tape = new tape_t(stce_10164);
    return tape;
  }
  static void remove(tape_t *&tape) {
    if (tape == 0)
      return;
    delete tape;
    tape = 0;
  }
  virtual stce_10079 stce_10074(const DCO_TAPE_INT stce_9806) {
    return stce_10383(stce_9806);
  }
  virtual void stce_10075(const DCO_TAPE_INT stce_9806, const stce_10079 &stce_10076) {
    stce_10383(stce_9806) += stce_10076;
  }
  typedef dco::helper::stce_10078<tape_t> callback_object_t;
  template <typename stce_10390>
    class stce_10391 {
public:
    typedef void (*stce_10392)(tape_t &stce_10393, const interpretation_settings &stce_9894, stce_10390 *stce_10394);
    typedef void (*stce_10395)(tape_t &stce_10393, stce_10390 *stce_10394);
    typedef void (*stce_10396)(stce_10390 *stce_10394);
  };
  class stce_10397 {
public:
    virtual void stce_10398(tape_t &stce_10393, const interpretation_settings &stce_9894, callback_object_t *stce_10394) = 0;
    virtual ~stce_10397() {};
  };
  template <typename stce_10390>
    class stce_10399 : public stce_10397 {
private:
    union {
      typename stce_10391<stce_10390>::stce_10396 stce_10400;
      typename stce_10391<stce_10390>::stce_10395 stce_10401;
      typename stce_10391<stce_10390>::stce_10392 stce_10402;
    } stce_10400;
    int stce_10403;
public:
    stce_10399(typename stce_10391<stce_10390>::stce_10396 stce_10404) : stce_10403(0) {
      stce_10400.stce_10400 = stce_10404;
    }
    stce_10399(typename stce_10391<stce_10390>::stce_10395 stce_10404) : stce_10403(1) {
      stce_10400.stce_10401 = stce_10404;
    }
    stce_10399(typename stce_10391<stce_10390>::stce_10392 stce_10404) : stce_10403(2) {
      stce_10400.stce_10402 = stce_10404;
    }
    void stce_10398(tape_t &stce_10393, const interpretation_settings &stce_9894, callback_object_t *stce_10394) {
      stce_10390 *stce_10405 = static_cast<stce_10390 *>(stce_10394);
      switch (stce_10403) {
      case 0:
        stce_10400.stce_10400(stce_10405);
        break;
      case 1:
        stce_10400.stce_10401(stce_10393, stce_10405);
        break;
      case 2:
        stce_10400.stce_10402(stce_10393, stce_9894, stce_10405);
        break;
      default:
        if (helper::stce_9817) {
          throw dco::exception::create<std::runtime_error>("unknown callback to run.");
        }
        break;
      }
    }
    ~stce_10399() {}
  };
  class stce_10364 {
    callback_object_t *stce_10394;
    stce_10397 *stce_10406;
    iterator_t stce_10407;
public:
    stce_10364(const iterator_t& stce_10152) :
      stce_10394(0), stce_10406(0), stce_10407(stce_10152){
    }
    bool operator<(const stce_10364& stce_9929) const
    {
      return stce_10407 < stce_9929.stce_10407;
    }
    void stce_10366() {
      delete stce_10394;
      if (stce_10406)
        delete stce_10406;
    }
    template<typename stce_10390, typename stce_10408>
      void stce_10409(stce_10408 stce_10404) {
      if (stce_10406)
        throw dco::exception::create<std::runtime_error>(
                                                         "currently not supported to insert external_adjoint_object_bases twice.");
      stce_10406 = new stce_10399<stce_10390>(stce_10404);
    }
    callback_object_t *&stce_10321() {
      return stce_10394;
    }
    callback_object_t * stce_10321() const {
      return stce_10394;
    }
    iterator_t &stce_10320() {
      return stce_10407;
    }
    iterator_t stce_10320() const {
      return stce_10407;
    }
    void stce_10398(tape_t &stce_10393, const interpretation_settings &stce_9894) {
      if (stce_10406)
        stce_10406->stce_10398(stce_10393, stce_9894, stce_10394);
    }
  };
  std::vector<stce_10364> stce_10365;
  size_t number_of_callbacks() {
    return stce_10365.size();
  }
  template <class stce_10410>
   void stce_10411(stce_10410 *stce_10394) {
       stce_10364 stce_10412(stce_10277::current());
       stce_10365.push_back(stce_10412);
       stce_10365.back().stce_10321() = stce_10394;
       stce_10394->set_tape(this);
       typedef typename stce_10410::active_t stce_10413;
       typedef dco::stce_10131::stce_10134<tape_t,stce_10413> stce_10414;
       stce_10414* stce_10130 = stce_10394;
       insert_callback(&dco::stce_10131::stce_10129<tape_t,stce_10413>, stce_10130);
   }
  template<class stce_10410, typename stce_10408>
    void insert_callback(stce_10408 stce_10406, stce_10410 *stce_10370) {
    if (stce_10365.back().stce_10321() == stce_10370) {
      stce_10365.back().template stce_10409<stce_10410>(stce_10406);
      stce_10365.back().stce_10320() = stce_10277::current();
      stce_10373< true>();
    } else {
      throw dco::exception::create<std::runtime_error>("please always insert most recently created external function object.");
    }
  }
  template<class stce_10410, typename FCN_PARAMETERS>
    stce_10410 *create_callback_object(const FCN_PARAMETERS &parameters) {
    const stce_9840<FCN_PARAMETERS> stce_10415(parameters);
    stce_10364 stce_10412(stce_10277::current());
    stce_10365.push_back(stce_10412);
    stce_10410 *stce_10394 = stce_10415.template create<stce_10410>();
    stce_10365.back().stce_10321() = stce_10394;
    stce_10394->set_tape(this);
    return stce_10394;
  }
  template<class stce_10410>
    stce_10410 *create_callback_object() {
    void *stce_9994;
    return create_callback_object<stce_10410>(stce_9994);
  }
  template <class stce_10410>
    inline stce_10410* create_ACM_object() {
    stce_10410* stce_9836 = create_callback_object<stce_10410>();
    typedef typename stce_10410::active_t stce_10413;
    typedef dco::ACM::baseclass<tape_t,stce_10413> stce_10416;
    stce_10416* stce_10130 = stce_9836;
    insert_callback(&dco::ACM::stce_10129<tape_t,stce_10413>, stce_10130);
    return stce_9836;
  }
  template <class stce_10410, typename FCN_PARAMETERS>
    inline stce_10410 *create_ACM_object(const FCN_PARAMETERS &parameters) {
    const stce_9840<FCN_PARAMETERS> stce_10415(parameters);
    stce_10364 stce_10412(stce_10277::current());
    stce_10365.push_back(stce_10412);
    stce_10410 *stce_10394 = stce_10415.template create<stce_10410>();
    stce_10365.back().stce_10321() = stce_10394;
    stce_10394->set_tape(this);
    typedef typename stce_10410::active_t stce_10413;
    typedef dco::ACM::baseclass<tape_t,stce_10413> stce_10416;
    stce_10416* stce_10130 = stce_10394;
    insert_callback(&dco::ACM::stce_10129<tape_t,stce_10413>, stce_10130);
    return stce_10394;
  }
    template <class stce_10410>
    void insert_ACM_object(stce_10410 *stce_10394) {
        stce_10364 stce_10412(stce_10277::current());
        stce_10365.push_back(stce_10412);
        stce_10365.back().stce_10321() = stce_10394;
        stce_10394->set_tape(this);
        typedef typename stce_10410::active_t stce_10413;
        typedef dco::ACM::baseclass<tape_t,stce_10413> stce_10416;
        stce_10416* stce_10130 = stce_10394;
        insert_callback(&dco::ACM::stce_10129<tape_t,stce_10413>, stce_10130);
    }
  void stce_10382(const iterator_t &stce_9859) {
    typedef typename std::vector<stce_10364>::iterator stce_10319;
    stce_10319 stce_10048 = stce_10365.begin();
    for (; stce_10048 != stce_10365.end() && stce_10048->stce_10320().index() < stce_9859.index(); ++stce_10048) {}
    for (stce_10319 stce_10158 = stce_10048; stce_10158 != stce_10365.end(); ++stce_10158)
      stce_10158->stce_10366();
    stce_10365.erase(stce_10048, stce_10365.end());
  }
public:
  void write_to_dot(const std::string &filename = "tape.dot")
  {
    write_to_dot(stce_10277::current(), stce_10277::start(), filename);
  }
  void write_to_dot(iterator_t stce_9858)
  {
    write_to_dot(stce_9858, stce_10277::start());
  }
  void write_to_dot(iterator_t stce_9815, iterator_t stce_9859,
                   const std::string &filename = "tape.dot") {
    stce_9844 stce_9845(filename);
    stce_10314(stce_9815, stce_9859 , stce_10284, stce_10365, stce_9845);
  }
  void write_to_csv(const std::string &filename = "tape.csv")
  {
    write_to_csv(stce_10277::current(), stce_10277::start(), filename);
  }
  void write_to_csv(iterator_t stce_9858)
  {
    write_to_csv(stce_9858, stce_10277::start());
  }
  void write_to_csv(iterator_t stce_9815, iterator_t stce_9859,
                   const std::string &filename = "tape.csv") {
    if(stce_10365.size() != 0)
      throw dco::exception::create<std::runtime_error>
        ("write_to_csv not implemented for callbacks.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 679);
    stce_10322(stce_9815, stce_9859, filename);
  }
    void stce_10417(const iterator_t& stce_9815, stce_10364& stce_10418,
                          const interpretation_settings &stce_10305)
  {
    if(stce_10418.stce_10320() != stce_9815)
      throw dco::exception::create<std::runtime_error>("invoking callback at wrong position", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 688);
    if (stce_10305.reset) {
      stce_10259(stce_9815);
      stce_10373< true>();
    }
    int stce_9833 = dco::helper::stce_9832<typename stce_9762::stce_9764>::stce_9833;
    for (DCO_TAPE_INT stce_10419 = 0; stce_10419 < stce_9833; stce_10419++) {
      stce_10303 = stce_10419;
      stce_10418.stce_10321()->stce_10081 = this;
      stce_10418.stce_10398(*this, stce_10305);
    }
  }
  template<class stce_10420>
    void stce_10417(stce_10420 &stce_10384, const iterator_t& stce_9815, stce_10364& stce_10418,
                          const interpretation_settings &stce_10305)
  {
    if(stce_10418.stce_10320() != stce_9815)
      throw dco::exception::create<std::runtime_error>("invoking callback at wrong position", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 712);
    if (stce_10305.reset) {
      stce_10259(stce_9815);
      stce_10373< true>();
    }
    int stce_9833 = dco::helper::stce_9832<typename stce_9762::stce_9764>::stce_9833;
    for (DCO_TAPE_INT stce_10419 = 0; stce_10419 < stce_9833; stce_10419++) {
      stce_10303 = stce_10419;
      stce_10418.stce_10321()->stce_10081 = &stce_10384;
      stce_10418.stce_10398(*this, stce_10305);
    }
  }
  struct stce_10421 {
    tape* stce_9948;
    bool stce_10359;
    stce_10421(tape* stce_9901) : stce_9948(stce_9901), stce_10359(stce_9901->sparse_interpret()) {}
    ~stce_10421() { stce_9948->sparse_interpret() = stce_10359; }
  };
  void stce_10306(iterator_t stce_9815, const iterator_t &stce_9859,
                                   const interpretation_settings &stce_10305) {
    stce_10421 stce_10422(this);
    if (stce_9815 > stce_10277::current()) {
      throw dco::exception::create<std::runtime_error>
        ("you try to use a tape position outside of the current tape.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 747);
    }
    if (stce_9859 > stce_9815) {
      throw dco::exception::create<std::runtime_error>
        ("adjoint interpretation: from < to.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 752);
    }
    stce_10360();
    typedef typename std::vector<stce_10364>::reverse_iterator stce_10319;
    for (stce_10319 stce_10158 = stce_10365.rbegin();
         stce_10158 != stce_10365.rend() && stce_9859 <= stce_10158->stce_10320();
         ++stce_10158) {
      const iterator_t& stce_10423 = stce_10158->stce_10320();
      if(stce_9815 < stce_10423) { continue; }
      stce_10277::stce_10177(stce_9815);
      stce_10361(stce_9815, stce_10423, stce_10284);
      stce_9815 = stce_10423;
      std::size_t stce_10424 = std::distance(stce_10365.rbegin(), stce_10158);
      if(stce_10305.stce_10059)
        stce_10284.stce_10297(stce_9815.index(), stce_10423.index());
      stce_10277::stce_10178();
      stce_10417(stce_9815, *stce_10158, stce_10305);
      stce_10158 = stce_10365.rbegin() + stce_10424;
      if (stce_10305.reset) {
        stce_10259(stce_9815);
        stce_10158->stce_10366();
        stce_10365.erase( --(stce_10158.base()) );
      }
    }
    stce_10277::stce_10177(stce_9815);
    stce_10361(stce_9815, stce_9859, stce_10284);
    if(stce_10305.stce_10059)
      stce_10284.stce_10297(stce_9815.index(), stce_9859.index());
    stce_10277::stce_10178();
  }
  template<class stce_10420>
    void stce_10310(stce_10420 &stce_10384,
                                       iterator_t stce_9815, const iterator_t &stce_9859,
                                       const interpretation_settings &stce_10305)
  {
    if (stce_9859 > stce_9815) {
      throw dco::exception::create<std::runtime_error>("adjoint interpretation: from < to.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 817);
    }
    if (stce_9815 > stce_10277::current()) {
      throw dco::exception::create<std::runtime_error>("you try to use a tape position outside of the current tape.", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape.hpp", 820);
    }
    typedef typename std::vector<stce_10364>::reverse_iterator stce_10319;
    for (stce_10319 stce_10158 = stce_10365.rbegin(); stce_10158 != stce_10365.rend() && stce_9859 <= stce_10158->stce_10320();
         ++stce_10158) {
      const iterator_t& stce_10423 = stce_10158->stce_10320();
      if(stce_9815 < stce_10423)
        continue;
      stce_10277::stce_10177(stce_9815);
      if(stce_10305.stce_10059)
        stce_10362(stce_9815, stce_10423, stce_10384);
      else {
        stce_10361(stce_9815, stce_10423, stce_10384);
      }
      stce_9815 = stce_10423;
      std::size_t stce_10424 = std::distance(stce_10365.rbegin(), stce_10158);
      stce_10417(stce_10384, stce_9815, *stce_10158, stce_10305);
      if (stce_9815 != stce_10277::start()) {
        --stce_9815;
      }
      stce_10158 = stce_10365.rbegin() + stce_10424;
      if (stce_10305.reset) {
        stce_10259(stce_9815);
        stce_10158->stce_10366();
        stce_10365.erase( --(stce_10158.base()) );
      }
    }
    stce_10277::stce_10177(stce_9815);
    if(stce_10305.stce_10059) { stce_10362(stce_9815, stce_9859, stce_10384); }
    else { stce_10361(stce_9815, stce_9859, stce_10384); }
    if (stce_10305.reset) { stce_10259(stce_9815); }
    stce_10277::stce_10178();
  }
};
}
}
       
namespace dco {
}
namespace dco {
  namespace stce_10425 {
    template<class stce_9771>struct stce_10426{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { return static_cast<stce_9771>(stce_9908 + stce_9909); } template<class stce_10427,class stce_10026>static inline const stce_9771 stce_10034(const stce_9771 &_value , const stce_10427 &stce_9908 , const stce_10026 &stce_9909 ) { (void) _value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(static_cast<stce_9771>(1.0)); } template<class stce_10427,class stce_10026>static inline const stce_9771 stce_10035(const stce_9771 &_value , const stce_10427 &stce_9908 , const stce_10026 &stce_9909 ) { (void) _value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(static_cast<stce_9771>(1.0)); } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << stce_9953 << "+" << stce_9954; return stce_9836.str(); } static std::string stce_9910() { return "static_cast<DCO_TAPE_REAL>(1.0)"; } static std::string stce_9912() { return "static_cast<DCO_TAPE_REAL>(1.0)"; } };
    template<class stce_9771>struct stce_10428{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { return static_cast<stce_9771>(stce_9908 - stce_9909); } template<class stce_10427,class stce_10026>static inline const stce_9771 stce_10034(const stce_9771 &_value , const stce_10427 &stce_9908 , const stce_10026 &stce_9909 ) { (void) _value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(static_cast<stce_9771>(1.0)); } template<class stce_10427,class stce_10026>static inline const stce_9771 stce_10035(const stce_9771 &_value , const stce_10427 &stce_9908 , const stce_10026 &stce_9909 ) { (void) _value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(static_cast<stce_9771>(-1.0)); } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << stce_9953 << "-" << stce_9954; return stce_9836.str(); } static std::string stce_9910() { return "static_cast<DCO_TAPE_REAL>(1.0)"; } static std::string stce_9912() { return "static_cast<DCO_TAPE_REAL>(-1.0)"; } };
    template<class stce_9771>struct stce_10429{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { return static_cast<stce_9771>(stce_9908 * stce_9909); } template<class stce_10427,class stce_10026>static inline const stce_9771 stce_10034(const stce_9771 &_value , const stce_10427 &stce_9908 , const stce_10026 &stce_9909 ) { (void) _value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(stce_9909); } template<class stce_10427,class stce_10026>static inline const stce_9771 stce_10035(const stce_9771 &_value , const stce_10427 &stce_9908 , const stce_10026 &stce_9909 ) { (void) _value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(stce_9908); } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << stce_9953 << "*" << stce_9954; return stce_9836.str(); } static std::string stce_9910() { return "arg2"; } static std::string stce_9912() { return "arg1"; } };
    template<class stce_9771>struct stce_10430{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { return static_cast<stce_9771>(stce_9908 / stce_9909); } template<class stce_10427,class stce_10026>static inline const stce_9771 stce_10034(const stce_9771 &_value , const stce_10427 &stce_9908 , const stce_10026 &stce_9909 ) { (void) _value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) / stce_9909); } template<class stce_10427,class stce_10026>static inline const stce_9771 stce_10035(const stce_9771 &_value , const stce_10427 &stce_9908 , const stce_10026 &stce_9909 ) { (void) _value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(-_value / stce_9909); } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << stce_9953 << "/" << stce_9954; return stce_9836.str(); } static std::string stce_9910() { return "static_cast<DCO_TAPE_REAL>(1.0) / arg2"; } static std::string stce_9912() { return "-_value / arg2"; } };
    template<class stce_9771>struct stce_10431{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { return stce_9908 + stce_9909; } template<class stce_10427>static inline const stce_9771 stce_10034(const stce_9771 &_value , const stce_10427 &stce_9908 , const stce_9771 &stce_9909 ) { (void)_value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(1.0); } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << stce_9953 << "+" << stce_9954; return stce_9836.str(); } static std::string stce_9910() { return "static_cast<DCO_TAPE_REAL>(1.0)"; } }; template<class stce_9771>struct stce_10432{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { return stce_9908 + stce_9909; } template<class stce_10026>static inline const stce_9771 stce_10035(const stce_9771 &_value , const stce_9771 &stce_9908 , const stce_10026 &stce_9909 ) { (void) _value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(1.0); } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << stce_9953 << "+" << stce_9954; return stce_9836.str(); } static std::string stce_9912() { return "static_cast<DCO_TAPE_REAL>(1.0)"; } };
    template<class stce_9771>struct stce_10433{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { return stce_9908 - stce_9909; } template<class stce_10427>static inline const stce_9771 stce_10034(const stce_9771 &_value , const stce_10427 &stce_9908 , const stce_9771 &stce_9909 ) { (void)_value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(1.0); } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << stce_9953 << "-" << stce_9954; return stce_9836.str(); } static std::string stce_9910() { return "static_cast<DCO_TAPE_REAL>(1.0)"; } }; template<class stce_9771>struct stce_10434{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { return stce_9908 - stce_9909; } template<class stce_10026>static inline const stce_9771 stce_10035(const stce_9771 &_value , const stce_9771 &stce_9908 , const stce_10026 &stce_9909 ) { (void) _value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(-1.0); } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << stce_9953 << "-" << stce_9954; return stce_9836.str(); } static std::string stce_9912() { return "static_cast<DCO_TAPE_REAL>(-1.0)"; } };
    template<class stce_9771>struct stce_10435{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { return stce_9908 * stce_9909; } template<class stce_10427>static inline const stce_9771 stce_10034(const stce_9771 &_value , const stce_10427 &stce_9908 , const stce_9771 &stce_9909 ) { (void)_value; (void)stce_9908; (void)stce_9909; return stce_9909; } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << stce_9953 << "*" << stce_9954; return stce_9836.str(); } static std::string stce_9910() { return "arg2"; } }; template<class stce_9771>struct stce_10436{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { return stce_9908 * stce_9909; } template<class stce_10026>static inline const stce_9771 stce_10035(const stce_9771 &_value , const stce_9771 &stce_9908 , const stce_10026 &stce_9909 ) { (void) _value; (void)stce_9908; (void)stce_9909; return stce_9908; } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << stce_9953 << "*" << stce_9954; return stce_9836.str(); } static std::string stce_9912() { return "arg1"; } };
    template<class stce_9771>struct stce_10437{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { return stce_9908 / stce_9909; } template<class stce_10427>static inline const stce_9771 stce_10034(const stce_9771 &_value , const stce_10427 &stce_9908 , const stce_9771 &stce_9909 ) { (void)_value; (void)stce_9908; (void)stce_9909; return static_cast<stce_9771>(1.0) / stce_9909; } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << stce_9953 << "/" << stce_9954; return stce_9836.str(); } static std::string stce_9910() { return "static_cast<DCO_TAPE_REAL>(1.0) / arg2"; } }; template<class stce_9771>struct stce_10438{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { return stce_9908 / stce_9909; } template<class stce_10026>static inline const stce_9771 stce_10035(const stce_9771 &_value , const stce_9771 &stce_9908 , const stce_10026 &stce_9909 ) { (void) _value; (void)stce_9908; (void)stce_9909; return -(stce_9908 / stce_9909) / stce_9909; } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << stce_9953 << "/" << stce_9954; return stce_9836.str(); } static std::string stce_9912() { return "-(arg1 / arg2) / arg2"; } };
    template<class stce_9771>struct stce_10439{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using std::sin; return sin(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using std::cos; (void)_value; return static_cast<stce_9771>(cos(stce_9839)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "sin" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("cos(x)"); return stce_9836; } };
    template<class stce_9771>struct stce_10440{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using std::cos; return cos(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using std::sin; (void)_value; return static_cast<stce_9771>(-sin(stce_9839)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "cos" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("-sin(x)"); return stce_9836; } };
    template<class stce_9771>struct stce_10441{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using std::tan; return tan(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using std::tan; (void)_value; return static_cast<stce_9771>((static_cast<stce_9771>(1.0) + (tan(stce_9839)*tan(stce_9839)))); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "tan" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("(static_cast<DCO_TAPE_REAL>(1.0) + (tan(x)*tan(x)))"); return stce_9836; } };
    template<class stce_9771>struct stce_10442{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using std::cosh; return cosh(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using std::sinh; (void)_value; return static_cast<stce_9771>(sinh(stce_9839)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "cosh" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("sinh(x)"); return stce_9836; } };
    template<class stce_9771>struct stce_10443{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using std::sinh; return sinh(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using std::cosh; (void)_value; return static_cast<stce_9771>(cosh(stce_9839)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "sinh" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("cosh(x)"); return stce_9836; } };
    template<class stce_9771>struct stce_10444{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using std::tanh; return tanh(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using std::tanh; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) - tanh(stce_9839)*tanh(stce_9839)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "tanh" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("static_cast<DCO_TAPE_REAL>(1.0) - tanh(x)*tanh(x)"); return stce_9836; } };
    template<class stce_9771>struct stce_10445{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using std::asin; return asin(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using std::sqrt; (void)_value; return static_cast<stce_9771>(1 / sqrt(static_cast<stce_9771>(1.0) - stce_9839*stce_9839)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "asin" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("1 / sqrt(static_cast<DCO_TAPE_REAL>(1.0) - x*x)"); return stce_9836; } };
    template<class stce_9771>struct stce_10446{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using std::acos; return acos(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using std::sqrt; (void)_value; return static_cast<stce_9771>(-1 / sqrt(static_cast<stce_9771>(1.0) - stce_9839*stce_9839)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "acos" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("-1 / sqrt(static_cast<DCO_TAPE_REAL>(1.0) - x*x)"); return stce_9836; } };
    template<class stce_9771>struct stce_10447{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using std::atan; return atan(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using std::atan; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) / (static_cast<stce_9771>(1.0) + stce_9839*stce_9839)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "atan" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("static_cast<DCO_TAPE_REAL>(1.0) / (static_cast<DCO_TAPE_REAL>(1.0) + x*x)"); return stce_9836; } };
    template<class stce_9771>struct stce_10448{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using std::exp; return exp(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using std::exp; (void)_value; return static_cast<stce_9771>(exp(stce_9839)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "exp" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("exp(x)"); return stce_9836; } };
    template<class stce_9771>struct stce_10449{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using std::sqrt; return sqrt(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using std::sqrt; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) / (static_cast<stce_9771>(2.0) * sqrt(stce_9839 + ((HAS_SQRT_AVOID_DIVZERO && stce_9839==0) ? DCO_SQRT_EPS : 0.0)))); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "sqrt" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("static_cast<DCO_TAPE_REAL>(1.0) / (static_cast<DCO_TAPE_REAL>(2.0) * sqrt(x + ((HAS_SQRT_AVOID_DIVZERO && x==0) ? DCO_SQRT_EPS : 0.0)))"); return stce_9836; } };
    template<class stce_9771>struct stce_10450{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using std::log; return log(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using std::log; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) / stce_9839); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "log" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("static_cast<DCO_TAPE_REAL>(1.0) / x"); return stce_9836; } };
    template<class stce_9771>struct stce_10451{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using ::erf; return erf(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using ::exp; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(2.0) / sqrt(static_cast<stce_9771>(3.14159265358979323846264338327950288)) * exp(-stce_9839 * stce_9839)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "erf" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("static_cast<DCO_TAPE_REAL>(2.0) / sqrt(static_cast<DCO_TAPE_REAL>(3.14159265358979323846264338327950288)) * exp(-x * x)"); return stce_9836; } };
    template<class stce_9771>struct stce_10452{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using ::erfc; return erfc(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using ::exp; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(-2.0) / sqrt(static_cast<stce_9771>(3.14159265358979323846264338327950288)) * exp(-stce_9839 * stce_9839)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "erfc" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("static_cast<DCO_TAPE_REAL>(-2.0) / sqrt(static_cast<DCO_TAPE_REAL>(3.14159265358979323846264338327950288)) * exp(-x * x)"); return stce_9836; } };
    template<class stce_9771>struct stce_10453{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using ::asinh; return asinh(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using ::sqrt; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.) / sqrt(static_cast<stce_9771>(1.) + (stce_9839*stce_9839))); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "asinh" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("static_cast<DCO_TAPE_REAL>(1.) / sqrt(static_cast<DCO_TAPE_REAL>(1.) + (x*x))"); return stce_9836; } };
    template<class stce_9771>struct stce_10454{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using ::acosh; return acosh(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using ::sqrt; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.) / sqrt((stce_9839*stce_9839) - static_cast<stce_9771>(1.))); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "acosh" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("static_cast<DCO_TAPE_REAL>(1.) / sqrt((x*x) - static_cast<DCO_TAPE_REAL>(1.))"); return stce_9836; } };
    template<class stce_9771>struct stce_10455{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using ::atanh; return atanh(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using ::atanh; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.) / (static_cast<stce_9771>(1.) - (stce_9839*stce_9839))); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "atanh" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("static_cast<DCO_TAPE_REAL>(1.) / (static_cast<DCO_TAPE_REAL>(1.) - (x*x))"); return stce_9836; } };
    template<class stce_9771>struct stce_10456{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using ::expm1; return expm1(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using ::exp; (void)_value; return static_cast<stce_9771>(exp(stce_9839)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "expm1" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("exp(x)"); return stce_9836; } };
    template<class stce_9771>struct stce_10457{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using ::log1p; return log1p(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using ::log1p; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) / (stce_9839 + 1)); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "log1p" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("static_cast<DCO_TAPE_REAL>(1.0) / (x + 1)"); return stce_9836; } };
    template<class stce_9771>struct stce_10458{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &arg) { using ::log10; return log10(arg); } template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9839) { using ::log10; (void)_value; return static_cast<stce_9771>(static_cast<stce_9771>(1.0) / (stce_9839*static_cast<stce_9771>(log(10)))); } static std::string stce_10032(const std::string& stce_9839) { std::stringstream stce_9836; stce_9836 << "log10" << "("<< stce_9839 << ")"; return stce_9836.str(); } static std::string stce_10033() { std::string stce_9836 = std::string("static_cast<DCO_TAPE_REAL>(1.0) / (x*static_cast<DCO_TAPE_REAL>(log(10)))"); return stce_9836; } };
    template<class stce_9771>struct stce_10459 {
      template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &stce_9908) {
        return -stce_9908;
      }
      template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9908 ) {
        (void)_value;
        (void)stce_9908;
        return static_cast<stce_9771>(-1.0);
      }
      static std::string stce_10032(const std::string& stce_9839) {
        return "-1 *" + stce_9839;
      }
      static std::string stce_10033() {
        return "-1";
      }
    };
    template<class stce_9771>struct stce_10460 {
      template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &stce_9908) {
        return stce_9908;
      }
      template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9908 ) {
        (void)_value;
        (void)stce_9908;
        return static_cast<stce_9771>(1.0);
      }
      static std::string stce_10032(const std::string& stce_9839) {
        return stce_9839;
      }
      static std::string stce_10033() {
        return "1";
      }
    };
    using ::fabs;
    template<class stce_9771>struct stce_10461 {
      template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &stce_9908) {
        return fabs(stce_9908);
      }
      template<class stce_9782>static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9908 ) {
        (void) _value;
        if (stce_9908 < 0) return static_cast<stce_9771>(-1.0);
        else return static_cast<stce_9771>(1.0);
      }
      static std::string stce_10032(const std::string& stce_9839) {
        return "fabs(" + stce_9839 + ")";
      }
      static std::string stce_10033() {
        return "((x >= 0) ? 1 : -1 )";
      }
    };
    using std::abs;
    template<class stce_9771> struct stce_10462 {
      template<class stce_9782> static inline const stce_9771 stce_10030(const stce_9782 &stce_9908) {
        return abs(stce_9908);
      }
      template<class stce_9782> static inline const stce_9771 stce_10031(const stce_9771 &_value , const stce_9782 &stce_9908 ) {
        (void) _value;
        if (stce_9908 < 0) return static_cast<stce_9771>(-1.0);
        else return static_cast<stce_9771>(1.0);
      }
      static std::string stce_10032(const std::string& stce_9839) {
        return "abs(" + stce_9839 + ")";
      }
      static std::string stce_10033() {
        return "((x >= 0) ? 1 : -1 )";
      }
    };
    using std::atan2; template<class stce_9771>struct stce_10463{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { (void) stce_9908; (void) stce_9909; return atan2(stce_9908,stce_9909); } template<class stce_10427,class stce_10026>static inline const stce_9771 stce_10034(stce_9771 const& _value , const stce_10427 &stce_9908, const stce_10026 &stce_9909) { (void) _value; (void) stce_9908; (void) stce_9909; return stce_9909 / (stce_9909 * stce_9909 + stce_9908 * stce_9908); } template<class stce_10427,class stce_10026>static inline const stce_9771 stce_10035(stce_9771 const& _value , const stce_10427 &stce_9908, const stce_10026 &stce_9909) { (void) _value; (void) stce_9908; (void) stce_9909; return -stce_9908 / (stce_9909 * stce_9909 + stce_9908 * stce_9908); } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << "atan2" << "(" << stce_9953 << "," << stce_9954 <<")"; return stce_9836.str(); } static std::string stce_9910() { return "arg2 / (arg2 * arg2 + arg1 * arg1)"; } static std::string stce_9912() { return "-arg1 / (arg2 * arg2 + arg1 * arg1)"; } }; template<class stce_9771>struct stce_10464{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &stce_9908, const stce_9771 &stce_9909) { (void) stce_9908; (void) stce_9909; return atan2(stce_9908,stce_9909); } template<class stce_9782>static inline const stce_9771 stce_10034(stce_9771 const& _value , const stce_9782 &stce_9908, const stce_9771 &stce_9909) { (void) _value; (void) stce_9908; (void) stce_9909; return stce_9909 / (stce_9909 * stce_9909 + stce_9908 * stce_9908); } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << "atan2" << "(" << stce_9953 << "," << stce_9954 <<")"; return stce_9836.str(); } static std::string stce_9910() { return "arg2 / (arg2 * arg2 + arg1 * arg1)"; } }; template<class stce_9771>struct stce_10465{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9771 &stce_9908, const stce_9782 &stce_9909) { (void) stce_9908; (void) stce_9909; return atan2(stce_9908,stce_9909); } template<class stce_9782>static inline const stce_9771 stce_10035(stce_9771 const& _value , const stce_9771 &stce_9908, const stce_9782 &stce_9909) { (void) _value; (void) stce_9908; (void) stce_9909; return -stce_9908 / (stce_9909 * stce_9909 + stce_9908 *stce_9908); } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << "atan2" << "(" << stce_9953 << "," << stce_9954 <<")"; return stce_9836.str(); } static std::string stce_9912() { return "-arg1 / (arg2 * arg2 + arg1 *arg1)"; } };
    template<class stce_9771>
    struct stce_10466 {
      template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) {
        using std::pow;
        return pow(stce_9908, stce_9909);
      }
      template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10034(stce_9771 const& _value , const stce_10427 &stce_9908, const stce_10026 &stce_9909) {
        using std::pow;
        (void)_value;
        return stce_9909 * pow(stce_9908, stce_9909 - static_cast<stce_9771>(1.0));
      }
      template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10035(stce_9771 const& _value, const stce_10427 &stce_9908, const stce_10026&) {
        using std::pow;
        using std::log;
        return log(stce_9908) * _value;
      }
      static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) {
        return "pow(" + stce_9953 + "," + stce_9954 + ")";
      }
      static std::string stce_9910() {
        return "arg2 * pow(arg1, arg2 - 1.0)";
      }
      static std::string stce_9912() {
        return "log(arg1) * pow(arg1, arg2)";
      }
    };
    template<class stce_9771>
    struct stce_10467 {
      template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &stce_9908, const stce_9771 &stce_9909) {
        using std::pow;
        return pow(stce_9908, stce_9909);
      }
      template<class stce_9782>static inline const stce_9771 stce_10034(stce_9771 const& _value , const stce_9782 &stce_9908, const stce_9771 &stce_9909) {
        using std::pow;
        (void) _value;
        return stce_9909 * pow(stce_9908, stce_9909 - static_cast<stce_9771>(1.0));
      }
      static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) {
        return "pow(" + stce_9953 + "," + stce_9954 + ")";
      }
      static std::string stce_9910() {
        return "arg2 * pow(arg1, arg2 - 1.0)";
      }
    };
    template<class stce_9771>
    struct stce_10468 {
      template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9771 &stce_9908, const stce_9782 &stce_9909) {
        using std::pow;
        return pow(stce_9908, stce_9909);
      }
      template<class stce_9782>static inline const stce_9771 stce_10035(stce_9771 const& _value , const stce_9771 &stce_9908, const stce_9782&) {
        using std::pow;
        using std::log;
        return log(stce_9908) * _value;
      }
      static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) {
        return "pow(" + stce_9953 + "," + stce_9954 + ")";
      }
      static std::string stce_9910() {
        return "log(arg1) * pow(arg1, arg2)";
      }
    };
    using ::hypot; template<class stce_9771>struct stce_10469{ template<class stce_10427, class stce_10026>static inline const stce_9771 stce_10030(const stce_10427 &stce_9908, const stce_10026 &stce_9909) { (void) stce_9908; (void) stce_9909; return hypot(stce_9908,stce_9909); } template<class stce_10427,class stce_10026>static inline const stce_9771 stce_10034(stce_9771 const& _value , const stce_10427 &stce_9908, const stce_10026 &stce_9909) { (void) _value; (void) stce_9908; (void) stce_9909; return stce_9908 / _value; } template<class stce_10427,class stce_10026>static inline const stce_9771 stce_10035(stce_9771 const& _value , const stce_10427 &stce_9908, const stce_10026 &stce_9909) { (void) _value; (void) stce_9908; (void) stce_9909; return stce_9909 / _value; } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << "hypot" << "(" << stce_9953 << "," << stce_9954 <<")"; return stce_9836.str(); } static std::string stce_9910() { return "arg1 / _value"; } static std::string stce_9912() { return "arg2 / _value"; } }; template<class stce_9771>struct stce_10470{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9782 &stce_9908, const stce_9771 &stce_9909) { (void) stce_9908; (void) stce_9909; return hypot(stce_9908,stce_9909); } template<class stce_9782>static inline const stce_9771 stce_10034(stce_9771 const& _value , const stce_9782 &stce_9908, const stce_9771 &stce_9909) { (void) _value; (void) stce_9908; (void) stce_9909; return stce_9908 / _value; } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << "hypot" << "(" << stce_9953 << "," << stce_9954 <<")"; return stce_9836.str(); } static std::string stce_9910() { return "arg1 / _value"; } }; template<class stce_9771>struct stce_10471{ template<class stce_9782>static inline const stce_9771 stce_10030(const stce_9771 &stce_9908, const stce_9782 &stce_9909) { (void) stce_9908; (void) stce_9909; return hypot(stce_9908,stce_9909); } template<class stce_9782>static inline const stce_9771 stce_10035(stce_9771 const& _value , const stce_9771 &stce_9908, const stce_9782 &stce_9909) { (void) _value; (void) stce_9908; (void) stce_9909; return stce_9909 / _value; } static std::string stce_10032(const std::string& stce_9953,const std::string& stce_9954) { std::stringstream stce_9836; stce_9836 << "hypot" << "(" << stce_9953 << "," << stce_9954 <<")"; return stce_9836.str(); } static std::string stce_9912() { return "arg2 / _value"; } };
  }
}
namespace dco {
namespace internal {
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10459<stce_9771> > operator -( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10459<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10459<stce_9771> > operator -( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10459<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10459<stce_9771> > operator -( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10459<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10459<stce_9771> > operator -( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10459<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10459<stce_9771> > operator -( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10459<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10460<stce_9771> > operator +( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10460<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10460<stce_9771> > operator +( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10460<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10460<stce_9771> > operator +( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10460<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10460<stce_9771> > operator +( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10460<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10460<stce_9771> > operator +( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10460<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10439<stce_9771> > sin( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10439<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10439<stce_9771> > sin( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10439<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10439<stce_9771> > sin( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10439<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10439<stce_9771> > sin( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10439<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10439<stce_9771> > sin( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10439<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10440<stce_9771> > cos( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10440<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10440<stce_9771> > cos( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10440<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10440<stce_9771> > cos( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10440<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10440<stce_9771> > cos( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10440<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10440<stce_9771> > cos( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10440<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10441<stce_9771> > tan( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10441<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10441<stce_9771> > tan( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10441<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10441<stce_9771> > tan( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10441<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10441<stce_9771> > tan( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10441<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10441<stce_9771> > tan( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10441<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10442<stce_9771> > cosh( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10442<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10442<stce_9771> > cosh( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10442<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10442<stce_9771> > cosh( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10442<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10442<stce_9771> > cosh( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10442<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10442<stce_9771> > cosh( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10442<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10443<stce_9771> > sinh( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10443<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10443<stce_9771> > sinh( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10443<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10443<stce_9771> > sinh( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10443<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10443<stce_9771> > sinh( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10443<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10443<stce_9771> > sinh( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10443<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10445<stce_9771> > asin( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10445<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10445<stce_9771> > asin( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10445<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10445<stce_9771> > asin( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10445<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10445<stce_9771> > asin( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10445<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10445<stce_9771> > asin( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10445<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10446<stce_9771> > acos( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10446<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10446<stce_9771> > acos( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10446<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10446<stce_9771> > acos( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10446<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10446<stce_9771> > acos( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10446<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10446<stce_9771> > acos( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10446<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10448<stce_9771> > exp( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10448<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10448<stce_9771> > exp( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10448<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10448<stce_9771> > exp( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10448<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10448<stce_9771> > exp( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10448<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10448<stce_9771> > exp( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10448<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10447<stce_9771> > atan( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10447<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10447<stce_9771> > atan( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10447<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10447<stce_9771> > atan( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10447<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10447<stce_9771> > atan( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10447<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10447<stce_9771> > atan( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10447<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10444<stce_9771> > tanh( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10444<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10444<stce_9771> > tanh( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10444<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10444<stce_9771> > tanh( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10444<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10444<stce_9771> > tanh( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10444<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10444<stce_9771> > tanh( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10444<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10449<stce_9771> > sqrt( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10449<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10449<stce_9771> > sqrt( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10449<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10449<stce_9771> > sqrt( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10449<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10449<stce_9771> > sqrt( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10449<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10449<stce_9771> > sqrt( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10449<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10450<stce_9771> > log( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10450<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10450<stce_9771> > log( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10450<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10450<stce_9771> > log( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10450<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10450<stce_9771> > log( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10450<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10450<stce_9771> > log( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10450<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10451<stce_9771> > erf( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10451<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10451<stce_9771> > erf( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10451<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10451<stce_9771> > erf( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10451<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10451<stce_9771> > erf( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10451<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10451<stce_9771> > erf( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10451<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10452<stce_9771> > erfc( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10452<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10452<stce_9771> > erfc( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10452<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10452<stce_9771> > erfc( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10452<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10452<stce_9771> > erfc( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10452<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10452<stce_9771> > erfc( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10452<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10456<stce_9771> > expm1( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10456<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10456<stce_9771> > expm1( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10456<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10456<stce_9771> > expm1( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10456<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10456<stce_9771> > expm1( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10456<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10456<stce_9771> > expm1( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10456<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10453<stce_9771> > asinh( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10453<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10453<stce_9771> > asinh( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10453<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10453<stce_9771> > asinh( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10453<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10453<stce_9771> > asinh( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10453<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10453<stce_9771> > asinh( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10453<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10454<stce_9771> > acosh( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10454<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10454<stce_9771> > acosh( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10454<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10454<stce_9771> > acosh( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10454<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10454<stce_9771> > acosh( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10454<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10454<stce_9771> > acosh( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10454<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10455<stce_9771> > atanh( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10455<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10455<stce_9771> > atanh( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10455<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10455<stce_9771> > atanh( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10455<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10455<stce_9771> > atanh( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10455<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10455<stce_9771> > atanh( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10455<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10457<stce_9771> > log1p( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10457<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10457<stce_9771> > log1p( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10457<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10457<stce_9771> > log1p( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10457<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10457<stce_9771> > log1p( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10457<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10457<stce_9771> > log1p( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10457<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10458<stce_9771> > log10( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10458<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10458<stce_9771> > log10( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10458<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10458<stce_9771> > log10( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10458<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10458<stce_9771> > log10( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10458<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10458<stce_9771> > log10( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10458<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10461<stce_9771> > fabs( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10461<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10461<stce_9771> > fabs( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10461<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10461<stce_9771> > fabs( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10461<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10461<stce_9771> > fabs( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10461<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10461<stce_9771> > fabs( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10461<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10462<stce_9771> > abs( const dco::internal::active_type<stce_9771, stce_9935> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::stce_10425::stce_10462<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10462<stce_9771> > abs( const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::stce_10425::stce_10462<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10462<stce_9771> > abs( const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::stce_10425::stce_10462<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10462<stce_9771> > abs( const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::stce_10425::stce_10462<stce_9771> >(stce_9953); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10462<stce_9771> > abs( const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953) { return dco::internal::stce_9776<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::stce_10425::stce_10462<stce_9771> >(stce_9953); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> > operator + (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10426<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10431<stce_9771> > operator + (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10431<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10432<stce_9771> > operator + (const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10432<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10431<stce_9771> > operator + (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10431<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10432<stce_9771> > operator + (const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10432<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10431<stce_9771> > operator + (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10431<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10432<stce_9771> > operator + (const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10432<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10431<stce_9771> > operator + (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10431<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10432<stce_9771> > operator + (const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10432<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10431<stce_9771> > operator + (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10431<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10432<stce_9771> > operator + (const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10432<stce_9771> >(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> > operator - (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10428<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10433<stce_9771> > operator - (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10433<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10434<stce_9771> > operator - (const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10434<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10433<stce_9771> > operator - (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10433<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10434<stce_9771> > operator - (const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10434<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10433<stce_9771> > operator - (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10433<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10434<stce_9771> > operator - (const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10434<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10433<stce_9771> > operator - (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10433<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10434<stce_9771> > operator - (const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10434<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10433<stce_9771> > operator - (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10433<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10434<stce_9771> > operator - (const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10434<stce_9771> >(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> > operator * (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10429<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10435<stce_9771> > operator * (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10435<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10436<stce_9771> > operator * (const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10436<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10435<stce_9771> > operator * (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10435<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10436<stce_9771> > operator * (const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10436<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10435<stce_9771> > operator * (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10435<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10436<stce_9771> > operator * (const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10436<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10435<stce_9771> > operator * (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10435<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10436<stce_9771> > operator * (const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10436<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10435<stce_9771> > operator * (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10435<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10436<stce_9771> > operator * (const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10436<stce_9771> >(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> > operator / (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10430<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10437<stce_9771> > operator / (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10437<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10438<stce_9771> > operator / (const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10438<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10437<stce_9771> > operator / (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10437<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10438<stce_9771> > operator / (const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10438<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10437<stce_9771> > operator / (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10437<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10438<stce_9771> > operator / (const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10438<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10437<stce_9771> > operator / (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10437<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10438<stce_9771> > operator / (const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10438<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10437<stce_9771> > operator / (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10437<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10438<stce_9771> > operator / (const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10438<stce_9771> >(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> > atan2 (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10463<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10464<stce_9771> > atan2 (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10464<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10465<stce_9771> > atan2 (const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10465<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10464<stce_9771> > atan2 (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10464<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10465<stce_9771> > atan2 (const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10465<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10464<stce_9771> > atan2 (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10464<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10465<stce_9771> > atan2 (const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10465<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10464<stce_9771> > atan2 (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10464<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10465<stce_9771> > atan2 (const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10465<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10464<stce_9771> > atan2 (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10464<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10465<stce_9771> > atan2 (const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10465<stce_9771> >(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> > pow (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10466<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10467<stce_9771> > pow (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10467<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10468<stce_9771> > pow (const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10468<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10467<stce_9771> > pow (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10467<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10468<stce_9771> > pow (const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10468<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10467<stce_9771> > pow (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10467<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10468<stce_9771> > pow (const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10468<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10467<stce_9771> > pow (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10467<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10468<stce_9771> > pow (const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10468<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10467<stce_9771> > pow (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10467<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10468<stce_9771> > pow (const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10468<stce_9771> >(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::active_type<stce_9771, stce_9935>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::active_type<stce_9771, stce_10476> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::active_type<stce_9771, stce_10476>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> > hypot (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473> &stce_9954) { return dco::internal::stce_9779<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>,dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>, dco::stce_10425::stce_10469<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10470<stce_9771> > hypot (const dco::internal::active_type<stce_9771, stce_9935> &stce_9953, const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10470<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10471<stce_9771> > hypot (const typename dco::internal::active_type<stce_9771, stce_9935>::VALUE_TYPE &stce_9953, const dco::internal::active_type<stce_9771, stce_9935> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::active_type<stce_9771, stce_9935>, dco::stce_10425::stce_10471<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10470<stce_9771> > hypot (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9953, const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10470<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10471<stce_9771> > hypot (const typename dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>, dco::stce_10425::stce_10471<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10470<stce_9771> > hypot (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10470<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10471<stce_9771> > hypot (const typename dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>, dco::stce_10425::stce_10471<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10470<stce_9771> > hypot (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9953, const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10470<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10471<stce_9771> > hypot (const typename dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>, dco::stce_10425::stce_10471<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10470<stce_9771> > hypot (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9953, const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9954) { return dco::internal::stce_9780<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10470<stce_9771> >(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10471<stce_9771> > hypot (const typename dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::VALUE_TYPE &stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9954) { return dco::internal::stce_9781<stce_9771, dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>, dco::stce_10425::stce_10471<stce_9771> >(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 > static inline bool operator == (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > static inline bool operator == (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator == (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > static inline bool operator == (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > static inline bool operator == (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > static inline bool operator == (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator == (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator == (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator == (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator == (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > static inline bool operator == (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator == (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator == (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator == (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator == (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > static inline bool operator == (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator == (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator == (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator == (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator == (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > static inline bool operator == (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator == (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator == (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator == (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator == (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9952(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > static inline bool operator == (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const typename dco::mode<dco::internal::active_type<stce_9771, stce_9935> >::passive_t& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9955(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > static inline bool operator == (const typename dco::mode<dco::internal::active_type<stce_9771, stce_9935> >::passive_t& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9956(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool operator == (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9955(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool operator == (const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9956(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool operator == (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9955(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool operator == (const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9956(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool operator == (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9955(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool operator == (const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9956(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool operator == (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9955(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool operator == (const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9956(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const stce_9782& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9955(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const stce_9782& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9956(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9955(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const stce_9782& stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9956(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9955(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const stce_9782& stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9956(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9955(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const stce_9782& stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9956(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9955(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator == (const stce_9782& stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9956(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 > static inline bool operator != (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > static inline bool operator != (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator != (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > static inline bool operator != (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > static inline bool operator != (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > static inline bool operator != (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator != (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator != (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator != (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator != (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > static inline bool operator != (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator != (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator != (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator != (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator != (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > static inline bool operator != (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator != (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator != (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator != (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator != (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > static inline bool operator != (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator != (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator != (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator != (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator != (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9957(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > static inline bool operator != (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const typename dco::mode<dco::internal::active_type<stce_9771, stce_9935> >::passive_t& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9958(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > static inline bool operator != (const typename dco::mode<dco::internal::active_type<stce_9771, stce_9935> >::passive_t& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9959(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool operator != (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9958(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool operator != (const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9959(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool operator != (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9958(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool operator != (const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9959(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool operator != (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9958(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool operator != (const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9959(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool operator != (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9958(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool operator != (const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9959(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const stce_9782& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9958(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const stce_9782& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9959(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9958(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const stce_9782& stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9959(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9958(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const stce_9782& stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9959(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9958(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const stce_9782& stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9959(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9958(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator != (const stce_9782& stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9959(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 > static inline bool operator < (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > static inline bool operator < (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator < (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > static inline bool operator < (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > static inline bool operator < (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > static inline bool operator < (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator < (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator < (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator < (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator < (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > static inline bool operator < (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator < (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator < (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator < (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator < (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > static inline bool operator < (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator < (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator < (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator < (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator < (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > static inline bool operator < (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator < (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator < (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator < (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator < (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9960(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > static inline bool operator < (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const typename dco::mode<dco::internal::active_type<stce_9771, stce_9935> >::passive_t& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9961(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > static inline bool operator < (const typename dco::mode<dco::internal::active_type<stce_9771, stce_9935> >::passive_t& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9962(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool operator < (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9961(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool operator < (const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9962(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool operator < (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9961(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool operator < (const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9962(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool operator < (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9961(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool operator < (const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9962(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool operator < (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9961(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool operator < (const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9962(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const stce_9782& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9961(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const stce_9782& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9962(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9961(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const stce_9782& stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9962(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9961(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const stce_9782& stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9962(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9961(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const stce_9782& stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9962(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9961(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator < (const stce_9782& stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9962(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 > static inline bool operator <= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > static inline bool operator <= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator <= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > static inline bool operator <= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > static inline bool operator <= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > static inline bool operator <= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > static inline bool operator <= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > static inline bool operator <= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > static inline bool operator <= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator <= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9963(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > static inline bool operator <= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const typename dco::mode<dco::internal::active_type<stce_9771, stce_9935> >::passive_t& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9964(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > static inline bool operator <= (const typename dco::mode<dco::internal::active_type<stce_9771, stce_9935> >::passive_t& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9965(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool operator <= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9964(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool operator <= (const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9965(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool operator <= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9964(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool operator <= (const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9965(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool operator <= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9964(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool operator <= (const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9965(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool operator <= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9964(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool operator <= (const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9965(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const stce_9782& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9964(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const stce_9782& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9965(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9964(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const stce_9782& stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9965(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9964(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const stce_9782& stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9965(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9964(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const stce_9782& stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9965(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9964(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator <= (const stce_9782& stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9965(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 > static inline bool operator > (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > static inline bool operator > (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator > (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > static inline bool operator > (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > static inline bool operator > (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > static inline bool operator > (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator > (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator > (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator > (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator > (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > static inline bool operator > (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator > (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator > (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator > (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator > (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > static inline bool operator > (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator > (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator > (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator > (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator > (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > static inline bool operator > (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator > (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator > (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator > (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator > (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9966(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > static inline bool operator > (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const typename dco::mode<dco::internal::active_type<stce_9771, stce_9935> >::passive_t& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9967(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > static inline bool operator > (const typename dco::mode<dco::internal::active_type<stce_9771, stce_9935> >::passive_t& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9968(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool operator > (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9967(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool operator > (const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9968(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool operator > (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9967(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool operator > (const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9968(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool operator > (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9967(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool operator > (const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9968(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool operator > (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9967(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool operator > (const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9968(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const stce_9782& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9967(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const stce_9782& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9968(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9967(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const stce_9782& stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9968(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9967(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const stce_9782& stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9968(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9967(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const stce_9782& stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9968(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9967(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator > (const stce_9782& stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9968(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 > static inline bool operator >= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10472, class stce_10473 > static inline bool operator >= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator >= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10474, class stce_10473 > static inline bool operator >= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, class stce_10475, class stce_10473 > static inline bool operator >= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10476 > static inline bool operator >= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10476 > static inline bool operator >= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10476 > static inline bool operator >= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10476 > static inline bool operator >= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::active_type<stce_9771, stce_10476>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10472, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9776<stce_9771, stce_10472, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10475, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9779<stce_9771, stce_10474, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10474, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9780<stce_9771, stce_10474, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, class stce_10475, class stce_10473 > static inline bool operator >= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const dco::internal::stce_9781<stce_9771, stce_10475, stce_10473>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9969(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > static inline bool operator >= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const typename dco::mode<dco::internal::active_type<stce_9771, stce_9935> >::passive_t& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9970(stce_9953,stce_9954); } template<class stce_9771, class stce_9935 > static inline bool operator >= (const typename dco::mode<dco::internal::active_type<stce_9771, stce_9935> >::passive_t& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9971(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool operator >= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9970(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool operator >= (const typename dco::mode<dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9971(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool operator >= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9970(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool operator >= (const typename dco::mode<dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9971(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool operator >= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9970(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool operator >= (const typename dco::mode<dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9971(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool operator >= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> >::passive_t& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9970(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool operator >= (const typename dco::mode<dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> >::passive_t& stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9971(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const dco::internal::active_type<stce_9771, stce_9935>& stce_9953, const stce_9782& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9970(stce_9953,stce_9954); } template<class stce_9771, class stce_9935, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const stce_9782& stce_9953, const dco::internal::active_type<stce_9771, stce_9935>& stce_9954) { return dco::internal::active_type<stce_9771, stce_9935>::data_t::stce_9989::stce_9971(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9970(stce_9953,stce_9954); } template<class stce_9771, class stce_9945, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const stce_9782& stce_9953, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9954) { return dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>::data_t::stce_9989::stce_9971(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9970(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const stce_9782& stce_9953, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>::data_t::stce_9989::stce_9971(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9970(stce_9953,stce_9954); } template<class stce_9771, class stce_9937, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const stce_9782& stce_9953, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9954) { return dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>::data_t::stce_9989::stce_9971(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9953, const stce_9782& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9970(stce_9953,stce_9954); } template<class stce_9771, class stce_9938, class stce_9939, typename stce_9782 > static inline typename dco::dco_type_constructable_from<stce_9782, bool>::type operator >= (const stce_9782& stce_9953, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9954) { return dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>::data_t::stce_9989::stce_9971(stce_9953,stce_9954); }
template<class stce_9771, class stce_9935 >
static inline void stce_10477(dco::internal::active_type<stce_9771, stce_9935> &stce_9839) {
    double stce_10057 = 0;
    get(stce_9839, stce_10057);
    stce_9839 = stce_10057;
}
template<class stce_9771, class stce_9935 >
    static inline std::istream &operator >> (std::istream &in, dco::internal::active_type<stce_9771, stce_9935> &stce_9839) {
      stce_9771 &stce_10057 = stce_9839._value();
      in >> stce_10057;
      return in;
    }
    using std::ceil; using std::floor; using std::isfinite; using std::isnan; using std::isnormal; using std::isinf; template<class stce_9771, class stce_9935 > static inline bool isnan(const dco::internal::active_type<stce_9771, stce_9935>& stce_9839) { return isnan(stce_9839._value()); } template<class stce_9771, class stce_9935 > static inline bool isnormal(const dco::internal::active_type<stce_9771, stce_9935>& stce_9839) { return isnormal(stce_9839._value()); } template<class stce_9771, class stce_9935 > static inline bool isinf(const dco::internal::active_type<stce_9771, stce_9935>& stce_9839) { return isinf(stce_9839._value()); } template<class stce_9771, class stce_9935 > static inline double ceil(const dco::internal::active_type<stce_9771, stce_9935>& stce_9839) { return ceil(stce_9839._value()); } template<class stce_9771, class stce_9935 > static inline double floor(const dco::internal::active_type<stce_9771, stce_9935>& stce_9839) { return floor(stce_9839._value()); } template<class stce_9771, class stce_9935 > static inline bool isfinite(const dco::internal::active_type<stce_9771, stce_9935>& stce_9839) { return isfinite(stce_9839._value()); }
    using std::ceil; using std::floor; using std::isfinite; using std::isnan; using std::isnormal; using std::isinf; template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool isnan(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9839) { return isnan(stce_9839._value()); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool isnormal(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9839) { return isnormal(stce_9839._value()); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool isinf(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9839) { return isinf(stce_9839._value()); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline double ceil(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9839) { return ceil(stce_9839._value()); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline double floor(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9839) { return floor(stce_9839._value()); } template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline bool isfinite(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9839) { return isfinite(stce_9839._value()); }
    using std::ceil; using std::floor; using std::isfinite; using std::isnan; using std::isnormal; using std::isinf; template<class stce_9771, class stce_9937, class stce_9939 > static inline bool isnan(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9839) { return isnan(stce_9839._value()); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool isnormal(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9839) { return isnormal(stce_9839._value()); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool isinf(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9839) { return isinf(stce_9839._value()); } template<class stce_9771, class stce_9937, class stce_9939 > static inline double ceil(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9839) { return ceil(stce_9839._value()); } template<class stce_9771, class stce_9937, class stce_9939 > static inline double floor(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9839) { return floor(stce_9839._value()); } template<class stce_9771, class stce_9937, class stce_9939 > static inline bool isfinite(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9839) { return isfinite(stce_9839._value()); }
    using std::ceil; using std::floor; using std::isfinite; using std::isnan; using std::isnormal; using std::isinf; template<class stce_9771, class stce_9938, class stce_9939 > static inline bool isnan(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9839) { return isnan(stce_9839._value()); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool isnormal(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9839) { return isnormal(stce_9839._value()); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool isinf(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9839) { return isinf(stce_9839._value()); } template<class stce_9771, class stce_9938, class stce_9939 > static inline double ceil(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9839) { return ceil(stce_9839._value()); } template<class stce_9771, class stce_9938, class stce_9939 > static inline double floor(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9839) { return floor(stce_9839._value()); } template<class stce_9771, class stce_9938, class stce_9939 > static inline bool isfinite(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9839) { return isfinite(stce_9839._value()); }
    using std::ceil; using std::floor; using std::isfinite; using std::isnan; using std::isnormal; using std::isinf; template<class stce_9771, class stce_9945, class stce_9939 > static inline bool isnan(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9839) { return isnan(stce_9839._value()); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool isnormal(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9839) { return isnormal(stce_9839._value()); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool isinf(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9839) { return isinf(stce_9839._value()); } template<class stce_9771, class stce_9945, class stce_9939 > static inline double ceil(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9839) { return ceil(stce_9839._value()); } template<class stce_9771, class stce_9945, class stce_9939 > static inline double floor(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9839) { return floor(stce_9839._value()); } template<class stce_9771, class stce_9945, class stce_9939 > static inline bool isfinite(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9839) { return isfinite(stce_9839._value()); }
    template<class stce_9771, class stce_9935 > static inline std::ostream& operator << (std::ostream& out, const dco::internal::active_type<stce_9771, stce_9935>& stce_9839) { out << stce_9839._value(); return out; }
    template<class stce_9771, class stce_9937, class stce_9938, class stce_9939 > static inline std::ostream& operator << (std::ostream& out, const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939>& stce_9839) { out << stce_9839._value(); return out; }
    template<class stce_9771, class stce_9937, class stce_9939 > static inline std::ostream& operator << (std::ostream& out, const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939>& stce_9839) { out << stce_9839._value(); return out; }
    template<class stce_9771, class stce_9938, class stce_9939 > static inline std::ostream& operator << (std::ostream& out, const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939>& stce_9839) { out << stce_9839._value(); return out; }
    template<class stce_9771, class stce_9945, class stce_9939 > static inline std::ostream& operator << (std::ostream& out, const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939>& stce_9839) { out << stce_9839._value(); return out; }
template <typename stce_9782, typename stce_9819>
  typename enable_if<mode<stce_9782>::is_dco_type || mode<stce_9819>::is_dco_type,
                     typename stce_9863<(mode<stce_9782>::order > mode<stce_9819>::order), typename mode<stce_9782>::type, typename mode<stce_9819>::type>::type
                                       >::type min(stce_9782 const& stce_9901, stce_9819 const& stce_10478) {
  if (passive_value(stce_9901) < passive_value(stce_10478)) return stce_9901;
  else return stce_10478;
}
template <typename stce_9782, typename stce_9819>
  typename enable_if<mode<stce_9782>::is_dco_type || mode<stce_9819>::is_dco_type,
                     typename stce_9863<(mode<stce_9782>::order > mode<stce_9819>::order), typename mode<stce_9782>::type, typename mode<stce_9819>::type>::type
  >::type max(stce_9782 const& stce_9901, stce_9819 const& stce_10478) {
  if (passive_value(stce_9901) > passive_value(stce_10478)) return stce_9901;
  else return stce_10478;
}
}
}
namespace dco {
  namespace instrument {
    extern bool dco_instrument;
    struct stce_10479 {
      struct iterator_t {
        DCO_TAPE_INT _progvarcounter;
        DCO_TAPE_INT _stackcounter;
        iterator_t(DCO_TAPE_INT stce_10480, DCO_TAPE_INT stce_10481) : _progvarcounter(stce_10480), _stackcounter(stce_10481) {}
        DCO_TAPE_INT stce_10482() {
          return _progvarcounter + _stackcounter * 2;
        }
        DCO_TAPE_INT &stce_10483() {
          return _progvarcounter;
        }
      };
      struct stce_10484 {
        DCO_TAPE_INT id;
        std::string stce_10485;
        DCO_TAPE_INT stce_10486;
        iterator_t stce_10158;
        DCO_TAPE_INT stce_10295;
        DCO_TAPE_INT stce_10487;
        DCO_TAPE_INT stce_10127;
        DCO_TAPE_INT stce_10128;
        DCO_TAPE_INT stce_10488;
        DCO_TAPE_INT stce_10489;
        stce_10484 *stce_10393;
        stce_10484(stce_10479 *tape, const std::string &name, stce_10484 *stce_10490) : id(++tape->stce_9853),
          stce_10485(name),
          stce_10486(-1),
          stce_10158(tape->get_position()),
          stce_10295(0),
          stce_10487(0),
          stce_10127(0),
          stce_10128(0),
          stce_10488(0),
          stce_10489(0),
          stce_10393(stce_10490) {
        }
      };
    public:
      struct callback_object_t {
        virtual ~callback_object_t() {
        }
        virtual void stce_10082(stce_9844 &stce_9845) {
          (void)stce_9845;
        }
      };
      struct interpretation_settings {
      };
      typedef void (*stce_10491)(stce_10479 &stce_10393, const interpretation_settings &stce_10305, callback_object_t *stce_10394);
      void stce_10492(stce_10491 stce_10493, callback_object_t *stce_10394) {
        (void) stce_10493;
        (void) stce_10394;
      }
      static stce_10479 *get_tape() {
        return NULL;
      }
      DCO_TAPE_INT tape_index() const {
        return 0;
      }
      DCO_TAPE_INT stce_9853;
      std::ofstream stce_10494;
      std::ofstream stce_10495;
      std::ofstream stce_10496;
      std::ofstream stce_10497;
      std::ofstream stce_10498;
      std::ofstream stce_10499;
      DCO_TAPE_INT _progvarcounter;
      DCO_TAPE_INT _stackcounter;
      stce_10484 *stce_10163;
      std::vector<stce_10484 *> stce_10500;
      std::stack<stce_10484 *> stce_10501;
      DCO_TAPE_INT stce_10502;
      stce_10479() : stce_9853(0),
        _progvarcounter(0),
        _stackcounter(0),
        stce_10163(NULL),
        stce_10502(0) { }
      ~stce_10479() {
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
      static void remove(stce_10479 *stce_9782) {
        delete stce_9782;
      }
      static stce_10479 *create(const DCO_TAPE_INT size = 0, const DCO_TAPE_INT stce_10503 = 0) {
        (void) size;
        (void) stce_10503;
        stce_10479 *stce_9836 = new stce_10479();
        stce_9836->_init_instrument();
        return stce_9836;
      }
      static stce_10479 *create(std::string stce_10504) {
        stce_10479 *stce_9836 = new stce_10479();
        stce_9836->_init_instrument(stce_10504);
        return stce_9836;
      }
      template<class stce_9771, class stce_9935>
      void register_variable(dco::internal::active_type<stce_9771, stce_9935> &stce_9839) {
        _progvarcounter++;
        stce_9839.tape_index = _progvarcounter;
        stce_9839.stce_10505 = stce_9853;
      }
      template<class stce_9771, class stce_9935>
      void register_variable(dco::internal::active_type<stce_9771, stce_9935> *stce_9839, DCO_TAPE_INT stce_10096) {
        for (DCO_TAPE_INT stce_9827 = 0; stce_9827 < stce_10096; stce_9827++)
          this->register_variable(stce_9839[stce_9827]);
      }
      template<class stce_9771, class stce_9935>
      void register_variable(std::vector<dco::internal::active_type<stce_9771, stce_9935> > &stce_9839) {
        register_variable(&(stce_9839[0]), stce_9839.size());
      }
      template<class stce_9771, class stce_9935>
      void register_output_variable(dco::internal::active_type<stce_9771, stce_9935> &stce_9839) {
        stce_9839 = 1 * stce_9839;
      }
      template<class stce_9771, class stce_9935>
      void register_output_variable(dco::internal::active_type<stce_9771, stce_9935> *stce_9839, const DCO_TAPE_INT stce_10096) {
        if (stce_10096 == 0) return;
        for (DCO_TAPE_INT stce_9827 = 0; stce_9827 < stce_10096; ++stce_9827)
          stce_9839[stce_9827] = 1 * stce_9839[stce_9827];
      }
      template<class stce_9771, class stce_9935>
      void register_output_variable(std::vector<dco::internal::active_type<stce_9771, stce_9935> > &stce_9839) {
        register_output_variable(&(stce_9839[0]), stce_9839.size());
      }
      void _init_instrument(const std::string &stce_10504 = "") {
        stce_9853 = 0;
        std::string stce_10506 = stce_10504;
        if (stce_10504.size() != 0) {
          stce_10506 += "_";
        }
        stce_10494.open((stce_10506 + "nodes.dot").c_str());
        stce_10495.open((stce_10506 + "edges.dot").c_str());
        stce_10494 << "digraph {" << std::endl;
        stce_10496.open((stce_10506 + "nodes.csv").c_str());
        stce_10497.open((stce_10506 + "edges.csv").c_str());
        stce_10498.open((stce_10506 + "ext_nodes.csv").c_str());
        stce_10499.open((stce_10506 + "ext_edges.csv").c_str());
        dco_instrument = true;
        std::string name = "ROOT";
        std::string stce_10393 = "caller";
        _dco_instrument_enter_function(name, stce_10393);
      }
      void _finish_instrument() {
        std::string name = "ROOT";
        std::string stce_10393 = "caller";
        _dco_instrument_leave_function(name, stce_10393);
        for (size_t stce_9827 = 0; stce_9827 < stce_10500.size(); ++stce_9827) {
          stce_10484 *stce_9815 = stce_10500[stce_9827];
          stce_10494 << stce_9815->id << "[shape=box,label=\"" << stce_9815->stce_10485 << "("
                    << stce_9815->stce_10295 << "," << stce_9815->stce_10127 << "," << stce_9815->stce_10128 << ")\"]" << std::endl;
          stce_10496 << stce_9815->stce_10488 << " "
                    << stce_9815->stce_10295 << " "
                    << stce_9815->stce_10489 << " "
                    << stce_9815->stce_10127 << " "
                    << stce_9815->stce_10128 << " "
                    << stce_9815->stce_10485 << std::endl;
          delete stce_9815;
        }
        stce_10500.clear();
        stce_10495 << "}" << std::endl;
        stce_10494.close();
        stce_10495.close();
        stce_10496.close();
        stce_10497.close();
        stce_10498.close();
        stce_10499.close();
        dco_instrument = false;
      }
      iterator_t get_position() {
        return iterator_t(_progvarcounter, _stackcounter);
      }
      void interpret_adjoint() { }
      void interpret_adjoint_to(iterator_t &stce_10158) {
        (void) stce_10158;
      }
      void reset() {
        std::cerr << "dco instrument -- possible error due to reset() call" << std::endl;
      }
      void reset_to(iterator_t &stce_10158) {
        std::cerr << "dco instrument -- possible error due to reset() call" << std::endl;
        (void) stce_10158;
      }
      void _dco_instrument_enter_function(const std::string &stce_10485 , const std::string &stce_10393 ) {
        (void) stce_10485;
        (void) stce_10393;
        stce_10484 *stce_9815 = new stce_10484(this, stce_10485, stce_10163);
        stce_10500.push_back(stce_9815);
        stce_10163 = stce_9815;
        if (! stce_10501.empty()) {
          stce_10484 *stce_10376 = stce_10501.top();
          DCO_TAPE_INT size = get_position().stce_10482() - stce_10376->stce_10158.stce_10482();
          stce_10495 << stce_10376->id << " -> " << stce_9815->id << "[label=\"" << size << "\"]" << std::endl;
          stce_10497 << stce_10376->id << " " << stce_9815->id << " " << size << std::endl;
          stce_10376->stce_10488++;
          stce_10376->stce_10295 += size;
          stce_10376->stce_10487 += size;
          stce_10376->stce_10158 = get_position();
        }
        stce_10501.push(stce_9815);
        stce_10502 = stce_9815->stce_10158.stce_10483();
        stce_9815->stce_10486 = stce_10502;
      }
      void _dco_instrument_leave_function(const std::string &stce_10485 , const std::string &stce_10393 ) {
        (void) stce_10485;
        (void) stce_10393;
        stce_10484 *stce_9815 = stce_10501.top();
        DCO_TAPE_INT size = get_position().stce_10482() - stce_9815->stce_10158.stce_10482();
        stce_9815->stce_10295 += size;
        stce_9815->stce_10487 += size;
        stce_10501.pop();
        if (!stce_10501.empty()) {
          stce_10484 *stce_10507 = stce_10501.top();
          stce_10507->stce_10158 = get_position();
          stce_10507->stce_10487 += stce_9815->stce_10487;
          stce_10502 = stce_10507->stce_10486;
          stce_10163 = stce_10507;
        }
      }
    };
  template<class stce_9771, stce_10479 *&stce_10508>
    struct stce_10509 {
    typedef internal::stce_9949 stce_9989;
    typedef stce_9803<stce_9771> mode_t;
  stce_9771 stce_10510;
  DCO_TAPE_INT tape_index;
  DCO_TAPE_INT _tape_index;
  mutable DCO_TAPE_INT stce_10505;
  mutable stce_10479::stce_10484 *stce_10393;
  stce_10509 &operator =(const stce_10509 &stce_9890) {
    tape_index = stce_9890.tape_index;
    stce_10505 = stce_9890.stce_10505;
    stce_10393 = stce_9890.stce_10393;
    _tape_index = stce_9890._tape_index;
    return *this;
  }
  stce_10509() :
      stce_10510(0.0), tape_index(0), _tape_index(0), stce_10505(0), stce_10393(0) {
    if (stce_10508)
      stce_10505 = stce_10508->stce_9853;
  }
  const stce_9771 &get_derivative() const {
    return stce_10510;
  }
  stce_9771 &get_derivative() {
    return stce_10510;
  }
  void clear() {
    tape_index = 0;
    stce_10505 = 0;
  }
  DCO_TAPE_INT stce_10511() const {
    return stce_9755 && tape_index == 0 ? 0 : 1;
  }
  template<class stce_9935 >
  static DCO_TAPE_INT stce_10512(const dco::internal::active_type<stce_9771, stce_9935> &stce_9839) {
    const stce_10509 &data = stce_9839;
    if (data.tape_index != 0) {
      if (data.stce_10393 != 0) {
        if (stce_10508->stce_10163 != data.stce_10393) {
          data.stce_10393->stce_10128++;
          data.stce_10393 = 0;
        }
      }
      if (data.tape_index <= stce_10508->stce_10502) {
        stce_10479::stce_10484 *stce_10513 = stce_10508->stce_10501.top();
        while (stce_10513) {
          if (data.stce_10505 < stce_10513->id) {
            stce_10513->stce_10127++;
            stce_10513 = stce_10513->stce_10393;
          } else {
            stce_10513 = NULL;
          }
        }
        data.stce_10505 = stce_10508->stce_10501.top()->id;
      }
    }
    return stce_9839.stce_10511();
  }
  template<class stce_9937, class stce_9938, class stce_9939 >
  static DCO_TAPE_INT stce_10512(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839) {
    return stce_10512(stce_9839.stce_9941) + stce_10512(stce_9839.stce_9943);
  }
  template<class stce_9945, class stce_9939 >
  static DCO_TAPE_INT stce_10512(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839) {
    return stce_10512(stce_9839.stce_9946);
  }
  template<class stce_9937, class stce_9939 >
  static DCO_TAPE_INT stce_10512(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839) {
    return stce_10512(stce_9839.stce_9941);
  }
  template<class stce_9938, class stce_9939 >
  static DCO_TAPE_INT stce_10512(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839) {
    return stce_10512(stce_9839.stce_9943);
  }
  template<class stce_9935 >
  static DCO_TAPE_INT stce_10514(const dco::internal::active_type<stce_9771, stce_9935> &stce_9839) {
    (void) stce_9839;
    return 0;
  }
  template<class stce_9937, class stce_9938, class stce_9939 >
  static DCO_TAPE_INT stce_10514(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839) {
    return stce_10514(stce_9839.stce_9941) + stce_10514(stce_9839.stce_9943) + 1;
  }
  template<class stce_9945, class stce_9939 >
  static DCO_TAPE_INT stce_10514(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839) {
    return stce_10514(stce_9839.stce_9946) + 1;
  }
  template<class stce_9937, class stce_9939 >
  static DCO_TAPE_INT stce_10514(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839) {
    return stce_10514(stce_9839.stce_9941) + 1;
  }
  template<class stce_9938, class stce_9939 >
  static DCO_TAPE_INT stce_10514(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839) {
    return stce_10514(stce_9839.stce_9943) + 1;
  }
  template<class stce_9888>
  DCO_ALWAYS_INLINE void stce_9889(const stce_9888 &stce_9890) {
    DCO_TAPE_INT stce_10515 = stce_10512(stce_9890);
    if (dco_instrument && (stce_10508->stce_10163 != NULL)) {
      DCO_TAPE_INT stce_10516 = stce_10514(stce_9890);
      stce_10508->stce_10163->stce_10489 += (stce_10516);
    }
    if (NULL != stce_10508 && !stce_10508->is_active()) {
      clear();
    }
    if (stce_10515 <= 0) {
      clear();
      return;
    }
    stce_10508->_stackcounter += (stce_10515 + 1);
    stce_10508->_progvarcounter++;
    tape_index = stce_10508->_progvarcounter;
    stce_10393 = stce_10508->stce_10163;
  }
};
  }
}
namespace dco {
  namespace internal {
    template<class stce_9771>
    struct stce_10517 {
    private:
      struct stce_10518 {
        DCO_TAPE_INT stce_10515;
        DCO_TAPE_INT stce_10519;
      };
      stce_10518 *stce_10520;
      DCO_TAPE_INT stce_10521;
      int *stce_10522;
      DCO_TAPE_INT stce_10523;
      stce_10517(const stce_10517<stce_9771> &) {
        throw dco::exception::create<std::runtime_error>("dco_error: Not implemented yet! Please report, why needed!", "/home/vagrant/dco_cpp_make_release/dco_cpp_dev/src/dco_tape_blob_pattern.hpp", 30);
      }
      stce_10517(DCO_TAPE_INT stce_10515, DCO_TAPE_INT stce_10480) {
        stce_10520 = new stce_10518[stce_9752(stce_10480)];
        stce_10521 = 0;
        stce_10522 = new int[size_t(stce_10515)];
        stce_10523 = 0;
      }
    public:
      struct stce_10524 {
        friend struct stce_10517<stce_9771>;
      private:
        DCO_TAPE_INT stce_10525;
        DCO_TAPE_INT stce_10526;
        stce_10524(DCO_TAPE_INT stce_10527, DCO_TAPE_INT stce_10528) : stce_10525(stce_10527), stce_10526(stce_10528) {}
      public:
        stce_10524() : stce_10525(0), stce_10526(0) {};
        DCO_TAPE_INT stce_10529() {
          return stce_10525;
        }
        DCO_TAPE_INT stce_10530() {
          if (stce_10526 < 0) return -stce_10526;
          return stce_10526;
        }
        bool stce_10531() {
          return stce_10526 < 0;
        }
      };
      struct iterator_t {
      private:
        DCO_TAPE_INT _progvarcounter;
      public:
        DCO_TAPE_INT stce_10480() const {
          return _progvarcounter;
        }
        iterator_t(const DCO_TAPE_INT stce_10532) : _progvarcounter(stce_10532) {}
      };
      iterator_t get_position() {
        return iterator_t(stce_10521);
      }
      typedef std::vector<stce_10524> sparse_jacobi;
      template<class stce_9782>
      sparse_jacobi *build_pattern(stce_9782 *stce_10529, DCO_TAPE_INT stce_10096, const iterator_t &stce_10533) {
        sparse_jacobi *stce_9836 = new sparse_jacobi();
        for (DCO_TAPE_INT stce_9827 = 0; stce_9827 < stce_10096; ++stce_9827) {
          stce_10534(stce_9827 + 1, stce_10529[stce_9827].tape_index, stce_10533, stce_9836);
        }
        return stce_9836;
      }
      void create_pattern_file(sparse_jacobi *stce_10535, std::string filename) {
        std::ofstream out(filename.c_str());
        for (size_t stce_9827 = 0; stce_9827 < stce_10535->size(); ++stce_9827) {
          stce_10524 &stce_10536 = (*stce_10535)[stce_9827];
          out << stce_10536.stce_10529() << " " << stce_10536.stce_10530() << " 1" << std::endl;
        }
        out.close();
      }
      void stce_10537(std::string filename) {
        std::ofstream out(filename.c_str());
        out << "digraph {" << std::endl;
        for (DCO_TAPE_INT stce_9827 = 1; stce_9827 <= stce_10521; ++stce_9827) {
          out << stce_9827 << std::endl;
        }
        for (DCO_TAPE_INT stce_9827 = 1; stce_9827 <= stce_10521; ++stce_9827) {
          int stce_10515 = stce_10520[stce_9827].stce_10515;
          int stce_10538 = stce_10520[stce_9827].stce_10519;
          for (int stce_10539 = 0; stce_10539 < stce_10515; ++stce_10539) {
            int stce_10048 = stce_10522[stce_10539 + stce_10538];
            out << stce_9827 << " -> " << stce_10048 << std::endl;
          }
        }
        out << "}" << std::endl;
        out.close();
      }
      void stce_10534(DCO_TAPE_INT stce_10527, DCO_TAPE_INT stce_9853, const iterator_t &stce_10533, sparse_jacobi *stce_10540) {
        if (stce_9853 == 0) return;
        if (stce_9853 <= stce_10533.stce_10480()) {
          stce_10540->push_back(stce_10524(stce_10527, stce_9853));
        } else {
          DCO_TAPE_INT stce_10515 = stce_10520[stce_9853].stce_10515;
          DCO_TAPE_INT stce_10538 = stce_10520[stce_9853].stce_10519;
          for (int stce_9827 = 0; stce_9827 < stce_10515; ++stce_9827) {
            stce_10534(stce_10527, stce_10522[stce_9827 + stce_10538], stce_10533, stce_10540);
          }
        }
      }
      static stce_10517 *create(DCO_TAPE_INT size, DCO_TAPE_INT stce_10503 = 0) {
        if (stce_10503 == 0) stce_10503 = size / 2;
        return new stce_10517(size, stce_10503);
      }
      template<class stce_9935 >
      void register_variable(dco::internal::active_type<stce_9771, stce_9935> &stce_10541) {
        stce_9935 &data = stce_10541.data();
        stce_10521++;
        data.tape_index = stce_10521;
      }
      template<stce_10517 *&global_tape>
      struct stce_10542 {
        typedef stce_9803<stce_9771> mode_t;
        DCO_TAPE_INT tape_index;
        stce_10542() : tape_index(0) {}
        void clear() {
          tape_index = 0;
        }
        template<class stce_9935 >
        static void interpret(const dco::internal::active_type<stce_9771, stce_9935> &stce_9839, bool stce_10531) {
          (void) stce_10531;
          global_tape->stce_10522[global_tape->stce_10523] = stce_9839.tape_index;
          global_tape->stce_10523++;
        }
        template<class stce_9937, class stce_9938, class stce_9939 >
        static void interpret(const dco::internal::stce_9779<stce_9771, stce_9937, stce_9938, stce_9939> &stce_9839, bool stce_10531) {
          interpret(stce_9839.stce_9941, stce_10531);
          interpret(stce_9839.stce_9943, stce_10531);
        }
        template<class stce_9945, class stce_9939 >
        static void interpret(const dco::internal::stce_9776<stce_9771, stce_9945, stce_9939> &stce_9839, bool stce_10531) {
          interpret(stce_9839.stce_9946, stce_10531);
        }
        template<class stce_9937, class stce_9939 >
        static void interpret(const dco::internal::stce_9780<stce_9771, stce_9937, stce_9939> &stce_9839, bool stce_10531) {
          interpret(stce_9839.stce_9941, stce_10531);
        }
        template<class stce_9938, class stce_9939 >
        static void interpret(const dco::internal::stce_9781<stce_9771, stce_9938, stce_9939> &stce_9839, bool stce_10531) {
          interpret(stce_9839.stce_9943, stce_10531);
        }
        template<class stce_9888, class stce_9794>
        static void stce_9889(const stce_9888 &stce_9890, stce_9794 &stce_10543) {
          const int stce_10519 = global_tape->stce_10523;
          interpret(stce_9890, false);
          stce_10542 &data = const_cast<stce_10542 &>(stce_10543);
          if (stce_10519 == global_tape->stce_10523) {
            data.clear();
          } else {
            global_tape->stce_10521++;
            data.tape_index = global_tape->stce_10521;
            global_tape->stce_10520[global_tape->stce_10521].stce_10515 = global_tape->stce_10523 - stce_10519;
            global_tape->stce_10520[global_tape->stce_10521].stce_10519 = stce_10519;
          }
        }
        template<class stce_9888>
        void stce_9889(const stce_9888 &stce_9890) {
          stce_9889(stce_9890, *this);
        }
      };
    };
  }
}
namespace dco {
  namespace instrument {
    typedef stce_10479 tape_t;
    extern tape_t *global_tape;
    typedef stce_10509<double, global_tape> data;
    typedef internal::active_type<double, data> type;
    typedef helper::stce_10102<type, tape_t> external_adjoint_object_t;
  }
  namespace p1t {
    typedef internal::stce_10517<double > tape_t;
    extern tape_t *global_tape;
    typedef tape_t::stce_10542<global_tape> data;
    typedef internal::active_type<double, data> type;
  }
  namespace p1f {
    typedef internal::stce_10054<32, double> data;
    typedef internal::active_type<double, data> type;
  }
  namespace p1f {
    inline void set(type &stce_9839, const bool stce_9753, const DCO_TAPE_INT stce_9806)
    {
      data &stce_10057 = stce_9839;
      stce_10057.stce_10056[stce_9806] = stce_9753;
    }
     inline void get(const type &stce_9839, bool &stce_9753, const DCO_TAPE_INT stce_9806)
    {
      stce_9753 = stce_9839.stce_10056[stce_9806];
    }
     inline void get(type &stce_9839, double &stce_9753)
    {
      stce_9753 = stce_9839._value();
    }
     inline std::string stce_10544(const type &stce_9839)
    {
      return stce_9839.stce_10056.to_string<char, std::char_traits<char>, std::allocator<char> >();
    }
  }
  namespace p1t {
     inline void set(type &stce_9839, const double stce_9753)
    {
      double &value = const_cast<double &>(stce_9839._value());
      value = stce_9753;
    }
     inline void get(type &stce_9839, double &stce_9753)
    {
      stce_9753 = stce_9839._value();
    }
  }
  namespace instrument {
  inline void set(type &stce_9839, const double stce_9753, const DCO_TAPE_INT what = 0)
    {
      stce_1000002();;
      if (what == 0) {
        double &value = const_cast<double &>(stce_9839._value());
        value = stce_9753;
      } else if (what == -1) {
      } else {
        throw exception::create<std::invalid_argument>("invalid kind of what (musst be 0 or -1)");
      }
    }
  inline void get(const type &stce_9839, double &stce_9753, const DCO_TAPE_INT what = 0)
    {
      stce_1000002();;
      if (what == 0)
        stce_9753 = stce_9839._value();
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
template<class stce_10413>
  struct jacobian_preaccumulator_t: stce_10073<typename dco::mode<stce_10413>::stce_10016> {
  typedef typename dco::mode<stce_10413>::stce_10016 stce_9764;
  typedef typename dco::mode<stce_10413>::tape_t tape_t;
  typedef typename tape_t::iterator_t stce_10545;
  typedef typename tape_t::stce_9920 stce_10546;
  tape_t* const stce_9948;
  stce_10545 stce_10547;
  DCO_TAPE_INT stce_10548;
  typedef stce_9793<stce_9764, tape_t> vector_t;
  vector_t stce_10284;
  bool stce_10549;
  void init() {}
  typedef stce_10143<memory_model::BLOB_TAPE, tape_t> stce_10277;
  typedef typename stce_10277::stce_9920 stce_10550;
  typedef typename stce_10277::iterator stce_10551;
  stce_10277* stce_10139;
  stce_10550* stce_9923;
  typedef typename stce_10413::data_t stce_10552;
  std::vector<stce_10552*> stce_10553;
  DCO_TAPE_INT stce_10303;
  jacobian_preaccumulator_t(tape_t * const tape) :
    stce_9948(tape),
    stce_10284(0, *tape, tape->stce_10164()),
    stce_10139(0),
    stce_9923(0),
    stce_10303(0)
  {
    if (stce_9948) {
      stce_10547 = stce_9948->get_position();
      stce_10548 = stce_10547.index();
      stce_10549 = true;
    } else {
      stce_10549 = false;
    }
  }
  ~jacobian_preaccumulator_t() {
    if (stce_10549)
      this->finish();
    if (stce_10139)
      delete stce_10139;
  }
  void start() {
    if (stce_9948) {
      stce_10547 = stce_9948->get_position();
      stce_10548 = stce_10547.index();
      stce_10549 = true;
    } else {
      stce_10549 = false;
    }
  }
  void stce_10554(stce_10413 &stce_9839) {
    if (dco::tape_index(stce_9839) && dco::tape_index(stce_9839) <= stce_10548) {
      throw dco::exception::create<std::runtime_error>("jacobian_preaccumulator_t: registering output, which is indeed an input (was alive already before starting point).");
    }
    stce_10553.push_back(&static_cast<stce_10552&>(stce_9839));
  }
  void register_output(stce_10413 &stce_9839)
  {
    if (stce_9948) { stce_10554(stce_9839); }
  }
  void register_output(std::vector<stce_10413> &stce_9839)
  {
    if (stce_9948) {
      stce_10553.reserve(stce_10553.size() + stce_9839.size());
      for (size_t stce_9827 = 0; stce_9827 < stce_9839.size(); ++stce_9827) {
        stce_10554(stce_9839[stce_9827]);
      }
    }
  }
  void finish() {
    if (!stce_10549) return;
    stce_10549 = false;
    if (!stce_9948) return;
    if (!stce_10553.size()) {
      stce_9948->reset_to(stce_10547);
      return;
    }
    const stce_10545 stce_10555 = stce_9948->get_position();
    const DCO_TAPE_INT stce_10556 = stce_10555.index() - stce_10547.index();
    stce_10284.resize(stce_10556 + 1);
    const DCO_TAPE_INT stce_10557 = stce_9948->size(stce_10555, stce_10547);
    const DCO_TAPE_INT stce_10558 = stce_10557 * (1 + static_cast<DCO_TAPE_INT>(stce_10553.size()));
    if (stce_10139) delete stce_10139;
    tape_options stce_10559;
    stce_10559.blob_size_in_byte() = stce_10558 * stce_10277::entry::size_of();
    stce_10139 = new stce_10277(stce_10559, *stce_9948);
    stce_10139->stce_10194().index() = stce_10555.index();
    stce_10284.stce_10300();
    for(typename std::vector<stce_10552*>::iterator stce_10158 = stce_10553.begin(); stce_10158 != stce_10553.end(); ++stce_10158) {
      stce_10552& stce_9815 = **stce_10158;
      if (stce_9815._tape_index() == 0) continue;
      stce_10284.stce_10297(stce_10556-1,0);
      stce_10043(stce_9815._tape_index()) = 1.0;
      stce_9923 = new stce_10550(stce_10139->stce_9924());
      stce_9948->stce_10310(*this,
                                           stce_10555, stce_10547,
                                           interpretation_settings(false, true));
      delete stce_9923;
    }
    stce_9948->reset_to(stce_10547);
    stce_10551 stce_10228 = stce_10139->current();
    for (typename std::vector<stce_10552*>::reverse_iterator stce_10158 = stce_10553.rbegin();
         stce_10158 != stce_10553.rend(); ++stce_10158) {
      stce_10546 stce_10560(stce_9948->stce_9924());
      for (;; --stce_10228) {
        stce_10560.insert(stce_10228->arg(), stce_10228->pval());
        if(stce_10228->stce_10151())
          break;
      }
      --stce_10228;
      (*stce_10158)->_tape_index() = stce_10560.index();
    }
    stce_10553.clear();
  }
  stce_9764& stce_10043(const DCO_TAPE_INT stce_9806) {
    return stce_10284[stce_9806 - stce_10548];
  }
  template<class stce_10309>
    void stce_10075(const DCO_TAPE_INT &stce_9806, const stce_10309 &stce_9981,
                            const stce_9764 &stce_10076) {
    stce_10075(stce_9806, stce_9981 * stce_10076);
  }
  void stce_10297(DCO_TAPE_INT stce_9858, DCO_TAPE_INT stce_9859) {
    stce_10284.stce_10297(stce_9858,stce_9859);
  }
  void stce_10289(DCO_TAPE_INT stce_9806) {
    stce_10284.stce_10289(stce_9806 - stce_10548);
  }
  stce_9764& operator[](const DCO_TAPE_INT stce_9806) {
    return stce_10284[stce_9806 - stce_10548];
  }
  template <typename iterator_t>
  void stce_10290(iterator_t& stce_9815) {
    stce_10075(stce_9815->arg(), stce_9815->pval()*(*this)[stce_9815.index()]);
  }
  virtual stce_9764 stce_10074(const DCO_TAPE_INT stce_9806) {
    return stce_10284[stce_9806 - stce_10548];
  }
  virtual void stce_10075(const DCO_TAPE_INT stce_9806, const stce_9764 &stce_10076) {
    if (stce_10076 == 0)
      return;
    if (stce_9806 <= stce_10548) {
      stce_9923->insert(stce_9806, stce_10076);
    } else {
      stce_10043(stce_9806) += stce_10076;
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
  typedef stce_9762<stce_9782, stce_9760, stce_9761> stce_10014;
public:
  static const memory_model::TYPE TAPE_MODEL = MEMORY_MODEL;
  typedef dco::internal::single_tape_data<ga1s> data_t;
  typedef typename stce_10014::value_t value_t;
  typedef typename stce_10014::value_t active_value_t;
  typedef typename stce_10014::stce_9764 derivative_t;
  typedef derivative_t stce_10016;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef internal::tape<TAPE_MODEL, stce_10014> tape_t;
  static DCOi_THREAD_LOCAL tape_t *global_tape;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10078<tape_t> callback_object_t;
  typedef dco::helper::stce_10084<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10102<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
template <class stce_9782, class stce_9760, class stce_9761, memory_model::TYPE MEMORY_MODEL>
  DCOi_THREAD_LOCAL typename ga1s<stce_9782,stce_9760,stce_9761,MEMORY_MODEL>::tape_t *ga1s<stce_9782,stce_9760,stce_9761,MEMORY_MODEL>::global_tape = 0;
}
       
namespace dco {
template <class stce_9782, class stce_9760 = stce_9782, class stce_9761 = stce_9760>
  class ga1s_mod {
  typedef stce_9762<stce_9782, stce_9760, stce_9761> stce_10014;
public:
  typedef dco::internal::single_tape_data<ga1s_mod> data_t;
  static const memory_model::TYPE TAPE_MODEL = memory_model::stce_9770;
  static const memory_model::TYPE stce_10561 = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10014::value_t value_t;
  typedef typename stce_10014::value_t active_value_t;
  typedef typename stce_10014::stce_9764 derivative_t;
  typedef derivative_t stce_10016;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef dco::internal::tape<TAPE_MODEL, stce_10014, stce_10561> tape_t;
  static DCOi_THREAD_LOCAL tape_t *global_tape;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10078<tape_t> callback_object_t;
  typedef dco::helper::stce_10084<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10102<type, tape_t> external_adjoint_object_t;
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
  typedef stce_9762<stce_9782, stce_9760, stce_9761> stce_10014;
public:
  typedef dco::internal::stce_9997<gtas> data_t;
  static const memory_model::TYPE TAPE_MODEL = memory_model::stce_9769;
  typedef typename stce_10014::value_t value_t;
  typedef typename stce_10014::value_t active_value_t;
  typedef typename stce_10014::stce_9764 derivative_t;
  typedef derivative_t stce_10016;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = true;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef internal::tape<TAPE_MODEL, stce_10014> tape_t;
  static tape_t *global_tape;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10078<tape_t> callback_object_t;
  typedef dco::helper::stce_10084<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10102<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
  typedef typename tape_t::logger logger;
  enum stce_10562 {
    check_branch = 1,
    check_fd = 2,
    check_identity = 4,
    stce_10563 = check_branch | check_fd | check_identity
  };
  static double& fd_eps() {
    assert(global_tape);
    return global_tape->stce_10246;
  }
  static double& ta_eps() {
    assert(global_tape);
    return global_tape->stce_10247;
  }
  static void init(const int stce_10388 = stce_10563) {
    assert(global_tape);
    global_tape->template stce_10262<data_t>();
    global_tape->stce_10243 = (stce_10388 & check_branch) != 0;
    global_tape->stce_10244 = (stce_10388 & check_fd) != 0;
    global_tape->stce_10006 = (stce_10388 & check_identity) != 0;
    global_tape->stce_10249 = 0;
    global_tape->stce_10250 = -1;
  }
  static DCO_TAPE_INT& throw_on_event() {
    assert(global_tape);
    return global_tape->stce_10250;
  }
  static bool& throw_on_warning() {
    assert(global_tape);
    return global_tape->stce_10245;
  }
};
template<typename stce_9782, typename stce_10026, typename stce_10564> typename gtas<stce_9782, stce_10026, stce_10564>::tape_t *gtas<stce_9782, stce_10026, stce_10564>::global_tape = 0;
}
       
namespace dco {
template <class stce_9782, class stce_9760 = stce_9782, class stce_9761 = stce_9760>
  class ga1sm {
  typedef stce_9762<stce_9782, stce_9760, stce_9761> stce_10014;
public:
  typedef dco::internal::stce_9992<ga1sm> data_t;
  static const memory_model::TYPE TAPE_MODEL = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10014::value_t value_t;
  typedef typename stce_10014::value_t active_value_t;
  typedef typename stce_10014::stce_9764 derivative_t;
  typedef derivative_t stce_10016;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef internal::tape<TAPE_MODEL, stce_10014> tape_t;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10078<tape_t> callback_object_t;
  typedef dco::helper::stce_10084<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10102<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
}
       
namespace dco {
template <class stce_9782, class stce_9760 = stce_9782, class stce_9761 = stce_9760>
  class ga1sm_mod {
  typedef stce_9762<stce_9782, stce_9760, stce_9761> stce_10014;
public:
  typedef dco::internal::stce_9992<ga1sm_mod> data_t;
  static const memory_model::TYPE TAPE_MODEL = memory_model::stce_9770;
  static const memory_model::TYPE stce_10561 = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10014::value_t value_t;
  typedef typename stce_10014::value_t active_value_t;
  typedef typename stce_10014::stce_9764 derivative_t;
  typedef derivative_t stce_10016;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef dco::internal::tape<TAPE_MODEL, stce_10014, stce_10561> tape_t;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10078<tape_t> callback_object_t;
  typedef dco::helper::stce_10084<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10102<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
}
namespace dco {
template <class stce_9782, const int stce_10565=1, class stce_9760=stce_9782, class stce_9761=stce_9760>
  class ga1vm {
  typedef stce_9762<stce_9782, stce_9760, dco::helper::stce_9822<stce_9761, stce_10565> > stce_10014;
public:
  typedef dco::internal::stce_9992<ga1vm> data_t;
  static const memory_model::TYPE TAPE_MODEL = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10014::value_t value_t;
  typedef typename stce_10014::value_t active_value_t;
  typedef typename stce_10014::stce_9764 derivative_t;
  typedef stce_9761 stce_10016;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef internal::tape<TAPE_MODEL, stce_10014> tape_t;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10078<tape_t> callback_object_t;
  typedef dco::helper::stce_10084<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10102<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
}
namespace dco {
template <class stce_9782, const int stce_10565 = 1, class stce_9760 = stce_9782, class stce_9761 = stce_9760>
  class ga1vm_mod {
  typedef stce_9762<stce_9782, stce_9760, dco::helper::stce_9822<stce_9761, stce_10565> > stce_10014;
public:
  typedef dco::internal::stce_9992<ga1vm_mod> data_t;
  static const memory_model::TYPE TAPE_MODEL = memory_model::stce_9770;
  static const memory_model::TYPE stce_10561 = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10014::value_t value_t;
  typedef typename stce_10014::value_t active_value_t;
  typedef typename stce_10014::stce_9764 derivative_t;
  typedef typename derivative_t::value_t stce_10016;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef dco::internal::tape<TAPE_MODEL, stce_10014, stce_10561> tape_t;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10078<tape_t> callback_object_t;
  typedef dco::helper::stce_10084<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10102<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
}
namespace dco {
template <class stce_9782, const int stce_10565 = 1, class stce_9760 = stce_9782, class stce_9761 = stce_9760>
  class ga1v {
  typedef stce_9762<stce_9782, stce_9760, dco::helper::stce_9822<stce_9761, stce_10565> > stce_10014;
public:
  typedef dco::internal::single_tape_data<ga1v> data_t;
  static const memory_model::TYPE TAPE_MODEL = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10014::value_t value_t;
  typedef typename stce_10014::value_t active_value_t;
  typedef typename stce_10014::stce_9764 derivative_t;
  typedef stce_9761 stce_10016;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef internal::tape<TAPE_MODEL, stce_10014> tape_t;
  static DCOi_THREAD_LOCAL tape_t *global_tape;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10078<tape_t> callback_object_t;
  typedef dco::helper::stce_10084<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10102<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
template <class stce_9782, const int stce_10565, class stce_9760, class stce_9761>
  DCOi_THREAD_LOCAL typename ga1v<stce_9782,stce_10565,stce_9760,stce_9761>::tape_t *ga1v<stce_9782,stce_10565,stce_9760,stce_9761>::global_tape;
}
namespace dco {
template <class stce_9782, const int stce_10565 = 1, class stce_9760 = stce_9782, class stce_9761 = stce_9760>
  class ga1v_mod {
  typedef stce_9762<stce_9782, stce_9760, dco::helper::stce_9822<stce_9761, stce_10565> > stce_10014;
public:
  typedef dco::internal::single_tape_data<ga1v_mod> data_t;
  static const memory_model::TYPE TAPE_MODEL = memory_model::stce_9770;
  static const memory_model::TYPE stce_10561 = DCO_TAPE_MEMORY_MODEL;
  typedef typename stce_10014::value_t value_t;
  typedef typename stce_10014::value_t active_value_t;
  typedef typename stce_10014::stce_9764 derivative_t;
  typedef typename derivative_t::value_t stce_10016;
  typedef typename mode<value_t>::passive_t passive_t;
  typedef dco::internal::active_type<value_t, data_t> type;
  typedef type active_t;
  static const bool is_dco_type = true;
  static const bool is_adjoint_type = true;
  static const bool is_tangent_type = false;
  static const int order = dco::mode<stce_9782>::order + 1;
  typedef dco::internal::tape<TAPE_MODEL, stce_10014, stce_10561> tape_t;
  static DCOi_THREAD_LOCAL tape_t *global_tape;
  typedef typename internal::local_gradient_t<type> local_gradient_t;
  typedef dco::helper::stce_10078<tape_t> callback_object_t;
  typedef dco::helper::stce_10084<type, tape_t> userdata_object_t;
  typedef dco::helper::stce_10102<type, tape_t> external_adjoint_object_t;
  typedef internal::jacobian_preaccumulator_t<type> jacobian_preaccumulator_t;
};
template <class stce_9782, const int stce_10565, class stce_9760, class stce_9761>
  DCOi_THREAD_LOCAL typename ga1v_mod<stce_9782,stce_10565,stce_9760,stce_9761>::tape_t *ga1v_mod<stce_9782,stce_10565,stce_9760,stce_9761>::global_tape;
}
namespace dco
{
  template <typename stce_10566>
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
  template <typename stce_10566> typename ginstrument<stce_10566>::tape_t *&ginstrument<stce_10566>::global_tape = instrument::global_tape;
}
       
namespace dco {
template <class stce_9782, typename stce_10567 = stce_9782>
  class gt1s {
public:
  typedef dco::internal::ts_data<gt1s> data_t;
  typedef stce_9782 value_t;
  typedef stce_9782 active_value_t;
  typedef stce_10567 derivative_t;
  typedef stce_10567 stce_10016;
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
          int stce_10568=1>
  class gt1v {
public:
    typedef dco::internal::stce_10050<gt1v> data_t;
    typedef stce_9782 value_t;
    typedef stce_9782 active_value_t;
    typedef dco::helper::stce_9822<stce_9782, stce_10568> derivative_t;
    typedef stce_9782 stce_10016;
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
  typedef ga1s<double> a1w;
  typedef double w_rtype;
}
namespace dco {
  namespace instrument {
    class instrumentation_data {
      const std::string stce_10197;
      const std::string stce_10569;
      bool stce_10570;
    public:
      instrumentation_data(std::string stce_10571) : stce_10197(stce_10571), stce_10569("caller"), stce_10570(false) {
        if (dco::instrument::dco_instrument && dco::instrument::global_tape) {
          dco::instrument::global_tape->_dco_instrument_enter_function(stce_10197, stce_10569);
          stce_10570 = true;
        }
      }
      ~instrumentation_data() {
        if (stce_10570) {
          dco::instrument::global_tape->_dco_instrument_leave_function(stce_10197, stce_10569);
        }
      }
    };
    inline void stce_10572(std::string stce_10571) {
      std::string stce_10569 = "caller";
      dco::instrument::global_tape->_dco_instrument_enter_function(stce_10571, stce_10569);
    }
    inline void stce_10573(std::string stce_10571) {
      std::string stce_10569 = "caller";
      dco::instrument::global_tape->_dco_instrument_leave_function(stce_10571, stce_10569);
    }
  }
}

