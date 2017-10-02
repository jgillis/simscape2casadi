struct NeSparsityPattern {
  int* mIr;
  int* mJc;
};

typedef int int32_T;
typedef double real_T;
typedef int size_t;
typedef int boolean_T;

struct Dummy {
  int dummy;
};

typedef struct Dummy NeAllocator;
typedef struct Dummy PmAllocator;
typedef struct Dummy NeDynamicSystem;
typedef struct Dummy NeDsMethodOutput;
typedef struct Dummy NeDynamicSystemInput;
typedef struct Dummy NeEquationData;
typedef struct Dummy NeVariableData;
typedef struct Dummy NeModeData;
typedef struct Dummy NeZCData;
typedef struct Dummy NeRange;
typedef struct Dummy NeAssertData;
typedef struct Dummy NeIntVector;
typedef struct Dummy NeSparsityPattern;
typedef struct Dummy NeRealVector;
typedef struct Dummy NeBoolVector;
typedef struct Dummy NeDsIoInfo;
typedef struct Dummy PmRealVector;
typedef struct Dummy NeObservableData;
typedef struct Dummy NeParameterData;
typedef struct Dummy PmSparsityPattern;
typedef struct Dummy PmIntVector;
typedef struct Dummy PmBoolVector;
typedef struct Dummy NeRtlEqInfo;


void  ne_rtl_call_method(void*, int, int, void*, void*, void*, void*, void*);
