#include "ruby.h"


#define RGET_CLASS(class_over, x) (rb_const_get(class_over, rb_intern(x)))
#define RGET_CLASS_TOP(x) (RGET_CLASS(rb_cObject, x))
#define RFCALL_10(name) (rb_funcall(self, rb_intern(name), 0))
#define RFCALL_10_ON(obj, name) (rb_funcall(obj, rb_intern(name), 0))
#define RFCALL_11(name, arg1) (rb_funcall(self, rb_intern(name), 1, arg1))
#define RFCALL_11_ON(obj, name, arg1) (rb_funcall(obj, rb_intern(name), 1, arg1))
#define RFCALL_12(name, arg1, arg2) (rb_funcall(self, rb_intern(name), 2, arg1, arg2))
#define RFCALL_12_ON(obj, name, arg1, arg2) (rb_funcall(obj, rb_intern(name), 2, arg1, arg2))
#define CR_ELEMENT_ACCESS(recvr, key) ( \
		RFCALL_11_ON(recvr, "[]", key) \
		)
#define CR_HKS(hash, cstr) (CR_ELEMENT_ACCESS(hash, ID2SYM(rb_intern(cstr))))
#define CR_HKS_SET(hash, cstr, value)(rb_funcall(hash, rb_intern("[]="), 2, ID2SYM(rb_intern(cstr)), value))
#define CR_RANGE_INC(start, end) (rb_funcall(RGET_CLASS_TOP("Range"), rb_intern("new"), 2, INT2FIX(start), INT2FIX(end)))
#define CR_RANGE_EXC(start, end) (rb_funcall(RGET_CLASS_TOP("Range"), rb_intern("new"), 3, INT2FIX(start), INT2FIX(end), Qtrue))

#define CR_SYM(cstr)(ID2SYM(rb_intern(cstr)))

/*Allocates an integer array on the heap
 with values of array. int_ptr should be
 an unallocated int*. array should not
 contain anything except ints.  */
#define CR_INT_ARY_R2C_STACK( array,int_ptr)\
 { int cr_internal_xaa12; \
	cr_internal_xaa12 = RARRAY_LEN(array);	\
	int_ptr =	ALLOCA_N(int, cr_internal_xaa12);\
	int cr_internal_xaa11;\
 	for (cr_internal_xaa11=0;\
		cr_internal_xaa11< cr_internal_xaa12;\
		cr_internal_xaa11++)\
			int_ptr[cr_internal_xaa11] = \
			FIX2INT(RARRAY_PTR(array)[cr_internal_xaa11]);\
		}

/* Get element of  a rank 1,2 or 3 tensor respectively */
#define CR_TELMT_R1(tensor, i)(rb_funcall(tensor, rb_intern("[]"), 1, INT2FIX(i)))
#define CR_TELMT_R2(tensor, i, j)(rb_funcall(tensor, rb_intern("[]"), 2, INT2FIX(i), INT2FIX(j)))
#define CR_TELMT_R3(tensor, i, j, k)(rb_funcall(tensor, rb_intern("[]"), 3, INT2FIX(i), INT2FIX(j), INT2FIX(k)))

// Set the value of a 4D tensor
#define CR_TELMT_R4_SET(tensor, i, j, k, l, value)(rb_funcall(tensor, rb_intern("[]="), 5, INT2FIX(i), INT2FIX(j), INT2FIX(k), INT2FIX(l), value))
 

/*Allocates an integer array on the heap
 with values of array. int_ptr should be
 an unallocated int*. array should not
 contain anything except ints; this is 
 checked in this macro. */
#define CR_INT_ARY_R2C_STACK_TCHECK( array,int_ptr)\
 { int cr_internal_xaa12; \
	cr_internal_xaa12 = RARRAY_LEN(array);	\
	int_ptr =	ALLOCA_N(int, cr_internal_xaa12);\
	int cr_internal_xaa11;\
 	for (cr_internal_xaa11=0;\
		cr_internal_xaa11< cr_internal_xaa12;\
		cr_internal_xaa11++)\
			int_ptr[cr_internal_xaa11] = \
			NUM2INT(RARRAY_PTR(array)[cr_internal_xaa11]);\
		}
					
#define CR_CHECK_CLASS(obj,cls){\
	if(!RTEST(rb_obj_is_kind_of(obj,cls))){\
		rb_raise(RGET_CLASS_TOP("TypeError"), "Expecting an instance of %s", StringValueCStr(cls));\
	}\
}


static VALUE cgsl; 
static VALUE ccode_runner; 
static VALUE ccode_runner_gs2; 
static VALUE ccode_runner_gs2_gsl_tensor_complexes; 
static VALUE ccode_runner_gs2_gsl_tensors; 
static VALUE cgsl_vector, cgsl_vector_complex;
/*void Init_code_runner_ext();*/
