#include "ruby.h"


#define RGET_CLASS(class_over, x) (rb_const_get(class_over, rb_intern(x)))
#define RGET_CLASS_TOP(x) (RGET_CLASS(rb_cObject, x))
#define RFCALL_10(name) (rb_funcall(self, rb_intern(name), 0))
#define RFCALL_10_ON(obj, name) (rb_funcall(obj, rb_intern(name), 0))
#define RFCALL_11(name, arg1) (rb_funcall(self, rb_intern(name), 1, arg1))
#define RFCALL_11_ON(obj, name, arg1) (rb_funcall(obj, rb_intern(name), 1, arg1))
#define CR_ELEMENT_ACCESS(recvr, key) ( \
		RFCALL_11_ON(recvr, "[]", key) \
		)
#define CR_HKS(hash, cstr) (CR_ELEMENT_ACCESS(hash, ID2SYM(rb_intern(cstr))))
#define CR_RANGE_INC(start, end) (rb_funcall(RGET_CLASS_TOP("Range"), rb_intern("new"), 2, INT2FIX(start), INT2FIX(end)))
#define CR_RANGE_EXC(start, end) (rb_funcall(RGET_CLASS_TOP("Range"), rb_intern("new"), 3, INT2FIX(start), INT2FIX(end), Qtrue))

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
					

static VALUE cgsl; 
static VALUE ccode_runner; 
static VALUE ccode_runner_gs2; 
static VALUE ccode_runner_gs2_gsl_tensor_complexes; 
static VALUE ccode_runner_gs2_gsl_tensors; 
static VALUE cgsl_vector, cgsl_vector_complex;
/*void Init_code_runner_ext();*/
