# distutils: language = c++

from libc.stdlib cimport malloc
cimport PyMain as pymain


def do_main(str arg1, str arg2):
	cdef char **c_argv

	args = [b'calling_from_cython'] + [bytes(arg1, encoding='utf8'), bytes(arg2, encoding='utf8')]

	c_argv = <char**>malloc(sizeof(char*) * len(args)) 
	for idx, s in enumerate(args):
		c_argv[idx] = s

	return pymain.main(len(args), c_argv)