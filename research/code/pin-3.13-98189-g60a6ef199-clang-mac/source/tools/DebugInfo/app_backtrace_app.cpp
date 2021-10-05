/*
 * Copyright 2002-2019 Intel Corporation.
 * 
 * This software is provided to you as Sample Source Code as defined in the accompanying
 * End User License Agreement for the Intel(R) Software Development Products ("Agreement")
 * section 1.L.
 * 
 * This software and the related documents are provided as is, with no express or implied
 * warranties, other than those that are expressly stated in the License.
 */

#include <stdio.h>

extern "C" void qux()
{
    printf("qux was called\n");
    fflush(stdout);
}

extern "C" void baz()
{
}

extern "C" void bar()
{
    baz();
    qux();
}

extern "C" void foo()
{
    bar();
    qux();
}


int main()
{
    foo();
    return 0;
}
