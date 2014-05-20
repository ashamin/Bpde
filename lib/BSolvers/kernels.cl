#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void explicitDerivative(__global double *H, __global double *V, __global double *mu,
        __global double *dx_l, __global double *dx_d, __global double *dx_u, 
        __global double *dy_l, __global double *dy_d, __global double *dy_u,
        __global int *I, __global int *J, __global double *Ha)
{
    int i = get_global_id(0);
    int j = get_global_id(1);
    int kT = i**J + j;
    int kH = j**I + i;
    Ha[kH] =(
             dx_l[kH]*H[kH-1] + dx_d[kH]*H[kH] + dx_u[kH]*H[kH+1] +
             dy_l[kT]*H[kH-*I] + dy_d[kT]*H[kH] + dy_u[kT]*H[kH+*I] +
             V[kH]
            )
            / mu[kH];
}

/// переписать, чтобы работало!!
__kernel void xTDMAKernel(__global double *a, __global double *b, __global double *c,
        __global double *x, __global double *d, __global int *n, __global int *step,
        __global double *loc_c, __global double *loc_d)
{
    int k = get_global_id(0);
    
    loc_c[k] = c[k] / b[k];  // c[0] - всегда 1  (-1)есть
    loc_d[k] = d[k] / b[k]; // d[0] - всегда 0 !!!

    for (int i = k+*step; i<k+*n**step; i+=*step){
        double tmp = b[i] - loc_c[i-*step]*a[i];
        loc_c[i] = c[i] / tmp;
        loc_d[i] = (d[i] - loc_d[i-*step]*a[i]) / tmp;
    }

    x[k + *n**step-*step] = loc_d[k + *n**step-*step];
    for (int i = k+(*n-2)**step; i>=0; i-=*step)
        x[i] = loc_d[i] - loc_c[i]*x[i+*step];
}