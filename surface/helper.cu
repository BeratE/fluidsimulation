#include <stdio.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__inline__ __device__
float3 operator+(float3 a, float3 b)
{
    return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
}

__inline__ __device__
float3 operator*(float s, float3 v)
{
    return make_float3(s*v.x, s*v.y, s*v.z);
}

__inline__ __device__
float3 lerp(float3 a, float3 b, float t)
{
    return (1.0-t) * a + t * b;
}

inline __device__
uint getGlobalIdx1d1d()
{
    return blockIdx.x * blockDim.x + threadIdx.x;
}
