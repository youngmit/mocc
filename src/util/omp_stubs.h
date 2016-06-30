/*
   Copyright 2016 Mitchell Young

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once

#ifdef _OPENMP
#error This file should not be included if OpenMP support is enabled
#endif

inline int omp_get_thread_num()
{
    return 0;
}

inline int omp_get_num_threads()
{
    return 1;
}

inline int omp_get_max_threads()
{
    return 1;
}

inline double omp_get_wtime()
{
    return 0.0;
}

inline void omp_set_num_threads( int i )
{
    return;
}

inline int omp_get_num_procs()
{
    return 1;
}
